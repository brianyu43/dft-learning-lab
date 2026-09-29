"""Independently compare ORCA log, Hessian, XYZ, and recorded candidate graph."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import numpy as np
from ase.io import read
from run_week2_xtb_screen import geometry_audit
from run_week2_frequency_resume import inspect_hessian
from week2_resume_provenance import validate_checkpoint

ROOT=Path(__file__).resolve().parents[1]
HARTREE_KCAL=627.509474
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

selection=json.loads((ROOT/'reports/week2-first-dft-selection.json').read_text())
out=ROOT/'structures/optimized/week2';out.mkdir(parents=True,exist_ok=True)
rows=[]
for row in selection['selected']:
    name=row['species'];job=f'w2-{name}-opt-freq';report=ROOT/f'reports/{job}.json'
    resumed=False;original=None
    resume_config=ROOT/'experiments/week2/frequency-resume.json'
    if name=='diflunisal--caffeine' and resume_config.exists():
        resume=json.loads(resume_config.read_text());newreport=ROOT/f"reports/{resume['job_id']}.json"
        if newreport.exists():
            original=json.loads(report.read_text());job=resume['job_id'];report=newreport;resumed=True
    if not report.exists():continue
    m=json.loads(report.read_text());run=Path(m['run_dir']);stem=Path(m['input']).stem
    if 'output_sha256' in m:
        assert sha(run/'output.log')==m['output_sha256']
    else:
        interrupted=validate_checkpoint(ROOT/'experiments/week2/pause-checkpoint.json')
        assert str(run.relative_to(ROOT))==interrupted['source_run_dir'] and m['stop_reason']
        assert sha(run/'output.log')==interrupted['optimization_log_sha256']
    assert sha(run/Path(m['input']).name)==m['source_sha256'][m['input']]
    assert sha(run/'runner_snapshot.py')==m['runner_sha256']
    r={'species':name,'job_id':job,'run_dir':str(run.relative_to(ROOT)),'engine_status':m['status'],
       'candidate_id':row['candidate_id'],'elapsed_seconds':m['elapsed_seconds'],
       'peak_sampled_process_tree_rss_bytes':m['peak_sampled_process_tree_rss_bytes'],
       'minimum_verified':False,'energy_Hartree':m['result'].get('energy')}
    hess=run/f'{stem}.hess';xyz=run/f'{stem}.xyz'
    if resumed:
        proof=validate_checkpoint(ROOT/'experiments/week2/pause-checkpoint.json')
        assert m['optimization_provenance']==proof
        assert sha(run/'resume_parser_snapshot.py')==m['result']['resume_parser_sha256']
        xyz=run/'saved-optimization.xyz';assert sha(xyz)==proof['geometry_sha256']
        assert sha(run/'provenance/optimization-output.log')==proof['optimization_log_sha256']
        r.update(elapsed_seconds=original['elapsed_seconds']+m['elapsed_seconds'],
                 timing_scope='Original serial optimization plus interrupted partial frequency work and resumed parallel Freq; not a continuous same-core Opt+Freq timing',
                 original_interrupted_wall_seconds=original['elapsed_seconds'],resumed_frequency_wall_seconds=m['elapsed_seconds'],
                 original_run_dir=original['run_dir'],frequency_MPI_processes=m['mpi_processes'],
                 optimization_provenance=proof,
                 peak_sampled_process_tree_rss_bytes=max(original['peak_sampled_process_tree_rss_bytes'],m['peak_sampled_process_tree_rss_bytes']))
    if hess.exists() and xyz.exists() and m['result']['frequency_analysis_complete']:
        atoms=read(xyz);assert len(atoms)==m['n_atoms'] and np.isfinite(atoms.positions).all()
        vib=inspect_hessian(hess,atoms,m['result']['frequencies_cm1'])
        frequencies=np.r_[np.zeros(6),vib]
        energy=float(re.search(r' E ([-\d.]+)',xyz.read_text().splitlines()[1]).group(1))
        energy_tolerance=1e-6 if resumed else 1e-10
        assert abs(energy-r['energy_Hartree'])<energy_tolerance
        audit=geometry_audit(row['candidate'],atoms)
        initial_xyz=ROOT/row['run_dir']/'xtbopt.xyz'
        assert sha(initial_xyz)==(original if resumed else m)['initial_structure_sha256']
        initial=read(initial_xyz);heavy=atoms.numbers>1
        p=initial.positions[heavy]-initial.positions[heavy].mean(axis=0)
        q=atoms.positions[heavy]-atoms.positions[heavy].mean(axis=0)
        u,_,vt=np.linalg.svd(p.T@q);rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
        rmsd=float(np.sqrt(np.mean(np.sum((p@rotation-q)**2,axis=1))))
        minimum=m['result']['minimum_verified'] is True and bool(np.all(frequencies[6:]>0)) and m['status']=='passed'
        target=out/f'{name}.xyz';shutil.copy2(xyz,target)
        text=(run/'output.log').read_text(errors='replace')
        timings={}
        for line in text.splitlines():
            match=re.match(r'^([A-Za-z][A-Za-z -]+?)\s+\.\.\.\s+([\d.]+) sec',line)
            if match:timings[match.group(1).strip()]=float(match.group(2))
        r.update(minimum_verified=minimum,vibrational_modes=len(frequencies)-6,
                 lowest_vibrational_frequency_cm1=float(min(frequencies[6:])),
                 imaginary_frequencies_cm1=frequencies[6:][frequencies[6:]<0].tolist(),
                 geometry_audit=audit,optimized_xyz=str(target.relative_to(ROOT)),optimized_xyz_sha256=sha(target),
                 input_xtb_to_DFT_heavy_atom_RMSD_angstrom=rmsd,
                 RMSD_scope='fixed atom mapping, proper rigid alignment of the whole system; no symmetry remapping or separate fragment alignment',
                 hessian_sha256=sha(hess),log_xyz_energy_matches=True,log_xyz_energy_tolerance_Hartree=energy_tolerance,
                 frequency_max_log_hessian_difference_cm1=float(max(abs(frequencies[6:]-m['result']['frequencies_cm1']))),
                 ORCA_reported_module_times_seconds=timings)
    rows.append(r)
result={'status':'partial','completed':len(rows),'expected':3,'checks':rows,'binding_energy_kcal_mol':None,
        'competitive_self_association_X_available':False,'cocrystal_formation_predicted':False,
        'scope':'one gas-phase dimer minimum from xTB ranking; electronic binding energy only, no thermal/solvent/crystal free energy; interrupted Opt and new Freq are independently linked when resumed'}
if len(rows)==3:
    good=all(r['minimum_verified'] and r.get('geometry_audit',{}).get('heavy_connectivity_preserved')
             and not r['geometry_audit'].get('unassigned_hydrogens')
             and not r['geometry_audit'].get('separated_geometry') for r in rows)
    result['status']='passed' if good else 'requires_geometry_or_frequency_review'
    if good:
        byname={r['species']:r for r in rows}
        a=byname['diflunisal'];b=byname['caffeine'];ab=byname['diflunisal--caffeine']
        result['binding_energy_kcal_mol']=(ab['energy_Hartree']-a['energy_Hartree']-b['energy_Hartree'])*HARTREE_KCAL
        result['binding_energy_definition']='B = E(AB optimized) - E(A optimized) - E(B optimized), same r2SCAN-3c model'
(ROOT/'reports/week2-drug-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
