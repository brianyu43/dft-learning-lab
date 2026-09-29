"""Freeze the original 17 starts and verify source coordinates without new DFT.

Preserves interrupted, failed, refined and frequency-only runs as separate
attempts. Historical minima are not certified for the new common grid.
"""
from __future__ import annotations

import collections
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

from run_week2_frequency_resume import inspect_hessian
from run_week2_xtb_screen import geometry_audit
from week2_resume_provenance import validate_checkpoint

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())


class Evidence:
    def __init__(self):
        self.hashes={}
        p='experiments/completion/preserved-original-evidence.json'
        self.hashes[p]=sha(ROOT/p)
        self.preserved=load(p)['files']

    def add(self,path,expected=None):
        p=(ROOT/path).resolve()
        if not p.is_relative_to(ROOT):raise ValueError('Evidence outside project')
        rel=str(p.relative_to(ROOT));actual=sha(p)
        if expected is not None and actual!=expected:raise ValueError(f'Hash mismatch: {rel}')
        if rel in self.preserved and actual!=self.preserved[rel]:raise ValueError(f'Changed preserved evidence: {rel}')
        if rel in self.hashes and actual!=self.hashes[rel]:raise ValueError(f'Changed during audit: {rel}')
        self.hashes[rel]=actual
        return actual


def attempt(config_path,candidate,evidence,*,verified_week2=None):
    cfg=load(config_path);evidence.add(config_path)
    report=f"reports/{cfg['job_id']}.json"
    if not (ROOT/report).exists():return None
    evidence.add(report);meta=load(report);run=Path(meta['run_dir']);result=meta['result']
    if meta['config_sha256']!=sha(ROOT/config_path):raise ValueError('Config/report mismatch')
    evidence.add(run/'job_config.json',meta['config_sha256'])
    evidence.add(run/'runner_snapshot.py',meta['runner_sha256'])
    for source,digest in meta['source_sha256'].items():
        dest=Path(meta['input']).name if source==meta['input'] else meta['assets'][source]
        evidence.add(run/dest,digest)
    proof=None
    if 'output_sha256' in meta:
        evidence.add(run/'output.log',meta['output_sha256'])
    else:
        proof=validate_checkpoint(ROOT/'experiments/week2/pause-checkpoint.json')
        if run!=(ROOT/proof['source_run_dir']).resolve():raise ValueError('Unexplained missing output digest')
        evidence.add(run/'output.log',proof['optimization_log_sha256'])
    if 'optimization_provenance' in meta:
        proof=validate_checkpoint(ROOT/'experiments/week2/pause-checkpoint.json')
        if meta['optimization_provenance']!=proof:raise ValueError('Resumed optimization provenance mismatch')
        for source,digest in load('experiments/week2/pause-checkpoint.json')['source_files_sha256'].items():evidence.add(source,digest)
        evidence.add(proof['checkpoint'],proof['checkpoint_sha256'])
        evidence.add(proof['optimized_geometry'],proof['geometry_sha256'])
        xyz=run/'saved-optimization.xyz'
        evidence.add(xyz,proof['geometry_sha256'])
    else:
        xyz=run/(Path(meta['input']).stem+'.xyz')
    atoms=read(xyz)
    if len(atoms)!=cfg['n_atoms'] or atoms.numbers.tolist()!=candidate['graph']['atomic_numbers']:
        raise ValueError('Candidate/optimized atom mapping changed')
    if not np.isfinite(atoms.positions).all():raise ValueError('Invalid coordinates')
    digest=evidence.add(xyz,result.get('optimized_xyz_sha256'))
    if verified_week2:
        evidence.add(verified_week2['optimized_xyz'],verified_week2['optimized_xyz_sha256'])
        if digest!=verified_week2['optimized_xyz_sha256']:raise ValueError('Week-two coordinate copy changed')
    audit=geometry_audit(candidate,atoms)
    if not audit['heavy_connectivity_preserved'] or audit['unassigned_hydrogens'] or audit.get('separated_geometry',False):
        raise ValueError('Review source geometry before using it')
    log=(run/'output.log').read_text()
    energies=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',log)
    if not energies or abs(float(energies[-1])-result['energy'])>1e-10:raise ValueError('Log/report energy mismatch')
    minimum=False;vib=None;hessian=None
    if result.get('frequency_analysis_complete'):
        hess=run/(Path(meta['input']).stem+'.hess')
        expected=result.get('hessian_sha256') or (verified_week2 or {}).get('hessian_sha256')
        if expected is None and str(hess.relative_to(ROOT)) not in evidence.preserved:
            raise ValueError('Hessian has no preexisting digest')
        hessian=evidence.add(hess,expected)
        vib=inspect_hessian(hess,atoms,result['frequencies_cm1'])
        minimum=bool(meta['status']=='passed' and result.get('minimum_verified') and np.all(vib>0))
        if minimum and ('ORCA TERMINATED NORMALLY' not in log or not result.get('scf_converged')):
            raise ValueError('Missing engine completion evidence')
        if minimum and not proof and 'THE OPTIMIZATION HAS CONVERGED' not in log:
            raise ValueError('No converged optimization provenance')
    if bool(result.get('minimum_verified'))!=minimum:raise ValueError('Minimum certificate could not be reproduced')
    return dict(job_id=meta['job_id'],config=config_path,report=report,status=meta['status'],
        source_xyz=str(xyz.relative_to(ROOT)),source_xyz_sha256=digest,
        archived_minimum_verified_under_own_settings=minimum,
        lowest_internal_frequency_cm1=float(min(vib)) if vib is not None else None,
        internal_frequency_count=len(vib) if vib is not None else 0,hessian_sha256=hessian,
        energy_Hartree=result['energy'],geometry_audit=audit,
        elapsed_seconds=meta['elapsed_seconds'],optimization_provenance=proof,
        eligible_for_new_common_grid_comparison=False)


def build():
    evidence=Evidence();rows=[]
    source_files=['reports/week2-first-dft-selection.json','reports/week2-drug-verification.json',
        'experiments/week3/drug/panel-plan.json','experiments/week3/drug/recovery-plan-r1.json',
        'reports/week3-recovery-r1/numerical-check.json','experiments/completion/geometry-precision-v3.json']
    for p in source_files:evidence.add(p)
    previous={r['candidate_id']:r for r in load(source_files[1])['checks']}
    for old in load(source_files[0])['selected']:
        cfgpath=f"experiments/week2/w2-{old['species']}-opt-freq.json"
        candidate=old['candidate'];attempts=[]
        resumed=old['species']=='diflunisal--caffeine'
        attempts.append(attempt(cfgpath,candidate,evidence,verified_week2=None if resumed else previous[old['candidate_id']]))
        if resumed:attempts.append(attempt('experiments/week2/frequency-resume.json',candidate,evidence,verified_week2=previous[old['candidate_id']]))
        rows.append(dict(candidate_id=old['candidate_id'],species=old['species'],candidate=candidate,
                         original_config=cfgpath,selection_role='primary',original_role='primary',
                         original_week=2,attempts=attempts))
    plan=load(source_files[2]);recovery=load(source_files[3])
    if recovery['parent_panel_plan_sha256']!=sha(ROOT/source_files[2]):raise ValueError('Recovery plan parent changed')
    replacements={r['parent_config']:r for r in recovery['replacements']}
    for cfgpath in [plan['first_job'],*plan['remaining_jobs']]:
        cfg=load(cfgpath);evidence.add(cfgpath,plan['input_hashes'][cfgpath])
        cp=next(p for p,dest in cfg['assets'].items() if dest=='candidate.json')
        candidate=load(cp);evidence.add(cp,cfg['candidate_sha256'])
        initial=next(p for p,dest in cfg['assets'].items() if dest=='initial.xyz')
        evidence.add(initial,cfg['initial_structure_sha256'])
        evidence.add(cfg['source_geometry'],cfg['source_geometry_sha256'])
        if sha(ROOT/initial)!=sha(ROOT/cfg['source_geometry']):raise ValueError('Changed initial/source copy')
        initial_audit=geometry_audit(candidate,read(ROOT/initial))
        if initial_audit!=cfg['initial_geometry_audit']:raise ValueError('Initial chemistry audit changed')
        found=attempt(cfgpath,candidate,evidence);attempts=[] if found is None else [found]
        if cfgpath in replacements:
            rep=replacements[cfgpath]
            evidence.add(rep['parent_report'],rep['parent_report_sha256'])
            evidence.add(rep['replacement_config'],rep['replacement_config_sha256'])
            attempts.append(attempt(rep['replacement_config'],candidate,evidence))
        rows.append(dict(candidate_id=cfg['candidate_id'],species=cfg['species'],candidate=candidate,
                         original_config=cfgpath,selection_role='primary' if cfg['role']=='primary' else 'alternative',
                         original_role=cfg['role'],original_week=3,attempts=attempts,
                         original_initial_xyz=initial,initial_geometry_audit=initial_audit))
    seeds={load(p)['candidate_id']:p for p in load(source_files[5])['jobs']}
    for row in rows:
        minima=[r for r in row['attempts'] if r['archived_minimum_verified_under_own_settings']]
        if minima:
            selected=minima[-1]
            row['source_geometry']=selected['source_xyz'];row['source_geometry_sha256']=selected['source_xyz_sha256']
            row['source_reason']='latest reviewed minimum within the original start lineage; not selected across starts by a formation label'
        else:
            row['source_geometry']=row['original_initial_xyz'];row['source_geometry_sha256']=sha(ROOT/row['source_geometry'])
            row['source_reason']='original frozen starting coordinates; no DFT attempt completed for this start'
        row['planned_representative_seed_config']=seeds.get(row['candidate_id'])
        row['new_independent_start']=False
        row['common_protocol_minimum_verified']=False
        evidence.add(row['source_geometry'],row['source_geometry_sha256'])
    counts=collections.Counter(r['species'] for r in rows)
    if len(rows)!=17 or len({r['candidate_id'] for r in rows})!=17 or len(counts)!=14 or max(counts.values())!=3:
        raise ValueError('Original 17 starts / 14 species / per-species cap changed')
    if collections.Counter(r['selection_role'] for r in rows)!=dict(primary=14,alternative=3):raise ValueError('Primary/alternative scope changed')
    for species in counts:
        if sum(r['species']==species and r['selection_role']=='primary' for r in rows)!=1:raise ValueError('One primary start per species required')
    ionic=next(r for r in rows if r['candidate_id']=='diflunisal--adenine-p36')
    if not ionic['initial_geometry_audit'].get('cross_fragment_proton_transfer'):raise ValueError('Direct ionic start lost its proton transfer')
    old_failure=load(source_files[4])
    if old_failure['status']=='passed':raise ValueError('Preserve the original failed numerical check')
    return dict(version='drug-start-lineage-r1',status='source_inventory_verified_no_new_DFT',
        original_start_count=len(rows),species_count=len(counts),starts_per_species=dict(counts),
        original_minimum_lineages=sum(any(a['archived_minimum_verified_under_own_settings'] for a in r['attempts']) for r in rows),
        starts_without_completed_DFT=sum(not r['attempts'] for r in rows),
        historical_attempt_count=sum(len(r['attempts']) for r in rows),
        new_independent_starts=0,rows=rows,input_hashes=evidence.hashes,
        script_sha256=sha(Path(__file__)),
        next_gate='representative_precision_passed plus raw proof before freezing the uniform full panel',
        full_panel_execution_enabled=False,rankings_computed=False,
        scope='Lineage and archived source validation; not common-grid minimum verification or a polymorph inventory')


def main():
    data=build();path=ROOT/'experiments/completion/drug-start-lineage-r1.json'
    output=json.dumps(data,ensure_ascii=False,indent=2)+'\n'
    if path.exists() and path.read_text()!=output:raise ValueError('Preserve existing lineage; version any change')
    if not path.exists():path.write_text(output)
    print(json.dumps({k:data[k] for k in ['status','original_start_count','species_count','original_minimum_lineages',
                                         'starts_without_completed_DFT','historical_attempt_count','new_independent_starts']},indent=2))
    print('source_files',len(data['input_hashes']),'sha256',sha(path))


if __name__=='__main__':main()
