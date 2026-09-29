"""Freeze local ORCA inputs without executing them; preserve source chemistry and evidence."""
import argparse
import json
from pathlib import Path

from ase.io import read

from prepare_week3_drug import LIBS, immutable
from research_run import ROOT, sha
from run_week2_xtb_screen import geometry_audit


GRIDS={
    'g3':dict(keywords='DefGrid3',block='',AngularGrid=6,IntAcc=4.959),
    'dense':dict(keywords='DefGrid3',block='%method\n  AngularGrid 7\n  IntAcc 6.0\nend\n',AngularGrid=7,IntAcc=6.),
    'strict':dict(keywords='DefGrid3',block='%method\n  AngularGrid 7\n  IntAcc 7.0\n  GridPruning Unpruned\n  HGridReduced false\nend\n',AngularGrid=7,IntAcc=7.),
}


def candidate_record(candidate_id):
    kind='dimers' if '--' in candidate_id else 'monomers'
    path=ROOT/f'reports/week2-xtb-{kind}.json'
    rows=json.loads(path.read_text())['results']
    row=next(r for r in rows if r['candidate_id']==candidate_id)
    return row['candidate'],path


def prepare_job(job_id,candidate_id,source,*,grid,optimization,role,source_evidence):
    source=Path(source).resolve()
    assert source.is_relative_to(ROOT)
    candidate,candidate_source=candidate_record(candidate_id)
    atoms=read(source)
    audit=geometry_audit(candidate,atoms)
    assert audit['heavy_connectivity_preserved'] and not audit['unassigned_hydrogens']
    assert not audit.get('separated_geometry',False)
    assert candidate['charge']==0 and candidate['multiplicity']==1
    xyz=ROOT/f'structures/completion/{job_id}.xyz'
    immutable(xyz,source.read_text())
    cp=ROOT/f'experiments/completion/candidates/{candidate_id}.json'
    immutable(cp,json.dumps(candidate,indent=2)+'\n')
    g=GRIDS[grid]
    task='VeryTightOpt Opt Freq' if optimization else 'EnGrad'
    text=f"! r2SCAN-3c VeryTightSCF {g['keywords']} {task}\n%maxcore 2000\n%pal nprocs 12 end\n{g['block']}\n* xyzfile 0 1 initial.xyz\n"
    inp=ROOT/f'inputs/orca/completion/{job_id}.inp';immutable(inp,text)
    assets={f'.tools/orca-mpi/lib/{lib}':lib for lib in LIBS}
    assets.update({str(xyz.relative_to(ROOT)):'initial.xyz',str(cp.relative_to(ROOT)):'candidate.json',
        'scripts/completion_orca.py':'completion_parser_snapshot.py',
        'scripts/week3_recovery.py':'provenance/gradient_check_snapshot.py',
        'scripts/run_week2_frequency_resume.py':'provenance/hessian_check_snapshot.py',
        'scripts/run_week2_xtb_screen.py':'provenance/geometry_audit_snapshot.py',
        'scripts/prepare_completion_drug.py':'provenance/preparer_snapshot.py',
        'scripts/orca_parallel.sh':'provenance/launcher_snapshot.sh',
        'environment/orca-mpi-explicit.txt':'provenance/orca-mpi-explicit.txt'})
    for i,path in enumerate(source_evidence):
        path=Path(path).resolve();assert path.is_relative_to(ROOT)
        assets[str(path.relative_to(ROOT))]=f'provenance/evidence-{i:02d}-{path.name}'
    cfg=dict(protocol_id='dft-completion-20260929-geometry-precision-v3',week=3,task_id='C02',
        job_id=job_id,engine='orca',executable='scripts/orca_parallel.sh',input=str(inp.relative_to(ROOT)),
        method='r2SCAN-3c',charge=0,multiplicity=1,solvent='vacuum',n_atoms=len(atoms),linear=False,
        requested=dict(optimization=optimization,frequencies=optimization,gradient=True),
        mpi_processes=12,omp_threads=1,timeout_seconds=14400 if optimization else 1800,rss_limit_gib=40,
        candidate_id=candidate_id,species=candidate['species'],role=role,grid_setting=grid,
        expected_scf_grid={k:g[k] for k in ['IntAcc','AngularGrid']},assets=assets,
        source_geometry=str(source.relative_to(ROOT)),source_geometry_sha256=sha(source),
        initial_structure_sha256=sha(xyz),candidate_sha256=sha(cp),candidate_source=str(candidate_source.relative_to(ROOT)),
        candidate_source_sha256=sha(candidate_source),completion_parser_sha256=sha(ROOT/'scripts/completion_orca.py'),
        expected_input_sha256=sha(inp),initial_geometry_audit=audit,
        new_independent_start=False,selection_uses_formation_labels=False,
        purpose='Numerical geometry refinement or fixed-geometry sensitivity; no extra independent pose or crystal claim')
    path=ROOT/f'experiments/completion/{job_id}.json';immutable(path,json.dumps(cfg,indent=2)+'\n')
    return str(path.relative_to(ROOT))


def freeze_plan(name,jobs,metadata):
    hashes={}
    for p in jobs:
        cfg=json.loads((ROOT/p).read_text())
        for q in [p,cfg['input'],*cfg['assets']]:hashes[q]=sha(ROOT/q)
    plan=dict(kind=name,jobs=jobs,input_hashes=hashes,automatic_retries=0,**metadata)
    path=ROOT/f'experiments/completion/{name}.json'
    immutable(path,json.dumps(plan,indent=2)+'\n')
    return path


def seed_precision():
    diagnostic=ROOT/'reports/completion/precision-diagnostic.json'
    data=json.loads(diagnostic.read_text())
    assert data['original_failed_gate_preserved'] and data['g3_vs_dense_contrast_within_target']
    oldplan=json.loads((ROOT/'experiments/completion/precision-diagnostic.json').read_text())
    mapping={'dif':'diflunisal-t02-s43','dif2':'diflunisal--diflunisal-p05'}
    jobs=[]
    for species,candidate in mapping.items():
        original=oldplan['fixed_sources'][species]
        source=ROOT/original['path'];assert sha(source)==original['sha256']
        jobs.append(prepare_job(f'w3-precision-v3-{species}-dense-opt-freq',candidate,source,
            grid='dense',optimization=True,role='representative_dense_grid_minimum',
            source_evidence=[diagnostic,ROOT/'docs/NUMERICAL_REVIEW_V2.ko.md']))
    return freeze_plan('geometry-precision-v3',jobs,dict(
        sources=oldplan['fixed_sources'],new_independent_start_count=0,
        parent_diagnostic_sha256=sha(diagnostic),
        follow_up='After verified dense minima, EnGrad at identical final coordinates with AngularGrid7/IntAcc7/unpruned/no H reduction; compare with dense gradients and self-binding energy.',
        acceptance=dict(max_self_binding_change_kcal_mol=.1,max_gradient_component_change_Hartree_Bohr=1e-4,
                        all_internal_frequencies_positive=True,geometry_and_chemistry_verified=True),
        scope='Representative geometry/grid decision before a uniform full panel; not proof for other species',
        panel_execution_enabled=False))


def seed_checks():
    from completion_batch import verify_report
    parent=ROOT/'experiments/completion/geometry-precision-v3.json'
    plan=json.loads(parent.read_text());jobs=[];parents=[]
    for path in plan['jobs']:
        meta=verify_report(path)
        assert meta['result']['minimum_verified'] and meta['result']['gradient_geometry_verified']
        source=Path(meta['run_dir'])/(Path(meta['input']).stem+'.xyz')
        assert sha(source)==meta['result']['optimized_xyz_sha256']
        report=ROOT/f"reports/{meta['job_id']}.json"
        species='dif2' if '--' in meta['candidate_id'] else 'dif'
        jobs.append(prepare_job(f'w3-precision-v3-{species}-strict-engrad',meta['candidate_id'],source,
            grid='strict',optimization=False,role='strict_grid_at_dense_minimum',source_evidence=[parent,report]))
        parents.append(dict(config=path,report=str(report.relative_to(ROOT)),report_sha256=sha(report),
                            geometry=str(source.relative_to(ROOT)),geometry_sha256=sha(source)))
    return freeze_plan('geometry-precision-v3-checks',jobs,dict(parents=parents,
        parent_plan_sha256=sha(parent),acceptance=plan['acceptance'],panel_execution_enabled=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['seed-precision','seed-checks'])
    args=parser.parse_args()
    print(seed_precision() if args.stage=='seed-precision' else seed_checks())
