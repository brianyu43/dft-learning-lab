"""Prepare the fixed 17-start vacuum panel only after raw seed checks pass.

Does not launch ORCA. Reuses two representative dense-grid jobs and adds 15
Opt/Freq jobs with unchanged start lineages. Other methods/solvents remain
separate downstream batches and are not silently enabled by this preparation.
"""
from __future__ import annotations

import json
from pathlib import Path

from ase.io import read

from completion_batch import verify_report
from prepare_completion_drug import GRIDS, LIBS
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_xtb_screen import geometry_audit


def verify_hashes(plan):
    for rel, expected in plan['input_hashes'].items():
        if sha(ROOT/rel)!=expected:raise ValueError(f'Changed frozen input: {rel}')


def prerequisites():
    lineage_path=ROOT/'experiments/completion/drug-start-lineage-r1.json'
    lineage=json.loads(lineage_path.read_text());verify_hashes(lineage)
    if lineage['original_start_count']!=17 or lineage['species_count']!=14:raise ValueError('Changed panel scope')
    seed_path=ROOT/'experiments/completion/geometry-precision-v3.json'
    check_path=ROOT/'experiments/completion/geometry-precision-v3-checks.json'
    if not check_path.exists():raise RuntimeError('Representative dense minima and their strict checks have not completed')
    seed=json.loads(seed_path.read_text());checks=json.loads(check_path.read_text())
    verify_hashes(seed);verify_hashes(checks)
    # The existing evaluator re-reads raw Hessians/coordinates/gradients and
    # validates the *contrast*. Re-evaluate before freezing its new digest.
    from evaluate_completion_precision import main as evaluate_raw_seed
    evaluate_raw_seed()
    review_path=ROOT/'reports/completion/geometry-precision-v3.json'
    review=json.loads(review_path.read_text())
    if review['status']!='representative_precision_passed':raise RuntimeError('Representative precision did not pass')
    if review['parent_plan_sha256']!=sha(seed_path) or review['check_plan_sha256']!=sha(check_path):
        raise ValueError('Seed review provenance mismatch')
    if review['acceptance']!=seed['acceptance']:raise ValueError('Changed acceptance policy')
    return lineage, seed, review_path


def new_job(row,review_path):
    candidate=row['candidate'];source=ROOT/row['source_geometry'];atoms=read(source)
    if sha(source)!=row['source_geometry_sha256']:raise ValueError('Changed lineage geometry')
    if atoms.numbers.tolist()!=candidate['graph']['atomic_numbers']:raise ValueError('Atom mapping mismatch')
    audit=geometry_audit(candidate,atoms)
    if not audit['heavy_connectivity_preserved'] or audit['unassigned_hydrogens'] or audit.get('separated_geometry',False):
        raise ValueError('Source chemistry needs review')
    job_id=f"w3-uniform-r4-{row['candidate_id']}-opt-freq"
    xyz=ROOT/f'structures/completion/uniform-r4/{job_id}.xyz';immutable(xyz,source.read_text())
    cp=ROOT/f"experiments/completion/uniform-r4/candidates/{row['candidate_id']}.json"
    immutable(cp,json.dumps(candidate,indent=2)+'\n')
    provenance=ROOT/f"experiments/completion/uniform-r4/lineages/{row['candidate_id']}.json"
    immutable(provenance,json.dumps(row,ensure_ascii=False,indent=2)+'\n')
    grid=GRIDS['dense']
    text=(f"! r2SCAN-3c VeryTightSCF {grid['keywords']} VeryTightOpt Opt Freq\n"
          f"%maxcore 2000\n%pal nprocs 12 end\n{grid['block']}\n* xyzfile 0 1 initial.xyz\n")
    inp=ROOT/f'inputs/orca/completion/uniform-r4/{job_id}.inp';immutable(inp,text)
    assets={f'.tools/orca-mpi/lib/{lib}':lib for lib in LIBS}
    assets.update({str(xyz.relative_to(ROOT)):'initial.xyz',str(cp.relative_to(ROOT)):'candidate.json',
        str(provenance.relative_to(ROOT)):'provenance/original-start-lineage.json',
        str(review_path.relative_to(ROOT)):'provenance/representative-precision-review.json',
        'scripts/completion_orca.py':'completion_parser_snapshot.py',
        'scripts/week3_recovery.py':'provenance/gradient_check_snapshot.py',
        'scripts/run_week2_frequency_resume.py':'provenance/hessian_check_snapshot.py',
        'scripts/run_week2_xtb_screen.py':'provenance/geometry_audit_snapshot.py',
        'scripts/prepare_completion_uniform_drug.py':'provenance/uniform_preparer_snapshot.py',
        'scripts/orca_parallel.sh':'provenance/launcher_snapshot.sh',
        'environment/orca-mpi-explicit.txt':'provenance/orca-mpi-explicit.txt'})
    cfg=dict(protocol_id='dft-completion-20260929-uniform-vacuum-r4',week=3,task_id='D01',
        job_id=job_id,engine='orca',executable='scripts/orca_parallel.sh',input=str(inp.relative_to(ROOT)),
        method='r2SCAN-3c',charge=0,multiplicity=1,solvent='vacuum',n_atoms=len(atoms),linear=False,
        requested=dict(optimization=True,frequencies=True,gradient=True),mpi_processes=12,omp_threads=1,
        timeout_seconds=14400,rss_limit_gib=40,candidate_id=row['candidate_id'],species=row['species'],
        role=row['original_role'],comparison_role=row['selection_role'],grid_setting='dense',
        expected_scf_grid={k:grid[k] for k in ('IntAcc','AngularGrid')},assets=assets,
        source_geometry=row['source_geometry'],source_geometry_sha256=sha(source),
        initial_structure_sha256=sha(xyz),candidate_sha256=sha(cp),
        completion_parser_sha256=sha(ROOT/'scripts/completion_orca.py'),expected_input_sha256=sha(inp),
        initial_geometry_audit=audit,original_config=row['original_config'],
        original_attempt_job_ids=[a['job_id'] for a in row['attempts']],
        first_DFT_attempt_for_existing_start=not bool(row['attempts']),new_independent_start=False,
        selection_uses_formation_labels=False,
        purpose='Common dense-grid vacuum Opt/Freq within the original selected start; no new pose or crystal claim')
    path=ROOT/f'experiments/completion/uniform-r4/{job_id}.json';immutable(path,json.dumps(cfg,indent=2)+'\n')
    return str(path.relative_to(ROOT))


def prepare():
    target=ROOT/'experiments/completion/uniform-vacuum-r4.json'
    if target.exists():
        verify_hashes(json.loads(target.read_text()))
        return target
    lineage,seed,review_path=prerequisites()
    jobs=[];mapping=[]
    for row in lineage['rows']:
        path=row['planned_representative_seed_config']
        if path:
            if path not in seed['jobs']:raise ValueError('Unknown representative job')
            meta=verify_report(path)
            if meta['candidate_id']!=row['candidate_id'] or meta['grid_setting']!='dense':raise ValueError('Seed mapping mismatch')
            if meta['method']!='r2SCAN-3c' or meta['solvent']!='vacuum':raise ValueError('Seed method/environment mismatch')
            reuse=True
        else:
            path=new_job(row,review_path);reuse=False
        jobs.append(path)
        mapping.append(dict(candidate_id=row['candidate_id'],species=row['species'],config=path,
                            selection_role=row['selection_role'],original_role=row['original_role'],
                            reused_representative=reuse,new_independent_start=False))
    if len(jobs)!=17 or sum(r['reused_representative'] for r in mapping)!=2:raise ValueError('Incomplete uniform panel')
    prereq_files=['experiments/completion/drug-start-lineage-r1.json',
                 'experiments/completion/geometry-precision-v3.json',
                 'experiments/completion/geometry-precision-v3-checks.json',
                 str(review_path.relative_to(ROOT)),
                 'scripts/evaluate_completion_precision.py','scripts/prepare_completion_drug.py']
    # Store these as explicit evidence as well as the ordinary job assets. The
    # batch verifies input_hashes; these are rechecked here at preparation.
    evidence={rel:sha(ROOT/rel) for rel in prereq_files}
    metadata=dict(original_start_count=17,new_independent_start_count=0,reused_seed_jobs=2,new_OptFreq_jobs=15,
        comparison_mapping=mapping,prerequisite_sha256=evidence,
        comparison_context=dict(method='r2SCAN-3c',environment='vacuum',
                                numerical_protocol_id='uniform-dense-r4',energy_quantity='electronic_energy'),
        execution_enabled_after_raw_representative_review=True,
        full_panel_numerical_precision_verified=False,
        follow_up='Strict-grid EnGrad at all accepted minima, full B/X contrast sensitivity, then separate method/solvent batches',
        max_per_job_hours=4,total_drug_budget_hours=80,
        scope='Original 17 starts / 14 species; representative pass is not whole-panel numerical certification')
    # Build the final freeze atomically, including prerequisite hashes in the
    # same map that the existing batch controller actually checks.
    hashes=dict(evidence)
    for path in jobs:
        cfg=json.loads((ROOT/path).read_text())
        for p in [path,cfg['input'],*cfg['assets']]:hashes[p]=sha(ROOT/p)
    plan=dict(kind='uniform-vacuum-r4',jobs=jobs,input_hashes=hashes,automatic_retries=0,**metadata)
    immutable(target,json.dumps(plan,indent=2)+'\n')
    return target


if __name__=='__main__':print(prepare())
