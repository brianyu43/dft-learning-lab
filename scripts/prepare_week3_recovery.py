"""One reviewed refinement of DIF-DIF p05, preserving the failed attempt."""
import json
from pathlib import Path

from ase.io import read
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_xtb_screen import geometry_audit


def prepare():
    parent_path = 'experiments/week3/drug/w3-diflunisal--diflunisal-p05-opt-freq.json'
    parent = json.loads((ROOT / parent_path).read_text())
    report_path = f"reports/{parent['job_id']}.json"
    report = json.loads((ROOT / report_path).read_text())
    run = Path(report['run_dir'])
    assert report['status'] == 'needs_review' and report['exit_code'] == 0 and not report['stop_reason']
    assert report['config_sha256'] == sha(ROOT / parent_path)
    assert report['result']['geometry_optimization_converged']
    assert report['result']['imaginary_frequencies_cm1'] == [-2.23]
    inspection_path = 'reports/week3-recovery-r1/mode-inspection.json'
    inspection = json.loads((ROOT / inspection_path).read_text())
    assert inspection['parent_report_sha256'] == sha(ROOT / report_path)
    source = run / (parent['job_id'] + '.xyz')
    assert sha(source) == inspection['geometry_sha256']
    assert sha(run / (parent['job_id'] + '.hess')) == inspection['hessian_sha256']
    candidate_path = f"experiments/week3/drug/candidates/{parent['candidate_id']}.json"
    candidate = json.loads((ROOT / candidate_path).read_text())
    audit = geometry_audit(candidate, read(source))
    assert audit['chemical_state_preserved_by_geometry'] and not audit['separated_geometry']
    job = 'w3-diflunisal--diflunisal-p05-refine-r1-opt-freq'
    initial = ROOT / f'structures/week3/drug_initial/{job}.xyz'
    immutable(initial, source.read_text())
    inp = ROOT / f'inputs/orca/week3/{job}.inp'
    immutable(inp, '! r2SCAN-3c VeryTightSCF VeryTightOpt DefGrid3 Opt Freq\n%maxcore 2000\n%pal nprocs 12 end\n\n* xyzfile 0 1 initial.xyz\n')
    cfg = dict(parent)
    cfg['assets'] = {k:v for k,v in parent['assets'].items() if v != 'initial.xyz'}
    cfg['assets'].update({str(initial.relative_to(ROOT)):'initial.xyz',
                          report_path:'provenance/failed-parent-report.json',
                          inspection_path:'provenance/mode-inspection.json',
                          str((run/'output.log').relative_to(ROOT)):'provenance/failed-parent-output.log',
                          str((run/(parent['job_id']+'.hess')).relative_to(ROOT)):'provenance/failed-parent.hess'})
    cfg.update(job_id=job,input=str(inp.relative_to(ROOT)),source_geometry=str(source.relative_to(ROOT)),
               source_geometry_sha256=sha(source),initial_structure_sha256=sha(initial),
               initial_geometry_audit=audit,source_stage='Undisplaced final geometry from reviewed non-minimum attempt',
               refinement_attempt=1,parent_job_id=parent['job_id'],
               parent_report_sha256=sha(ROOT/report_path),
               numerical_settings='VeryTightSCF VeryTightOpt DefGrid3',
               numerical_comparability_check_required=True,
               purpose='One reviewed tighter numerical refinement of the same starting candidate; not an additional pose')
    path = ROOT / f'experiments/week3/drug/{job}.json'
    immutable(path,json.dumps(cfg,indent=2)+'\n')
    # A new plan supplements, rather than mutates, the original 14-start plan.
    plan = dict(schema=1,parent_panel_plan_sha256=sha(ROOT/'experiments/week3/drug/panel-plan.json'),
                automatic_retries=0,reviewed_refinement_attempts=1,new_independent_start_count=0,
                inspection=inspection_path,inspection_sha256=sha(ROOT/inspection_path),
                common_grid_policy=dict(max_absolute_energy_delta_kcal_mol=.1,
                    max_gradient_component_delta_Hartree_Bohr=1e-4,
                    ranking_energy='Same-geometry original TightSCF/default-grid single point after refined minimum validation',
                    sensitivity_scope='Local numerical sensitivity only, not a bound for all candidates or model error'),
                replacements=[dict(parent_config=parent_path,parent_config_sha256=sha(ROOT/parent_path),
                    parent_report=report_path,parent_report_sha256=sha(ROOT/report_path),
                    replacement_config=str(path.relative_to(ROOT)),replacement_config_sha256=sha(path),
                    common_grid_job_id='w3-diflunisal--diflunisal-p05-refine-r1-common-grid',
                    common_grid_check_report='reports/week3-recovery-r1/numerical-check.json')])
    plan_path = ROOT/'experiments/week3/drug/recovery-plan-r1.json'
    immutable(plan_path,json.dumps(plan,indent=2)+'\n')
    print(str(path.relative_to(ROOT)))
    print(str(plan_path.relative_to(ROOT)))


if __name__ == '__main__':
    prepare()
