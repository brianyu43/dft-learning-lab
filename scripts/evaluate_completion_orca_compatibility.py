"""Read-only review of the two water runtime checks, separate from pilot scores.

Output-format checks are provisional until exercised on the actual new jobs.
This never calls the historical parser, which would overwrite result.json.
"""
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read

from completion_batch import verify_report
from completion_drug_panel import require, verify_hashes
from prepare_completion_orca_compatibility import PLAN
from prepare_week3_drug import immutable
from research_run import ROOT, sha, now
from run_week2_frequency_resume import inspect_hessian
from run_week2_xtb_screen import geometry_audit
from week3_recovery import gradient_in_frame

NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'
HARTREE_KCAL = 627.509474


def numbers(pattern, text):
    values = [float(x.replace('D', 'E').replace('d', 'e'))
              for x in re.findall(pattern, text, re.I | re.M)]
    require(values and np.isfinite(values).all(), 'Missing or nonfinite output: ' + pattern)
    return values


def method_evidence(log, method):
    require('****END OF INPUT****' in log, 'Missing end of input marker')
    runtime = log.split('****END OF INPUT****', 1)[1]
    token = {'r2SCAN-3c': r'R2SCAN-?3C', 'wB97X-3c': r'WB97X-?3C'}[method]
    recognized = re.findall(r'^The\s+' + token + r'\s+composite method is recognized\s*$',
                            runtime, re.I | re.M)
    require(recognized, 'No runtime composite-method recognition (credits/input do not count)')
    basis = re.findall(r'^Your calculation utilizes the basis:\s*(\S+)', log, re.M)
    expected = 'def2-mTZVPP' if method == 'r2SCAN-3c' else 'vDZP'
    require(basis and basis[0].lower() == expected.lower(), 'Unexpected orbital basis')
    require(re.search(r'\bDFTD4\b', runtime), 'No actual D4 output')
    return dict(method=method, orbital_basis=basis[0], runtime_recognition_count=len(recognized),
                detailed_parameter_and_ecp_review_required=True)


def smd_energy_evidence(log):
    """Check the final SCF/CDS/composite sum without adding CDS twice."""
    require('****END OF INPUT****' in log, 'Missing end of input marker')
    runtime = log.split('****END OF INPUT****', 1)[1]
    solvents = re.findall(r'^\s*Solvent\s+\.\.\.\s*(\S+)', runtime, re.M | re.I)
    require(solvents and {s.upper() for s in solvents} == {'ETHANOL'}, 'Actual solvent is not ethanol')
    require('SMD-CDS' in runtime and 'CPCM SOLVATION MODEL' in runtime, 'Missing SMD/CPCM output')
    final = list(re.finditer(r'FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')', runtime))
    require(final, 'Missing final energy')
    prefix = runtime[:final[-1].start()]
    require('TOTAL SCF ENERGY' in prefix, 'Missing SCF energy decomposition')
    block = prefix.rsplit('TOTAL SCF ENERGY', 1)[1]
    total = numbers(r'^\s*Total Energy\s*:\s*(' + NUMBER + r')\s+Eh', block)[-1]
    cds = numbers(r'^\s*SMD CDS \(Gcds\)\s*:\s*(' + NUMBER + r')\s+Eh', block)[-1]
    dielectric = numbers(r'^\s*CPCM Dielectric\s*:\s*(' + NUMBER + r')\s+Eh', block)[-1]
    corrected = numbers(r'Total Energy after SMD CDS correction\s*=\s*(' + NUMBER + r')\s+Eh', prefix)[-1]
    cds_kcal = numbers(r'SMD CDS free energy correction energy\s*:\s*(' + NUMBER + r')\s+Kcal/mol', prefix)[-1]
    dispersion = numbers(r'^\s*Dispersion correction\s+(' + NUMBER + ')\s*$', block)[-1]
    gcp = numbers(r'^\s*gCP correction\s+(' + NUMBER + ')\s*$', block)[-1]
    final_energy = float(final[-1].group(1).replace('D', 'E'))
    require(abs(total - corrected) < 2e-8, 'CDS-corrected SCF total disagrees')
    require(abs(cds - cds_kcal / HARTREE_KCAL) < 2e-8, 'CDS unit conversion disagrees')
    require(abs(final_energy - (total + dispersion + gcp)) < 2e-8,
            'Composite energy sum disagrees; CDS must not be added twice')
    return dict(energy_quantity='orca_smd_total_energy_including_cds',
                total_scf_including_cds_Hartree=total, cds_Hartree=cds,
                dielectric_Hartree=dielectric, dispersion_Hartree=dispersion,
                gcp_Hartree=gcp, final_energy_Hartree=final_energy,
                cds_added_again=False, molecular_formation_free_energy=False)


def audit_job(config_path):
    cfg = json.loads((ROOT / config_path).read_text())
    meta = verify_report(config_path)
    require(meta['exit_code'] == 0 and not meta.get('stop_reason'), 'Incomplete engine run')
    require(all(meta.get(key) == value for key, value in cfg.items()), 'Config/report mismatch')
    run = Path(meta['run_dir'])
    require(sha(run / 'job_config.json') == meta['config_sha256'], 'Archived config mismatch')
    require(sha(run / 'runner_snapshot.py') == meta['runner_sha256'], 'Runner snapshot mismatch')
    log = (run / 'output.log').read_text()
    require('ORCA TERMINATED NORMALLY' in log and 'SCF CONVERGED' in log
            and 'SCF NOT CONVERGED' not in log, 'Raw termination/SCF failure')
    require(re.search(r'Program Version\s+6\.1\.1\b', log), 'Unreviewed ORCA version')
    require(set(numbers(r'running with\s+(\d+)\s+parallel MPI-processes', log)) == {12}, 'MPI count mismatch')
    runtime = log.split('****END OF INPUT****', 1)[-1]
    require(set(numbers(r'Total Charge\s+Charge\s+\.+\s*(' + NUMBER + ')', runtime)) == {0}, 'Charge mismatch')
    require(set(numbers(r'Multiplicity\s+Mult\s+\.+\s*(' + NUMBER + ')', runtime)) == {1}, 'Multiplicity mismatch')
    require(numbers(r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s*(' + NUMBER + ')', runtime)[0] == 6., 'SCF radial grid mismatch')
    require(numbers(r'Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.\.\.\s*(\d+)', runtime)[0] == 7, 'SCF angular grid mismatch')
    energy = numbers(r'FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')', runtime)[-1]
    require(meta['result']['energy_unit'] == 'Hartree' and abs(energy - meta['result']['energy']) < 1e-9,
            'Raw/report energy disagreement')
    minimum = cfg['requested']['optimization']
    stem = Path(cfg['input']).stem
    geometry = run / (stem + '.xyz') if minimum else run / 'initial.xyz'
    initial, atoms = read(run / 'initial.xyz'), read(geometry)
    require(initial.numbers.tolist() == atoms.numbers.tolist() == [8, 1, 1]
            and np.isfinite(atoms.positions).all(), 'Water atom order/geometry mismatch')
    require(sha(run / 'initial.xyz') == cfg['initial_structure_sha256'], 'Initial geometry mismatch')
    require(sha(run / 'candidate.json') == cfg['candidate_sha256'], 'Candidate mismatch')
    require(geometry_audit(json.loads((run / 'candidate.json').read_text()), atoms)
            ['chemical_state_preserved_by_geometry'], 'Water connectivity changed')
    _, gradient, gradient_path, gradient_hash = gradient_in_frame(meta, atoms)
    require(gradient_hash == meta['result']['engrad_sha256'], 'Gradient digest mismatch')
    require(sha(geometry) == meta['result']['geometry_sha256'], 'Final geometry digest mismatch')
    paths = [ROOT / config_path, ROOT / f"reports/{cfg['job_id']}.json", geometry,
             run / 'output.log', run / 'runner_snapshot.py', run / 'job_config.json', Path(gradient_path)]
    frequencies = None
    if minimum:
        require('THE OPTIMIZATION HAS CONVERGED' in runtime, 'Optimization incomplete')
        hessian = run / (stem + '.hess')
        require(sha(hessian) == meta['result']['hessian_sha256'], 'Hessian digest mismatch')
        block = runtime.rsplit('VIBRATIONAL FREQUENCIES', 1)[-1].split('NORMAL MODES', 1)[0]
        modes = [(int(i), float(f)) for i, f in re.findall(r'^\s*(\d+):\s+([-+]?\d+\.\d+)\s+cm', block, re.M)]
        require([i for i, _ in modes] == list(range(9)) and all(abs(f) <= .001 for _, f in modes[:6]),
                'Incomplete water frequencies')
        frequencies = inspect_hessian(hessian, atoms, [f for _, f in modes[6:]])
        require(len(frequencies) == 3 and np.all(frequencies > 0), 'Nonpositive internal mode')
        paths.append(hessian)
    evidence = method_evidence(log, cfg['method'])
    smd = smd_energy_evidence(log) if cfg['solvent'] != 'vacuum' else None
    if smd is None:
        require('CPCM SOLVATION MODEL' not in runtime and 'SMD CDS (Gcds)' not in runtime,
                'Vacuum result unexpectedly contains solvent')
    return dict(job_id=cfg['job_id'], method_evidence=evidence, solvent=cfg['solvent'],
                energy_quantity=cfg['energy_quantity'], energy_Hartree=energy,
                smd_energy_evidence=smd, geometry_sha256=sha(geometry),
                max_gradient_component_Hartree_Bohr=float(np.abs(gradient).max()),
                internal_frequencies_cm1=frequencies.tolist() if frequencies is not None else None,
                elapsed_seconds=meta['elapsed_seconds'],
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                contributes_to_pilot_scores=False)


def evaluate():
    plan = json.loads((ROOT / PLAN).read_text()); verify_hashes(plan)
    configs = [json.loads((ROOT / path).read_text()) for path in plan['jobs']]
    require(len(configs) == 2 and len({cfg['job_id'] for cfg in configs}) == 2, 'Changed compatibility batch')
    missing = [cfg['job_id'] for cfg in configs if not (ROOT / f"reports/{cfg['job_id']}.json").exists()]
    require(not missing, 'Compatibility runs not completed: ' + ', '.join(missing))
    rows = [audit_job(path) for path in plan['jobs']]
    record = dict(checked_at=now(), status='automated_checks_passed_detailed_runtime_review_pending',
                  plan_sha256=sha(ROOT / PLAN), evaluator_sha256=sha(Path(__file__)), rows=rows,
                  manual_review=plan['required_post_run_review'],
                  runtime_compatibility_verified=False, full_panel_execution_enabled=False,
                  method_sensitivity_complete=False, solvent_sensitivity_complete=False,
                  scope='Two water jobs only; inspect actual ECP/D4 parameters and analytic Hessian output before release')
    target = ROOT / 'reports/completion/orca-compatibility-automated-r1.json'
    immutable(target, json.dumps(record, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(evaluate())
