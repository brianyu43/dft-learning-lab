"""Evaluate the prepared finite-displacement and real D3(BJ) compatibility jobs."""
import json
from pathlib import Path
import sys

import numpy as np

from completion_qe_bridge import load_scf, require
from research_run import ROOT, now, save, sha


def main():
    path = ROOT / 'experiments/completion/qe-force-and-method-checks-r1.json'
    plan = json.loads(path.read_text())
    for source, digest in plan['input_hashes'].items():
        require(sha(ROOT / source) == digest, f'Frozen source changed: {source}')
    settings = plan['force_check']
    center = load_scf(ROOT / settings['parent_report'])
    runs = []
    for config_path in plan['jobs']:
        config = json.loads((ROOT / config_path).read_text())
        report_path = ROOT / 'reports' / (config['job_id'] + '.json')
        require(report_path.exists(), f'Calculation has not finished: {config["job_id"]}')
        meta = json.loads(report_path.read_text())
        require(meta['config_sha256'] == sha(ROOT / config_path), 'Executed configuration differs from plan')
        expected = load_scf(ROOT / config['parent_report'])['atoms'].copy()
        if 'displacement_angstrom' in config:
            expected.positions[config['displaced_atom_index'], config['displaced_cartesian_axis']] += config['displacement_angstrom']
        result = load_scf(report_path, expected_method=config['expected_bridge_method'], expected_atoms=expected)
        runs.append(dict(config=config, result=result, report=report_path))
    signs = {row['config'].get('displacement_angstrom'): row for row in runs
             if 'displacement_angstrom' in row['config']}
    h = settings['step_angstrom']
    require(set(signs) == {-h, h}, 'Incomplete or mismatched displacement pair')
    for row in signs.values():
        audit = row['result']['audit']
        for key in ('method', 'ecutwfc_Ry', 'ecutrho_Ry', 'k_grid', 'k_offset', 'valence_electrons'):
            require(audit[key] == center['audit'][key], f'Derivative settings differ: {key}')
    plus = signs[h]['result']['energy_eV']
    minus = signs[-h]['result']['energy_eV']
    derivative_force = -(plus - minus) / (2 * h)
    reference_force = float(center['forces'][settings['atom_index'], settings['cartesian_axis']])
    difference = abs(derivative_force - reference_force)
    limit = settings['absolute_force_difference_limit_eV_angstrom']
    compatible = [row for row in runs if row['config']['expected_bridge_method'] == 'PBE-D3BJ']
    require(len(compatible) == 1, 'Expected one D3(BJ) compatibility result')
    passed = np.isfinite(difference) and difference <= limit
    output = dict(status='passed' if passed else 'needs_review', checked_at=now(),
                  plan_sha256=sha(path), evaluator_sha256=sha(Path(__file__)),
                  results=[row['result']['audit'] for row in runs], center=center['audit'],
                  finite_difference=dict(step_angstrom=h, atom_index=settings['atom_index'],
                                         cartesian_axis=settings['cartesian_axis'],
                                         plus_energy_eV=plus, minus_energy_eV=minus,
                                         force_from_energy_eV_angstrom=derivative_force,
                                         force_from_XML_eV_angstrom=reference_force,
                                         absolute_difference_eV_angstrom=difference,
                                         threshold_eV_angstrom=limit, passed=bool(passed)),
                  pbe_d3bj_options_verified_in_actual_output=True,
                  production_alloy_or_constrained_relaxation_completed=False,
                  limits=['One MoS2 chalcogen direction; finite step and SCF error remain in the derivative',
                          'PBE-D3BJ compatibility is not validation of CI/CO energy rankings',
                          'A complete constrained QE optimization still requires its own trajectory audit'])
    save(ROOT / 'reports/completion/qe-force-and-method-checks-r1.json', output)
    print(json.dumps({k: output[k] for k in ['status', 'finite_difference', 'pbe_d3bj_options_verified_in_actual_output']}, indent=2))
    if not passed:
        sys.exit(1)


if __name__ == '__main__':
    main()
