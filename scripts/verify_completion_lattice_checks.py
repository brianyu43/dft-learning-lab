"""Audit completed monolayer brackets/precision from raw QE, without engine calls."""
import argparse
import json
from pathlib import Path

import numpy as np

from completion_qe_bridge import load_scf, require
from prepare_completion_alloy_probe import _raw_paths
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from verify_completion_equilibrium import verify as verify_equilibrium
from week3_materials_pipeline import centered
from week3_qe import load


SETTINGS = {'selected': (80, 6, 20), 'middle': (100, 9, 25), 'reference': (120, 12, 30)}


def equilibrium(material):
    audit = verify_equilibrium(material)
    result = load(audit['job_id'])
    sources = dict(audit['source_sha256'])
    sources[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    return result, audit, sources


def point(job, expected_atoms, cut, k, eq_audit, sources):
    """Require the reviewed settings and the exact geometry implied by equilibrium."""
    report = ROOT / f'reports/{job}.json'
    require(report.is_file(), f'Unfinished point: {job}; no completed comparison written')
    meta = json.loads(report.read_text())
    require(meta['status'] == 'passed' and meta['exit_code'] == 0 and meta['stop_reason'] is None,
            f'Point did not terminate successfully: {job}')
    cfg_path = ROOT / f'experiments/week3/materials/{job}.json'
    cfg = json.loads(cfg_path.read_text()); run = Path(meta['run_dir'])
    require(sha(cfg_path) == meta['config_sha256'] == sha(run / 'job_config.json'), 'Configuration digest differs')
    require(all(meta.get(key) == value for key, value in cfg.items()), 'Config/report values differ')
    require(meta['executable_sha256'] == eq_audit['source_sha256'].get(meta['executable'])
            == sha(ROOT / meta['executable']), 'Executable differs from the equilibrium')
    require(meta['thread_environment'] == dict(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1',
            MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1'), 'Thread environment changed')
    eq_meta = json.loads((ROOT / f"reports/{eq_audit['job_id']}.json").read_text())
    expected_assets = eq_meta['assets']
    require(meta['assets'] == expected_assets and meta['charge'] == 0, 'Pseudopotential set or charge changed')
    eq_run = Path(eq_meta['run_dir'])
    for dest in expected_assets.values():
        key = str((eq_run / dest).relative_to(ROOT))
        require(sha(run / dest) == eq_audit['source_sha256'][key], 'Pseudo bytes differ from equilibrium')
    result = load_scf(report, expected_method='C09', expected_atoms=expected_atoms)
    audit = result['audit']
    require(audit['ecutwfc_Ry'] == cut and audit['ecutrho_Ry'] == 4 * cut
            and audit['k_grid'] == [k, k, 1] and audit['conv_thr_Ry'] == 1e-10
            and audit['openmp_threads_in_log'] == 4, 'Point numerical settings differ')
    for p in _raw_paths(meta, report) + [cfg_path, run / 'job_config.json', ROOT / meta['geometry_source']]:
        sources[str(p.relative_to(ROOT))] = sha(p)
    record = dict(audit, elapsed_seconds=meta['elapsed_seconds'], run_ended_at=meta['ended_at'],
                  lattice_a_angstrom=float(result['atoms'].cell.lengths()[0]),
                  vacuum_gap_angstrom=float(result['atoms'].cell[2, 2] - np.ptp(result['atoms'].positions[:, 2])),
                  forces_eV_angstrom=result['forces'].tolist())
    return result, record


def brackets(material):
    eq, audit, sources = equilibrium(material)
    rows = []
    for direction, label in ((-1, 'minus'), (1, 'plus')):
        atoms = eq['atoms'].copy(); a = float(atoms.cell.lengths()[0])
        cell = atoms.cell.array.copy(); cell[:2] *= (a + direction * .005) / a
        atoms.set_cell(cell, scale_atoms=True)
        require(np.allclose(atoms.positions[:, 2], eq['atoms'].positions[:, 2], atol=1e-10, rtol=0),
                'Lattice probe altered internal z coordinates')
        result, record = point(f'w3-{material}-lattice-probe-{label}', atoms, 120, 12, audit, sources)
        delta = result['energy_eV'] - eq['energy_eV']
        record.update(delta_a_angstrom=direction * .005, delta_energy_eV=delta, passed=delta > 0)
        rows.append(record)
    passed = all(row['passed'] for row in rows)
    return dict(status='passed' if passed else 'needs_review', material=material,
                phase='fixed_internal_coordinates_lattice_brackets', equilibrium=audit['summary'],
                points=rows, source_sha256=sources, both_energy_differences_positive=passed,
                full_materials_baseline_complete=False,
                scope='Two local lattice displacements only; no full phonon or composition precision claim')


def precision(material):
    eq, audit, sources = equilibrium(material)
    results, records = {}, {}
    for label, (cut, k, vacuum) in SETTINGS.items():
        atoms = centered(eq['atoms'], vacuum)
        result, record = point(f'w3-{material}-precision-{label}', atoms, cut, k, audit, sources)
        require(abs(record['vacuum_gap_angstrom'] - vacuum) <= 1e-7, 'Incorrect precision vacuum gap')
        require(np.allclose(result['atoms'].positions - result['atoms'].positions[0],
                            eq['atoms'].positions - eq['atoms'].positions[0], atol=1e-7, rtol=0),
                'Precision comparison changed internal geometry')
        results[label], records[label] = result, record
    reference = results['reference']
    for label, result in results.items():
        de = abs(result['energy_eV'] - reference['energy_eV']) * 1000 / len(result['atoms'])
        df = float(np.linalg.norm(result['forces'] - reference['forces'], axis=1).max())
        records[label].update(delta_energy_meV_atom=de, max_delta_force_eV_angstrom=df,
                              passed=de <= 1 and df <= .005)
    selected = next(label for label in SETTINGS if records[label]['passed'])
    return dict(status='comparisons_verified', material=material, phase='combined_numerical_sensitivity',
                equilibrium=audit['summary'], comparisons=records,
                least_cost_passing_for_this_composition=selected,
                selection_for_both_compositions_verified=False, source_sha256=sources,
                full_materials_baseline_complete=False,
                scope='Combined cutoff/k/vacuum change relative to 120/12/30; no isolated parameter error or infinite-basis claim')


def selection():
    """Independently verify both species and the pipeline's common-setting decision."""
    proofs = {name: dict(brackets=brackets(name), precision=precision(name)) for name in ('mos2', 'mose2')}
    require(all(p['brackets']['status'] == 'passed' for p in proofs.values()), 'Lattice bracket failed')
    chosen = next(label for label in SETTINGS
                  if all(p['precision']['comparisons'][label]['passed'] for p in proofs.values()))
    cut, k, vacuum = SETTINGS[chosen]
    reference = json.loads((ROOT / 'reports/week3-materials.json').read_text())
    require(all(key in reference for key in ('lattice_brackets', 'monolayer_combined_sensitivity',
                'setting_selection', 'common_lattice_angstrom', 'strain_fraction')), 'Pipeline selection is unfinished')
    saved = reference['setting_selection']
    require((saved['label'], saved['ecutwfc_Ry'], saved['k'], saved['vacuum_angstrom'])
            == (chosen, cut, k, vacuum), 'Pipeline selected a different common setting')
    a_common = float(np.mean([p['precision']['equilibrium']['lattice_a_angstrom'] for p in proofs.values()]))
    require(abs(reference['common_lattice_angstrom'] - a_common) <= 1e-8, 'Common lattice mismatch')
    sources = {}
    for name, proof in proofs.items():
        eq = proof['precision']['equilibrium']
        require(reference['equilibria'][name] == eq, 'Pipeline equilibrium provenance differs')
        require(abs(reference['strain_fraction'][name] - (a_common / eq['lattice_a_angstrom'] - 1)) <= 1e-10,
                'Pipeline strain differs')
        saved_rows = reference['lattice_brackets'][name]
        require(len(saved_rows) == 2, 'Incomplete stored bracket')
        for saved_row, raw in zip(saved_rows, proof['brackets']['points']):
            require(saved_row['job_id'] == raw['job_id'] and saved_row['delta_a_angstrom'] == raw['delta_a_angstrom']
                    and abs(saved_row['delta_energy_eV'] - raw['delta_energy_eV']) <= 1e-8,
                    'Stored bracket differs from raw output')
        for label, raw in proof['precision']['comparisons'].items():
            saved_row = reference['monolayer_combined_sensitivity'][name][label]
            require(saved_row['job_id'] == raw['job_id'] and saved_row['passed'] == raw['passed']
                    and abs(saved_row['delta_energy_meV_atom'] - raw['delta_energy_meV_atom']) <= 1e-8
                    and abs(saved_row['max_delta_force_eV_angstrom'] - raw['max_delta_force_eV_angstrom']) <= 1e-10,
                    'Stored precision comparison differs from raw output')
        for value in proof.values():
            for path, digest in value['source_sha256'].items():
                require(path not in sources or sources[path] == digest, 'Source changed during audit')
                sources[path] = digest
    # This report continues to evolve. Preserve the inspected sections, not a
    # digest that would incorrectly invalidate the audit when later stages finish.
    snapshot_keys = ('equilibria', 'lattice_brackets', 'monolayer_combined_sensitivity',
                     'setting_selection', 'common_lattice_angstrom', 'strain_fraction')
    return dict(status='monolayer_lattice_and_precision_verified', common_setting_label=chosen,
                common_lattice_angstrom=a_common, proofs=proofs, source_sha256=sources,
                pipeline_selection_snapshot={key: reference[key] for key in snapshot_keys},
                full_materials_baseline_complete=False,
                scope='Both monolayers only; common-cell relaxation and heterobilayer precision remain separate gates')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('material', choices=('mos2', 'mose2', 'both'))
    parser.add_argument('phase', choices=('brackets', 'precision', 'selection'))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.phase == 'selection':
        require(args.material == 'both', 'Common selection requires both compositions')
        result = selection()
    else:
        require(args.material != 'both', 'Individual comparison requires one composition')
        result = {'brackets': brackets, 'precision': precision}[args.phase](args.material)
    if args.write:
        path = ROOT / f'reports/completion/lattice-check-audits-r1/{args.material}-{args.phase}.json'
        immutable(path, json.dumps(result, indent=2) + '\n'); print(path)
    else:
        print(json.dumps(result, indent=2))
    if result['status'] == 'needs_review':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
