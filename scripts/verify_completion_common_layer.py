"""Audit a completed strained monolayer without starting or modifying QE."""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from ase.io import read
from ase.io.espresso import read_fortran_namelist

from completion_qe_bridge import require, NUMBER, HARTREE_EV, BOHR_ANGSTROM
from prepare_completion_alloy_probe import verify_optimized, _raw_paths
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from week3_materials_pipeline import summarize
from week3_qe import load


def verify(material):
    require(material in ('mos2', 'mose2'), 'Unreviewed material')
    job = f'w3-{material}-common-lattice-relax'
    report = ROOT / f'reports/{job}.json'
    require(report.is_file(), 'Terminal report is not available; common-layer completion is unproven')
    selection_path = ROOT / 'reports/completion/lattice-check-audits-r1/both-selection.json'
    selection = json.loads(selection_path.read_text())
    require(selection['status'] == 'monolayer_lattice_and_precision_verified', 'Joint selection did not pass')
    sources = dict(selection['source_sha256'])
    for p, digest in sources.items():
        require(sha(ROOT / p) == digest, f'Changed selection evidence: {p}')
    snapshot = selection['pipeline_selection_snapshot']
    settings = snapshot['setting_selection']
    a_common = selection['common_lattice_angstrom']
    eq = verify_optimized(f'w3-{material}-lattice-vc', snapshot['equilibria'][material])
    expected = eq['atoms'].copy()
    cell = expected.cell.array.copy()
    cell[:2] *= a_common / expected.cell.lengths()[0]
    expected.set_cell(cell, scale_atoms=True)
    expected.center(vacuum=settings['vacuum_angstrom'] / 2, axis=2)

    meta = json.loads(report.read_text())
    run = Path(meta['run_dir'])
    config_path = ROOT / f'experiments/week3/materials/{job}.json'
    config = json.loads(config_path.read_text())
    require(meta['status'] == 'passed' and meta['exit_code'] == 0 and meta['stop_reason'] is None,
            'Common-layer run did not complete normally')
    require(sha(config_path) == meta['config_sha256'] == sha(run / 'job_config.json')
            and all(meta.get(k) == v for k, v in config.items()), 'Config or report changed')
    require(meta['calculation'] == 'relax' and meta['role'] == 'strained_monolayer_internal_relaxation'
            and meta['common_lattice_angstrom'] == a_common and meta['omp_threads'] == 4
            and meta['ecutwfc_Ry'] == settings['ecutwfc_Ry']
            and meta['ecutrho_Ry'] == 4 * settings['ecutwfc_Ry']
            and meta['k_grid'] == [settings['k'], settings['k'], 1]
            and meta['initial_vacuum_gap_angstrom'] == settings['vacuum_angstrom'], 'Strained-layer protocol changed')
    require(sha(ROOT / meta['executable']) == meta['executable_sha256'], 'QE executable changed')
    initial = read(run / Path(meta['input']).name, format='espresso-in')
    require(initial.get_chemical_symbols() == expected.get_chemical_symbols()
            and np.allclose(initial.cell, expected.cell, rtol=0, atol=1e-7)
            and np.allclose(initial.positions, expected.positions, rtol=0, atol=1e-7),
            'Input does not follow audited equilibrium, mean lattice and selected vacuum')
    with (run / Path(meta['input']).name).open() as handle:
        inp, _ = read_fortran_namelist(handle)
    require(inp['control']['calculation'] == 'relax' and inp['ions']['ion_dynamics'] == 'bfgs'
            and inp['control']['etot_conv_thr'] == 1e-6 and inp['control']['forc_conv_thr'] == 1e-4
            and inp['electrons']['conv_thr'] == 1e-10, 'Relaxation convergence settings changed')
    value = load(job)
    value = verify_optimized(job, summarize(value))
    atoms = value['atoms']
    require(np.allclose(atoms.cell, initial.cell, rtol=0, atol=1e-7)
            and value['nelec'] == eq['nelec'] and value['max_force_eV_angstrom'] <= .005,
            'Final fixed cell, electron count or force failed')
    log = (run / 'output.log').read_text()
    require('JOB DONE.' in log and re.search(r'Program PWSCF v\.7\.6\b', log), 'Missing QE normal termination')
    require(re.findall(r'Serial multi-threaded version, running on\s+(\d+) processor cores', log) == ['4']
            and meta['thread_environment'] == dict(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1',
                MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1'), 'Actual thread configuration differs')
    force_rows = re.findall(r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*(' + NUMBER + r')\s+('
                           + NUMBER + r')\s+(' + NUMBER + ')', log)
    require(len(force_rows) >= len(atoms) and 'Ry/au' in log, 'Missing force output')
    forces = np.array([[float(x.replace('D', 'E').replace('d', 'e')) for x in row]
                       for row in force_rows[-len(atoms):]]) * HARTREE_EV / (2 * BOHR_ANGSTROM)
    force_error = float(np.abs(forces - value['forces']).max())
    require(force_error <= .51e-8 * HARTREE_EV / (2 * BOHR_ANGSTROM), 'XML/text force conversion differs')
    xml = next((run / 'scratch').glob('*.save/data-file-schema.xml'))
    tree = ET.parse(xml).getroot()
    paths = _raw_paths(meta, report) + [selection_path, config_path, run / 'job_config.json',
        ROOT / meta['geometry_source'], ROOT / meta['executable'], Path(__file__),
        ROOT / 'scripts/prepare_completion_alloy_probe.py', ROOT / 'scripts/completion_qe_bridge.py',
        ROOT / 'scripts/week3_qe.py', ROOT / 'scripts/week3_materials_pipeline.py']
    sources.update({str(p.relative_to(ROOT)): sha(p) for p in paths})
    return dict(status='common_lattice_internal_relaxation_verified', material=material, job_id=job,
        summary=summarize(value), run_ended_at=meta['ended_at'], elapsed_seconds=meta['elapsed_seconds'],
        selected_settings=settings, common_lattice_angstrom=a_common,
        imposed_strain_fraction=a_common / eq['atoms'].cell.lengths()[0] - 1,
        cell_angstrom=atoms.cell.array.tolist(), atomic_positions_angstrom=atoms.positions.tolist(),
        forces_eV_angstrom=value['forces'].tolist(),
        stress_kbar=(value['stress_eV_angstrom3'] * 1602.176634).tolist(),
        initial_layer_thickness_angstrom=float(np.ptp(initial.positions[:, 2])),
        relaxed_layer_thickness_angstrom=float(np.ptp(atoms.positions[:, 2])),
        maximum_atom_displacement_angstrom=float(np.linalg.norm(atoms.positions - initial.positions, axis=1).max()),
        final_scf_error_Ry=2 * float(tree.findtext('output/convergence_info/scf_conv/scf_error')),
        completed_scf_cycles=log.count('convergence has been achieved'),
        xml_text_force_max_component_difference_eV_angstrom=force_error,
        peak_sampled_process_tree_rss_bytes=meta['peak_sampled_process_tree_rss_bytes'],
        runtime_omp_threads=4, source_sha256=sources,
        heterobilayer_baseline_complete=False, phonon_stability_verified=False,
        scope='Internal relaxation of one strained isolated layer at the imposed common cell. '
              'Nonzero in-plane stress is expected; this is not a zero-stress bilayer equilibrium or phonon test.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('material', choices=('mos2', 'mose2'))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = verify(args.material)
    if args.write:
        path = ROOT / f'reports/completion/common-layer-audits-r1/{args.material}.json'
        immutable(path, json.dumps(result, indent=2) + '\n')
        print(path)
    else:
        print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__': main()
