"""Read-only terminal audit of either monolayer equilibrium; no engine calls."""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np

from prepare_completion_alloy_probe import verify_optimized, _raw_paths
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from week3_materials_pipeline import summarize
from week3_qe import load


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(material):
    require(material in ('mos2', 'mose2'), 'Only the two frozen equilibrium jobs are supported')
    job = f'w3-{material}-lattice-vc'
    report = ROOT / f'reports/{job}.json'
    require(report.is_file(), 'Terminal report is not available; no equilibrium completion claim')
    meta = json.loads(report.read_text())
    require(meta['status'] == 'passed' and meta['exit_code'] == 0 and meta['stop_reason'] is None,
            'Job did not terminate successfully')
    config_path = ROOT / f'experiments/week3/materials/{job}.json'
    config = json.loads(config_path.read_text()); run = Path(meta['run_dir'])
    require(sha(config_path) == meta['config_sha256'] == sha(run / 'job_config.json'), 'Configuration digest differs')
    require(all(meta.get(key) == value for key, value in config.items()), 'Config/report values differ')
    require(meta['calculation'] == 'vc-relax' and meta['ecutwfc_Ry'] == 120 and meta['ecutrho_Ry'] == 480
            and meta['k_grid'] == [12, 12, 1] and meta['omp_threads'] == 4 and meta['conv_thr_Ry'] == 1e-10,
            'Equilibrium protocol changed')
    require(sha(ROOT / meta['executable']) == meta['executable_sha256'], 'QE executable differs')
    result = load(job); summary = summarize(result)
    result = verify_optimized(job, summary)  # Closed XML, XC, final cutoffs, k-grid and SCF/Opt convergence.
    expected_symbols = ['Mo', 'S', 'S'] if material == 'mos2' else ['Mo', 'Se', 'Se']
    require(result['atoms'].get_chemical_symbols() == expected_symbols, 'Wrong equilibrium composition')
    require(result['nelec'] == (26 if material == 'mos2' else 46), 'Unexpected valence electron count')
    log = (run / 'output.log').read_text()
    require('JOB DONE.' in log and re.search(r'Program PWSCF v\.7\.6\b', log), 'Missing normal QE 7.6 runtime completion')
    threads = re.findall(r'Serial multi-threaded version, running on\s+(\d+) processor cores', log)
    require(threads == ['4'], 'Actual OpenMP thread count is not four')
    require(meta['thread_environment'] == dict(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1',
            MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1'), 'Thread environment differs')
    xml = next((run / 'scratch').glob('*.save/data-file-schema.xml'))
    root = ET.parse(xml).getroot()
    paths = _raw_paths(meta, report) + [config_path, run / 'job_config.json', ROOT / meta['geometry_source'],
        ROOT / meta['executable'], Path(__file__), ROOT / 'scripts/prepare_completion_alloy_probe.py',
        ROOT / 'scripts/completion_qe_bridge.py', ROOT / 'scripts/week3_qe.py',
        ROOT / 'scripts/week3_materials_pipeline.py']
    return dict(status='equilibrium_verified_brackets_and_precision_pending', material=material, job_id=job,
        run_ended_at=meta['ended_at'], summary=summary,
        cell_angstrom=result['atoms'].cell.array.tolist(), cell_angles_degrees=result['atoms'].cell.angles().tolist(),
        atomic_positions_angstrom=result['atoms'].positions.tolist(),
        max_in_plane_stress_kbar=float(np.abs(result['stress_eV_angstrom3'][:2, :2]).max()) * 1602.176634,
        final_scf_error_Ry=2 * float(root.findtext('output/convergence_info/scf_conv/scf_error')),
        ionic_converged=root.findtext('output/convergence_info/opt_conv/convergence_achieved'),
        completed_scf_cycles=log.count('convergence has been achieved'),
        runtime_omp_threads=4, elapsed_seconds=meta['elapsed_seconds'],
        memory={key: meta[key] for key in ('rss_limit_gib', 'peak_sampled_process_tree_rss_bytes',
            'memory_samples', 'memory_access_denied', 'system_swap_increase_peak_bytes', 'memory_measurement')},
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        lattice_energy_brackets_verified=False, composition_numerical_precision_verified=False,
        phonon_stability_verified=False, full_materials_baseline_complete=False,
        scope='Converged equilibrium at the frozen C09 settings; local lattice brackets and higher-level validations remain')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('material', choices=('mos2', 'mose2'))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args(); result = verify(args.material)
    if args.write:
        path = ROOT / f'reports/completion/equilibrium-audits-r1/{args.material}.json'
        immutable(path, json.dumps(result, indent=2) + '\n'); print(path)
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
