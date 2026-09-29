"""Prepare/evaluate bounded QE OpenMP comparisons without launching an engine."""
import argparse
import json
from pathlib import Path
import statistics

import numpy as np

from completion_qe_bridge import load_scf, require
from prepare_week3_drug import immutable
from research_run import ROOT, sha

PARENT = 'reports/w2-bilayer-selected-r0-d6p3-pair.json'
SOURCE_CONFIG = 'experiments/week2/w2-bilayer-selected-r0-d6p3-pair.json'
PLAN = 'experiments/completion/qe-thread-benchmark-r1.json'
OUTPUT = 'reports/completion/qe-thread-benchmark-r1.json'
ORDER = [(1, 4), (1, 8), (1, 12), (2, 12), (2, 8), (2, 4)]
POLICY = dict(max_energy_difference_eV=1e-5, max_force_vector_difference_eV_angstrom=1e-5,
              minimum_speedup_in_both_repeats=1.1)


def verify_hashes(hashes):
    for p, digest in hashes.items():
        require((ROOT / p).is_file() and sha(ROOT / p) == digest, f'Changed dependency: {p}')


def prepare():
    target = ROOT / PLAN
    if target.exists():
        saved = json.loads(target.read_text())
        verify_hashes(saved['input_hashes'])
        require(saved['policy'] == POLICY and [(r['repeat'], r['threads']) for r in saved['mapping']] == ORDER,
                'Saved benchmark policy changed')
        return target
    original = load_scf(ROOT / PARENT)
    meta = json.loads((ROOT / PARENT).read_text())
    cfg = json.loads((ROOT / SOURCE_CONFIG).read_text())
    require(meta['config_sha256'] == sha(ROOT / SOURCE_CONFIG) and meta['exit_code'] == 0
            and not meta.get('stop_reason'), 'Original run/config is incomplete')
    require(sha(ROOT / cfg['executable']) == meta['executable_sha256'], 'QE executable changed')
    for source, digest in meta['source_sha256'].items():
        require(sha(ROOT / source) == digest, f'Original input changed: {source}')
    run = Path(meta['run_dir'])
    sources = [PARENT, SOURCE_CONFIG, cfg['executable'], cfg['input'], *cfg['assets'],
               'scripts/completion_qe_thread_benchmark.py', 'scripts/completion_qe_bridge.py',
               'scripts/research_run.py', 'scripts/completion_batch.py', 'scripts/prepare_week3_drug.py',
               'reports/week2-summary.json', str((run / 'output.log').relative_to(ROOT))]
    sources += [str(p.relative_to(ROOT)) for p in (run / 'scratch').glob('*.save/data-file-schema.xml')]
    jobs, mapping = [], []
    for repeat, threads in ORDER:
        job_id = f'w3-qe-omp-r1-n{threads}-r{repeat}'
        inp = f'inputs/qe/completion/{job_id}.in'
        immutable(ROOT / inp, (ROOT / cfg['input']).read_text())
        new = {**cfg, 'protocol_id': 'dft-pilot-completion-qe-threads-r1', 'week': 3,
               'job_id': job_id, 'task_id': 'QE-OpenMP-performance-validation', 'input': inp,
               'omp_threads': threads, 'timeout_seconds': 900, 'rss_limit_gib': 40,
               'parent_report': PARENT, 'parent_report_sha256': sha(ROOT / PARENT),
               'benchmark_repeat': repeat, 'purpose': 'Thread-count comparison at identical archived six-atom coordinates',
               'accuracy_scope': 'Parallel numerical equivalence and timing at one MoS2 bilayer; no new scientific sample'}
        config = f'experiments/completion/qe-threads/{job_id}.json'
        immutable(ROOT / config, json.dumps(new, indent=2) + '\n')
        require(sha(ROOT / inp) == meta['source_sha256'][cfg['input']], 'Benchmark input bytes changed')
        jobs.append(config); sources.extend([config, inp])
        mapping.append(dict(repeat=repeat, threads=threads, config=config, job_id=job_id))
    plan = dict(kind='qe_openmp_thread_benchmark_r1', jobs=jobs, mapping=mapping, policy=POLICY,
        parent_report=PARENT, input_hashes={p: sha(ROOT / p) for p in sorted(set(sources))},
        scientific_atom_count=len(original['atoms']), maximum_total_seconds=5400,
        automatic_execution=False, automatic_retries=0, numerical_precision_verified=False,
        execution_precondition='After the current materials controller exits, with the single-engine lock and remaining pilot budget',
        input_prefix_policy='Byte-identical QE input, including prefix; scratch is isolated inside each unique run directory',
        selection_policy='Only consider 8 or 12 threads if both repeat speedups over corresponding 4-thread runs are at least 1.1; choose fastest eligible median',
        scope='Six fresh SCFs, not six independent structures; two timing repeats do not establish universal scaling')
    verify_hashes(plan['input_hashes'])
    immutable(target, json.dumps(plan, indent=2) + '\n')
    return target


def evaluate():
    require((ROOT / PLAN).is_file(), 'Benchmark plan is missing')
    plan = json.loads((ROOT / PLAN).read_text()); verify_hashes(plan['input_hashes'])
    require(plan['policy'] == POLICY and [(r['repeat'], r['threads']) for r in plan['mapping']] == ORDER,
            'Benchmark policy or order changed')
    require(plan['jobs'] == [r['config'] for r in plan['mapping']], 'Job coverage changed')
    for row in plan['mapping']:
        require((ROOT / f"reports/{row['job_id']}.json").exists(), 'Benchmark incomplete; no result report written')
    original = load_scf(ROOT / PARENT)
    original_meta = json.loads((ROOT / PARENT).read_text())
    rows, evidence = [], dict(plan['input_hashes'])
    evidence[PLAN] = sha(ROOT / PLAN)
    for row in plan['mapping']:
        report = ROOT / f"reports/{row['job_id']}.json"
        m = json.loads(report.read_text()); cfg = json.loads((ROOT / row['config']).read_text())
        run = Path(m['run_dir'])
        require(m['exit_code'] == 0 and not m.get('stop_reason')
                and m['config_sha256'] == sha(ROOT / row['config']) == sha(run / 'job_config.json')
                and all(m.get(k) == v for k, v in cfg.items()), 'Benchmark run/config is incomplete or changed')
        require(m['omp_threads'] == row['threads'] and m['benchmark_repeat'] == row['repeat'], 'Thread/repeat mapping changed')
        require(m['executable_sha256'] == original_meta['executable_sha256'], 'Benchmark executable changed')
        require(m['source_sha256'][cfg['input']] == original_meta['source_sha256'][original_meta['input']], 'SCF input changed')
        require(all(m['source_sha256'][p] == original_meta['source_sha256'][p] for p in cfg['assets']), 'Pseudopotential changed')
        require(m['thread_environment'] == dict(OMP_NUM_THREADS=str(row['threads']), OPENBLAS_NUM_THREADS='1',
                MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1'), 'Nested thread settings changed')
        value = load_scf(report, expected_atoms=original['atoms'])
        de = abs(value['energy_eV'] - original['energy_eV'])
        df = float(np.linalg.norm(value['forces'] - original['forces'], axis=1).max())
        seconds = float(m['elapsed_seconds'])
        require(np.isfinite(seconds) and seconds > 0, 'Invalid measured duration')
        passed = de <= POLICY['max_energy_difference_eV'] and df <= POLICY['max_force_vector_difference_eV_angstrom']
        rows.append(dict(**row, elapsed_seconds=seconds, delta_energy_eV=de,
            max_delta_force_vector_eV_angstrom=df, numerical_equivalence_passed=passed, raw_audit=value['audit']))
        files = [report, run / 'job_config.json', run / 'runner_snapshot.py', run / 'output.log',
                 run / Path(cfg['input']).name, *(run / d for d in cfg['assets'].values()),
                 *list((run / 'scratch').glob('*.save/data-file-schema.xml'))]
        evidence.update({str(p.relative_to(ROOT)): sha(p) for p in files})
    times = {(r['repeat'], r['threads']): r['elapsed_seconds'] for r in rows}
    timing = {str(n): dict(median_seconds=statistics.median(times[r, n] for r in (1, 2)),
        speedup_vs_4_by_repeat=[times[r, 4] / times[r, n] for r in (1, 2)]) for n in (4, 8, 12)}
    eligible = [4] + [n for n in (8, 12) if min(timing[str(n)]['speedup_vs_4_by_repeat']) >= POLICY['minimum_speedup_in_both_repeats']]
    passed = all(r['numerical_equivalence_passed'] for r in rows)
    result = dict(status='parallel_equivalence_and_timing_reviewed' if passed else 'stopped_requires_review',
        rows=rows, timing=timing, policy=POLICY, all_numerical_checks_passed=passed,
        recommended_threads_for_this_case=min(eligible, key=lambda n: timing[str(n)]['median_seconds']) if passed else None,
        production_settings_changed=False, broader_scaling_verified=False, scientific_project_complete=False,
        source_sha256=evidence)
    verify_hashes(evidence)
    immutable(ROOT / OUTPUT, json.dumps(result, indent=2) + '\n')
    return ROOT / OUTPUT


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('prepare', 'evaluate'))
    args = parser.parse_args()
    print(prepare() if args.stage == 'prepare' else evaluate())
