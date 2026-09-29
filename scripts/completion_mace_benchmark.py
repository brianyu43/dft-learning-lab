"""Bounded CPU/MPS comparison of one frozen MACE model; no DFT replacement.

Preparation uses an audited archived geometry. Execution takes both existing
engine locks and is explicitly separate from the active QE/ORCA controller.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import statistics
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / 'experiments/completion/mace-device-benchmark-r1.json'
CHECKPOINT = 'data/completion/materials-models/checkpoints/mace-mp-0b3-medium.model'
CHECKPOINT_SHA = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
BASELINE = 'reports/w2-bilayer-selected-r0-d6p3-pair.json'
VARIANTS = {'cpu64_4': ('cpu', 'float64', 4), 'cpu32_4': ('cpu', 'float32', 4),
            'cpu64_8': ('cpu', 'float64', 8), 'cpu32_8': ('cpu', 'float32', 8),
            'mps32_4': ('mps', 'float32', 4)}
POLICY = dict(energy_meV_atom=1., force_eV_angstrom=.005, gamma_meV_angstrom2=.1,
              minimum_speedup=1.1, warmups=2, repeats=5, worker_timeout_seconds=180,
              maximum_workers=10, rss_limit_gib=40)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, allow_nan=False) + '\n'
    with path.open('x') as stream:
        stream.write(text)


def verify_plan(plan):
    require(plan['policy'] == POLICY and plan['variants'] == {k: list(v) for k, v in VARIANTS.items()},
            'Benchmark policy changed')
    for path, digest in plan['source_sha256'].items():
        require(sha(ROOT / path) == digest, f'Changed benchmark dependency: {path}')
    require(plan['source_sha256'][CHECKPOINT] == CHECKPOINT_SHA, 'Unreviewed checkpoint')
    for case in plan['cases'].values():
        pair = case['parts']['pair']
        indices = case['layer_indices']
        require(sorted(indices[0] + indices[1]) == list(range(len(pair['numbers']))), 'Invalid partition')
        for label, ids in zip(('bottom', 'top'), indices):
            part = case['parts'][label]
            require(part['cell_angstrom'] == pair['cell_angstrom'] and part['pbc'] == pair['pbc'],
                    'Fragment cell or periodicity differs')
            require(part['numbers'] == [pair['numbers'][i] for i in ids]
                    and part['positions_angstrom'] == [pair['positions_angstrom'][i] for i in ids],
                    'Fragment is not an exact subset')
        area = float(np.linalg.norm(np.cross(*np.array(pair['cell_angstrom'])[:2])))
        require(area > 0 and abs(area - case['area_angstrom2']) < 1e-10, 'Invalid area')


def prepare():
    if PLAN.exists():
        verify_plan(json.loads(PLAN.read_text()))
        return PLAN
    from completion_qe_bridge import load_scf
    source = load_scf(ROOT / BASELINE)
    pair = source['atoms'].copy()
    require(pair.numbers.tolist() == [42, 16, 16, 42, 16, 16], 'Unexpected archived atom order')
    # Explicit performance fixture, not an optimized MoSe2 layer or DFT result.
    pair.numbers[4:6] = 34
    pair.set_array('layer_id', np.array([0, 0, 0, 1, 1, 1]))
    cases = {}
    for name, atoms in [('n6', pair), ('n24', pair.repeat((2, 2, 1)))]:
        ids = [np.flatnonzero(atoms.arrays['layer_id'] == i).tolist() for i in (0, 1)]
        if name == 'n24':
            # One inner-face S/Se exchange exercises the same total composition.
            a = max((i for i in ids[0] if atoms.numbers[i] == 16), key=lambda i: atoms.positions[i, 2])
            b = min((i for i in ids[1] if atoms.numbers[i] == 34), key=lambda i: atoms.positions[i, 2])
            atoms.numbers[a], atoms.numbers[b] = 34, 16
        parts = {}
        for part, item in [('pair', atoms), ('bottom', atoms[ids[0]]), ('top', atoms[ids[1]])]:
            parts[part] = dict(numbers=item.numbers.tolist(), positions_angstrom=item.positions.tolist(),
                               cell_angstrom=item.cell.array.tolist(), pbc=item.pbc.tolist())
        cases[name] = dict(parts=parts, layer_indices=ids,
                          area_angstrom2=float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))))
    meta = json.loads((ROOT / BASELINE).read_text())
    run = Path(meta['run_dir'])
    paths = [ROOT / BASELINE, run / 'output.log', run / Path(meta['input']).name,
             run / 'runner_snapshot.py', *list((run / 'scratch').glob('*.save/data-file-schema.xml')),
             *[run / p for p in meta['assets'].values()], Path(__file__), ROOT / CHECKPOINT,
             ROOT / 'scripts/completion_qe_bridge.py', ROOT / 'scripts/research_run.py',
             ROOT / 'environment/mps-requirements.lock', ROOT / 'reports/completion/mace-mp0b3-checkpoint-metadata-r1.json']
    plan = dict(status='prepared_not_executed', policy=POLICY, variants={k: list(v) for k, v in VARIANTS.items()},
        model='MACE-MP-0b3 medium', dispersion_correction='none', cases=cases,
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        source_geometry_audit=source['audit'],
        scope='Backend and dtype benchmark on explicitly unrelaxed performance fixtures. '
              'No C09 accuracy, D3 calibration, minimum, or scientific model suitability is established.')
    verify_plan(plan)
    dump_new(PLAN, plan)
    return PLAN


def deltas(plan, reference, target):
    """Independent observable checks, including cancellation in pair-minus-layers."""
    require(set(reference) == set(plan['cases']) == set(target), 'Missing or extra benchmark case')
    rows = {}
    for name, case in plan['cases'].items():
        metrics = {}
        energies = []
        for data in (reference[name], target[name]):
            require(set(data) == set(case['parts']), 'Incomplete pair or layer results')
        for part, geom in case['parts'].items():
            a, b = reference[name][part], target[name][part]
            ea, eb = float(a['energy_eV']), float(b['energy_eV'])
            fa, fb = np.array(a['forces_eV_angstrom']), np.array(b['forces_eV_angstrom'])
            n = len(geom['numbers'])
            require(fa.shape == fb.shape == (n, 3) and np.isfinite([ea, eb]).all()
                    and np.isfinite(fa).all() and np.isfinite(fb).all(), 'Invalid energy/force data')
            de, df = abs(eb - ea) * 1000 / n, float(np.linalg.norm(fb - fa, axis=1).max())
            metrics[part] = dict(delta_energy_meV_atom=de, max_delta_force_eV_angstrom=df,
                                passed=de <= POLICY['energy_meV_atom'] and df <= POLICY['force_eV_angstrom'])
        for result in (reference[name], target[name]):
            energies.append((result['pair']['energy_eV'] - result['bottom']['energy_eV']
                             - result['top']['energy_eV']) * 1000 / case['area_angstrom2'])
        dg = abs(energies[1] - energies[0])
        rows[name] = dict(parts=metrics, reference_gamma_meV_angstrom2=energies[0],
            target_gamma_meV_angstrom2=energies[1], delta_gamma_meV_angstrom2=dg,
            passed=all(r['passed'] for r in metrics.values()) and dg <= POLICY['gamma_meV_angstrom2'])
    return rows


def worker(plan, variant, output):
    started = time.monotonic()
    # Only the supervisor passes these already-held open file descriptions.
    # A standalone worker command must never bypass the one-engine rule.
    inherited = os.environ.get('MACE_BENCHMARK_LOCK_FDS', '').split(',')
    require(len(inherited) == 2 and all(x.isdigit() for x in inherited), 'Worker requires supervisor locks')
    for fd, name in zip(map(int, inherited), ('.completion-controller.lock', '.research-run.lock')):
        actual, expected = os.fstat(fd), (ROOT / 'experiments' / name).stat()
        require((actual.st_dev, actual.st_ino) == (expected.st_dev, expected.st_ino), 'Incorrect inherited lock')
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    import importlib.metadata
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    device, dtype, threads = VARIANTS[variant]
    require(os.environ.get('PYTORCH_ENABLE_MPS_FALLBACK') == '0', 'MPS fallback must be disabled')
    require(os.environ.get('PYTORCH_MPS_FAST_MATH', '0') == '0', 'Unreviewed MPS fast math')
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    if device == 'mps':
        require(torch.backends.mps.is_available(), 'MPS unavailable; do not substitute CPU')
        torch.mps.set_per_process_memory_fraction(.5)
    # Cast on CPU: loading stored float64 directly onto MPS is unsupported.
    model = torch.load(ROOT / CHECKPOINT, map_location='cpu', weights_only=False)
    model = model.double() if dtype == 'float64' else model.float()
    calc = MACECalculator(models=model, device=device, default_dtype=dtype,
                          energy_units_to_eV=1., length_units_to_A=1.)
    observed = set()

    def capture(_module, _inputs, out):
        for key in ('energy', 'forces'):
            value = out[key]
            require(value.device.type == device and str(value.dtype) == 'torch.' + dtype,
                    f'Unexpected forward output placement: {key}')
            observed.add((key, str(value.device), str(value.dtype)))

    hook = calc.models[0].register_forward_hook(capture)
    result = dict(variant=variant, device=device, dtype=dtype, threads=torch.get_num_threads(),
        actual_forward_output_placement=[], fallback_enabled=False, results={},
        packages={p: importlib.metadata.version(p) for p in ('torch', 'mace-torch', 'numpy', 'ase')},
        load_and_import_seconds=time.monotonic() - started, plan_sha256=sha(PLAN),
        model_sha256=CHECKPOINT_SHA, phase='single_structure_ASE_calculator_including_neighbor_graph_and_stress')
    for name, case in plan['cases'].items():
        result['results'][name] = {}
        for part, geom in case['parts'].items():
            atoms = Atoms(numbers=geom['numbers'], positions=geom['positions_angstrom'],
                          cell=geom['cell_angstrom'], pbc=geom['pbc'])
            samples = []
            warmup_seconds = []
            for i in range(POLICY['warmups'] + POLICY['repeats']):
                if device == 'mps': torch.mps.synchronize()
                tick = time.monotonic()
                # Direct calculate forces fresh inference; ASE's property cache is not timed.
                calc.calculate(atoms, properties=['energy', 'forces'])
                if device == 'mps': torch.mps.synchronize()
                elapsed = time.monotonic() - tick
                energy = float(calc.results['energy'])
                forces = np.array(calc.results['forces'], dtype=float)
                require(np.isfinite(energy) and forces.shape == (len(atoms), 3)
                        and np.isfinite(forces).all(), 'Nonfinite model result')
                row = dict(energy_eV=energy, forces_eV_angstrom=forces.tolist(), elapsed_seconds=elapsed)
                if i < POLICY['warmups']: warmup_seconds.append(elapsed)
                else: samples.append(row)
            result['results'][name][part] = dict(**{k: samples[-1][k] for k in ('energy_eV', 'forces_eV_angstrom')},
                median_seconds=statistics.median(s['elapsed_seconds'] for s in samples),
                samples=samples, warmup_seconds=warmup_seconds)
    hook.remove()
    require(len(observed) == 2, 'Missing device proof for energy and forces')
    result.update(status='completed', actual_forward_output_placement=sorted(observed),
                  elapsed_seconds_including_import=time.monotonic() - started,
                  allocated_gpu_bytes=torch.mps.current_allocated_memory() if device == 'mps' else None)
    dump_new(output, result)


def run(plan):
    """Locks precede directory creation/import; no overlap and no automatic retry."""
    import psutil
    with ExitStack() as stack:
        lock_fds = []
        for name in ('.completion-controller.lock', '.research-run.lock'):
            lock = stack.enter_context((ROOT / 'experiments' / name).open('a'))
            try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError: raise RuntimeError('Another research engine/controller is active; no benchmark launched')
            lock_fds.append(lock.fileno())
        def terminate(*_):
            raise KeyboardInterrupt('Benchmark terminated')
        previous = signal.signal(signal.SIGTERM, terminate)
        stack.callback(signal.signal, signal.SIGTERM, previous)
        out = ROOT / 'runs/completion-mace-device-benchmark-r1'
        out.mkdir()  # Preserve any partial/failed attempt for review.
        shutil.copy2(PLAN, out / 'plan.json')
        shutil.copy2(__file__, out / 'runner_snapshot.py')
        state = dict(status='running', started_at=datetime.now(timezone.utc).isoformat(),
                     steps=[], plan_sha256=sha(PLAN), automatic_retries=0)
        wall_start = time.monotonic()
        try:
            for repeat, order in [(1, list(VARIANTS)), (2, list(reversed(VARIANTS)))]:
                for variant in order:
                    require(sha(PLAN) == state['plan_sha256'], 'Benchmark plan changed during execution')
                    verify_plan(plan)
                    target = out / f'{variant}-r{repeat}.json'
                    env = dict(os.environ, OMP_NUM_THREADS=str(VARIANTS[variant][2]),
                        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
                        PYTORCH_ENABLE_MPS_FALLBACK='0', PYTORCH_MPS_FAST_MATH='0',
                        MACE_BENCHMARK_LOCK_FDS=','.join(map(str, lock_fds)))
                    command = [str(ROOT / '.tools/mps-env/bin/python'), str(Path(__file__)),
                               'worker', '--variant', variant, '--output', str(target)]
                    tick, peak = time.monotonic(), 0
                    with target.with_suffix('.log').open('x') as log:
                        proc = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log,
                                                stderr=subprocess.STDOUT, pass_fds=lock_fds)
                        try:
                            while proc.poll() is None:
                                try: peak = max(peak, psutil.Process(proc.pid).memory_info().rss)
                                except psutil.NoSuchProcess: break
                                require(time.monotonic() - tick < POLICY['worker_timeout_seconds'], 'Worker timeout')
                                require(peak < POLICY['rss_limit_gib'] * 1024**3, 'Worker RSS limit')
                                time.sleep(.25)
                            require(proc.wait() == 0 and target.exists(), f'Failed worker: {variant} repeat {repeat}')
                        finally:
                            if proc.poll() is None:
                                proc.terminate()
                                try: proc.wait(timeout=10)
                                except subprocess.TimeoutExpired: proc.kill(); proc.wait()
                    state['steps'].append(dict(variant=variant, repeat=repeat,
                        output=str(target.relative_to(ROOT)), output_sha256=sha(target),
                        elapsed_seconds=time.monotonic() - tick, peak_worker_rss_bytes=peak))
                    (out / 'progress.json').write_text(json.dumps(state, indent=2) + '\n')
            state['status'] = 'completed'
        except BaseException as exc:
            state.update(status='stopped_requires_review', reason=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            state.update(ended_at=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - wall_start)
            dump_new(out / 'manifest.json', state)
    return out


def evaluate(plan, run_dir):
    state = json.loads((run_dir / 'manifest.json').read_text())
    require(state['status'] == 'completed' and state['plan_sha256'] == sha(PLAN)
            and len(state['steps']) == POLICY['maximum_workers'], 'Incomplete benchmark')
    workers = {}
    for row in state['steps']:
        path = ROOT / row['output']
        require(sha(path) == row['output_sha256'], 'Changed worker results')
        value = json.loads(path.read_text())
        key = (row['repeat'], row['variant'])
        device, dtype, threads = VARIANTS[row['variant']]
        require(key not in workers and value['status'] == 'completed' and value['variant'] == row['variant']
                and value['plan_sha256'] == sha(PLAN) and value['model_sha256'] == CHECKPOINT_SHA
                and value['threads'] == threads and value['device'] == device and value['dtype'] == dtype
                and value['fallback_enabled'] is False, 'Worker context mismatch')
        placement = value['actual_forward_output_placement']
        require(len(placement) == 2 and {x[0] for x in placement} == {'energy', 'forces'}
                and all(x[1].split(':')[0] == device and x[2] == 'torch.' + dtype for x in placement),
                'Missing forward tensor device/dtype evidence')
        for parts in value['results'].values():
            for result in parts.values():
                samples = result['samples']
                require(len(samples) == POLICY['repeats'] and len(result['warmup_seconds']) == POLICY['warmups'],
                        'Incorrect repeat count')
                times = [s['elapsed_seconds'] for s in samples]
                require(all(np.isfinite(t) and t > 0 for t in times)
                        and result['median_seconds'] == statistics.median(times), 'Invalid timing summary')
        workers[key] = value
    require(set(workers) == {(r, v) for r in (1, 2) for v in VARIANTS}, 'Missing or duplicate worker')
    comparisons = []
    for repeat in (1, 2):
        base = workers[(repeat, 'cpu64_4')]['results']
        for variant in VARIANTS:
            target = workers[(repeat, variant)]['results']
            # Every measured sample must pass, not just the most favorable one.
            sample_checks = []
            for index in range(POLICY['repeats']):
                sample = {n: {p: x['samples'][index] for p, x in parts.items()} for n, parts in target.items()}
                sample_checks.append(deltas(plan, base, sample))
            comparisons.append(dict(repeat=repeat, variant=variant, versus='cpu64_4',
                samples=sample_checks, passed=all(c['passed'] for s in sample_checks for c in s.values())))
        comparisons.append(dict(repeat=repeat, variant='mps32_4', versus='cpu32_4',
            final_sample_only=True, cases=deltas(plan, workers[(repeat, 'cpu32_4')]['results'],
                                               workers[(repeat, 'mps32_4')]['results'])))
    choices = {}
    for name in plan['cases']:
        times = {v: [sum(workers[(r, v)]['results'][name][p]['median_seconds'] for p in ('pair', 'bottom', 'top'))
                     for r in (1, 2)] for v in VARIANTS}
        eligible = [v for v in VARIANTS if all(c['passed'] for c in comparisons
                    if c['variant'] == v and c['versus'] == 'cpu64_4')]
        require('cpu64_4' in eligible, 'CPU reference is internally unstable')
        fastest_cpu = min((v for v in eligible if v.startswith('cpu')), key=lambda v: statistics.median(times[v]))
        speedups = [times[fastest_cpu][i] / times['mps32_4'][i] for i in (0, 1)]
        use_gpu = 'mps32_4' in eligible and min(speedups) >= POLICY['minimum_speedup']
        choices[name] = dict(cpu_candidate=fastest_cpu, seconds_per_pair_plus_two_layers=times,
            gpu_speedup_over_fastest_eligible_cpu_by_repeat=speedups,
            recommended_backend='mps32_4' if use_gpu else fastest_cpu,
            gpu_precision_passed='mps32_4' in eligible, gpu_speed_threshold_passed=min(speedups) >= POLICY['minimum_speedup'])
    return dict(status='device_comparison_complete', comparisons=comparisons, recommendations=choices,
        actual_model_inference_executed=True, dft_acceleration_verified=False,
        scientific_model_suitability_verified=False, settings_automatically_changed=False,
        total_benchmark_seconds=state['elapsed_seconds'], plan_sha256=sha(PLAN),
        manifest_sha256=sha(run_dir / 'manifest.json'), scope=plan['scope'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'run', 'evaluate', 'worker'])
    parser.add_argument('--variant', choices=VARIANTS)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.action == 'prepare': print(prepare()); return
    plan = json.loads(PLAN.read_text()); verify_plan(plan)
    if args.action == 'run': print(run(plan))
    elif args.action == 'evaluate':
        require(args.output is not None and not args.output.exists(), 'Use a new output path')
        result = evaluate(plan, ROOT / 'runs/completion-mace-device-benchmark-r1')
        dump_new(args.output, result)
        print(json.dumps(result['recommendations'], indent=2))
    else:
        require(args.variant is not None and args.output is not None, 'Worker requires variant/output')
        worker(plan, args.variant, args.output)


if __name__ == '__main__': main()
