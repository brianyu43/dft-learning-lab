"""Bounded local research jobs. All input paths and metadata come from reviewed JSON.

Unlike smoke.py, optimization and Hessian completion are evaluated independently.
No automatic restarts. --parse-run reparses existing outputs without running an engine.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import time

import psutil

ROOT = Path(__file__).resolve().parents[1]
NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def now():
    return datetime.now(timezone.utc).isoformat()


def swap_used():
    try:
        return psutil.swap_memory().used
    except (OSError, psutil.AccessDenied):
        return None


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def save(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(path)


def numeric_matches(pattern, text):
    return [float(x.replace('D', 'E').replace('d', 'e')) for x in re.findall(pattern, text)]


def parse_output(run, meta):
    text = (run / 'output.log').read_text(errors='replace')
    engine = meta['engine']
    wanted = meta['requested']
    normal = False
    if engine == 'orca':
        energies = numeric_matches(r'FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')', text)
        normal = 'ORCA TERMINATED NORMALLY' in text
        scf = 'SCF CONVERGED' in text and 'SCF NOT CONVERGED' not in text
        opt = 'THE OPTIMIZATION HAS CONVERGED' in text
        version = re.search(r'Program Version\s+(\S+)', text)
        unit = 'Hartree'
        # ORCA prints the projected translations/rotations first, followed by 3N-6 vibrations.
        block = text.rsplit('VIBRATIONAL FREQUENCIES', 1)[-1] if 'VIBRATIONAL FREQUENCIES' in text else ''
        block = block.split('NORMAL MODES', 1)[0]
        modes = [(int(i), float(f)) for i, f in re.findall(r'^\s*(\d+):\s+([-+]?\d+\.\d+)\s+cm', block, re.M)]
        natoms = meta.get('n_atoms')
        expected = (3 * natoms - (5 if meta.get('linear') else 6)) if natoms else None
        external = 5 if meta.get('linear') else 6
        prefix_valid = len(modes) == 3 * natoms and all(abs(f) <= 0.001 for _, f in modes[:external]) if natoms else False
        vib = [f for _, f in modes[external:]] if prefix_valid else []
        freq_complete = bool(expected and len(vib) == expected and normal)
        minimum = bool(opt and freq_complete and all(f > 0 for f in vib)) if wanted.get('frequencies') else None
        extra = {'frequencies_cm1': vib, 'all_reported_modes_cm1': modes,
                 'frequency_count_expected': expected, 'frequency_count_observed': len(vib),
                 'projected_external_modes_verified': prefix_valid,
                 'imaginary_frequencies_cm1': [f for f in vib if f < 0],
                 'lowest_vibrational_frequency_cm1': min(vib) if vib else None,
                 'minimum_verified': minimum, 'frequency_analysis_complete': freq_complete,
                 'optimization_cycles': [int(x) for x in re.findall(r'GEOMETRY OPTIMIZATION CYCLE\s+(\d+)', text)]}
    elif engine == 'qe':
        energies = numeric_matches(r'!\s+total energy\s+=\s+(' + NUMBER + ')\s+Ry', text)
        normal = 'JOB DONE.' in text
        scf = 'convergence has been achieved' in text and 'convergence NOT achieved' not in text
        opt = 'End of BFGS Geometry Optimization' in text or bool(re.search(r'bfgs converged in', text, re.I))
        version = re.search(r'Program PWSCF v\.([^\s]+)', text)
        unit = 'Ry'
        forces = numeric_matches(r'Total force\s*=\s*(' + NUMBER + ')', text)
        extra = {'minimum_verified': None, 'frequency_analysis_complete': False,
                 'total_force_Ry_Bohr': forces[-1] if forces else None,
                 'accuracy_scope': meta.get('accuracy_scope', 'week1 baseline; precision convergence not established')}
    else:
        raise ValueError(engine)
    finite = bool(energies and math.isfinite(energies[-1]))
    result = {'engine_version': version.group(1) if version else None,
              'normal_termination': normal, 'scf_converged': scf,
              'geometry_optimization_converged': opt if wanted.get('optimization') else None,
              'energy': energies[-1] if finite else None, 'energy_unit': unit,
              'energy_history': energies, **extra}
    checks = [meta.get('exit_code') == 0, not meta.get('stop_reason'), normal, scf, finite]
    if wanted.get('optimization'):
        checks.append(opt)
    if wanted.get('frequencies'):
        checks += [extra['frequency_analysis_complete'], extra['minimum_verified'] is True]
    result['status'] = 'passed' if all(checks) else 'needs_review'
    save(run / 'result.json', result)
    return result


def stop_group(proc):
    if proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
    except ProcessLookupError:
        proc.wait()


def run_job(config_path):
    config = json.loads(config_path.read_text())
    assert config['engine'] in ('orca', 'qe')
    assert config['timeout_seconds'] > 0
    executable = (ROOT / config['executable']).resolve()
    if not executable.is_file():
        raise FileNotFoundError(executable)
    input_path = (ROOT / config['input']).resolve()
    job_id = config['job_id']
    if not re.fullmatch(r'[a-zA-Z0-9_-]+', job_id):
        raise ValueError('invalid job_id')
    (ROOT / 'experiments').mkdir(exist_ok=True)
    with (ROOT / 'experiments/.research-run.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        run = ROOT / 'runs' / f'{stamp}-{job_id}'
        run.mkdir(parents=True)
        shutil.copy2(input_path, run / input_path.name)
        source_hashes = {config['input']: sha(input_path)}
        for source, destination in config.get('assets', {}).items():
            target = (run / destination).resolve()
            if not target.is_relative_to(run):
                raise ValueError('asset outside run')
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / source, target)
            source_hashes[source] = sha(ROOT / source)
        (run / 'scratch').mkdir()
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS=str(config.get('omp_threads', 4)), OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        env['PATH'] = str(executable.parent) + os.pathsep + env.get('PATH', '')
        args = [input_path.name] if config['engine'] == 'orca' else ['-in', input_path.name]
        command = [str(executable), *args]
        meta = {**config, 'started_at': now(), 'status': 'running', 'run_dir': str(run),
                'config_sha256': sha(config_path), 'source_sha256': source_hashes,
                'executable_sha256': sha(executable), 'command': command,
                'runner_sha256': sha(Path(__file__)), 'stop_reason': None,
                'thread_environment': {k: env[k] for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}}
        save(run / 'manifest.json', meta)
        print(json.dumps({'event': 'started', 'job_id': job_id, 'run_dir': str(run)}), flush=True)
        start = time.monotonic()
        peak = 0
        samples = 0
        denied = 0
        initial_swap = swap_used()
        peak_swap_delta = 0 if initial_swap is not None else None
        last_status = 0
        with (run / 'output.log').open('w') as out, (run / 'stderr.log').open('w') as err:
            proc = subprocess.Popen(command, cwd=run, env=env, stdout=out, stderr=err, start_new_session=True)
            meta['pid'] = proc.pid
            try:
                while proc.poll() is None:
                    elapsed = time.monotonic() - start
                    try:
                        parent = psutil.Process(proc.pid)
                        members = [parent, *parent.children(recursive=True)]
                        rss = 0
                        for member in members:
                            try:
                                rss += member.memory_info().rss
                            except psutil.NoSuchProcess:
                                pass
                        samples += 1
                        peak = max(peak, rss)
                        if rss > config.get('rss_limit_gib', 40) * 1024**3:
                            meta['stop_reason'] = 'process_tree_rss_limit'
                    except (psutil.AccessDenied, OSError):
                        denied += 1
                    except psutil.NoSuchProcess:
                        pass
                    current_swap = swap_used()
                    if initial_swap is not None and current_swap is not None:
                        peak_swap_delta = max(peak_swap_delta, current_swap - initial_swap)
                    if elapsed > config['timeout_seconds']:
                        meta['stop_reason'] = 'timeout'
                    if meta['stop_reason']:
                        stop_group(proc)
                        break
                    if elapsed - last_status >= 10:
                        save(run / 'progress.json', {'updated_at': now(), 'elapsed_seconds': round(elapsed, 2),
                             'peak_sampled_process_tree_rss_bytes': peak if samples else None,
                             'memory_samples': samples, 'memory_access_denied': denied})
                        last_status = elapsed
                    time.sleep(0.5)
            except (KeyboardInterrupt, BaseException):
                meta['stop_reason'] = 'interrupted_or_runner_error'
                stop_group(proc)
                raise
            finally:
                stop_group(proc)
                rss_os = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
                meta.update(ended_at=now(), elapsed_seconds=round(time.monotonic()-start, 3),
                            exit_code=proc.returncode, peak_sampled_process_tree_rss_bytes=peak if samples else None,
                            memory_samples=samples, memory_access_denied=denied,
                            max_single_child_rss_bytes=int(rss_os),
                            system_swap_increase_peak_bytes=peak_swap_delta,
                            memory_measurement='0.5s sampled summed process-tree RSS; may double-count shared pages and miss peaks. OS ru_maxrss is maximum single-child RSS in macOS bytes. Swap is system-wide, not attributable solely to this job.')
                meta['status'] = 'stopped' if meta['stop_reason'] else 'completed_pending_parse'
                save(run / 'manifest.json', meta)
        result = parse_output(run, meta)
        meta.update(status=result['status'], output_sha256=sha(run / 'output.log'), result=result)
        save(run / 'manifest.json', meta)
        save(ROOT / 'reports' / f'{job_id}.json', meta)
        print(json.dumps({'event': 'finished', 'job_id': job_id, 'run_dir': str(run),
                         'elapsed_seconds': meta['elapsed_seconds'], 'result': result}), flush=True)
        return 0 if result['status'] == 'passed' else 1


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--config', type=Path)
    group.add_argument('--parse-run', type=Path)
    args = parser.parse_args()
    if args.parse_run:
        meta = json.loads((args.parse_run / 'manifest.json').read_text())
        result = parse_output(args.parse_run, meta)
        meta.update(result=result, status=result['status'], reparsed_at=now(), parser_sha256=sha(Path(__file__)))
        save(args.parse_run / 'manifest.json', meta)
        save(ROOT / 'reports' / f"{meta['job_id']}.json", meta)
        print(json.dumps(result, indent=2))
        return 0 if result['status'] == 'passed' else 1
    return run_job(args.config.resolve())


if __name__ == '__main__':
    raise SystemExit(main())
