"""Bounded local MPI launcher diagnostics; does not run a chemistry engine."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--loopback', action='store_true')
args = parser.parse_args()
mpi = ROOT / '.tools/orca-mpi'
env = os.environ.copy()
env['PATH'] = str(mpi / 'bin') + os.pathsep + env['PATH']
env['DYLD_LIBRARY_PATH'] = str(mpi / 'lib')
commands = [[str(mpi / 'bin/mpirun'), '--version']]
opts = ['--host', 'localhost:18', '--mca', 'oob_tcp_if_include', 'lo0',
        '--mca', 'btl_tcp_if_include', 'lo0'] if args.loopback else []
commands.append([str(mpi / 'bin/mpirun'), *opts, '--bind-to', 'none', '-np', '2', '/usr/bin/true'])
checks = []
for command in commands:
    started = time.monotonic()
    proc = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, start_new_session=True)
    timeout = False
    try:
        out, err = proc.communicate(timeout=20)
    except subprocess.TimeoutExpired:
        timeout = True
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate()
    checks.append(dict(command=command, exit_code=proc.returncode, timeout=timeout,
                       elapsed_seconds=time.monotonic()-started, stdout=out, stderr=err))
result = dict(observed_at=datetime.now(timezone.utc).isoformat(),
              status='passed' if all(c['exit_code']==0 for c in checks) else 'requires_review',
              scope='MPI launcher only; no ORCA energy/gradient validation',
              openmpi_version='4.1.1', loopback=args.loopback, checks=checks)
suffix = 'loopback' if args.loopback else 'default'
(ROOT/f'reports/week2-mpi-launcher-{suffix}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['status']=='passed' else 1)
