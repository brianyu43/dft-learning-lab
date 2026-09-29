"""Bounded on-demand QE SCFs for the fixed-registry optimizer.

The caller must supply a real eligibility review of reference, force calibration,
24-atom cost and trajectory budget. This component does not decide those gates.
It has no CLI and does not launch anything on import or construction.
"""
from __future__ import annotations

from copy import deepcopy
import fcntl
import hashlib
import io
import json
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

import numpy as np
import psutil
from ase.io import read, write
from ase.io.espresso import read_fortran_namelist

from completion_batch import spent
from completion_qe_bridge import load_scf, require
from completion_qe_relaxation import settings_fingerprint
from prepare_week3_drug import immutable
from research_run import ROOT, now, save, sha


def geometry_key(atoms):
    fields = dict(numbers=atoms.numbers.tolist(), positions=atoms.positions.tolist(),
                  cell=atoms.cell.array.tolist(), pbc=atoms.pbc.tolist(),
                  charges=atoms.get_initial_charges().tolist(),
                  magmoms=atoms.get_initial_magnetic_moments().tolist(), masses=atoms.get_masses().tolist())
    return hashlib.sha256(json.dumps(fields, sort_keys=True, allow_nan=False).encode()).hexdigest()


def render_from_reference(meta, atoms, job_id):
    """Keep reference namelists/pseudopotentials; change only run paths/positions."""
    require(re.fullmatch(r'[a-zA-Z0-9_-]+', job_id), 'Invalid job name')
    inp = Path(meta['run_dir']) / Path(meta['input']).name
    with inp.open() as handle:
        namelist, _ = read_fortran_namelist(handle)
    require(namelist['control']['calculation'] == 'scf', 'Reference must be SCF')
    require(namelist['control'].get('restart_mode', 'from_scratch') == 'from_scratch'
            and namelist['electrons'].get('startingpot', 'atomic') != 'file'
            and namelist['electrons'].get('startingwfc', 'atomic+random') != 'file',
            'Reference depends on a previous wavefunction or restart')
    data = deepcopy(namelist)
    data['control'].update(calculation='scf', prefix=job_id.replace('-', '_'),
                           pseudo_dir='./pseudo', outdir='./scratch', tprnfor=True, tstress=True)
    pseudos = {Path(dest).stem: Path(dest).name for dest in meta['assets'].values()
               if dest.startswith('pseudo/') and dest.endswith('.upf')}
    require(set(atoms.get_chemical_symbols()) == set(pseudos), 'Reference pseudo/species coverage differs')
    detached = atoms.copy(); detached.calc = None; detached.set_constraint()
    output = io.StringIO()
    write(output, detached, format='espresso-in', input_data=data, pseudopotentials=pseudos,
          kpts=tuple(meta['k_grid']), koffset=(0, 0, 0))
    text = output.getvalue()
    with io.StringIO(text) as handle:
        restored, _ = read_fortran_namelist(handle)
    for group in ('system', 'electrons'):
        require(dict(restored[group]) == dict(namelist[group]), 'Physical/numerical settings changed during serialization')
    recovered = read(io.StringIO(text), format='espresso-in')
    require(np.array_equal(recovered.numbers, atoms.numbers)
            and np.allclose(recovered.positions, atoms.positions, rtol=0, atol=1e-9)
            and np.allclose(recovered.cell, atoms.cell, rtol=0, atol=1e-9), 'Serialized SCF geometry differs')
    return text


class BoundedQESCFProvider:
    """Callable report supplier, only usable inside its exclusive context.

    validate_eligibility is called before opening the context and before each
    new SCF. It must recheck the scientific prerequisites in the calling plan.
    A completed provider context is not a converged optimization claim.
    """
    def __init__(self, reference_report, reference_atoms, campaign_id, *, validate_eligibility,
                 maximum_new_scfs, per_scf_timeout_seconds, total_seconds_limit, method='C09'):
        require(callable(validate_eligibility), 'An explicit eligibility validator is required')
        require(re.fullmatch(r'[a-zA-Z0-9_-]{1,64}', campaign_id), 'Invalid campaign name')
        require(type(maximum_new_scfs) is int and 0 < maximum_new_scfs <= 100, 'Invalid SCF count limit')
        require(type(per_scf_timeout_seconds) is int and 0 < per_scf_timeout_seconds <= 14400,
                'Invalid per-SCF time limit')
        require(type(total_seconds_limit) is int and per_scf_timeout_seconds <= total_seconds_limit <= 14400,
                'Invalid per-trajectory time limit')
        self.reference_report = Path(reference_report).resolve()
        require(self.reference_report.is_relative_to(ROOT), 'Reference lies outside the project')
        self.reference = reference_atoms.copy(); self.reference.calc = None; self.reference.set_constraint()
        initial = load_scf(self.reference_report, expected_method=method, expected_atoms=self.reference)
        self.meta = json.loads(self.reference_report.read_text())
        run = Path(self.meta['run_dir'])
        require(self.meta['exit_code'] == 0 and not self.meta.get('stop_reason'), 'Reference run incomplete')
        require(sha(run / 'job_config.json') == self.meta['config_sha256'], 'Reference config changed')
        self.config = json.loads((run / 'job_config.json').read_text())
        require(all(self.meta.get(k) == v for k, v in self.config.items()), 'Reference config/report mismatch')
        require(self.config['omp_threads'] == 4 and self.config['rss_limit_gib'] <= 40,
                'Unreviewed reference resource settings')
        require(sha(ROOT / self.config['executable']) == self.meta['executable_sha256'], 'QE executable changed')
        self.fingerprint = settings_fingerprint(self.reference_report)
        files = [self.reference_report, run / 'job_config.json', run / Path(self.meta['input']).name,
                 ROOT / self.config['executable'], ROOT / 'reports/week2-summary.json']
        files += [ROOT / p for p in self.config['assets']]
        files += [ROOT / 'scripts' / p for p in ('completion_qe_provider.py', 'completion_qe_relaxation.py',
                  'completion_qe_bridge.py', 'completion_batch.py', 'research_run.py', 'prepare_week3_drug.py')]
        for source, destination in self.config['assets'].items():
            require(sha(ROOT / source) == self.meta['source_sha256'][source], 'Current reference asset changed')
        self.frozen = {str(p.relative_to(ROOT)): sha(p) for p in files}
        self.cache = {geometry_key(self.reference): (self.reference_report, initial['audit']['report_sha256'])}
        self.campaign_id = campaign_id
        self.validate_eligibility = validate_eligibility
        self.maximum_new_scfs = maximum_new_scfs
        self.per_scf_timeout_seconds = per_scf_timeout_seconds
        self.total_seconds_limit = total_seconds_limit
        self.method = method
        self.directory = ROOT / 'runs' / ('completion-qe-relaxation-' + campaign_id)
        self.lock = None; self.proc = None; self.started = None; self.old_sigterm = None
        self.state = dict(campaign_id=campaign_id, status='not_started', steps=[],
                          source_sha256=self.frozen, maximum_new_scfs=maximum_new_scfs,
                          per_scf_timeout_seconds=per_scf_timeout_seconds, total_seconds_limit=total_seconds_limit,
                          automatic_retries=0, optimization_converged=None)

    def _verify_frozen(self):
        for path, expected in self.frozen.items():
            require(sha(ROOT / path) == expected, 'Changed frozen provider dependency: ' + path)

    def _save(self):
        self.state['updated_at'] = now()
        save(self.directory / 'state.json', self.state)

    def __enter__(self):
        require(self.lock is None and self.started is None, 'Provider contexts cannot be reused')
        self.validate_eligibility(); self._verify_frozen()
        lock = (ROOT / 'experiments/.completion-controller.lock').open('a')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            require(not self.directory.exists(), 'Existing provider attempt requires a new reviewed campaign')
            self.directory.mkdir()
            self.old_sigterm = signal.signal(signal.SIGTERM, self._terminated)
            self.lock = lock; self.started = time.monotonic()
            self.state.update(status='running', started_at=now()); self._save()
        except BaseException:
            lock.close()
            raise
        return self

    @staticmethod
    def _terminated(signum, frame):
        raise KeyboardInterrupt('QE provider terminated')

    def _stop_child(self):
        if self.proc is None or self.proc.poll() is not None:
            return
        try:
            members = psutil.Process(self.proc.pid).children(recursive=True)
        except psutil.Error:
            members = []
        try:
            self.proc.terminate()  # research_run handles SIGTERM and stops its QE group.
        except ProcessLookupError:
            return
        try:
            self.proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            for member in members:
                try: member.kill()
                except psutil.Error: pass
            self.proc.kill(); self.proc.wait(timeout=5)

    def __exit__(self, exc_type, exc, traceback):
        try:
            self._stop_child()
            incomplete = any(step['status'] != 'raw_validated' for step in self.state['steps'])
            self.state.update(status='closed' if exc is None and not incomplete else 'stopped_requires_review',
                              ended_at=now(), reason=None if exc is None else f'{type(exc).__name__}: {exc}')
            self._save()
        finally:
            signal.signal(signal.SIGTERM, self.old_sigterm)
            self.lock.close(); self.lock = None

    def _budget(self):
        require(self.lock is not None, 'Provider must hold its exclusive context')
        require(len(self.state['steps']) < self.maximum_new_scfs, 'SCF count budget exhausted')
        remaining = self.total_seconds_limit - (time.monotonic() - self.started)
        require(remaining >= self.per_scf_timeout_seconds, 'Trajectory time budget cannot cover another SCF')
        original = json.loads((ROOT / 'reports/week2-summary.json').read_text())
        pilot_remaining = original['remaining_budget_hours']['semiconductor'] * 3600 - spent()['qe']
        require(pilot_remaining >= self.per_scf_timeout_seconds, 'Remaining semiconductor pilot budget is insufficient')
        return dict(trajectory_remaining_seconds=remaining, pilot_remaining_seconds=pilot_remaining)

    def _launch(self, config_path):
        with (self.directory / 'controller.log').open('a') as log:
            command = [sys.executable, 'scripts/research_run.py', '--config', str(config_path)]
            self.proc = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            self.state['runner_pid'] = self.proc.pid; self._save()
            try:
                code = self.proc.wait(timeout=self.per_scf_timeout_seconds + 60)
            except BaseException:
                self._stop_child()
                raise
            finally:
                if self.proc.poll() is not None: self.state['runner_pid'] = None
        return code

    def __call__(self, atoms):
        require(self.lock is not None, 'Provider must hold its exclusive context')
        require(all(step['status'] == 'raw_validated' for step in self.state['steps']),
                'A prior SCF attempt is unresolved; no retry inside this campaign')
        self._verify_frozen()
        require(np.array_equal(atoms.numbers, self.reference.numbers)
                and np.array_equal(atoms.pbc, self.reference.pbc), 'Atom identities or periodicity changed')
        require(np.allclose(atoms.cell, self.reference.cell, rtol=0, atol=1e-10), 'Fixed cell changed')
        require(np.isfinite(atoms.positions).all(), 'Nonfinite requested coordinates')
        for getter in ('get_initial_charges', 'get_initial_magnetic_moments', 'get_masses'):
            require(np.array_equal(getattr(atoms, getter)(), getattr(self.reference, getter)()), 'Electronic state or masses changed')
        key = geometry_key(atoms)
        if key in self.cache:
            path, expected = self.cache[key]
            require(sha(path) == expected, 'Cached report changed')
            load_scf(path, expected_method=self.method, expected_atoms=atoms)
            require(settings_fingerprint(path) == self.fingerprint, 'Cached settings changed')
            return path
        self.validate_eligibility()
        budget = self._budget()
        job_id = f'{self.campaign_id}-scf-{len(self.state["steps"]):03d}'
        require(not (ROOT / f'reports/{job_id}.json').exists()
                and not list((ROOT / 'runs').glob('*-' + job_id)), 'Existing SCF attempt must not be retried')
        geometry = f'structures/completion/qe-relaxation/{job_id}.extxyz'
        detached = atoms.copy(); detached.calc = None; detached.set_constraint()
        stream = io.StringIO(); write(stream, detached, format='extxyz')
        inp = f'inputs/qe/completion/{job_id}.in'
        text = render_from_reference(self.meta, atoms, job_id)
        immutable(ROOT / geometry, stream.getvalue()); immutable(ROOT / inp, text)
        cfg = deepcopy(self.config)
        cfg.pop('computational_identity_sha256', None)
        cfg['assets'] = {source: dest for source, dest in cfg['assets'].items() if dest != 'initial.extxyz'}
        cfg['assets'].update({geometry: 'initial.extxyz',
                             str(self.reference_report.relative_to(ROOT)): 'provenance/reference-report.json',
                             'scripts/completion_qe_provider.py': 'provenance/provider_snapshot.py'})
        cfg.update(protocol_id='completion-constrained-qe-scf-r1', week=3, task_id='S06-SCF', job_id=job_id,
                   input=inp, geometry_source=geometry, requested=dict(optimization=False, frequencies=False),
                   timeout_seconds=self.per_scf_timeout_seconds, rss_limit_gib=40, omp_threads=4,
                   role='SCF_at_requested_constrained_trajectory_coordinates',
                   requested_geometry_sha256=key, initial_structure_sha256=sha(ROOT / geometry),
                   initial_cell_angstrom=atoms.cell.array.tolist(),
                   initial_vacuum_gap_angstrom=float(atoms.cell[2, 2] - np.ptp(atoms.positions[:, 2])),
                   reference_report_sha256=sha(self.reference_report),
                   accuracy_scope='One fixed-geometry SCF; optimization convergence requires a whole-trajectory audit')
        config_path = ROOT / f'experiments/completion/qe-relaxation/{job_id}.json'
        immutable(config_path, json.dumps(cfg, indent=2) + '\n')
        step = dict(job_id=job_id, config=str(config_path.relative_to(ROOT)), config_sha256=sha(config_path),
                    geometry_sha256=key, status='prepared', budget_at_launch=budget, started_at=now())
        self.state['steps'].append(step); self._save()
        code = self._launch(config_path)
        step.update(exit_code=code, ended_at=now(), status='completed_pending_raw_validation'); self._save()
        require(code == 0, 'SCF process failed; no automatic retry')
        path = ROOT / f'reports/{job_id}.json'
        meta = json.loads(path.read_text())
        require(meta['config_sha256'] == step['config_sha256'] and meta['exit_code'] == 0
                and not meta.get('stop_reason'), 'Executed SCF config/status mismatch')
        result = load_scf(path, expected_method=self.method, expected_atoms=atoms)
        require(settings_fingerprint(path) == self.fingerprint, 'SCF physical/numerical settings changed')
        step.update(status='raw_validated', audit=result['audit']); self._save()
        self.cache[key] = (path, sha(path))
        return path
