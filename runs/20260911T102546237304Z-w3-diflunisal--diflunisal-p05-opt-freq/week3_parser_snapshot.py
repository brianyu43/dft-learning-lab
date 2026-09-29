"""Week-three ORCA verifier v2: monomers have no separated_geometry field."""
import argparse
import json
from pathlib import Path
import signal

import numpy as np
from ase.io import read

import research_run
from run_week2_xtb_screen import geometry_audit
from run_week2_frequency_resume import inspect_hessian

ROOT = Path(__file__).resolve().parents[1]
BASE_PARSE = research_run.parse_output


def checked_parse(run, meta):
    result = BASE_PARSE(run, meta)
    result['week3_verification_passed'] = False
    try:
        candidate = json.loads((run / 'candidate.json').read_text())
        assert research_run.sha(run / 'candidate.json') == meta['candidate_sha256']
        assert research_run.sha(run / 'initial.xyz') == meta['initial_structure_sha256']
        assert research_run.sha(run / 'week3_parser_snapshot.py') == meta['week3_parser_sha256']
        assert research_run.sha(run / Path(meta['input']).name) == meta['source_sha256'][meta['input']]
        atoms = read(run / (Path(meta['input']).stem + '.xyz'))
        initial = read(run / 'initial.xyz')
        assert len(atoms) == len(initial) == meta['n_atoms']
        assert atoms.get_chemical_symbols() == initial.get_chemical_symbols()
        assert np.isfinite(atoms.positions).all()
        audit = geometry_audit(candidate, atoms)
        initial_audit = geometry_audit(candidate, initial)
        assert audit['heavy_connectivity_preserved']
        assert not audit['unassigned_hydrogens'] and not audit.get('separated_geometry', False)
        result.update(geometry_audit=audit, initial_geometry_audit=initial_audit,
                      optimized_xyz_sha256=research_run.sha(run / (Path(meta['input']).stem + '.xyz')))
        if meta['requested'].get('frequencies'):
            hess = run / (Path(meta['input']).stem + '.hess')
            frequencies = inspect_hessian(hess, atoms, result['frequencies_cm1'])
            assert np.all(frequencies > 0), 'Nonpositive internal vibration requires review'
            result.update(hessian_sha256=research_run.sha(hess), hessian_geometry_matches=True,
                          full_precision_frequencies_cm1=frequencies.tolist())
        assert result['status'] == 'passed', 'Engine or minimum validation incomplete'
        result['week3_verification_passed'] = True
    except (AssertionError, ValueError, KeyError, FileNotFoundError, IndexError) as exc:
        result.update(status='needs_review', minimum_verified=False,
                      week3_validation_error=f'{type(exc).__name__}: {exc}')
    result['week3_parser_sha256'] = meta['week3_parser_sha256']
    research_run.save(run / 'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    assert config['week'] == 3 and config['engine'] == 'orca'
    assert config['mpi_processes'] == 12 and config['omp_threads'] == 1
    assert config['week3_parser_sha256'] == research_run.sha(Path(__file__))
    assert not (ROOT / f"reports/{config['job_id']}.json").exists(), 'Review existing attempt before retry'
    research_run.parse_output = checked_parse
    return research_run.run_job(args.config)


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt('terminated')))
    raise SystemExit(main())
