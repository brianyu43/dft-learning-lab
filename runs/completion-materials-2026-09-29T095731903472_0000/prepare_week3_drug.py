"""Prepare immutable week-three jobs; preparation does not launch an engine."""
import argparse
import json
from pathlib import Path

from ase.io import read
from research_run import sha
from run_week2_xtb_screen import geometry_audit

ROOT = Path(__file__).resolve().parents[1]
LIBS = ['libmpi.40.dylib', 'libopen-rte.40.dylib', 'libopen-pal.40.dylib',
        'libopen-orted-mpir.dylib', 'libz.1.dylib']


def immutable(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_text() == text, f'Preserve and review existing file: {path}'
    else:
        path.write_text(text)


def prepare(candidate_id, direct=False, role='primary'):
    kind = 'dimers' if '--' in candidate_id else 'monomers'
    rows = json.loads((ROOT / f'reports/week2-xtb-{kind}.json').read_text())['results']
    row = next(r for r in rows if r['candidate_id'] == candidate_id)
    candidate = row['candidate']
    source = ROOT / candidate['path'] if direct else ROOT / row['run_dir'] / 'xtbopt.xyz'
    if not direct:
        assert row['result']['status'] == 'passed'
    atoms = read(source)
    audit = geometry_audit(candidate, atoms)
    assert audit['heavy_connectivity_preserved'] and not audit['unassigned_hydrogens'] and not audit.get('separated_geometry', False)
    suffix = '-direct-ion' if direct else ''
    job = f'w3-{candidate_id}{suffix}-opt-freq'
    cfgpath = ROOT / f'experiments/week3/drug/{job}.json'
    candidate_path = ROOT / f'experiments/week3/drug/candidates/{candidate_id}.json'
    immutable(candidate_path, json.dumps(candidate, indent=2) + '\n')
    initial = ROOT / f'structures/week3/drug_initial/{candidate_id}{suffix}.xyz'
    immutable(initial, source.read_text())
    inp = ROOT / f'inputs/orca/week3/{job}.inp'
    text = '! r2SCAN-3c TightSCF TightOpt Opt Freq\n%maxcore 2000\n%pal nprocs 12 end\n\n* xyzfile 0 1 initial.xyz\n'
    immutable(inp, text)
    assets = {f'.tools/orca-mpi/lib/{lib}': lib for lib in LIBS}
    assets.update({str(initial.relative_to(ROOT)): 'initial.xyz',
                   str(candidate_path.relative_to(ROOT)): 'candidate.json',
                   'scripts/run_week3_drug_v2.py': 'week3_parser_snapshot.py',
                   'scripts/run_week2_xtb_screen.py': 'provenance/geometry_audit_snapshot.py',
                   'scripts/run_week2_frequency_resume.py': 'provenance/hessian_check_snapshot.py',
                   'scripts/orca_parallel.sh': 'provenance/parallel_launcher_snapshot.sh',
                   'environment/orca-mpi-explicit.txt': 'provenance/orca-mpi-explicit.txt'})
    config = dict(protocol_id='dft-pilot-20260909-v1-week3', week=3, task_id='DR03', job_id=job,
                  engine='orca', executable='scripts/orca_parallel.sh', input=str(inp.relative_to(ROOT)),
                  method='r2SCAN-3c', charge=0, multiplicity=1, solvent='vacuum', n_atoms=len(atoms), linear=False,
                  requested={'optimization': True, 'frequencies': True}, mpi_processes=12, omp_threads=1,
                  timeout_seconds=14400, rss_limit_gib=40, assets=assets,
                  candidate_id=candidate_id, species=row['species'], role=role,
                  source_geometry=str(source.relative_to(ROOT)), source_geometry_sha256=sha(source),
                  initial_structure_sha256=sha(initial), candidate_sha256=sha(candidate_path),
                  week3_parser_sha256=sha(ROOT / 'scripts/run_week3_drug_v2.py'),
                  state_scope='Starting chemical graph and actual hydrogen ownership recorded; isolated 1:1 dimer or monomer, not crystal stoichiometry',
                  source_stage='direct proton-transfer initial geometry' if direct else 'GFN2-xTB optimized candidate',
                  selection_uses_formation_labels=False, initial_geometry_audit=audit,
                  purpose='Week-three vacuum contact/control comparison; no crystal formation or polymorph claim')
    immutable(cfgpath, json.dumps(config, indent=2) + '\n')
    return str(cfgpath.relative_to(ROOT))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate_id')
    parser.add_argument('--direct', action='store_true')
    parser.add_argument('--role', default='primary')
    args = parser.parse_args()
    print(prepare(args.candidate_id, args.direct, args.role))
