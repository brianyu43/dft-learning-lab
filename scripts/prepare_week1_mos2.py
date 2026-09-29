"""Freeze a first fixed-cell monolayer relaxation; this is not a convergence study."""
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
PSEUDOS = 'data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'
table = json.loads((ROOT / 'data/literature/week1/nc-sr-04_pbe_standard.json').read_text())
headers = {}
for element in ('Mo', 'S', 'Se'):
    content = (ROOT / PSEUDOS / f'{element}.upf').read_text()
    header = ET.fromstring(re.search(r'<PP_HEADER.*?/>', content, re.S).group(0)).attrib
    assert header['element'].strip() == element
    assert header['functional'] == 'PBE' and header['pseudo_type'] == 'NC'
    assert header['relativistic'] == 'scalar' and header['has_so'] == 'F'
    headers[element] = header
atoms = read(ROOT / 'structures/materials/MoS2_initial.cif')
assert atoms.get_chemical_formula() == 'MoS2' and len(atoms) == 3
# Official hints are Hartree, whereas QE uses Ry. Mo normal hint 40 Ha => 80 Ry.
cutoff = 2 * max(table[e]['hn'] for e in ('Mo', 'S'))
input_file = ROOT / 'inputs/qe/week1_mos2_relax.in'
write(input_file, atoms, format='espresso-in',
      pseudopotentials={'Mo': 'Mo.upf', 'S': 'S.upf'}, kpts=(6, 6, 1), koffset=(0, 0, 0),
      input_data={
          'control': {'calculation': 'relax', 'prefix': 'week1_mos2', 'pseudo_dir': './pseudo',
                      'outdir': './scratch', 'tprnfor': True, 'tstress': True, 'nstep': 80,
                      'etot_conv_thr': 1e-6, 'forc_conv_thr': 3e-4},
          'system': {'ecutwfc': cutoff, 'ecutrho': 4*cutoff, 'input_dft': 'vdw-df-c09',
                     'occupations': 'fixed', 'nosym': False, 'nspin': 1},
          'electrons': {'conv_thr': 1e-8, 'electron_maxstep': 150, 'mixing_beta': 0.4},
          'ions': {'ion_dynamics': 'bfgs'}})
config = {
    'protocol_id': 'dft-pilot-20260909-v1', 'week': 1, 'job_id': 'week1-mos2-relax',
    'task_id': 'SR01-baseline', 'engine': 'qe', 'executable': '.tools/qe-7.6/bin/pw.x',
    'input': str(input_file.relative_to(ROOT)),
    'assets': {f'{PSEUDOS}/{e}.upf': f'pseudo/{e}.upf' for e in ('Mo', 'S')},
    'method': 'vdW-DF-C09 with PseudoDojo 0.4 standard PBE-generated scalar NC pseudopotentials',
    'purpose': 'First MoS2 monolayer internal-coordinate relaxation at a fixed initial cell',
    'accuracy_scope': 'Numerical baseline only. Cutoff, k-grid, vacuum, lattice and phonon convergence not established. PBE-generated pseudos with C09 valence XC are an explicit approximation; exact literature pseudo revision is unverified.',
    'n_atoms': 3, 'charge': 0, 'SOC': False, 'lattice_a_angstrom': 3.18,
    'cell_c_angstrom': float(atoms.cell[2,2]), 'initial_vacuum_gap_angstrom': 20.0,
    'ecutwfc_Ry': cutoff, 'ecutrho_Ry': 4*cutoff, 'k_grid': [6,6,1], 'k_offset': [0,0,0],
    'cutoff_provenance': {'table': 'data/literature/week1/nc-sr-04_pbe_standard.json',
                          'hints_unit': 'Hartree', 'Mo_normal': table['Mo']['hn'], 'S_normal': table['S']['hn']},
    'pseudopotential_headers': headers,
    'requested': {'optimization': True, 'frequencies': False},
    'timeout_seconds': 43200, 'rss_limit_gib': 40, 'omp_threads': 4}
(ROOT / 'experiments/week1/mos2-relax.json').write_text(json.dumps(config, indent=2)+'\n')
print(f'Prepared MoS2: {cutoff:g}/{4*cutoff:g} Ry, k=6x6x1, fixed a=3.18 A, C09')
