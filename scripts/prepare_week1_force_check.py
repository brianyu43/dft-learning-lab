"""One targeted stricter SCF after QE's force-correction warning, not a full sweep."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[1]
BOHR_QE_76 = 0.529177210903  # Installed QE Modules/constants.f90, Angstrom.
meta = json.loads((ROOT / 'reports/week1-mos2-relax.json').read_text())
assert meta['status'] == 'passed'
parent = Path(meta['run_dir'])
source = parent / 'scratch/week1_mos2.save/data-file-schema.xml'
structure = ET.parse(source).getroot().find('output/atomic_structure')
cell = np.array([np.fromstring(row.text, sep=' ') for row in structure.find('cell')]) * BOHR_QE_76
sites = list(structure.find('atomic_positions'))
positions = np.array([np.fromstring(row.text, sep=' ') for row in sites]) * BOHR_QE_76
atoms = Atoms([row.attrib['name'] for row in sites], positions=positions, cell=cell, pbc=True)
write(ROOT / 'inputs/qe/week1_mos2_force_check.in', atoms, format='espresso-in',
      pseudopotentials={'Mo': 'Mo.upf', 'S': 'S.upf'}, kpts=(6,6,1), koffset=(0,0,0),
      input_data={'control': {'calculation': 'scf', 'prefix': 'week1_mos2_check',
                              'pseudo_dir': './pseudo', 'outdir': './scratch', 'tprnfor': True},
                  'system': {'ecutwfc': 80, 'ecutrho': 320, 'input_dft': 'vdw-df-c09',
                             'occupations': 'fixed', 'nspin': 1},
                  'electrons': {'conv_thr': 1e-10, 'mixing_beta': 0.4, 'electron_maxstep': 150}})
config = {k: meta[k] for k in ('protocol_id','week','engine','executable','assets','method','n_atoms','charge','SOC','omp_threads','rss_limit_gib')}
config.update(job_id='week1-mos2-force-check', task_id='SR01-force-validation',
              input='inputs/qe/week1_mos2_force_check.in', timeout_seconds=1800,
              parent_run_id=parent.name,
              purpose='Targeted force/energy sensitivity check at unchanged final relaxed geometry after SCF-correction warning',
              requested={'optimization': False, 'frequencies': False},
              accuracy_scope='Only conv_thr 1e-8 to 1e-10 Ry at one fixed geometry; cutoff/k/vacuum/lattice/phonon convergence remain unverified',
              parent_geometry_source=str(source.relative_to(ROOT)),
              ecutwfc_Ry=80, ecutrho_Ry=320, k_grid=[6,6,1], conv_thr_Ry=1e-10)
(ROOT / 'experiments/week1/mos2-force-check.json').write_text(json.dumps(config, indent=2)+'\n')
print('Prepared one fixed-geometry SCF force check, conv_thr=1e-10 Ry')
