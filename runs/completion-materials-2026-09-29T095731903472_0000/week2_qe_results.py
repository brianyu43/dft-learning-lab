"""QE XML in explicit QE 7.6 atomic units, independent of ASE reader constants."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BOHR=0.529177210903
HARTREE=4.3597447222071e-18/1.602176634e-19

def load_job(job):
    meta=json.loads((ROOT/f'reports/{job}.json').read_text())
    assert meta['status']=='passed'
    run=Path(meta['run_dir'])
    assert hashlib.sha256((run/'output.log').read_bytes()).hexdigest()==meta['output_sha256']
    assert hashlib.sha256((run/Path(meta['input']).name).read_bytes()).hexdigest()==meta['source_sha256'][meta['input']]
    assert hashlib.sha256((run/'runner_snapshot.py').read_bytes()).hexdigest()==meta['runner_sha256']
    for source,dest in meta['assets'].items():
        assert hashlib.sha256((run/dest).read_bytes()).hexdigest()==meta['source_sha256'][source]
    paths=list((run/'scratch').glob('*.save/data-file-schema.xml'));assert len(paths)==1
    output=ET.parse(paths[0]).getroot().find('output');s=output.find('atomic_structure')
    symbols=[a.attrib['name'] for a in s.find('atomic_positions')]
    pos=np.array([np.fromstring(a.text,sep=' ') for a in s.find('atomic_positions')])*BOHR
    cell=np.array([np.fromstring(a.text,sep=' ') for a in s.find('cell')])*BOHR
    forces=np.fromstring(output.find('forces').text,sep=' ').reshape(-1,3)*HARTREE/BOHR
    energy_ha=float(output.find('total_energy/etot').text)
    actual_gap=float(cell[2,2]-np.ptp(pos[:,2]))
    assert abs(actual_gap-meta['vacuum_gap_angstrom'])<1e-8
    nelec=float(output.find('band_structure/nelec').text)
    assert abs(nelec-sum({'Mo':14,'S':6}[s] for s in symbols)+meta['charge'])<1e-8
    assert len(symbols)==meta['n_atoms'] and np.isfinite(pos).all() and np.isfinite(forces).all()
    assert abs(2*energy_ha-meta['result']['energy'])<1e-7
    return {'job_id':job,'meta':meta,'symbols':symbols,'positions':pos,'cell':cell,'forces':forces,
            'energy_eV':energy_ha*HARTREE,'max_force_eV_angstrom':float(np.linalg.norm(forces,axis=1).max()),
            'area_angstrom2':float(np.linalg.norm(np.cross(cell[0],cell[1]))),
            'actual_vacuum_gap_angstrom':actual_gap,'number_of_electrons':nelec,
            'xml_sha256':hashlib.sha256(paths[0].read_bytes()).hexdigest()}
