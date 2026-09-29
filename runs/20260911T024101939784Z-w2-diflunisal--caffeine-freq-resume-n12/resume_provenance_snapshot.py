"""Validate the completed optimization independently of its interrupted Hessian."""
import hashlib
import json
from pathlib import Path
import re

from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_checkpoint(path):
    checkpoint=json.loads(path.read_text())
    assert checkpoint['optimization_converged'] is True
    for source,expected in checkpoint['source_files_sha256'].items():
        assert sha(ROOT/source)==expected, f'Changed checkpoint evidence: {source}'
    xyz=ROOT/checkpoint['optimized_geometry']
    assert sha(xyz)==checkpoint['geometry_sha256']
    run=ROOT/checkpoint['source_run_dir']
    original=json.loads((run/'job_config.json').read_text())
    assert original['method']=='r2SCAN-3c' and original['charge']==0 and original['multiplicity']==1
    assert original['requested']['optimization'] is True and original['n_atoms']==50
    source_xyz=run/(Path(original['input']).stem+'.xyz')
    assert sha(source_xyz)==sha(xyz)
    atoms=read(xyz);assert len(atoms)==50
    log=(run/'output.log').read_text()
    assert 'THE OPTIMIZATION HAS CONVERGED' in log
    assert max(map(int,re.findall(r'GEOMETRY OPTIMIZATION CYCLE\s+(\d+)',log)))==87
    energy=float(re.search(r' E ([-\d.]+)',xyz.read_text().splitlines()[1]).group(1))
    final=float(re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',log)[-1])
    assert abs(energy-final)<1e-10
    return dict(checkpoint=str(path.relative_to(ROOT)),checkpoint_sha256=sha(path),
                optimized_geometry=str(xyz.relative_to(ROOT)),geometry_sha256=sha(xyz),
                optimization_log=str((run/'output.log').relative_to(ROOT)),optimization_log_sha256=sha(run/'output.log'),
                source_run_dir=checkpoint['source_run_dir'],optimization_converged=True,
                optimization_cycles=87,optimization_energy_Hartree=energy,
                method='r2SCAN-3c',charge=0,multiplicity=1,n_atoms=50,
                frequency_or_minimum_verified=False)
