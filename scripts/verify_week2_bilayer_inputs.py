"""Read-only geometric and settings validation before or after the bilayer runs."""
import json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.io.espresso import read_fortran_namelist

ROOT=Path(__file__).resolve().parents[1]
batch=json.loads((ROOT/'experiments/week2/bilayer-batch.json').read_text());atoms={}
for path in batch['jobs']:
    c=json.loads((ROOT/path).read_text());a=read(ROOT/c['input'],format='espresso-in');atoms[c['job_id']]=a
    assert len(a)==c['n_atoms']
    assert abs(a.cell[2,2]-np.ptp(a.positions[:,2])-c['vacuum_gap_angstrom'])<1e-8
    with (ROOT/c['input']).open() as handle:nml,extra=read_fortran_namelist(handle)
    assert nml['control']['calculation']=='scf' and nml['system']['nspin']==1
    assert nml['system'].get('tot_charge',0)==c['charge']==0
    assert nml['system']['ecutwfc']==c['ecutwfc_Ry'] and nml['system']['ecutrho']==c['ecutrho_Ry']
    assert nml['electrons']['conv_thr']==c['conv_thr_Ry']
    kline=extra.index('K_POINTS automatic');assert list(map(int,extra[kline+1].split()))==c['k_grid']+[0,0,0]
comparisons=json.loads((ROOT/'experiments/week2/bilayer-comparisons.json').read_text())['configurations']
for c in comparisons:
    pair=atoms[c['jobs']['pair']]
    for key,slice_ in [('bottom',slice(0,3)),('top',slice(3,6))]:
        a=atoms[c['jobs'][key]]
        np.testing.assert_allclose(a.cell,pair.cell,atol=1e-10,rtol=0)
        np.testing.assert_allclose(a.positions,pair.positions[slice_],atol=1e-10,rtol=0)
    assert abs(pair.positions[3,2]-pair.positions[0,2]-c['Mo_plane_distance_angstrom'])<1e-8
    shift=c['registry_fraction'][0]*pair.cell[0]+c['registry_fraction'][1]*pair.cell[1]
    shift[2]+=c['Mo_plane_distance_angstrom']
    np.testing.assert_allclose(pair.positions[3:]-pair.positions[:3],np.tile(shift,(3,1)),atol=1e-8,rtol=0)
    np.testing.assert_allclose(pair.positions[1,:2],pair.positions[2,:2],atol=1e-8,rtol=0)
    assert abs(pair.positions[1,2]+pair.positions[2,2]-2*pair.positions[0,2])<1e-8
result={'status':'passed','unique_inputs':len(batch['jobs']),'matched_comparisons':len(comparisons),
        'charge_and_spin':'0, nonspinpolarized','scope':'input validation; not an engine result',
        'geometry_policy':'fixed same-orientation layers; matched frozen single-layer references',
        'registry_and_horizontal_mirror_symmetry_verified':True,
        'vacuum_metadata':'actual empty gap per system; isolated references retain the bilayer cell and have a larger empty gap'}
(ROOT/'reports/week2-bilayer-input-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
