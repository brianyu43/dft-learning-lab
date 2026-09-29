"""Prepare finite pure-bilayer controls after the monolayer scan is assessed."""
import argparse
import hashlib
import json
import numpy as np
from ase import Atoms
from ase.io import read
from prepare_week2_materials import ROOT,OUT,write_job

parser=argparse.ArgumentParser();parser.add_argument('--cutoff',type=int,required=True)
parser.add_argument('--k',type=int,required=True);parser.add_argument('--vacuum',type=int,required=True)
args=parser.parse_args()
analysis=json.loads((ROOT/'reports/week2-monolayer-convergence.json').read_text())
assert analysis['status']=='initial_axes_passed_combined_setting_pending','Review monolayer axes before bilayer preparation'
combined=json.loads((ROOT/'reports/week2-monolayer-combined.json').read_text())
assert combined['status']=='passed','Combined monolayer setting must pass before bilayer use'
assert combined['selected_setting']=={'cutoff':args.cutoff,'k_grid':args.k,'vacuum':args.vacuum}
base=read(ROOT/'structures/optimized/week1/MoS2_fixed_cell_relaxed.cif')
base.positions[:,2]-=base.positions[0,2]
height=np.ptp(base.positions[:,2]);a1,a2=base.cell.array[:2]
registry=[(0,0),(1/3,2/3),(2/3,1/3)]
jobs=[];records=[];seen={}

def make_bilayer(r,d,v):
    cell=base.cell.array.copy();cell[2,2]=d+height+v
    bottom=base.positions.copy();top=base.positions.copy()
    bottom[:,2]+=cell[2,2]/2-d/2;top[:,2]+=cell[2,2]/2+d/2
    top+=r[0]*a1+r[1]*a2
    atoms=Atoms(['Mo','S','S']*2,positions=np.vstack([bottom,top]),cell=cell,pbc=True)
    assert abs(atoms.positions[3,2]-atoms.positions[0,2]-d)<1e-10
    assert abs(cell[2,2]-np.ptp(atoms.positions[:,2])-v)<1e-10
    return atoms

def add_atoms(atoms,cut,k,v,label,**metadata):
    # Canonical scientific identity includes coordinates and all physical settings.
    # Reuse below occurs only when the same full existing input file is referenced again.
    key=hashlib.sha256(json.dumps({'symbols':atoms.get_chemical_symbols(),'positions':atoms.positions.tolist(),
                                  'cell':atoms.cell.array.tolist(),'cutoff':cut,'k':k,'vac':v},sort_keys=True).encode()).hexdigest()
    if key in seen:return seen[key]
    job=f'w2-bilayer-{label}'
    actual_gap=float(atoms.cell[2,2]-np.ptp(atoms.positions[:,2]))
    path=write_job(job,atoms,cut,k,actual_gap,configuration_id=label,
                   reference_bilayer_vacuum_gap_angstrom=v,**metadata)
    jobs.append(path);seen[key]=job
    return job

settings=[('selected',args.cutoff,args.k,args.vacuum),('strict',120,12,30)]
for setting,cut,k,v in settings:
    geometries=[(0,6.3)] if setting=='strict' else [(0,6.0),(0,6.3),(0,6.6),(1,6.3),(2,6.3)]
    for rid,d in geometries:
        tag=f'{setting}-r{rid}-d{d:.1f}'.replace('.','p');atoms=make_bilayer(registry[rid],d,v)
        ref={}
        for part,geometry in [('pair',atoms),('bottom',atoms[:3]),('top',atoms[3:])]:
            ref[part]=add_atoms(geometry,cut,k,v,f'{tag}-{part}',
                                calculation_role=part,registry_fraction=list(registry[rid]),Mo_plane_distance_angstrom=d,
                                orientation='R-type: identical in-plane orientation, no twist',comparison_setting=setting)
        records.append({'configuration':tag,'setting':setting,'registry_id':rid,'registry_fraction':list(registry[rid]),
                        'Mo_plane_distance_angstrom':d,'jobs':ref,'area_angstrom2':float(np.linalg.norm(np.cross(a1,a2)))})
(OUT/'bilayer-batch.json').write_text(json.dumps({'batch_id':'week2-bilayer-control','budget_seconds':14400,'jobs':jobs},indent=2)+'\n')
(OUT/'bilayer-comparisons.json').write_text(json.dumps({'configurations':records,'deduplication':'Reference exactly the existing immutable job input when scientific identity matches',
                                                    'equilibrium_in_plane_lattice_verified':False,'reference_pseudo_revision_matched':False},indent=2)+'\n')
print(json.dumps({'jobs':len(jobs),'comparisons':len(records)}))
