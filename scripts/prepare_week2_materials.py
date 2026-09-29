"""Generate seven orthogonal monolayer sensitivity calculations from week1 geometry."""
import json
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/week2'
OUT.mkdir(parents=True,exist_ok=True)
(ROOT/'inputs/qe/week2').mkdir(parents=True,exist_ok=True)
(ROOT/'structures/week2/materials').mkdir(parents=True,exist_ok=True)
PROTOCOL='dft-pilot-20260909-v1-week2'
PSEUDOS='data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'
base=read(ROOT/'structures/optimized/week1/MoS2_fixed_cell_relaxed.cif')

def write_job(job_id,atoms,cut,k,vac,**extra):
    prefix=job_id.replace('-','_')
    inp=ROOT/f'inputs/qe/week2/{job_id}.in'
    write(inp,atoms,format='espresso-in',pseudopotentials={'Mo':'Mo.upf','S':'S.upf'},kpts=(k,k,1),koffset=(0,0,0),
          input_data={'control':{'calculation':'scf','prefix':prefix,'pseudo_dir':'./pseudo','outdir':'./scratch','tprnfor':True},
                      'system':{'ecutwfc':cut,'ecutrho':4*cut,'input_dft':'vdw-df-c09','occupations':'fixed','nspin':1},
                      'electrons':{'conv_thr':1e-10,'electron_maxstep':200,'mixing_beta':0.4}})
    cfg={'protocol_id':PROTOCOL,'week':2,'job_id':job_id,'task_id':'SR01','engine':'qe','executable':'.tools/qe-7.6/bin/pw.x',
         'input':str(inp.relative_to(ROOT)),'assets':{f'{PSEUDOS}/{e}.upf':f'pseudo/{e}.upf' for e in ('Mo','S')},
         'method':'vdW-DF-C09; PseudoDojo nc-sr-04 standard PBE scalar NC',
         'purpose':'week2 numerical sensitivity at fixed geometry','requested':{'optimization':False,'frequencies':False},
         'n_atoms':len(atoms),'charge':0,'SOC':False,'ecutwfc_Ry':cut,'ecutrho_Ry':4*cut,'k_grid':[k,k,1],
         'vacuum_gap_angstrom':vac,'conv_thr_Ry':1e-10,'lattice_a_angstrom':float(atoms.cell.lengths()[0]),
         'timeout_seconds':43200,'rss_limit_gib':40,'omp_threads':4,
         'accuracy_scope':'Numerical sensitivity only at this composition/geometry; exact literature pseudo revision unverified',**extra}
    path=OUT/f'{job_id}.json';path.write_text(json.dumps(cfg,indent=2)+'\n')
    write(ROOT/f'structures/week2/materials/{job_id}.extxyz',atoms)
    return str(path.relative_to(ROOT))

def mono(vac):
    atoms=base.copy();atoms.center(vacuum=vac/2,axis=2)
    assert abs(atoms.cell[2,2]-np.ptp(atoms.positions[:,2])-vac)<1e-10
    return atoms

if __name__=='__main__':
    rows=[(80,12,30,'cutoff'),(100,12,30,'cutoff'),(120,12,30,'anchor'),(120,6,30,'k_grid'),(120,9,30,'k_grid'),(120,12,20,'vacuum'),(120,12,25,'vacuum')]
    jobs=[write_job(f'w2-mono-c{c}-k{k}-v{v}',mono(v),c,k,v,scan_axis=axis,configuration_id='MoS2-monolayer') for c,k,v,axis in rows]
    (OUT/'materials-initial-batch.json').write_text(json.dumps({'batch_id':'week2-mono-convergence','budget_seconds':14400,'jobs':jobs},indent=2)+'\n')
    print('Prepared',len(jobs),'monolayer SCFs')
