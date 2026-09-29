"""Freeze the first three DFT jobs from valid xTB results, without outcome labels."""
import hashlib
import json
from pathlib import Path
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
mono=json.loads((ROOT/'reports/week2-xtb-monomers.json').read_text())
dim=json.loads((ROOT/'reports/week2-xtb-dimers.json').read_text())
assert mono['status']=='completed'
inventory=json.loads((ROOT/'experiments/week2/chemistry-candidates.json').read_text())
expected={r['candidate_id'] for r in inventory['dimers'] if r['species']=='diflunisal--caffeine'}
observed={r['candidate_id'] for r in dim['results'] if r['species']=='diflunisal--caffeine'}
assert expected==observed, 'Complete every DIF-CAF start before selecting the first DFT geometry'

def usable(row,require_same_state=False):
    r=row['result'];a=r.get('geometry_audit') or {}
    return (r['status']=='passed' and a.get('heavy_connectivity_preserved') and not a.get('unassigned_hydrogens')
            and not a.get('separated_geometry') and (a.get('chemical_state_preserved_by_geometry') or not require_same_state))

selected=[]
for species,state in [('diflunisal','diflunisal-t02'),('caffeine','caffeine-t00')]:
    pool=[r for r in mono['results'] if r['species']==species and r['candidate']['state_id']==state and usable(r,True)]
    assert pool,f'No valid first monomer for {species}'
    selected.append(min(pool,key=lambda x:x['result']['energy_Hartree']))
pool=[r for r in dim['results'] if r['species']=='diflunisal--caffeine' and usable(r)]
assert pool
selected.append(min(pool,key=lambda x:x['result']['energy_Hartree']))
out=ROOT/'inputs/orca/week2';out.mkdir(parents=True,exist_ok=True)
jobs=[]
for row in selected:
    name=row['species'];job=f'w2-{name}-opt-freq';source=ROOT/row['run_dir']/'xtbopt.xyz'
    atoms=read(source);xyz=source.read_text().splitlines();n=int(xyz[0]);assert n==row['candidate']['n_atoms']
    inp=out/f'{job}.inp'
    inp.write_text('! r2SCAN-3c TightSCF TightOpt Opt Freq\n%maxcore 2000\n\n* xyz 0 1\n'+'\n'.join(xyz[2:2+n])+'\n*\n')
    config={'protocol_id':'dft-pilot-20260909-v1-week2','week':2,'job_id':job,'task_id':'DR01',
            'engine':'orca','executable':'.tools/orca-6.1.1/orca','input':str(inp.relative_to(ROOT)),
            'method':'r2SCAN-3c','purpose':'First cost and minimum verification; one selected starting geometry, not four-case ranking',
            'charge':0,'multiplicity':1,'n_atoms':n,'linear':False,'solvent':'vacuum',
            'requested':{'optimization':True,'frequencies':True},'timeout_seconds':43200,'rss_limit_gib':40,'omp_threads':1,
            'selected_candidate_id':row['candidate_id'],'parent_run_id':row['run_dir'],
            'selection_rule':'lowest valid GFN2-xTB geometry; fixed chemical state for monomers; retained composition/connectivity and no separation for dimer',
            'selection_uses_outcome_labels':False,'selection_label_scope':'geometry ranking within preselected chemical species only',
            'first_cost_pair_is_known_positive_by_design':True,
            'initial_structure_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'initial_geometry_audit':row['result']['geometry_audit']}
    cfg=ROOT/f'experiments/week2/{job}.json';cfg.write_text(json.dumps(config,indent=2)+'\n');jobs.append(str(cfg.relative_to(ROOT)))
(ROOT/'experiments/week2/drug-first-dft-batch.json').write_text(json.dumps({'batch_id':'week2-drug-first-dft','budget_seconds':21600,'jobs':jobs},indent=2)+'\n')
(ROOT/'reports/week2-first-dft-selection.json').write_text(json.dumps({'selected':selected,'uses_outcome_labels_for_pose_ranking':False,
                                                                   'first_cost_pair_is_known_positive_by_design':True,
                                                                   'all_diflunisal_caffeine_starts_attempted':len(expected)},indent=2)+'\n')
print(json.dumps({'prepared_DFT_jobs':len(jobs),'candidates':[r['candidate_id'] for r in selected]}))
