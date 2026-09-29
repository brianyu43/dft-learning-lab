"""Inspect the initial scan without calling finite reference settings exact."""
import json
from pathlib import Path
import numpy as np
from week2_qe_results import ROOT,load_job

batch=json.loads((ROOT/'experiments/week2/materials-initial-batch.json').read_text())
rows=[]
for config in batch['jobs']:
    cfg=json.loads((ROOT/config).read_text());p=ROOT/f"reports/{cfg['job_id']}.json"
    if not p.exists():continue
    r=load_job(cfg['job_id']);c=r['cell'];pos=r['positions']
    assert abs(c[0,0]-3.18)<1e-9
    assert abs(c[2,2]-np.ptp(pos[:,2])-cfg['vacuum_gap_angstrom'])<1e-8
    r.update(cutoff_Ry=cfg['ecutwfc_Ry'],k=cfg['k_grid'][0],vacuum_angstrom=cfg['vacuum_gap_angstrom'])
    rows.append(r)
result={'status':'partial','completed':len(rows),'expected':7,'axes':{},'thresholds':{'energy_meV_atom':1.0,'force_eV_angstrom':0.005},
        'reference_is_exact':False,'monolayer_geometry':'week1 internally relaxed at fixed a=3.18 A'}
reference=next((r for r in rows if (r['cutoff_Ry'],r['k'],r['vacuum_angstrom'])==(120,12,30)),None)
if reference:
    ref_relative=reference['positions']-reference['positions'][0]
    for r in rows:
        np.testing.assert_allclose(r['positions']-r['positions'][0],ref_relative,atol=1e-8,rtol=0)
        r['energy_delta_reference_meV_atom']=(r['energy_eV']-reference['energy_eV'])*1000/3
        r['max_force_delta_reference_eV_angstrom']=float(np.linalg.norm(r['forces']-reference['forces'],axis=1).max())
        r['within_reference_thresholds']=abs(r['energy_delta_reference_meV_atom'])<=1 and r['max_force_delta_reference_eV_angstrom']<=0.005
    for axis,field,condition in [('cutoff','cutoff_Ry',lambda r:r['k']==12 and r['vacuum_angstrom']==30),
                                 ('k_grid','k',lambda r:r['cutoff_Ry']==120 and r['vacuum_angstrom']==30),
                                 ('vacuum','vacuum_angstrom',lambda r:r['cutoff_Ry']==120 and r['k']==12)]:
        selected=sorted([r for r in rows if condition(r)],key=lambda r:r[field])
        if len(selected)==3:
            low,high=selected[-2:]
            de=abs(low['energy_eV']-high['energy_eV'])*1000/3
            df=float(np.linalg.norm(low['forces']-high['forces'],axis=1).max())
            passing=[r[field] for r in selected if r['within_reference_thresholds']]
            result['axes'][axis]={'tested_values':[r[field] for r in selected],'last_interval_energy_meV_atom':de,
                                  'last_interval_force_eV_angstrom':df,'last_interval_passed':de<=1 and df<=0.005,
                                  'smallest_tested_value_within_reference_tolerance':min(passing) if passing else None}
    if len(rows)==7:
        passed=all(a['last_interval_passed'] for a in result['axes'].values())
        result['status']='initial_axes_passed_combined_setting_pending' if passed else 'extension_or_precision_review_required'
        if passed:result['candidate_combined_setting']={axis:a['smallest_tested_value_within_reference_tolerance'] for axis,a in result['axes'].items()}
result['rows']=[{k:v for k,v in r.items() if k not in ('meta','positions','cell','forces')} for r in rows]
(ROOT/'reports/week2-monolayer-convergence.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
