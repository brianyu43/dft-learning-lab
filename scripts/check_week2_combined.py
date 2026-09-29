"""Prepare or assess the combined monolayer setting, without starting an engine."""
import argparse
import json
import numpy as np
from prepare_week2_materials import ROOT, OUT, mono, write_job
from week2_qe_results import load_job

p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');args=p.parse_args()
scan=json.loads((ROOT/'reports/week2-monolayer-convergence.json').read_text())
assert scan['status']=='initial_axes_passed_combined_setting_pending'
c=scan['candidate_combined_setting'];cut,k,v=c['cutoff'],c['k_grid'],c['vacuum']
job=f'w2-mono-combined-c{cut}-k{k}-v{v}'
if args.prepare:
    cfg=write_job(job,mono(v),cut,k,v,scan_axis='combined_setting',configuration_id='MoS2-monolayer')
    (OUT/'materials-combined-batch.json').write_text(json.dumps({'batch_id':'week2-mono-combined','budget_seconds':3600,'jobs':[cfg]},indent=2)+'\n')
    print(job)
else:
    r=load_job(job);ref=load_job('w2-mono-c120-k12-v30')
    np.testing.assert_allclose(r['positions']-r['positions'][0],ref['positions']-ref['positions'][0],atol=1e-8,rtol=0)
    np.testing.assert_allclose(r['cell'][:2],ref['cell'][:2],atol=1e-9,rtol=0)
    de=(r['energy_eV']-ref['energy_eV'])*1000/3
    df=float(np.linalg.norm(r['forces']-ref['forces'],axis=1).max())
    result={'job_id':job,'selected_setting':c,'reference_setting':{'cutoff':120,'k_grid':12,'vacuum':30},
            'energy_delta_meV_atom':de,'max_force_delta_eV_angstrom':df,
            'status':'passed' if abs(de)<=1 and df<=0.005 else 'combined_setting_requires_refinement',
            'reference_is_exact':False,'scope':'fixed-geometry MoS2 monolayer only; interaction convergence tested separately'}
    (ROOT/'reports/week2-monolayer-combined.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
