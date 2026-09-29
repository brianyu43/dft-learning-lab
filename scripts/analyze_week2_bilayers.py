"""Check matched frozen-layer references before reporting interaction energies."""
import json
import numpy as np
from week2_qe_results import ROOT, load_job

inventory=json.loads((ROOT/'experiments/week2/bilayer-comparisons.json').read_text())
rows=[]
for cfg in inventory['configurations']:
    if not all((ROOT/f'reports/{j}.json').is_file() for j in cfg['jobs'].values()):continue
    r={part:load_job(job) for part,job in cfg['jobs'].items()}
    pair=r['pair'];bottom=r['bottom'];top=r['top']
    for part,slice_ in [(bottom,slice(0,3)),(top,slice(3,6))]:
        np.testing.assert_allclose(pair['cell'],part['cell'],rtol=0,atol=1e-9)
        np.testing.assert_allclose(pair['positions'][slice_],part['positions'],rtol=0,atol=1e-8)
        assert pair['symbols'][slice_]==part['symbols']
        for key in ('ecutwfc_Ry','ecutrho_Ry','k_grid','conv_thr_Ry','method','charge','SOC'):
            assert pair['meta'][key]==part['meta'][key],key
        for source in pair['meta']['assets']:
            assert pair['meta']['source_sha256'][source]==part['meta']['source_sha256'][source]
    area=pair['area_angstrom2'];assert abs(area-cfg['area_angstrom2'])<1e-8
    eint=pair['energy_eV']-bottom['energy_eV']-top['energy_eV']
    # Rigidly moving the upper layer along z changes the energy by -sum(Fz).
    upper_force=float(pair['forces'][3:,2].sum())
    rows.append({**cfg,'energy_eV':{p:q['energy_eV'] for p,q in r.items()},
                 'interaction_energy_eV_cell':eint,'gamma_meV_angstrom2':1000*eint/area,
                 'upper_layer_force_z_eV_angstrom':upper_force,'area_angstrom2':area,
                 'matched_frozen_references_verified':True,'xml_sha256':{p:q['xml_sha256'] for p,q in r.items()}})
result={'status':'partial','completed_comparisons':len(rows),'expected_comparisons':len(inventory['configurations']),
        'rows':rows,'equilibrium_in_plane_lattice_verified':False,'phonon_stability_verified':False,
        'reference_is_exact':False,'heterobilayer_or_Se_accuracy_verified':False}
if len(rows)==len(inventory['configurations']):
    low=next(r for r in rows if r['setting']=='selected' and r['registry_id']==0 and r['Mo_plane_distance_angstrom']==6.3)
    high=next(r for r in rows if r['setting']=='strict')
    for r in [high]:
        # Vacuum varies but relative internal geometry and area must be identical.
        a=load_job(low['jobs']['pair']);b=load_job(r['jobs']['pair'])
        np.testing.assert_allclose(a['positions']-a['positions'][0],b['positions']-b['positions'][0],rtol=0,atol=1e-8)
    precision=abs(low['gamma_meV_angstrom2']-high['gamma_meV_angstrom2'])
    shifted=[r for r in rows if r['setting']=='selected' and r['registry_id'] in (1,2)]
    symmetry=abs(shifted[0]['gamma_meV_angstrom2']-shifted[1]['gamma_meV_angstrom2'])
    curve=sorted([r for r in rows if r['setting']=='selected' and r['registry_id']==0],key=lambda r:r['Mo_plane_distance_angstrom'])
    imin=min(range(len(curve)),key=lambda i:curve[i]['gamma_meV_angstrom2'])
    bracketed=0<imin<len(curve)-1
    result.update(status='completed_checks_passed' if precision<=0.1 and symmetry<=0.1 else 'numerical_precision_review_required',
                  representative_setting_difference_meV_angstrom2=precision,interaction_precision_passed=precision<=0.1,
                  equivalent_registry_difference_meV_angstrom2=symmetry,registry_symmetry_passed=symmetry<=0.1,
                  lowest_sampled_r0_distance_angstrom=curve[imin]['Mo_plane_distance_angstrom'],
                  r0_minimum_bracketed_by_samples=bracketed,equilibrium_layer_distance_verified=False,
                  precision_scope='one r0, D=6.3 A geometry; not a global bound for all registries or compositions',
                  distance_scope='three rigid-layer samples, no interpolation or relaxation; boundary minimum requires extension')
(ROOT/'reports/week2-bilayer-analysis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
