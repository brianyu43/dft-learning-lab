"""Cluster optimized starting structures, retaining failures and changed-state records."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from ase.io import read
from run_week2_xtb_screen import geometry_audit

ROOT=Path(__file__).resolve().parents[1]
state_smiles={s['state_id']:s['smiles'] for m in json.loads((ROOT/'reports/week1-molecule-identities.json').read_text())['molecules'] for s in m['tautomer_inventory']}

def make_mol(candidate,positions):
    states=candidate.get('state_ids',[candidate.get('state_id')])
    ms=[Chem.AddHs(Chem.MolFromSmiles(state_smiles[s])) for s in states]
    mol=ms[0]
    for m in ms[1:]:mol=Chem.CombineMols(mol,m)
    assert mol.GetNumAtoms()==len(positions)
    conf=Chem.Conformer(len(positions))
    for i,p in enumerate(positions):conf.SetAtomPosition(i,p)
    mol.AddConformer(conf)
    return mol

def contacts(candidate,positions):
    if candidate['kind']!='dimer':return []
    z=candidate['graph']['atomic_numbers'];n=candidate['fragment_atom_counts'][0];heavy=[i for i,x in enumerate(z) if x>1]
    result=set()
    for h,element in enumerate(z):
        if element!=1:continue
        d=min(heavy,key=lambda j:np.linalg.norm(positions[h]-positions[j]))
        if z[d] not in (7,8):continue
        for a in heavy:
            if (a<n)==(d<n) or z[a] not in (7,8):continue
            dh=positions[d]-positions[h];ah=positions[a]-positions[h]
            angle=np.degrees(np.arccos(np.clip(np.dot(dh,ah)/np.linalg.norm(dh)/np.linalg.norm(ah),-1,1)))
            if np.linalg.norm(ah)<2.6 and angle>120:result.add((d,a))
    return sorted(result)

summaries={};all_records=[];transitions=[]
for phase in ('monomers','dimers'):
    batch=json.loads((ROOT/f'reports/week2-xtb-{phase}.json').read_text())
    assert batch['status']=='completed'
    species=sorted({r['species'] for r in batch['results']})
    for name in species:
        source=[r for r in batch['results'] if r['species']==name]
        valid=[];excluded=[]
        for r in source:
            a=r['result'].get('geometry_audit') or {}
            if a.get('status')=='geometry_checked':
                c=r['candidate'];initial=geometry_audit(c,read(ROOT/r['run_dir']/'input.xyz'))
                initial_owners={int(h):owner for h,owner in c['graph']['hydrogen_owners'].items()}
                final_owners=initial_owners.copy()
                for change in initial['hydrogen_owner_changes']:initial_owners[change['hydrogen']]=change['to']
                for change in a['hydrogen_owner_changes']:final_owners[change['hydrogen']]=change['to']
                net=[{'hydrogen':h,'input_owner':initial_owners[h],'output_owner':final_owners[h]}
                     for h in initial_owners if initial_owners[h]!=final_owners[h]]
                transitions.append({'candidate_id':r['candidate_id'],'species':name,'engine_status':r['result']['status'],
                                    'initial_motif':c.get('motif'),'net_hydrogen_owner_changes_from_input':net,
                                    'input_changes_relative_to_state_graph':initial['hydrogen_owner_changes'],
                                    'output_changes_relative_to_state_graph':a['hydrogen_owner_changes'],
                                    'output_hydrogens_match_state_graph':not a['hydrogen_owner_changes'] and not a['unassigned_hydrogens']})
            if r['result']['status']!='passed' or not a.get('heavy_connectivity_preserved') or a.get('unassigned_hydrogens') or a.get('separated_geometry'):
                excluded.append({'candidate_id':r['candidate_id'],'engine_status':r['result']['status'],'audit':a});continue
            valid.append(r)
        clusters=[]
        for row in sorted(valid,key=lambda r:r['result']['energy_Hartree']):
            c=row['candidate'];a=row['result']['geometry_audit'];positions=read(ROOT/row['run_dir']/'xtbopt.xyz').positions
            mol=make_mol(c,positions);heavy=Chem.RemoveHs(mol)
            change_signature=tuple((x['from'],x['to']) for x in a['hydrogen_owner_changes'])
            state_signature=tuple(c.get('state_ids',[c.get('state_id')]))
            contact=contacts(c,positions)
            key=(state_signature,change_signature,tuple(contact))
            match=None;rmsd=None
            for i,cluster in enumerate(clusters):
                if cluster['key']!=key:continue
                rms=rdMolAlign.GetBestRMS(Chem.Mol(heavy),cluster['mol'],maxMatches=10000)
                if rms<0.75:match=i;rmsd=rms;break
            if match is None:
                match=len(clusters);clusters.append({'key':key,'mol':heavy,'representative':row['candidate_id'],
                                                     'energy_Hartree':row['result']['energy_Hartree'],'members':[],
                                                     'contacts':contact,'hydrogen_owner_changes':a['hydrogen_owner_changes'],
                                                     'initial_motif':c.get('motif')})
            clusters[match]['members'].append({'candidate_id':row['candidate_id'],'heavy_RMSD_angstrom':rmsd})
            all_records.append({'candidate_id':row['candidate_id'],'species':name,'cluster_id':match,
                                'energy_Hartree':row['result']['energy_Hartree'],'changed_hydrogen_owners':bool(a['hydrogen_owner_changes'])})
        output=[]
        for i,c in enumerate(clusters):
            output.append({k:v for k,v in c.items() if k not in ('key','mol')}|{'cluster_id':i,'state_ids':list(c['key'][0]),
                          'relative_energy_kcal_mol':(c['energy_Hartree']-clusters[0]['energy_Hartree'])*627.509474 if clusters else None})
        shortlist=[]
        if output:
            shortlist=[output[0]['representative']];seen={output[0]['initial_motif']}
            for c in output[1:]:
                if c['initial_motif'] not in seen and len(shortlist)<3:shortlist.append(c['representative']);seen.add(c['initial_motif'])
            for c in output:
                if len(shortlist)<3 and c['representative'] not in shortlist:shortlist.append(c['representative'])
        summaries[name]={'started':len(source),'engine_passed':sum(r['result']['status']=='passed' for r in source),
                         'geometry_usable':len(valid),'excluded_but_retained':excluded,'cluster_count':len(clusters),
                         'optimized_geometry_clusters':output,'provisional_max3_shortlist':shortlist,
                         'shortlist_requires_state_review':True}
result={'status':'analyzed','method':'GFN2-xTB optimized geometry clusters, not Hessian-verified minima',
        'clustering':'Symmetry-aware RDKit heavy-atom RMSD <0.75 A within explicit starting states, changed-proton signature, and geometric X-H...Y contacts; index signatures may conservatively over-split symmetry-related contacts',
        'selection_uses_outcome_labels':False,'scope':'Candidate preparation and cost calibration only; no cocrystal classifier performance',
        'hydrogen_assignment_scope':'Nearest-heavy-atom assignment only. Raw geometry_audit changes compare with the original neutral state graph, which differs from the intentionally proton-shifted input for 12 candidates. net_hydrogen_owner_changes_from_input compares actual input and output coordinates; it is not an electronic charge or full trajectory assignment.',
        'hydrogen_endpoint_transitions':transitions,'species':summaries,'records':all_records}
pt=[r for r in transitions if r['initial_motif']=='proton_transfer_start']
pt_converged=[r for r in pt if r['engine_status']=='passed']
result['proton_transfer_state_review']={
    'initial_geometries':len(pt),'converged':len(pt_converged),
    'converged_outputs_matching_neutral_template_H_owners':sum(r['output_hydrogens_match_state_graph'] for r in pt_converged),
    'DFT_ionic_basin_tested':False,
    'next_stage_action':'Retain original proton-transfer XYZ inputs for a direct DFT state check within the per-species cap. An xTB-relaxed shortlist alone may omit a DFT ionic basin.',
    'scope':'gas-phase GFN2 geometry endpoints; not an electronic charge assignment or a test of the reported crystal protonation state'}
(ROOT/'reports/week2-xtb-analysis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({name:{k:s[k] for k in ('started','engine_passed','geometry_usable','cluster_count')} for name,s in summaries.items()},indent=2))
