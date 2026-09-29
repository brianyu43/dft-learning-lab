"""Deterministic, label-independent monomer and dimer starting geometry construction.

Geometry heuristics label proton ownership and clashes, not formal electronic charge.
The inventory is distinct from a DFT shortlist and from confirmed local minima.
"""
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'structures/week2/drug_initial'
OUT.mkdir(parents=True,exist_ok=True)
TABLE=Chem.GetPeriodicTable()

def coords(mol):
    return np.array(mol.GetConformer().GetPositions())

def xyz(mol,positions):
    return str(len(positions))+'\nGenerated starting geometry; not a quantum result\n'+'\n'.join(
        f'{a.GetSymbol():2s} {p[0]: .12f} {p[1]: .12f} {p[2]: .12f}' for a,p in zip(mol.GetAtoms(),positions))+'\n'

def graph_info(mol):
    return {'atomic_numbers':[a.GetAtomicNum() for a in mol.GetAtoms()],
            'bonds':[[b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondTypeAsDouble()] for b in mol.GetBonds()],
            'hydrogen_owners':{str(a.GetIdx()):a.GetNeighbors()[0].GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==1},
            'formula':rdMolDescriptors.CalcMolFormula(mol)}

def embed(smiles,seed):
    m=Chem.AddHs(Chem.MolFromSmiles(smiles))
    p=AllChem.ETKDGv3();p.randomSeed=seed;p.numThreads=1
    if AllChem.EmbedMolecule(m,p):return None
    if not AllChem.MMFFHasAllMoleculeParams(m):return None
    if AllChem.MMFFOptimizeMolecule(m,maxIters=1500):return None
    return m

def align(v,target):
    v=v/np.linalg.norm(v);target=np.array(target,dtype=float);target/=np.linalg.norm(target)
    axis=np.cross(v,target);s=np.linalg.norm(axis);c=np.clip(np.dot(v,target),-1,1)
    if s<1e-12:
        if c>0:return np.eye(3)
        axis=np.cross(v,[1,0,0] if abs(v[0])<0.8 else [0,1,0]);axis/=np.linalg.norm(axis)
        return Rotation.from_rotvec(np.pi*axis).as_matrix()
    return Rotation.from_rotvec(np.arctan2(s,c)*axis/s).as_matrix()

def donors(m):
    return [(a.GetIdx(),h.GetIdx()) for a in m.GetAtoms() if a.GetAtomicNum() in (7,8)
            for h in a.GetNeighbors() if h.GetAtomicNum()==1]

def acceptors(m):
    result=[]
    for a in m.GetAtoms():
        if a.GetAtomicNum()==8 and not any(n.GetAtomicNum()==1 for n in a.GetNeighbors()):result.append(a.GetIdx())
        if a.GetAtomicNum()==7 and a.GetIsAromatic() and not any(n.GetAtomicNum()==1 for n in a.GetNeighbors()):
            # Aromatic nitrogens bonded to carbonyl carbons are amide-like, not ring imine acceptors.
            amide=any(any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(n).GetAtomicNum()==8 for b in n.GetBonds()) for n in a.GetNeighbors() if n.GetAtomicNum()==6)
            if not amide:result.append(a.GetIdx())
    return result

def no_clash(a,b,pa,pb,transferred=None):
    for i,aa in enumerate(a.GetAtoms()):
        for j,bb in enumerate(b.GetAtoms()):
            if transferred==(i,j):continue
            if aa.GetAtomicNum()>1 and bb.GetAtomicNum()>1:
                threshold=0.72*(TABLE.GetRvdw(aa.GetAtomicNum())+TABLE.GetRvdw(bb.GetAtomicNum()))
            else:threshold=0.95*(TABLE.GetRcovalent(aa.GetAtomicNum())+TABLE.GetRcovalent(bb.GetAtomicNum()))
            if np.linalg.norm(pa[i]-pb[j])<threshold:return False
    return True

def hbond_positions(a,b,donor,acceptor,angle,transfer=False):
    heavy,h=donor;pa=coords(a);pb=coords(b)
    pa=(pa-pa[heavy])@align(pa[h]-pa[heavy],[1,0,0]).T
    neigh=[n.GetIdx() for n in b.GetAtomWithIdx(acceptor).GetNeighbors() if n.GetAtomicNum()>1]
    outward=pb[acceptor]-np.mean(pb[neigh],axis=0)
    if np.linalg.norm(outward)<1e-8:return None
    pb=(pb-pb[acceptor])@align(outward,[-1,0,0]).T
    pb=pb@Rotation.from_rotvec(np.radians(angle)*np.array([1,0,0])).as_matrix().T
    pb+=np.array([2.85,0,0])
    if transfer:pa[h]=pb[acceptor]+np.array([-1.04,0,0])
    if not no_clash(a,b,pa,pb,(h,acceptor) if transfer else None):return None
    return pa,pb

def planar(m):
    p=coords(m);idx=[a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()]
    if len(idx)<3:idx=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]
    center=p[idx].mean(axis=0);_,_,vt=np.linalg.svd(p[idx]-center)
    return (p-center)@align(vt[-1],[0,0,1]).T

def stack_positions(a,b,angle,offset,height):
    pa=planar(a);pb=planar(b)@Rotation.from_rotvec(np.radians(angle)*np.array([0,0,1])).as_matrix().T
    pb+=np.array([*offset,height])
    if not no_clash(a,b,pa,pb):return None
    return pa,pb

def save_candidate(identifier,m,positions,extra):
    path=OUT/f'{identifier}.xyz';path.write_text(xyz(m,positions))
    return {'candidate_id':identifier,'path':str(path.relative_to(ROOT)),
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'charge':0,'multiplicity':1,
            'n_atoms':m.GetNumAtoms(),'graph':graph_info(m),**extra}

def main():
    identities=json.loads((ROOT/'reports/week1-molecule-identities.json').read_text())['molecules']
    monomers=[];by_state={};failures=[]
    for record in identities:
        for state in record['tautomer_inventory']:
            for seed in (17,43,89):
                m=embed(state['smiles'],seed)
                if m is None:
                    failures.append({'state_id':state['state_id'],'seed':seed,'stage':'RDKit MMFF embedding'});continue
                by_state[(state['state_id'],seed)]=m
                monomers.append(save_candidate(f"{state['state_id']}-s{seed}",m,coords(m),
                                {'kind':'monomer','species':record['name'],'state_id':state['state_id'],
                                 'state_smiles':state['smiles'],'seed':seed,'state_role':'pubchem_representation' if state['is_pubchem_representation'] else 'alternative_tautomer'}))
    # Two neutral N-H placements retained for purines. Alternative keto/enol/amino-imino states remain in monomer screening.
    use={'diflunisal':['diflunisal-t02'],'adenine':['adenine-t06','adenine-t07'],
         'caffeine':['caffeine-t00'],'theobromine':['theobromine-t00'],
         'xanthine':['xanthine-t00','xanthine-t01']}
    pairs=[('diflunisal',b) for b in ('adenine','caffeine','theobromine','xanthine')]+[(a,a) for a in use]
    dimers=[];pair_summary=[]
    for an,bn in pairs:
        species=f'{an}--{bn}';pools={'hydrogen_bond':[],'stacked_contact':[],'proton_transfer_start':[]}
        combinations=list(itertools.product(use[an],use[bn]))
        if an==bn:combinations=[p for p in combinations if p[0]<=p[1]]
        for sa,sb in combinations:
            a,b=by_state[(sa,17)],by_state[(sb,43)]
            for d,acceptor,angle in itertools.product(donors(a),acceptors(b),(0,60,120,180,240,300)):
                pos=hbond_positions(a,b,d,acceptor,angle)
                if pos:pools['hydrogen_bond'].append((a,b,*pos,sa,sb,{'donor_atom':d[0],'hydrogen_atom':d[1],'acceptor_atom':len(coords(a))+acceptor,'torsion_degrees':angle}))
                if an=='diflunisal' and bn=='adenine':
                    # Acid OH donor only; keep total composition/charge while moving one existing H.
                    heavy=a.GetAtomWithIdx(d[0]);is_acid=any(n.GetAtomicNum()==6 and any(z.GetOtherAtom(n).GetAtomicNum()==8 and z.GetBondTypeAsDouble()==2 for z in n.GetBonds()) for n in heavy.GetNeighbors())
                    if is_acid:
                        pos=hbond_positions(a,b,d,acceptor,angle,True)
                        if pos:pools['proton_transfer_start'].append((a,b,*pos,sa,sb,{'hydrogen_atom':d[1],'from_atom':d[0],'to_atom':a.GetNumAtoms()+acceptor,'torsion_degrees':angle}))
            if an!=bn:
                for d,acceptor,angle in itertools.product(donors(b),acceptors(a),(0,60,120,180,240,300)):
                    pos=hbond_positions(b,a,d,acceptor,angle)
                    if pos:
                        pb,pa=pos
                        pools['hydrogen_bond'].append((a,b,pa,pb,sa,sb,{'donor_atom':a.GetNumAtoms()+d[0],'hydrogen_atom':a.GetNumAtoms()+d[1],'acceptor_atom':acceptor,'torsion_degrees':angle,'direction':'B-donor-to-A'}))
            for angle,offset,height in itertools.product((0,60,120,180,240,300),((0,0),(1,0),(-1,0),(0,1),(0,-1)),(3.4,3.7)):
                pos=stack_positions(a,b,angle,offset,height)
                if pos:pools['stacked_contact'].append((a,b,*pos,sa,sb,{'rotation_degrees':angle,'lateral_offset_angstrom':list(offset),'initial_plane_separation_angstrom':height}))
        selected=[]
        quotas={'hydrogen_bond':24,'stacked_contact':12 if pools['proton_transfer_start'] else 24,'proton_transfer_start':12}
        for motif,quota in quotas.items():
            pool=pools[motif]
            # Evenly cover state/order/angle combinations deterministically rather than taking the first state only.
            if pool:
                indices=np.unique(np.linspace(0,len(pool)-1,min(quota,len(pool)),dtype=int))
                selected.extend((motif,pool[i]) for i in indices)
        # No conventional donors in caffeine homodimer: use the remaining contact coverage.
        if len(selected)<48:
            used={(x[0],str(x[1][-1]),x[1][4],x[1][5]) for x in selected}
            for entry in pools['stacked_contact']:
                key=('stacked_contact',str(entry[-1]),entry[4],entry[5])
                if key not in used:selected.append(('stacked_contact',entry));used.add(key)
                if len(selected)==48:break
        for index,(motif,entry) in enumerate(selected):
            a,b,pa,pb,sa,sb,detail=entry;m=Chem.CombineMols(a,b)
            dimers.append(save_candidate(f'{species}-p{index:02d}',m,np.vstack([pa,pb]),
                          {'kind':'dimer','species':species,'state_ids':[sa,sb],'fragment_atom_counts':[a.GetNumAtoms(),b.GetNumAtoms()],
                           'motif':motif,'motif_detail':detail,'seeds':[17,43],
                           'state_role':'proton-transfer geometry candidate' if motif=='proton_transfer_start' else 'neutral starting pair'}))
        pair_summary.append({'species':species,'candidate_count':len(selected),
                             'motifs':{k:sum(x[0]==k for x in selected) for k in pools}})
    result={'protocol_id':'dft-pilot-20260909-v1-week2','algorithm':'ETKDGv3/MMFF; directed X-H...acceptor and stacked contacts',
            'uses_outcome_labels':False,'DFT_selected':False,'embedding_failures':failures,
            'monomers':monomers,'dimers':dimers,'pairs':pair_summary,
            'chemical_review':'Canonical acid/phenol diflunisal and diketo caffeine used for first timing. Purine 7H/9H candidates retained; proton-transfer contacts are candidates, not measured ionic charge.'}
    (ROOT/'experiments/week2/chemistry-candidates.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'monomers':len(monomers),'dimers':len(dimers),'embedding_failures':len(failures),'pairs':pair_summary}))

if __name__=='__main__':main()
