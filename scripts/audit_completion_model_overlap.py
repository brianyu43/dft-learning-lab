"""Audit public CSV exposure before selecting an independent cocrystal panel."""
import csv
from collections import defaultdict
from functools import lru_cache
import json

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from research_run import ROOT,sha,save,now

UNCHARGER=rdMolStandardize.Uncharger()
TAUTOMERS=rdMolStandardize.TautomerEnumerator()


@lru_cache(maxsize=None)
def identities(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError(f'Invalid SMILES: {smiles}')
    exact=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
    norm=UNCHARGER.uncharge(mol);Chem.RemoveStereochemistry(norm)
    neutral=Chem.MolToSmiles(norm,canonical=True,isomericSmiles=False)
    tautomer=canonical_tautomer(neutral)
    return {'exact':exact,'neutral_no_stereo':neutral,'tautomer_family':tautomer}


@lru_cache(maxsize=None)
def canonical_tautomer(neutral):
    return Chem.MolToSmiles(TAUTOMERS.Canonicalize(Chem.MolFromSmiles(neutral)),isomericSmiles=False)


def main():
    base=ROOT/'data/completion/deep-cocrystal'
    source=json.loads((base/'source.json').read_text())
    for row in source['files']:assert sha(ROOT/row['path'])==row['sha256']
    split_indices={};summary={};errors=[];api_indices={}
    for split in ['train','val','test_canonical','test_augmented']:
        path=base/f'data/{split}.csv'
        indices={k:defaultdict(set) for k in ['exact','neutral_no_stereo','tautomer_family']}
        molecules={k:set() for k in indices};count=0
        with path.open() as f:
            for number,row in enumerate(csv.DictReader(f),start=2):
                count+=1
                try:a,b=identities(row['SMILES1']),identities(row['SMILES2'])
                except ValueError as exc:
                    errors.append(dict(split=split,line=number,error=str(exc)));continue
                assert row['cocrystallization'] in ['0','1']
                for level in indices:
                    pair=tuple(sorted((a[level],b[level])))
                    indices[level][pair].add(int(row['cocrystallization']))
                    molecules[level].update(pair)
        split_indices[split]=indices;api_indices[split]=molecules
        summary[split]=dict(raw_rows=count,unique_pairs={k:len(v) for k,v in indices.items()},
            label_conflict_pairs={k:sum(len(labels)>1 for labels in v.values()) for k,v in indices.items()},
            source_sha256=sha(path))
    identity_source=json.loads((ROOT/'data/completion/identities/source.json').read_text())
    molecules={}
    for row in identity_source['rows']:
        assert row['status']=='downloaded' and sha(ROOT/row['path'])==row['sha256']
        p=row['properties'];mol=Chem.MolFromSmiles(p['SMILES']);assert mol is not None
        assert rdMolDescriptors.CalcMolFormula(mol)==p['MolecularFormula']
        assert Chem.MolToInchiKey(mol)==p['InChIKey']
        assert Chem.GetFormalCharge(mol)==0 and len(Chem.GetMolFrags(mol))==1
        molecules[row['key']]=dict(name=row['name'],CID=p['CID'],SMILES=p['SMILES'],formula=p['MolecularFormula'],
            total_atoms=Chem.AddHs(mol).GetNumAtoms(),heavy_atoms=mol.GetNumHeavyAtoms(),identities=identities(p['SMILES']))
    pairs=[]
    for api in ['PARA','ASPI','CARB']:
        for cf in sorted(set(molecules)-{'PARA','ASPI','CARB'}):
            matches={}
            for split,indices in split_indices.items():
                matches[split]={k:sorted(indices[k].get(tuple(sorted([molecules[api]['identities'][k],molecules[cf]['identities'][k]])),set())) for k in indices}
            exposed=any(matches[s][k] for s in ['train','val'] for k in matches[s])
            pairs.append(dict(api=api,coformer=cf,matches=matches,
                train_or_validation_match=exposed,
                held_out_from_downloaded_train_val=not exposed and not errors,
                seen_API_in_train=any(molecules[api]['identities'][k] in api_indices['train'][k] for k in api_indices['train']),
                final_model_independence_verified=False))
    out=ROOT/'reports/completion/model-overlap-2022.json'
    result=dict(updated_at=now(),repository=source['repo'],commit=source['commit'],split_summary=summary,
        invalid_rows=errors,molecules=molecules,pairs=pairs,
        matched_train_or_val_pair_count=sum(p['train_or_validation_match'] for p in pairs),
        no_match_in_downloaded_train_val_count=sum(p['held_out_from_downloaded_train_val'] for p in pairs),
        scope='Public CSV exposure audit only. Final pretrained-weight training membership is not established by CSV split names.',
        identity_rules={'exact':'canonical isomeric SMILES, entire structure',
            'neutral_no_stereo':'RDKit Uncharger then stereo removal; entire structure, no largest-fragment selection',
            'tautomer_family':'RDKit canonical tautomer of neutral/stereo-free structure; conservative leakage flag, not equality of calculated states'},
        policy='A train/validation match excludes the pair from a held-out model benchmark. API-seen and pair-seen are distinct. No-match alone does not prove final-model or unpublished-training independence.',
        RDKit_version=Chem.rdBase.rdkitVersion)
    save(out,result)
    print(json.dumps({k:result[k] for k in ['split_summary','matched_train_or_val_pair_count','no_match_in_downloaded_train_val_count','invalid_rows']},indent=2))


if __name__=='__main__':main()
