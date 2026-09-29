"""Transcribe visually checked SI tables as condition-specific records, not universal labels."""
import csv
import io
import json
from pathlib import Path

from prepare_week3_drug import immutable
from research_run import ROOT,sha

COFORMERS=['BITZ','CEBG','DUPK','ESAL','NICC','NICM','OXAL','SUCA','TELZ','PYRD']
METHODS=[('contact_preparation',None),('slurry','n-heptane'),('slurry','pyridine'),
         ('slurry','diethyl ether'),('LAG','n-heptane'),('LAG','diethyl ether'),('cosublimation',None)]
LABELS={'x':'physical_mixture','yes':'multicomponent_form_observed','na':'not_attempted','inc':'inconclusive'}
# Each literal row follows the printed seven columns in Tables S1, S2, S3.
RAW={
 'PARA':[
  'x x na x x x x','x x na x x x x','x x na x x x x','x x na x x x x',
  'na x na x x x x','x x na x x x x','na yes na yes yes yes yes',
  'x x na x x x x','na x na x x x x','na na yes na na na na'],
 'ASPI':[
  'x x na x x x x','x x na x x x x','x x na x x x x','na x na x x x x',
  'na x na x x x x','x x na x x x x','na x na x x x x',
  'na x na x x x x','x x na x x x x','na na x na na na na'],
 'CARB':[
  'x x na x x x x','yes yes na yes yes yes x','x inc na x x x x','na yes na yes yes yes x',
  'na x na x x x x','yes yes na yes yes yes x','na yes na yes yes yes x',
  'na yes na yes yes yes x','x yes na yes yes yes x','na na x na na na na'],
}


def main():
    source=ROOT/'data/completion/literature/cg2c00433_si_001.pdf'
    assert sha(source)=='9b6358abe407a5a17517dc2f7d0096d2dfffe7dcee41fabed283aa162d2e6eaa'
    rows=[]
    for api,table_rows in RAW.items():
        assert len(table_rows)==10
        for coformer,line in zip(COFORMERS,table_rows):
            labels=line.split();assert len(labels)==7
            for column,((method,solvent),label) in enumerate(zip(METHODS,labels)):
                rows.append(dict(record_id=f'cg2c00433-{api}-{coformer}-c{column+1}',api=api,coformer=coformer,
                    method=method,solvent=solvent,outcome=LABELS[label],raw_code={'na':'n. a.','inc':'Inconclusive'}.get(label,label),
                    input_molar_ratio=None if method=='contact_preparation' else '1:1',
                    product_stoichiometry=None,independent_trial_count=None,
                    source_table={'PARA':'S1','ASPI':'S2','CARB':'S3'}[api],
                    source_pdf_pages={'PARA':[5,6],'ASPI':[12],'CARB':[17]}[api],
                    endpoint='solvate' if coformer=='PYRD' else 'cocrystal_screen',
                    scope='Reported table cell; repeat counts and precise per-cell conditions are not supplied by this table'))
    assert len(rows)==210 and len({r['record_id'] for r in rows})==210
    # Table S4 retains the ambiguity of CARB/DUPK, even when later conditions
    # returned a physical mixture. Its rows may overlap S3; do not count as replicates.
    extra=[]
    groups={
      'grinding':['dry','dichloromethane','dichloroethane','diethyl ether','diisopropyl ether','acetone','methanol','n-butanol','ethyl acetate','n-heptane'],
      'slurry':['diethyl ether','diisopropyl ether','n-butanol','n-heptane'],
      'solvent_evaporation':['acetone','methanol','ethanol','n-butanol','ethyl acetate','acetonitrile'],
      'cooling_crystallization':['acetone','methanol','ethanol','n-butanol','ethyl acetate','acetonitrile'],
    }
    for method,solvents in groups.items():
        for solvent in solvents:
            for ratio in ['1:1','2:1','1:2']:
                attempted=method=='grinding' or ratio=='1:1'
                extra_peaks=ratio=='1:1' and ((method=='grinding' and solvent=='dichloromethane') or (method=='slurry' and solvent=='n-heptane'))
                extra.append(dict(api='CARB',coformer='DUPK',method=method,solvent=None if solvent=='dry' else solvent,
                    input_molar_ratio=ratio,outcome='physical_mixture' if attempted else 'not_attempted',
                    raw_code='x, x^a' if extra_peaks else 'x' if attempted else 'n. a.',
                    qualifier='Two additional peaks, otherwise physical mixture' if extra_peaks else None,
                    source_table='S4',source_pdf_pages=[23],independent_trial_count=None))
    assert len(extra)==78
    data=dict(source_doi='10.1021/acs.cgd.2c00433',supplement_doi='10.1021/acs.cgd.2c00433.s001',
        source=str(source.relative_to(ROOT)),source_sha256=sha(source),
        review='SI pages 5,6,12,17,23 extracted and rendered; column labels and table footnotes visually checked.',
        source_method_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9337750/#sec2',
        method_conditions={
            'slurry':dict(input_ratio='1:1',mixture_mass_mg=200,solvent_volume_uL=[300,600],temperature_C=[10,30],duration='up to 4 weeks',readout='PXRD',per_pair_exact_duration_reported=False),
            'LAG':dict(input_ratio='1:1',mixture_mass_mg=100,solvent_amount='1-2 drops',frequency_Hz=30,duration_minutes=15,readout=['PXRD','IR'],per_pair_temperature_reported=False),
            'cosublimation':dict(input_ratio='1:1',pressure_mbar=200,temperature='elevated, pair-specific values not extracted'),
            'contact_preparation':dict(input_ratio=None,temperature='sequential melting, pair-specific thermal history not fully tabulated'),
        },
        records=rows,additional_ambiguous_pair_records=extra,
        pair_count=30,table_cell_count=len(rows),additional_table_cell_count=len(extra),
        tables_are_not_independent_replicate_counts=True,negative_label_scope='Physical mixture under the reported attempt, not impossibility under every condition',
        status='source_pool_curated_not_selected_for_independent_calculation',
        outcomes_known_before_future_score_calculation=True)
    path=ROOT/'data/completion/cocrystal-screen-2022.json'
    immutable(path,json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    fields=['record_id','api','coformer','method','solvent','raw_code','outcome','input_molar_ratio','source_table','endpoint']
    buf=io.StringIO();writer=csv.DictWriter(buf,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    immutable(path.with_suffix('.csv'),buf.getvalue())
    summary={api:{out:sum(r['api']==api and r['outcome']==out for r in rows) for out in LABELS.values()} for api in RAW}
    print(json.dumps(dict(path=str(path),sha256=sha(path),summary=summary),indent=2))


if __name__=='__main__':main()
