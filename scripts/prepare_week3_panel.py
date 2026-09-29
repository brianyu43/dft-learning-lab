"""Freeze a small control panel before DFT outcomes are known; no engine launches."""
import json
from pathlib import Path

from prepare_week3_drug import ROOT, immutable, prepare
from research_run import sha

monos = json.loads((ROOT/'reports/week2-xtb-monomers.json').read_text())['results']
dimers = json.loads((ROOT/'reports/week2-xtb-dimers.json').read_text())['results']


def usable(row):
    r = row['result']; a = r.get('geometry_audit', {})
    return r['status']=='passed' and a.get('heavy_connectivity_preserved') and not a.get('unassigned_hydrogens') and not a.get('separated_geometry')


jobs = []
selections = []
for species, state in [('adenine','adenine-t06'), ('theobromine','theobromine-t00'), ('xanthine','xanthine-t00')]:
    row = min((r for r in monos if r['species']==species and r['candidate']['state_id']==state
               and usable(r) and r['result']['geometry_audit']['chemical_state_preserved_by_geometry']),
              key=lambda r: r['result']['energy_Hartree'])
    cfg = prepare(row['candidate_id'], role='primary')
    jobs.append(cfg); selections.append(dict(candidate_id=row['candidate_id'], role='primary', species=species))

for species in ['diflunisal--diflunisal','caffeine--caffeine','diflunisal--adenine','adenine--adenine',
                'diflunisal--theobromine','theobromine--theobromine','diflunisal--xanthine','xanthine--xanthine']:
    row = min((r for r in dimers if r['species']==species and usable(r)), key=lambda r:r['result']['energy_Hartree'])
    cfg = prepare(row['candidate_id'], role='primary')
    jobs.append(cfg); selections.append(dict(candidate_id=row['candidate_id'],role='primary',species=species))

# Paper SI p15/Fig.S32: CAF N9 is our atom32; ADE N1 is our atom28.
# These mappings follow the molecular graph, not atom labels in an arbitrary XYZ.
caf = [r for r in dimers if r['species']=='diflunisal--caffeine' and usable(r)
       and r['candidate'].get('motif_detail',{}).get('donor_atom')==2
       and r['candidate'].get('motif_detail',{}).get('acceptor_atom')==32]
row = min(caf,key=lambda r:r['result']['energy_Hartree']); assert row['candidate_id']=='diflunisal--caffeine-p06'
cfg = prepare(row['candidate_id'], role='literature_contact_COOH_to_N9')
jobs.insert(0,cfg); selections.append(dict(candidate_id=row['candidate_id'],role='literature_contact_COOH_to_N9',species=row['species']))
ion = next(r for r in dimers if r['candidate_id']=='diflunisal--adenine-p36')
assert ion['candidate']['state_ids']==['diflunisal-t02','adenine-t06']
assert ion['candidate']['motif_detail']['to_atom']==28
cfg = prepare(ion['candidate_id'],direct=True,role='literature_N1_proton_transfer_probe')
jobs.append(cfg);selections.append(dict(candidate_id=ion['candidate_id'],direct_initial=True,role='literature_N1_proton_transfer_probe',species=ion['species']))

pilot = 'experiments/week3/drug/w3-diflunisal--caffeine-p14-opt-freq.json'
plan = dict(status='prepared_pending_first_cost_review', execution_enabled=False,
            protocol_id='dft-pilot-20260909-v1-week3', first_job=pilot, remaining_jobs=jobs,
            new_start_count_including_first=1+len(jobs), selections=selections,
            total_drug_budget_hours=80, remaining_after_week2_hours=75.88367425,
            first_cost_review='Use the first uninterrupted 12-process Opt/Freq actual cost before starting the remaining finite panel',
            parallelism=dict(ORCA_MPI=12, ORCA_OpenMP=1, QE_OpenMP=4, concurrent_engines=1),
            selection_uses_DFT_outcomes=False, first_cases_known_by_design=True,
            characterization_informed_candidates=['diflunisal--caffeine-p06','diflunisal--adenine-p36-direct-ion'],
            characterization_scope='SI spectroscopic hypotheses, not measured atomic crystal coordinates; source-informed retrospective comparison',
            candidate_cap='At most three DFT starts per species including week-two starts. Additional alternatives require budget/state review.',
            duplicate_policy='This bounded panel uses combined Opt/Freq for auditable full-job timing. Record any duplicate Hessian cost; do not count converged equivalent structures as independent minima.',
            comparison_policy='Report primary-only B/X and all-available-pose B/X separately with start counts; uneven pose coverage is not a fair exhaustive ranking.',
            missing_states='Primary panel alone does not cover every ADE/XAN tautomer or every binding motif; follow-up remains conditional.',
            frequency_check='Every reported new minimum requires complete positive internal modes and matching XYZ/Hessian geometry',
            input_hashes={p:sha(ROOT/p) for p in [pilot,*jobs]},
            journal_si_sha256=sha(ROOT/'data/week3/literature/ANIE-64-e202507835-s001.pdf'))
assert plan['new_start_count_including_first']==14 and len(jobs)==13
immutable(ROOT/'experiments/week3/drug/panel-plan.json',json.dumps(plan,indent=2)+'\n')
print(json.dumps({'new_starts':14,'remaining_prepared':13,'execution_enabled':False}))
