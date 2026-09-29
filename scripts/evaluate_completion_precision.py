"""Recheck raw coordinates/Hessians/gradients before a numerical protocol decision."""
import json
from pathlib import Path

import numpy as np
from ase.io import read

from completion_batch import verify_report
from research_run import ROOT,sha,save,now
from run_week2_frequency_resume import inspect_hessian
from run_week2_xtb_screen import geometry_audit
from week3_recovery import gradient_in_frame

HARTREE_KCAL=627.509474


def main():
    parent_path=ROOT/'experiments/completion/geometry-precision-v3.json'
    check_path=ROOT/'experiments/completion/geometry-precision-v3-checks.json'
    parent=json.loads(parent_path.read_text());checks=json.loads(check_path.read_text())
    assert checks['parent_plan_sha256']==sha(parent_path)
    policy=parent['acceptance'];assert checks['acceptance']==policy
    measured={};rows=[]
    for proof,check_config in zip(checks['parents'],checks['jobs'],strict=True):
        assert sha(ROOT/proof['report'])==proof['report_sha256']
        dense=verify_report(proof['config']);strict=verify_report(check_config)
        assert dense['candidate_id']==strict['candidate_id']
        assert dense['method']==strict['method']=='r2SCAN-3c'
        assert dense['solvent']==strict['solvent']=='vacuum'
        for key in ['charge','multiplicity','n_atoms']:assert dense[key]==strict[key]
        atoms=read(ROOT/proof['geometry'])
        assert sha(ROOT/proof['geometry'])==proof['geometry_sha256']==dense['result']['optimized_xyz_sha256']
        assert sha(Path(strict['run_dir'])/'initial.xyz')==proof['geometry_sha256']
        path=Path(dense['run_dir'])/(Path(dense['input']).stem+'.hess')
        assert sha(path)==dense['result']['hessian_sha256']
        frequencies=inspect_hessian(path,atoms,dense['result']['frequencies_cm1'])
        assert np.all(frequencies>0)
        candidate=json.loads((Path(dense['run_dir'])/'candidate.json').read_text())
        audit=geometry_audit(candidate,atoms)
        assert audit['chemical_state_preserved_by_geometry']
        energies={};gradients={}
        for label,meta in [('dense',dense),('strict',strict)]:
            energy,gradient,_,digest=gradient_in_frame(meta,atoms)
            assert digest==meta['result']['engrad_sha256']
            energies[label]=energy;gradients[label]=gradient
        change=float(np.max(np.abs(gradients['strict']-gradients['dense'])))
        species='dif2' if '--' in dense['candidate_id'] else 'dif'
        assert species not in measured
        measured[species]=energies
        rows.append(dict(species=species,candidate_id=dense['candidate_id'],energies_Hartree=energies,
            max_gradient_component_change_Hartree_Bohr=change,
            gradient_within_target=change<=policy['max_gradient_component_change_Hartree_Bohr'],
            minimum_frequency_cm1=float(min(frequencies)),n_frequencies=len(frequencies),
            geometry_sha256=proof['geometry_sha256'],minimum_verified_from_raw_hessian=True,
            source_reports=[proof['report'],f"reports/{strict['job_id']}.json"]))
    assert set(measured)=={'dif','dif2'}
    binding={label:(measured['dif2'][label]-2*measured['dif'][label])*HARTREE_KCAL for label in ['dense','strict']}
    delta=binding['strict']-binding['dense']
    energy_pass=abs(delta)<=policy['max_self_binding_change_kcal_mol']
    passed=energy_pass and all(r['gradient_within_target'] for r in rows)
    record=dict(status='representative_precision_passed' if passed else 'stopped_requires_review',
        updated_at=now(),parent_plan_sha256=sha(parent_path),check_plan_sha256=sha(check_path),
        rows=rows,self_binding_kcal_mol=binding,self_binding_change_kcal_mol=delta,
        self_binding_within_target=energy_pass,acceptance=policy,
        full_panel_numerical_precision_verified=False,full_panel_execution_enabled=False,
        next_action='Freeze common full-panel settings and finish all comparisons' if passed else 'Review the remaining numerical or geometry sensitivity without relaxing thresholds',
        scope='Only the saved DIF and DIF-DIF starting basin; no global-minimum or cocrystal-formation conclusion')
    save(ROOT/'reports/completion/geometry-precision-v3.json',record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
