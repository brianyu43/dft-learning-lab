"""Verify that all 505 unique xTB attempts retain matching inputs and outputs."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
inventory=json.loads((ROOT/'experiments/week2/chemistry-candidates.json').read_text())
records=[]
for phase in ('monomers','dimers'):
    data=json.loads((ROOT/f'reports/week2-xtb-{phase}.json').read_text())
    assert data['status']=='completed'
    assert {r['candidate_id'] for r in data['results']}=={r['candidate_id'] for r in inventory[phase]}
    assert len(data['results'])==len({r['candidate_id'] for r in data['results']})
    assert data['inventory_sha256']==sha(ROOT/'experiments/week2/chemistry-candidates.json')
    for row in data['results']:
        run=ROOT/row['run_dir'];meta=json.loads((run/'manifest.json').read_text())
        assert sha(run/'output.log')==meta['output_sha256']
        assert sha(run/'input.xyz')==meta['candidate_sha256']==row['candidate']['sha256']
        assert sha(ROOT/row['candidate']['path'])==meta['candidate_sha256']
        assert json.loads((run/'candidate.json').read_text())==row['candidate']
        assert meta['executable_sha256']==data['executable_sha256']
        assert meta['result']==row['result']
        if (run/'xtbopt.xyz').exists():assert sha(run/'xtbopt.xyz')==meta['optimized_xyz_sha256']
        if 'runner_snapshot' in meta:assert sha(ROOT/meta['runner_snapshot'])==meta['runner_sha256']
    records.append({'phase':phase,'unique_attempts':len(data['results']),'hash_checks':'passed',
                    'geometries_with_optimizer_convergence':sum(r['result']['status']=='passed' for r in data['results'])})
analysis=json.loads((ROOT/'reports/week2-xtb-analysis.json').read_text())
points=[r for r in analysis['hydrogen_endpoint_transitions'] if r['initial_motif']=='proton_transfer_start']
result={'status':'passed','records':records,'proton_transfer_starts':len(points),
        'proton_transfer_starts_converged':sum(r['engine_status']=='passed' for r in points),
        'converged_PT_starts_matching_neutral_template_H_owners':sum(r['engine_status']=='passed' and r['output_hydrogens_match_state_graph'] for r in points),
        'assignment_scope':'nearest-heavy-atom geometry; not formal electronic charge'}
(ROOT/'reports/week2-xtb-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
