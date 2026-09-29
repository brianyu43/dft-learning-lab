"""Check frozen week-two evidence and exact preservation of the week-one snapshot."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
index=json.loads((ROOT/'reports/week2-evidence-sha256.json').read_text())
failures=[p for p,h in index['files'].items() if sha(ROOT/p)!=h]
visual_qa_path=ROOT/'reports/week2-figure-visual-qa.json';visual_qa_count=0
if visual_qa_path.exists():
    visual_qa=json.loads(visual_qa_path.read_text())['png_sha256'];visual_qa_count=len(visual_qa)
    failures.extend('visual_qa:'+p for p,h in visual_qa.items() if sha(ROOT/p)!=h)
parallel_index=ROOT/'reports/week2-parallel-evidence-sha256.json'
parallel_count=0
if parallel_index.exists():
    parallel=json.loads(parallel_index.read_text())['files'];parallel_count=len(parallel)
    failures.extend('parallel:'+p for p,h in parallel.items() if sha(ROOT/p)!=h)
pause_preservation=ROOT/'reports/week2-resume-preservation.json';pause_count=0
if pause_preservation.exists():
    record=json.loads(pause_preservation.read_text())
    old=json.loads((ROOT/record['archived_original_index']).read_text())['files'];pause_count=len(old)
    failures.extend('pause_archive:'+p for p,h in old.items() if sha(ROOT/record['files'][p])!=h)
asset_count=0
for manifest in [*(ROOT/'runs').glob('*-w2-*/manifest.json'),*(ROOT/'runs').glob('*-p2-*/manifest.json')]:
    meta=json.loads(manifest.read_text());run=manifest.parent
    if sha(run/'job_config.json')!=meta['config_sha256']:failures.append('config:'+str(run))
    for source,destination in meta.get('assets',{}).items():
        asset_count+=1
        if sha(run/destination)!=meta['source_sha256'][source]:failures.append('asset:'+str(run/destination))
week1=json.loads((ROOT/'reports/week1-evidence-sha256.json').read_text())['files']
archives={'README.md':'docs/archive/README-before-week2.md',
          'experiments/protocol-v1.json':'experiments/archive/protocol-v1-before-week2.json'}
preserved={};unchanged=0
for path,expected in week1.items():
    if sha(ROOT/path)==expected:unchanged+=1;continue
    archive=archives.get(path)
    if archive and sha(ROOT/archive)==expected:preserved[path]=archive
    else:failures.append('week1:'+path)
result={'status':'passed' if not failures else 'failed','week2_hashed_files_checked':len(index['files']),
        'visually_checked_PNGs_matching_final_files':visual_qa_count,
        'parallel_validation_files_checked':parallel_count,
        'paused_files_preserved_in_archives':pause_count,'copied_DFT_assets_checked':asset_count,
        'week1_original_files_unchanged':unchanged,'week1_live_metadata_preserved_in_exact_archives':preserved,
        'mismatches':failures,'scope':'byte-level integrity, separate from scientific interpretation and numerical validation'}
(ROOT/'reports/week2-evidence-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if not failures else 1)
