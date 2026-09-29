"""Fingerprint pre-campaign scientific evidence; never relabel or rewrite old runs."""
import argparse
import json
from pathlib import Path

from prepare_week3_drug import immutable
from research_run import ROOT,sha,now

TARGET=ROOT/'experiments/completion/preserved-original-evidence.json'


def freeze():
    if TARGET.exists():return verify()
    paths=set();runs=[]
    for manifest in (ROOT/'runs').glob('*/manifest.json'):
        meta=json.loads(manifest.read_text())
        if meta.get('engine') not in ['orca','qe']:continue
        if not meta.get('started_at','').startswith(('2026-09-09','2026-09-10','2026-09-11')):continue
        assert meta.get('status')!='running',f'Inspect original incomplete run before freezing: {manifest}'
        run=manifest.parent;runs.append(dict(path=str(run.relative_to(ROOT)),status=meta['status'],job_id=meta.get('job_id'),purpose=meta.get('purpose')))
        for name,digest in meta.get('input_sha256',{}).items():
            assert sha(run/name)==digest
            paths.add(run/name)
        for source,digest in meta.get('source_sha256',{}).items():
            destination=Path(meta['input']).name if source==meta['input'] else meta['assets'][source]
            copied=run/destination
            # Engine libraries/checkpoints are large and are covered by their
            # original manifest hashes, not duplicated in this evidence index.
            if copied.suffix in ['.dylib','.gbw']:continue
            assert sha(copied)==digest,f'Original copied asset changed: {copied}'
            paths.add(copied)
        if meta.get('output_sha256'):assert sha(run/'output.log')==meta['output_sha256']
        for p in run.iterdir():
            if p.is_file() and (p.suffix in ['.json','.log','.xyz','.hess','.engrad','.inp','.in','.py','.sh']):paths.add(p)
        paths.update(run.rglob('data-file-schema.xml'))
        if meta.get('job_id'):
            rp=ROOT/f"reports/{meta['job_id']}.json"
            if rp.exists():paths.add(rp)
    for sub in ['reports/week3-recovery-r1','experiments/week3/drug','inputs/orca/week3','structures/week3/drug_initial']:
        for p in (ROOT/sub).rglob('*'):
            if p.is_file():paths.add(p)
    for p in (ROOT/'reports').glob('week[12]*'):
        if p.is_file():paths.add(p)
    paths.update(ROOT/p for p in ['reports/week3-summary.json','reports/week3-run-state.json'])
    record=dict(created_at=now(),scope='Preserved pre-2026-09-29 ORCA/QE scientific records and first three week drug provenance; checkpoint binaries excluded explicitly.',
        original_run_count=len(runs),runs=runs,files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)},
        note='Fingerprint taken after historical outputs and asset digests were checked. This index does not imply scientific acceptance of failed/interrupted calculations.')
    immutable(TARGET,json.dumps(record,indent=2)+'\n')
    return verify()


def verify():
    record=json.loads(TARGET.read_text())
    for p,digest in record['files'].items():assert sha(ROOT/p)==digest,f'Preserved evidence changed: {p}'
    return dict(status='preserved',original_run_count=record['original_run_count'],file_count=len(record['files']),manifest_sha256=sha(TARGET))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');args=parser.parse_args()
    print(json.dumps(freeze() if args.freeze else verify()))
