"""ORCA validation for the completion campaign; local minima and fixed-geometry SP are distinct."""
import argparse
import json
from pathlib import Path
import re
import signal

import numpy as np
from ase.io import read

import research_run
from research_run import ROOT, sha, save
from run_week2_xtb_screen import geometry_audit
from run_week2_frequency_resume import inspect_hessian
from week3_recovery import gradient_in_frame

BASE_PARSE=research_run.parse_output


def checked_parse(run,meta):
    result=BASE_PARSE(run,meta)
    result['scientific_verification_passed']=False
    try:
        assert sha(run/Path(meta['input']).name)==meta['expected_input_sha256']
        for source,digest in meta['source_sha256'].items():
            dest=Path(meta['input']).name if source==meta['input'] else meta['assets'][source]
            assert sha(run/dest)==digest
        assert sha(run/'initial.xyz')==meta['initial_structure_sha256']
        assert sha(run/'candidate.json')==meta['candidate_sha256']
        assert sha(run/'completion_parser_snapshot.py')==meta['completion_parser_sha256']
        text=(run/'output.log').read_text()
        ranks=[int(x) for x in re.findall(r'running with\s+(\d+)\s+parallel MPI-processes',text,re.I)]
        assert ranks and all(x==meta['mpi_processes'] for x in ranks), 'Unverified MPI process count'
        if 'expected_scf_grid' in meta:
            accuracies=[float(x) for x in re.findall(r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s+([\d.]+)',text)]
            angular=[int(x) for x in re.findall(r'Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.\.\.\s+(\d+)',text)]
            # The first grid construction describes the ground-state XC grid;
            # later response/Hessian grids have their own settings.
            expected=meta['expected_scf_grid']
            assert accuracies and abs(accuracies[0]-expected['IntAcc'])<.0006
            assert angular and angular[0]==expected['AngularGrid']
            result['applied_scf_grid_verified']={'IntAcc':accuracies[0],'AngularGrid':angular[0]}
        candidate=json.loads((run/'candidate.json').read_text())
        initial=read(run/'initial.xyz')
        geometry=run/(Path(meta['input']).stem+'.xyz') if meta['requested'].get('optimization') else run/'initial.xyz'
        atoms=read(geometry)
        assert len(atoms)==len(initial)==meta['n_atoms']
        assert atoms.get_chemical_symbols()==initial.get_chemical_symbols()
        assert np.isfinite(atoms.positions).all()
        audit=geometry_audit(candidate,atoms)
        assert audit['heavy_connectivity_preserved']
        assert not audit['unassigned_hydrogens'] and not audit.get('separated_geometry',False)
        result.update(geometry_audit=audit,geometry_sha256=sha(geometry),
                      optimized_xyz_sha256=sha(geometry) if meta['requested'].get('optimization') else None)
        # Independent .engrad coordinate/energy check for both SP and Opt+Freq.
        check_meta={**meta,'run_dir':str(run),'result':result}
        energy,gradient,path,digest=gradient_in_frame(check_meta,atoms)
        result.update(engrad_sha256=digest,gradient_Hartree_Bohr=gradient.tolist(),
                      gradient_geometry_verified=True,max_gradient_component_Hartree_Bohr=float(np.abs(gradient).max()))
        if meta['requested'].get('frequencies'):
            hess=run/(Path(meta['input']).stem+'.hess')
            frequencies=inspect_hessian(hess,atoms,result['frequencies_cm1'])
            result.update(full_precision_frequencies_cm1=frequencies.tolist(),hessian_sha256=sha(hess),hessian_geometry_matches=True)
            assert np.all(frequencies>0), 'A nonpositive internal mode requires review'
        assert result['status']=='passed', 'Engine/minimum verification incomplete'
        result['observed_mpi_process_counts']=sorted(set(ranks))
        result['scientific_verification_passed']=True
    except (AssertionError,ValueError,KeyError,FileNotFoundError,IndexError,StopIteration) as exc:
        result.update(status='needs_review',minimum_verified=False,
                      completion_validation_error=f'{type(exc).__name__}: {exc}')
    save(run/'result.json',result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    cfg=json.loads(args.config.read_text())
    assert cfg['engine']=='orca' and cfg['mpi_processes']==12 and cfg['omp_threads']==1
    assert cfg['completion_parser_sha256']==sha(Path(__file__))
    assert sha(ROOT/cfg['input'])==cfg['expected_input_sha256']
    assert cfg['rss_limit_gib']<=40 and cfg['timeout_seconds']<=14400
    assert not (ROOT/f"reports/{cfg['job_id']}.json").exists(), 'Review an existing attempt; never overwrite'
    research_run.parse_output=checked_parse
    return research_run.run_job(args.config)


if __name__=='__main__':
    signal.signal(signal.SIGTERM,lambda *_:(_ for _ in ()).throw(KeyboardInterrupt('terminated')))
    raise SystemExit(main())
