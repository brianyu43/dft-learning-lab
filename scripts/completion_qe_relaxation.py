"""Fixed-cell ASE optimization using separately verified QE SCF evaluations.

The provider must schedule each requested geometry under the existing engine
lock and resource budget, then return its report path. This module never
launches QE itself or makes the conditional 24-atom work automatically eligible.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.io.espresso import read_fortran_namelist
from ase.optimize import BFGS

from completion_alloy_geometry import FixMoLayerCentroids, centroids, validate_geometry
from completion_qe_bridge import load_scf, require
from research_run import sha


def settings_fingerprint(report_path):
    """Compare physical/numerical settings, excluding run names and coordinates."""
    meta=json.loads(Path(report_path).read_text())
    inp=Path(meta['run_dir'])/Path(meta['input']).name
    with inp.open() as stream:namelist,_=read_fortran_namelist(stream)
    pseudos={dest:meta['source_sha256'][source] for source,dest in meta['assets'].items() if dest.startswith('pseudo/')}
    return dict(system=dict(namelist['system']),electrons=dict(namelist['electrons']),
                k_grid=meta['k_grid'],pseudopotential_sha256=pseudos,
                executable_sha256=meta['executable_sha256'])


class VerifiedQESCFCalculator(Calculator):
    """Cache only at identical ASE state; verify every newly supplied report."""
    implemented_properties=['energy','forces']

    def __init__(self, initial_atoms, reference_report, obtain_scf_report, method='C09'):
        super().__init__()
        self.reference=initial_atoms.copy();self.reference.calc=None;self.reference.set_constraint()
        self.reference_report=Path(reference_report).resolve()
        self.method=method;self.obtain_scf_report=obtain_scf_report
        initial=load_scf(self.reference_report,expected_method=method,expected_atoms=self.reference)
        self.reference_report_sha256=initial['audit']['report_sha256']
        self.fingerprint=settings_fingerprint(self.reference_report)
        self.evaluation_history=[]

    def calculate(self,atoms=None,properties=('energy','forces'),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        self.results={}
        require(np.array_equal(self.atoms.numbers,self.reference.numbers),'Fixed atom identities/order changed')
        require(np.array_equal(self.atoms.pbc,self.reference.pbc),'Fixed periodicity changed')
        require(np.allclose(self.atoms.cell,self.reference.cell,rtol=0,atol=1e-10),'Fixed cell changed')
        require(np.isfinite(self.atoms.positions).all(),'Non-finite requested coordinates')
        require(np.array_equal(self.atoms.get_initial_charges(),self.reference.get_initial_charges())
                and np.array_equal(self.atoms.get_initial_magnetic_moments(),self.reference.get_initial_magnetic_moments()),
                'Charge/spin annotations changed')
        if np.array_equal(self.atoms.positions,self.reference.positions):
            report=self.reference_report
            require(sha(report)==self.reference_report_sha256,'Reference report changed')
        else:
            query=self.atoms.copy();query.calc=None;query.set_constraint()
            report=Path(self.obtain_scf_report(query)).resolve()
        # Never accept a report merely because the backend says it succeeded.
        result=load_scf(report,expected_method=self.method,expected_atoms=self.atoms)
        require(settings_fingerprint(report)==self.fingerprint,'SCF settings changed during optimization')
        self.results=dict(energy=result['energy_eV'],forces=result['forces'].copy())
        self.evaluation_history.append(dict(report=str(report),audit=result['audit'],
            positions_angstrom=self.atoms.positions.tolist(),cell_angstrom=self.atoms.cell.array.tolist(),
            atomic_numbers=self.atoms.numbers.tolist()))

    def revalidate_history(self):
        """Re-read immutable output evidence before accepting a final trajectory."""
        for point in self.evaluation_history:
            require(sha(point['report'])==point['audit']['report_sha256'],'A trajectory report changed')
            atoms=self.reference.copy()
            atoms.set_positions(point['positions_angstrom'],apply_constraint=False)
            result=load_scf(point['report'],expected_method=self.method,expected_atoms=atoms)
            for key in ('input_sha256','output_sha256','xml_sha256','energy_eV'):
                require(result['audit'][key]==point['audit'][key],'Trajectory SCF evidence changed')
            require(settings_fingerprint(point['report'])==self.fingerprint,'Trajectory settings changed')
        return len(self.evaluation_history)


def minimize_fixed_registry(initial_atoms, calculator, *, fmax, max_steps, maxstep,
                            initial_curvature=70., on_step=None, allow_analytic_test=False):
    """Hold both four-Mo centroids while relaxing the allowed internal motion.

    The source SCF at the initial geometry is required for a real QE run.
    Production input eligibility, budget, calibration and controller lifecycle
    are checked by its caller, not inferred from a successful optimizer return.
    """
    require(np.isfinite(fmax) and 0<fmax<=.005,'Unreviewed projected force threshold')
    require(isinstance(max_steps,int) and not isinstance(max_steps,bool) and 0<max_steps<=100,'Invalid step budget')
    require(np.isfinite(maxstep) and 0<maxstep<=.1,'Invalid maximum BFGS step')
    require(np.isfinite(initial_curvature) and initial_curvature>0,'Invalid initial BFGS curvature')
    validate_geometry(initial_atoms)
    actual_qe=isinstance(calculator,VerifiedQESCFCalculator)
    require(actual_qe or allow_analytic_test,'Only verified QE or explicitly marked analytic tests are supported')
    if actual_qe:
        require(not initial_atoms.info.get('fixture_only',False),'Fixture coordinates cannot enter production relaxation')
    atoms=initial_atoms.copy();atoms.set_constraint()
    cell=atoms.cell.array.copy();numbers=atoms.numbers.copy();layers=atoms.arrays['layer_id'].copy()
    origin=centroids(atoms);constraint=FixMoLayerCentroids.from_atoms(atoms)
    atoms.set_constraint(constraint);atoms.calc=calculator
    trace=[]
    optimizer=BFGS(atoms,logfile=None,maxstep=maxstep,alpha=initial_curvature)

    def observe():
        require(np.array_equal(atoms.numbers,numbers) and np.array_equal(atoms.arrays['layer_id'],layers),
                'Atom/layer identity changed during optimization')
        drift=float(np.abs(centroids(atoms)-origin).max())
        require(drift<=1e-8,'Fixed Mo centroid drift')
        require(np.allclose(atoms.cell,cell,rtol=0,atol=1e-10),'Cell changed during fixed-cell optimization')
        raw=atoms.get_forces(apply_constraint=False);projected=atoms.get_forces()
        require(np.isfinite(raw).all() and np.isfinite(projected).all(),'Non-finite force')
        expected=raw.copy();constraint.adjust_forces(atoms,expected)
        require(np.allclose(projected,expected,rtol=0,atol=1e-12),'Constraint force projection mismatch')
        record=dict(step=int(optimizer.nsteps),energy_eV=float(atoms.get_potential_energy()),
                    raw_max_force_eV_angstrom=float(np.linalg.norm(raw,axis=1).max()),
                    projected_max_force_eV_angstrom=float(np.linalg.norm(projected,axis=1).max()),
                    reaction_max_force_eV_angstrom=float(np.linalg.norm(raw-projected,axis=1).max()),
                    max_Mo_centroid_drift_angstrom=drift,positions_angstrom=atoms.positions.tolist(),
                    scf_report=calculator.evaluation_history[-1]['report'] if actual_qe else None)
        trace.append(record)
        if on_step is not None:on_step(record)

    optimizer.attach(observe,interval=1)
    converged=optimizer.run(fmax=fmax,steps=max_steps)
    require(bool(converged)==(trace[-1]['projected_max_force_eV_angstrom']<fmax),'Optimizer/force convergence mismatch')
    rechecked=calculator.revalidate_history() if actual_qe else 0
    return dict(atoms=atoms,trace=trace,converged=bool(converged),
                evidence_kind='verified_QE_SCF' if actual_qe else 'analytic_test_only',
                source_SCF_evaluations_rechecked=rechecked,
                final_projected_force_threshold_eV_angstrom=fmax,
                initial_BFGS_curvature_eV_angstrom2=initial_curvature,
                Mo_reference_centroids_angstrom=origin.tolist(),fixed_cell_angstrom=cell.tolist(),
                constraint=constraint.todict(),
                scope='Fixed cell and two Mo centroids; no unconstrained minimum, phonon or project-completion claim')
