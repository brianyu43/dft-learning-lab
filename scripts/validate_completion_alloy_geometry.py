"""Save reproducible, explicitly illustrative geometry and constraint checks."""
import io
import json
import os
from pathlib import Path
import sys
import unittest

import ase
from ase.io import read, write
import numpy as np
import spglib

from completion_alloy_geometry import build_family, decoration_equivalence, validate_geometry
from prepare_week3_drug import immutable
from research_run import ROOT, now, save, sha
from week3_qe import monolayer


def illustration(family, output):
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    colors = {'Mo': '#4B6475', 'S': '#E5B643', 'Se': '#AA4C73'}
    labels = {'C0': 'C0: no exchange', 'CI': 'CI: facing surfaces', 'CO': 'CO: outer surfaces'}
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 6.7), sharex=True, sharey=True)
    for ax, (label, atoms) in zip(axes, family.items()):
        z = atoms.positions[:, 2] - atoms.positions[:, 2].mean()
        horizontal = atoms.positions[:, 0] + 0.37 * atoms.positions[:, 1]
        for element in colors:
            mask = np.array(atoms.get_chemical_symbols()) == element
            ax.scatter(horizontal[mask], z[mask], s=95 if element == 'Mo' else 115,
                       color=colors[element], edgecolors='white', linewidths=0.8, zorder=3)
        swapped = atoms.info['exchange_indices']
        if swapped:
            ax.scatter(horizontal[swapped], z[swapped], s=230, facecolors='none', edgecolors='#222222',
                       linewidths=1.2, zorder=4)
            for i in swapped:
                ax.annotate(f'#{i}', (horizontal[i], z[i]), xytext=(6, 8), textcoords='offset points', fontsize=9)
        for layer in (0, 1):
            mask = (atoms.arrays['layer_id'] == layer) & (atoms.numbers == 42)
            ax.axhline(z[mask].mean(), color='#9AADB8', linestyle='--', linewidth=0.8, zorder=1)
        ax.set_title(labels[label], fontsize=13, pad=16)
        ax.set_xlabel('Oblique side projection: x + 0.37y (Å)')
        ax.grid(axis='y', alpha=0.12)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(-6.2, 6.2)
    axes[0].set_ylabel('z relative to cell contents (Å)')
    handles = [Line2D([], [], marker='o', linestyle='', color=c, label=e, markersize=9) for e, c in colors.items()]
    handles.append(Line2D([], [], marker='o', linestyle='', markerfacecolor='none', color='#222222',
                          label='Exchanged site (zero-based index)', markersize=12))
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.13), ncol=4, frameon=False)
    fig.suptitle('24-atom comparison: equal total composition, different exchange faces', fontsize=15, y=0.98)
    fig.text(0.5, 0.06, 'ILLUSTRATIVE INPUT CHECK — a = 3.24 Å, D = 6.70 Å, r0; not optimized material results.\n'
             'All panels: Mo8S8Se8. CI and CO have the same composition in each layer. Cell vacuum omitted from view.',
             ha='center', fontsize=10, color='#3A4650')
    fig.subplots_adjust(left=0.07, right=0.99, top=0.87, bottom=0.28, wspace=0.12)
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix('.svg'))
    plt.close(fig)


def main():
    sys.path.insert(0, str(ROOT / 'tests'))
    from test_completion_alloy_geometry import constraint_optimizer_trial
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromName('test_completion_alloy_geometry')
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    report_dir = ROOT / 'reports/completion'
    report_dir.mkdir(exist_ok=True)
    (report_dir / 'alloy-geometry-tests.log').write_text(stream.getvalue())
    if not result.wasSuccessful():
        raise RuntimeError('Geometry tests failed; no successful fixture validation will be written')

    bottom, top = monolayer('S', 3.24), monolayer('Se', 3.24)
    structures, equivalences = [], []
    for rid in range(3):
        for distance in (6.45, 6.7, 6.95):
            family = build_family(bottom, top, rid, distance, 25)
            for label, atoms in family.items():
                atoms.info.update(fixture_only=True, scientific_equilibrium_verified=False)
                path = ROOT / f'structures/completion/geometry-fixtures/{label}-r{rid}-d{distance:.2f}.extxyz'
                text = io.StringIO()
                write(text, atoms, format='extxyz')
                immutable(path, text.getvalue())
                recovered = read(path)
                check = validate_geometry(recovered)
                np.testing.assert_allclose(recovered.positions, atoms.positions, rtol=0, atol=1e-8)
                structures.append(dict(path=str(path.relative_to(ROOT)), sha256=sha(path), **check))
            checks = {str(tol): decoration_equivalence(family['CI'], family['CO'], tol)
                      for tol in (1e-5, 1e-4, 1e-3)}
            if any(check['equivalent'] for check in checks.values()):
                raise RuntimeError('CI and CO cannot be counted as independent configurations')
            equivalences.append(dict(registry_id=rid, distance_angstrom=distance, checks=checks))
    trial = constraint_optimizer_trial()
    if not trial['converged'] or trial['maximum_centroid_drift_angstrom'] > 1e-12:
        raise RuntimeError('Centroid constraint trial failed')
    illustration(build_family(bottom, top, 0, 6.7, 25), report_dir / 'alloy-geometry-fixture.png')
    sources = ['scripts/completion_alloy_geometry.py', 'scripts/validate_completion_alloy_geometry.py',
               'tests/test_completion_alloy_geometry.py', 'scripts/week3_materials_pipeline.py',
               'scripts/week3_qe.py', 'structures/optimized/week1/MoS2_fixed_cell_relaxed.cif']
    output = dict(status='geometry_and_constraint_fixture_checks_passed', checked_at=now(),
                  scientific_alloy_calculations_completed=False, production_execution_enabled=False,
                  scope='Illustrative a=3.24 A and D=6.7 A with +/-0.25 A; not a fitted equilibrium',
                  ase_version=ase.__version__, spglib_version=spglib.__version__,
                  unit_tests_run=result.testsRun, unit_test_failures=len(result.failures),
                  unit_test_errors=len(result.errors), structures=structures,
                  ci_co_symmetry_checks=equivalences, constraint_trial=trial,
                  source_sha256={p: sha(ROOT / p) for p in sources},
                  production_preconditions=[
                      'Verified common-lattice monolayers and D0 from the completed nine-point reference',
                      'Repeat composition, geometry and CI/CO symmetry checks on actual production structures',
                      'Review measured compute budget before any 24-atom engine job',
                      'Verify input/output settings and separated layers in the actual QE workflow',
                      'Before constrained material relaxation, validate QE/ASE force and unit transfer and cost'],
                  references=['https://spglib.readthedocs.io/en/v2.7.0/api/autodoc/spglib.spg.html',
                              'https://docs.ase-lib.org/ase/constraints.html'])
    save(report_dir / 'alloy-geometry-fixture-validation.json', output)
    print(json.dumps({k: output[k] for k in ['status', 'unit_tests_run', 'constraint_trial']}, indent=2))


if __name__ == '__main__':
    main()
