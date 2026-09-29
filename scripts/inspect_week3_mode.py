"""Inspect an archived ORCA mode without changing its geometry or result."""
import json
import os
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from research_run import ROOT, sha, save, now
from run_week2_frequency_resume import inspect_hessian, BOHR_ANGSTROM


def read_modes(path):
    lines = path.read_text().splitlines()
    at = lines.index('$normal_modes') + 1
    n, m = map(int, lines[at].split())
    modes = np.full((n, m), np.nan)
    at += 1
    while at < len(lines) and not lines[at].startswith('$'):
        if not lines[at].strip() or lines[at].lstrip().startswith('#'):
            at += 1
            continue
        columns = list(map(int, lines[at].split()))
        at += 1
        for row in range(n):
            values = lines[at].split()
            assert int(values[0]) == row and len(values) == len(columns) + 1
            modes[row, columns] = list(map(float, values[1:]))
            at += 1
    assert np.isfinite(modes).all()
    at = lines.index('$atoms')
    count = int(lines[at + 1])
    records = [line.split() for line in lines[at + 2:at + 2 + count]]
    atoms = Atoms([r[0] for r in records],
                  positions=np.array([[float(v) for v in r[2:5]] for r in records]) * BOHR_ANGSTROM,
                  masses=[float(r[1]) for r in records])
    return atoms, modes


def rigid_fraction(atoms, vector, indices):
    masses = atoms.get_masses()[indices]
    p = atoms.positions[indices]
    p = p - np.average(p, axis=0, weights=masses)
    # Cartesian velocity = translation + omega cross r.
    matrix = np.zeros((3 * len(indices), 6))
    for i, r in enumerate(p):
        matrix[3*i:3*i+3, :3] = np.eye(3)
        matrix[3*i:3*i+3, 3:] = np.stack([np.cross(e, r) for e in np.eye(3)], axis=1)
    weights = np.repeat(np.sqrt(masses), 3)
    v = vector[indices].reshape(-1)
    params = np.linalg.lstsq(matrix * weights[:, None], v * weights, rcond=None)[0]
    residual = np.sum(((v - matrix @ params) * weights)**2)
    total = np.sum((v * weights)**2)
    return dict(fraction_explained=float(1 - residual / total),
                total_mass_weighted_amplitude=float(total),
                residual_mass_weighted_amplitude=float(residual))


def main():
    job = 'w3-diflunisal--diflunisal-p05-opt-freq'
    report = ROOT / f'reports/{job}.json'
    meta = json.loads(report.read_text())
    run = Path(meta['run_dir'])
    hess = run / f'{job}.hess'
    final = read(run / f'{job}.xyz')
    freq = inspect_hessian(hess, final, meta['result']['frequencies_cm1'])
    assert len(freq) == 150 and np.flatnonzero(freq < 0).tolist() == [0]
    atoms, modes = read_modes(hess)
    assert modes.shape == (156, 156)
    vector = modes[:, 6].reshape(-1, 3)
    assert abs(np.linalg.norm(vector) - 1) < 1e-8
    mass = atoms.get_masses()
    cm = np.average(atoms.positions, axis=0, weights=mass)
    # ORCA normal_modes already are Cartesian vectors after inverse-sqrt mass conversion.
    # Independently check removal of whole-system translation and rotation.
    translation = np.sum(mass[:, None] * vector, axis=0)
    rotation = np.sum(np.cross(atoms.positions - cm, mass[:, None] * vector), axis=0)
    assert np.linalg.norm(translation) < 1e-7 and np.linalg.norm(rotation) < 1e-6
    parts = [rigid_fraction(atoms, vector, np.arange(a, b)) for a, b in [(0, 26), (26, 52)]]
    fraction = 1 - sum(p['residual_mass_weighted_amplitude'] for p in parts) / sum(p['total_mass_weighted_amplitude'] for p in parts)
    displacement = vector * (0.10 / np.linalg.norm(vector, axis=1).max())
    plus, minus = atoms.copy(), atoms.copy()
    plus.positions += displacement
    minus.positions -= displacement
    out = ROOT / 'reports/week3-recovery-r1'
    out.mkdir(exist_ok=True)
    write(out / 'mode6-visualization.xyz', [minus, atoms, plus], format='xyz')
    torsions = []
    for offset in [0, 26]:
        ids = [4+offset, 5+offset, 6+offset, 7+offset]
        center = atoms.get_dihedral(*ids)
        delta = lambda value: (value - center + 180) % 360 - 180
        torsions.append(dict(indices_zero_based=ids, center_degrees=center,
                             plus_delta_degrees=delta(plus.get_dihedral(*ids)),
                             minus_delta_degrees=delta(minus.get_dihedral(*ids))))
    result = dict(checked_at=now(),parent_job_id=job,parent_report_sha256=sha(report),
                  hessian_sha256=sha(hess),geometry_sha256=sha(run / f'{job}.xyz'),
                  hessian_geometry_verified=True,mode_index_zero_based=6,
                  frequency_cm1=float(freq[0]),internal_mode_count=len(freq),negative_mode_count=1,
                  whole_system_translation_residual=translation.tolist(),
                  whole_system_rotation_residual=rotation.tolist(),
                  per_fragment_rigid_motion=parts,combined_rigid_motion_fraction=fraction,
                  biphenyl_dihedrals_at_visualization_amplitude=torsions,
                  visualization_max_atom_displacement_angstrom=.10,
                  visualization_scope='Arbitrary amplified displacement, not thermal amplitude or calculated alternative minimum',
                  interpretation='Predominantly relative fragment motion' if fraction > .8 else 'Mixed intermolecular and internal motion',
                  cause='Unresolved: tighter optimization and numerical integration required; this mode inspection alone does not prove numerical noise',
                  retry_starts_from='Undisplaced final optimized XYZ; no extra independent starting pose')
    save(out / 'mode-inspection.json', result)
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    candidate = json.loads((run / 'candidate.json').read_text())
    p = atoms.positions - atoms.positions.mean(axis=0)
    _, _, vt = np.linalg.svd(p[atoms.numbers > 1])
    p, arrows = p @ vt.T, displacement @ vt.T
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), layout='constrained')
    colors = {1:'#bec4cc', 6:'#4b5868', 8:'#d04444', 9:'#278e68'}
    for ax, dims in zip(axes, [(0, 1), (0, 2)]):
        q, v = p[:, dims], arrows[:, dims] * 5
        for i,j,_ in candidate['graph']['bonds']:
            ax.plot(*q[[i,j]].T, color='#8c97a4', lw=1.5, zorder=1)
        ax.scatter(*q.T, c=[colors[z] for z in atoms.numbers], s=[22 if z==1 else 70 for z in atoms.numbers], zorder=2)
        ax.quiver(q[:,0],q[:,1],v[:,0],v[:,1],angles='xy',scale_units='xy',scale=1,color='#c46d0b',width=.003)
        ax.set_aspect('equal'); ax.axis('off'); ax.margins(.1)
    fig.suptitle(f'DIF-DIF mode 6: {freq[0]:.3f} cm$^{{-1}}$\nTwo orthogonal views; arrows amplified for inspection; no new optimization', fontsize=13)
    fig.savefig(out / 'mode6-review.png', dpi=180)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
