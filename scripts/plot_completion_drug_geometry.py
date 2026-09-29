"""Two common-frame projections of three archived, verified DIF-CAF geometries."""
import json
import os
from pathlib import Path

import numpy as np
from ase.io import read

from completion_drug_geometry import align
from research_run import ROOT, sha

os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = ROOT / 'reports/figures/completion/drug-poses-r1'


def main():
    review_path = ROOT / 'reports/completion/drug-archived-geometry-r1.json'
    lineage_path = ROOT / 'experiments/completion/drug-start-lineage-r1.json'
    review = json.loads(review_path.read_text()); lineage = json.loads(lineage_path.read_text())
    assert sha(lineage_path) == review['lineage_sha256']
    assert sha(ROOT / 'scripts/completion_drug_geometry.py') == review['script_sha256']
    by_id = {r['candidate_id']: r for r in lineage['rows']}
    reviewed = {r['candidate_id']: r for r in review['rows']}
    source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in (review_path, lineage_path, Path(__file__))}
    ref_row = by_id['diflunisal-t02-s43']
    ref_path = ROOT / ref_row['source_geometry']
    assert sha(ref_path) == ref_row['source_geometry_sha256']
    source_hashes[str(ref_path.relative_to(ROOT))] = sha(ref_path)
    ref = read(ref_path); heavy = np.flatnonzero(ref.numbers > 1)
    center = ref.positions[heavy].mean(axis=0)
    _, _, vt = np.linalg.svd(ref.positions[heavy] - center)
    x, y = vt[:2].copy()
    if np.dot(x, ref.positions[9] - ref.positions[1]) < 0: x *= -1
    if np.dot(y, ref.positions[17] - ref.positions[16]) < 0: y *= -1
    frame = np.column_stack([x, y, np.cross(x, y)])
    assert abs(np.linalg.det(frame) - 1.) < 1e-12
    labels = {
        'p36': 'Phenol OH to carbonyl O39',
        'p14': 'Phenol OH to carbonyl O29',
        'p06': 'Carboxylic OH to N9 (atom 32)',
    }
    rows = []
    for short, label in labels.items():
        cid = 'diflunisal--caffeine-' + short
        row, geometry = by_id[cid], reviewed[cid]
        path = ROOT / geometry['source_geometry']
        assert geometry['archived_minimum_verified'] and sha(path) == geometry['source_geometry_sha256']
        source_hashes[str(path.relative_to(ROOT))] = sha(path)
        report = ROOT / geometry['minimum_report']
        assert sha(report) == lineage['input_hashes'][str(report.relative_to(ROOT))]
        source_hashes[str(report.relative_to(ROOT))] = sha(report)
        atoms = read(path)
        assert atoms.numbers.tolist() == row['candidate']['graph']['atomic_numbers']
        assert row['candidate']['fragment_atom_counts'] == [26, 24]
        residual, rotation, pc, qc = align(atoms.positions[heavy], ref.positions[heavy])
        positions = ((atoms.positions - pc) @ rotation + qc - center) @ frame
        old_dist = np.linalg.norm(atoms.positions[:, None, :] - atoms.positions[None, :, :], axis=2)
        new_dist = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)
        delta = float(np.max(np.abs(old_dist - new_dist)))
        assert delta < 1e-10  # A single rigid transform preserves the full complex.
        contact, = geometry['final']['contacts']
        assert abs(np.linalg.norm(positions[contact['hydrogen']] - positions[contact['acceptor']])
                   - contact['H_acceptor_distance_angstrom']) < 1e-10
        rows.append(dict(candidate_id=cid, short=short, description=label,
                         atomic_numbers=atoms.numbers.tolist(), bonds=row['candidate']['graph']['bonds'],
                         displayed_positions_angstrom=positions.tolist(), contact=contact,
                         reference_DIF_heavy_alignment_RMSD_angstrom=residual,
                         maximum_pair_distance_change_angstrom=delta,
                         rigid_transform=dict(rotation_to_DIF=rotation.tolist(), source_center=pc.tolist(),
                                              target_center=qc.tolist(), display_frame=frame.tolist(), display_origin=center.tolist()),
                         source_geometry=str(path.relative_to(ROOT)), source_geometry_sha256=sha(path)))

    colors = {1: '#bfc5cf', 7: '#2769bb', 8: '#d73d49', 9: '#14916f'}
    carbon_dif, carbon_caf, contact_color = '#303b4b', '#c78a25', '#8542ad'
    all_positions = np.concatenate([np.array(row['displayed_positions_angstrom']) for row in rows])
    limits = [[float(all_positions[:, axis].min() - 1.1), float(all_positions[:, axis].max() + 1.1)] for axis in range(3)]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none', 'svg.hashsalt': 'dft-drug-poses-r1'})
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.7))
    for col, row in enumerate(rows):
        p = np.array(row['displayed_positions_angstrom']); z = np.array(row['atomic_numbers']); c = row['contact']
        for view, dims in enumerate(((0, 1), (0, 2))):
            ax = axes[view, col]; q = p[:, dims]
            for i, j, _ in row['bonds']:
                ax.plot(*q[[i, j]].T, color='#aab2be', lw=1.2, zorder=1)
            for i in range(len(z)):
                color = (carbon_dif if i < 26 else carbon_caf) if z[i] == 6 else colors[int(z[i])]
                ax.scatter(*q[i], s=19 if z[i] == 1 else 57, c=color, edgecolors='white', linewidths=.45, zorder=3)
            h, acceptor, donor = (c[k] for k in ('hydrogen', 'acceptor', 'donor'))
            ax.plot(*q[[h, acceptor]].T, color=contact_color, lw=2., ls='--', zorder=4)
            if view == 0:
                for index, offset in ((donor, (-21, 15)), (acceptor, (21, -18))):
                    symbol = 'N' if z[index] == 7 else 'O'
                    ax.annotate(f'{symbol}{index}', q[index], xytext=offset, textcoords='offset points', ha='center',
                                fontsize=8.5, color=contact_color, zorder=5,
                                arrowprops=dict(arrowstyle='-', color=contact_color, lw=.6),
                                bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=1.5))
                ax.set_title(f"{row['short']}  |  {row['description']}\n"
                             f"H···A {c['H_acceptor_distance_angstrom']:.3f} Å  /  D–H···A {c['DHA_angle_degrees']:.1f}°",
                             fontsize=11, pad=10)
            ax.set_xlim(limits[0]); ax.set_ylim(limits[dims[1]]); ax.set_aspect('equal'); ax.axis('off')
            ax.text(.01, .98, 'XY projection' if view == 0 else 'XZ projection', transform=ax.transAxes,
                    va='top', color='#637082', fontsize=8.5)
            # The scale bar is in the displayed Cartesian coordinate plane.
            sx, sy = limits[0][0] + .7, limits[dims[1]][0] + .45
            ax.plot([sx, sx + 2.], [sy, sy], color='#465363', lw=2)
            ax.text(sx + 1., sy + .15, '2 Å', ha='center', va='bottom', fontsize=8, color='#465363')
    legend = [Line2D([], [], marker='o', ls='', color=color, markersize=6, label=label)
              for label, color in [('C in DIF', carbon_dif), ('C in CAF', carbon_caf),
                                   ('N', colors[7]), ('O', colors[8]), ('F', colors[9]), ('H', colors[1])]]
    legend.append(Line2D([], [], color=contact_color, ls='--', label='Geometric H···A contact'))
    fig.suptitle('Three archived diflunisal–caffeine molecular-pair minima', fontsize=17, y=.98, weight='bold')
    fig.text(.5, .935, 'r2SCAN-3c / vacuum · original-run geometries · common-grid refinement remains pending',
             ha='center', fontsize=10.5, color='#4b5c70')
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, .071), ncol=7, frameon=False, fontsize=9)
    fig.text(.5, .041, 'All pairs share one DIF reference frame. Each entire complex is rotated rigidly; both projections use common coordinate limits.',
             ha='center', fontsize=8.5, color='#4b5c70')
    fig.text(.5, .020, 'Connectivity follows preserved molecular templates. Contacts are geometric descriptors; these are not measured crystal structures or polymorphs.',
             ha='center', fontsize=8.5, color='#4b5c70')
    fig.subplots_adjust(left=.028, right=.98, bottom=.145, top=.855, wspace=.065, hspace=.04)
    OUT.mkdir(parents=True, exist_ok=True)
    png, svg = OUT / 'archived-dif-caf-poses.png', OUT / 'archived-dif-caf-poses.svg'
    fig.savefig(png, dpi=190, facecolor='white'); fig.savefig(svg, metadata={'Date': None}); plt.close(fig)
    record = dict(source_sha256=source_hashes, rows=rows, common_coordinate_limits_angstrom=limits,
                  projections=[[0, 1], [0, 2]], generated_from_actual_archived_coordinates=True,
                  estimated_coordinates=False, cocrystal_polymorphs=False,
                  artifacts={str(path.relative_to(ROOT)): sha(path) for path in (png, svg)},
                  scope='Archived molecular geometry illustration; no method/solvent ranking or common-grid acceptance')
    (OUT / 'figure-data.json').write_text(json.dumps(record, indent=2) + '\n')
    print(png)


if __name__ == '__main__':
    main()
