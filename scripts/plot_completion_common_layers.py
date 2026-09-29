"""Static research figure of audited initial/final common-cell layer thickness."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ase.io import read

from prepare_week3_drug import immutable
from research_run import ROOT, sha
from verify_completion_common_layer import verify


def plot():
    folder = ROOT / 'reports/figures/completion/common-layers-r1'
    paths = {kind: folder / f'layer-thickness.{kind}' for kind in ('png', 'svg', 'json')}
    qa_path = ROOT / 'reports/completion/common-layer-figure-r1.json'
    if qa_path.exists():
        qa = json.loads(qa_path.read_text())
        for section in ('source_sha256', 'output_sha256'):
            for p, digest in qa[section].items():
                if sha(ROOT / p) != digest:
                    raise ValueError(f'Changed figure evidence: {p}')
        return paths
    if any(p.exists() for p in paths.values()):
        raise ValueError('Partial figure already exists; preserve and review it')
    rows, sources = [], {}
    for material in ('mos2', 'mose2'):
        audit = verify(material)
        saved = ROOT / f'reports/completion/common-layer-audits-r1/{material}.json'
        if json.loads(saved.read_text()) != audit:
            raise ValueError('Saved audit differs from actual raw common-layer results')
        meta = json.loads((ROOT / f"reports/{audit['job_id']}.json").read_text())
        initial = read(Path(meta['run_dir']) / Path(meta['input']).name, format='espresso-in')
        z_initial = initial.positions[:, 2]
        z_final = np.array(audit['atomic_positions_angstrom'])[:, 2]
        before, after = float(max(z_initial) - min(z_initial)), float(max(z_final) - min(z_final))
        if abs(before - audit['initial_layer_thickness_angstrom']) > 1e-12 or abs(after - audit['relaxed_layer_thickness_angstrom']) > 1e-12:
            raise ValueError('Thickness does not reproduce from coordinates')
        rows.append(dict(material=material, input_thickness_angstrom=before,
            relaxed_thickness_angstrom=after, delta_thickness_angstrom=after - before,
            imposed_strain_percent=100 * audit['imposed_strain_fraction'],
            common_lattice_angstrom=audit['common_lattice_angstrom'],
            max_final_force_eV_angstrom=audit['summary']['max_force_eV_angstrom'],
            source_job_id=audit['job_id']))
        sources.update(audit['source_sha256'])
        sources[str(saved.relative_to(ROOT))] = sha(saved)
    sources[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    if rows[0]['common_lattice_angstrom'] != rows[1]['common_lattice_angstrom']:
        raise ValueError('The two layers do not share the imposed lattice')
    data = dict(rows=rows, units='angstrom', y_axis_range_angstrom=[3.02, 3.46],
        thickness_definition='z_max - z_min of each isolated three-atom layer (outer S/Se planes)',
        comparison='Before versus after internal relaxation at the same fixed imposed cell for each material',
        independent_material_count=2, relaxed_geometries_per_material=1,
        uncertainty_intervals_estimated=False, elastic_fit_performed=False,
        phonon_stability_verified=False, source_sha256=sources,
        matplotlib_version=matplotlib.__version__)
    ink, blue, gray = '#263442', '#23648a', '#6b747d'
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'svg.hashsalt': 'common-layer-thickness-r1', 'text.color': ink,
                         'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink}):
        fig, ax = plt.subplots(figsize=(7.8, 5.4))
        fig.subplots_adjust(left=.14, right=.95, bottom=.27, top=.72)
        for i, row in enumerate(rows):
            before, after = row['input_thickness_angstrom'], row['relaxed_thickness_angstrom']
            ax.plot([i - .11, i + .11], [before, after], color=gray, linewidth=1.2, zorder=2)
            ax.scatter(i - .11, before, marker='o', s=64, facecolor='white', edgecolor=gray, linewidth=1.5, zorder=3)
            ax.scatter(i + .11, after, marker='D', s=56, color=blue, zorder=3)
            direction = 1 if before > after else -1
            ax.annotate(f'{before:.4f}', (i - .11, before), xytext=(-8, 12 * direction),
                        textcoords='offset points', ha='right', va='center', color=gray, fontsize=10)
            ax.annotate(f'{after:.4f}', (i + .11, after), xytext=(8, -12 * direction),
                        textcoords='offset points', ha='left', va='center', color=blue, fontsize=10)
            ax.text(i + .38, (before + after) / 2, f"Δt = {after - before:+.4f} Å", va='center', fontsize=10)
        ax.set_xlim(-.6, 1.95)
        ax.set_ylim(*data['y_axis_range_angstrom'])
        ax.set_xticks([0, 1], [r'MoS$_2$' + f"\nImposed strain {rows[0]['imposed_strain_percent']:+.3f}%",
                              r'MoSe$_2$' + f"\nImposed strain {rows[1]['imposed_strain_percent']:+.3f}%"])
        ax.tick_params(axis='x', pad=9, length=0)
        ax.set_yticks(np.arange(3.05, 3.451, .1))
        ax.set_ylabel('Layer thickness t (Å)')
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e2e7eb', linewidth=.6)
        ax.set_axisbelow(True)
        fig.text(.14, .93, 'Layer thickness after internal relaxation', fontsize=17, fontweight='bold')
        fig.text(.14, .87, f"Fixed common lattice a = {rows[0]['common_lattice_angstrom']:.6f} Å · C09 · 100/400 Ry · 9×9×1 k",
                 fontsize=10)
        handles = [Line2D([], [], marker='o', linestyle='none', markerfacecolor='white', markeredgecolor=gray,
                          markersize=7, label='Input at imposed cell'),
                   Line2D([], [], marker='D', linestyle='none', color=blue, markersize=6, label='Relaxed at same cell')]
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.13, .835), frameon=False, ncol=2, fontsize=10)
        fig.text(.14, .065, 'Thickness = separation of the outer S/Se planes. Focused y-axis.\n'
                 'One imposed strain per material; no elastic fit or phonon-stability claim.',
                 fontsize=9, color='#4b5967', linespacing=1.5)
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(paths['png'], dpi=240, metadata={'Software': 'matplotlib'})
        fig.savefig(paths['svg'], metadata={'Date': None, 'Creator': 'matplotlib'})
        plt.close(fig)
    immutable(paths['json'], json.dumps(data, indent=2) + '\n')
    qa = dict(status='generated_from_two_raw_verified_common_layers', plotted_rows=len(rows),
        source_sha256=sources, output_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths.values()},
        displayed_values_recomputed_from_coordinates=True, shared_y_scale=True,
        legend_has_marker_shape_and_fill_distinctions=True, uncertainty_intervals_estimated=False,
        full_materials_baseline_complete=False, visual_inspection_recorded_separately=True)
    immutable(qa_path, json.dumps(qa, indent=2) + '\n')
    return paths


if __name__ == '__main__':
    print(json.dumps({k: str(v) for k, v in plot().items()}, indent=2))
