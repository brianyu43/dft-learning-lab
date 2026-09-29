"""Plot the three actual lattice points after revalidating their raw QE evidence."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from prepare_week3_drug import immutable
from research_run import ROOT, sha
from verify_completion_lattice_checks import brackets


def plot(material):
    audit = brackets(material)
    saved = ROOT / f'reports/completion/lattice-check-audits-r1/{material}-brackets.json'
    if saved.exists():
        if json.loads(saved.read_text()) != audit:
            raise ValueError('Saved bracket audit differs from current raw evidence')
    else:
        immutable(saved, json.dumps(audit, indent=2) + '\n')
    folder = ROOT / 'reports/figures/completion/lattice-brackets-r1'
    paths = {kind: folder / f'{material}.{kind}' for kind in ('png', 'svg', 'json')}
    qa_path = ROOT / f'reports/completion/lattice-figure-r1-{material}.json'
    if qa_path.exists():
        qa = json.loads(qa_path.read_text())
        for section in ('source_sha256', 'output_sha256'):
            for path, digest in qa[section].items():
                if sha(ROOT / path) != digest:
                    raise ValueError(f'Changed figure evidence: {path}')
        return paths
    if any(path.exists() for path in paths.values()):
        raise ValueError('Partial figure output exists; inspect it before writing')
    points = [dict(delta_a_angstrom=-.005, delta_energy_meV=1000 * audit['points'][0]['delta_energy_eV']),
              dict(delta_a_angstrom=0., delta_energy_meV=0.),
              dict(delta_a_angstrom=.005, delta_energy_meV=1000 * audit['points'][1]['delta_energy_eV'])]
    sources = dict(audit['source_sha256'])
    for path in (saved, Path(__file__)):
        sources[str(path.relative_to(ROOT))] = sha(path)
    data = dict(material=material, equilibrium_a_angstrom=audit['equilibrium']['lattice_a_angstrom'],
                points=points, energy_reference='Same-composition final equilibrium total energy',
                energy_denominator='One three-atom cell', bracket_status=audit['status'],
                curve_fit_performed=False, interpolation_is_only_a_visual_guide=True,
                source_sha256=sources, matplotlib_version=matplotlib.__version__)
    label = {'mos2': r'MoS$_2$', 'mose2': r'MoSe$_2$'}[material]
    color = '#146b5b' if audit['status'] == 'passed' else '#9a5620'
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'svg.hashsalt': 'lattice-brackets-r1-' + material}):
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        fig.subplots_adjust(left=.15, right=.94, bottom=.22, top=.78)
        x = [p['delta_a_angstrom'] for p in points]
        y = [p['delta_energy_meV'] for p in points]
        ax.plot(x, y, '--', color=color, linewidth=1, alpha=.65)
        ax.scatter(x, y, s=65, color=color, zorder=3)
        ax.axhline(0, color='#aab4bf', linewidth=.7, zorder=0)
        for xx, yy in zip(x, y):
            ax.annotate(f'{yy:.6f}', (xx, yy), xytext=(0, 12), textcoords='offset points',
                        ha='center', fontsize=10, color='#263442')
        low, high = min(y), max(y)
        span = max(high - low, .01)
        ax.set_ylim(low - .15 * span, high + .35 * span)
        ax.set_xlim(-.0065, .0065)
        ax.set_xticks(x, ['−0.005', '0', '+0.005'])
        ax.set_xlabel(r'Change in lattice length $a$ (Å)')
        ax.set_ylabel('Energy change (meV / 3-atom cell)')
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e4e9ef', linewidth=.6)
        ax.set_axisbelow(True)
        fig.text(.15, .92, label + ' · Local lattice check', fontsize=17, fontweight='bold', color='#172b3a')
        fig.text(.15, .85, f"Equilibrium a = {data['equilibrium_a_angstrom']:.9f} Å · C09 · 120/480 Ry · 12×12×1 k", fontsize=10)
        fig.text(.15, .065, 'Three computed points; dashed segments are guides, with no fitted energy curve.\n'
                 'Internal z coordinates and cell height are fixed. Phonon stability is untested.',
                 fontsize=9, color='#4b5967', linespacing=1.5)
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(paths['png'], dpi=240, metadata={'Software': 'matplotlib'})
        fig.savefig(paths['svg'], metadata={'Date': None, 'Creator': 'matplotlib'})
        plt.close(fig)
    immutable(paths['json'], json.dumps(data, indent=2) + '\n')
    qa = dict(status='generated_from_verified_three_point_bracket', material=material,
              bracket_status=audit['status'], source_sha256=sources,
              output_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths.values()},
              curve_fit_performed=False, full_materials_baseline_complete=False,
              visual_inspection_recorded_separately=True)
    immutable(qa_path, json.dumps(qa, indent=2) + '\n')
    return paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('material', choices=('mos2', 'mose2'))
    args = parser.parse_args()
    print(json.dumps({key: str(path) for key, path in plot(args.material).items()}, indent=2))
