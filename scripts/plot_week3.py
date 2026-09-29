"""Scientific figures from completed records only; no inferred error bars."""
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    out=ROOT/'reports/figures/week3';out.mkdir(parents=True,exist_ok=True)
    p=ROOT/'reports/week3-drug-comparison.json'
    if p.exists():
        d=json.loads(p.read_text());a=d['scores']['primary_only'];b=d['scores']['all_available_poses']
        if all(row['X_kcal_mol'] is not None and row['B_kcal_mol'] is not None for row in a):
            fig,axes=plt.subplots(1,2,figsize=(10,4.3),layout='constrained');x=np.arange(4)
            for ax,key,title in zip(axes,['B_kcal_mol','X_kcal_mol'],['Heterodimer binding B','Self-association contrast X']):
                ax.bar(x-.17,[r[key] for r in a],width=.34,label='Primary starts',color='#3875a3')
                ax.bar(x+.17,[r[key] for r in b],width=.34,label='All available poses',color='#d58d37')
                ax.set_xticks(x,[r['coformer'].capitalize() for r in a],rotation=20)
                ax.axhline(0,color='black',lw=.7);ax.set_ylabel('Electronic energy (kcal/mol)');ax.set_title(title);ax.legend(fontsize=8)
                ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
            fig.suptitle('Vacuum r2SCAN-3c; unequal added-pose coverage\nIsolated dimers, not crystal formation free energies',fontsize=11)
            fig.savefig(out/'drug-controls.png',dpi=180);plt.close(fig)
    p=ROOT/'reports/week3-materials.json'
    if p.exists():
        d=json.loads(p.read_text());grid=d.get('nine_point_grid',[])
        if len(grid)==9:
            fig,ax=plt.subplots(figsize=(7,4.8),layout='constrained')
            for rid in range(3):
                rows=sorted((r for r in grid if r['registry_id']==rid),key=lambda r:r['Mo_plane_distance_angstrom'])
                ax.plot([r['Mo_plane_distance_angstrom'] for r in rows],[r['gamma_meV_angstrom2'] for r in rows],'-o',label=f'r{rid}')
            ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.2);ax.legend(title='Registry')
            ax.set_xlabel('Mo-plane separation (Angstrom)');ax.set_ylabel('Layer interaction (meV / Angstrom$^2$)')
            ax.set_title('MoS2 / MoSe2: fixed layers at a common in-plane lattice\nC09, no twist; three sampled registries')
            fig.savefig(out/'heterobilayer-reference.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
