"""Static scientific figures from saved, verified week-two data only."""
import json
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=ROOT/'reports/figures/week2';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
files=[]
def save(fig,name):
    fig.savefig(OUT/f'{name}.png',bbox_inches='tight');fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
    files.append(str((OUT/f'{name}.png').relative_to(ROOT)))

scan=json.loads((ROOT/'reports/week2-monolayer-convergence.json').read_text())
if scan['completed']==7:
    fig,axes=plt.subplots(2,3,figsize=(11.5,6.5))
    for column,(field,label,condition) in enumerate(zip(['cutoff_Ry','k','vacuum_angstrom'],['Cutoff (Ry)','k x k x 1 mesh','Vacuum gap (Angstrom)'],
                                      [lambda r:r['k']==12 and r['vacuum_angstrom']==30,
                                       lambda r:r['cutoff_Ry']==120 and r['vacuum_angstrom']==30,
                                       lambda r:r['cutoff_Ry']==120 and r['k']==12])):
        ax=axes[0,column];force_ax=axes[1,column]
        rows=sorted([r for r in scan['rows'] if condition(r)],key=lambda r:r[field])
        ax.axhspan(-1,1,color='#daf0e5',label='Pilot tolerance: +/-1 meV/atom')
        ax.plot([r[field] for r in rows],[r['energy_delta_reference_meV_atom'] for r in rows],'o-',color='#176b88')
        ax.axhline(0,color='#888888',lw=.7);ax.set_xlabel(label);ax.set_xticks([r[field] for r in rows]);ax.grid(alpha=.2)
        force_ax.axhspan(0,.005,color='#daf0e5')
        force_ax.plot([r[field] for r in rows],[r['max_force_delta_reference_eV_angstrom'] for r in rows],'o-',color='#176b88')
        force_ax.set_xlabel(label);force_ax.set_xticks([r[field] for r in rows]);force_ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Energy change\n(meV/atom)')
    axes[1,0].set_ylabel('Max force change\n(eV/Angstrom)')
    fig.suptitle('MoS2 monolayer: fixed a = 3.18 Angstrom, C09 / PseudoDojo PBE NC',fontsize=12)
    fig.text(.5,.005,'Reference: 120 Ry, 12 x 12 x 1, 30 Angstrom. Green band is an operational tolerance, not a statistical interval.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.05,1,.93]);save(fig,'monolayer-convergence')

path=ROOT/'reports/week2-bilayer-analysis.json'
if path.exists():
    data=json.loads(path.read_text());rows=data['rows']
    if rows:
        fig,axes=plt.subplots(1,2,figsize=(10,4))
        curve=sorted([r for r in rows if r['setting']=='selected' and r['registry_id']==0],key=lambda r:r['Mo_plane_distance_angstrom'])
        axes[0].plot([r['Mo_plane_distance_angstrom'] for r in curve],[r['gamma_meV_angstrom2'] for r in curve],'o-',color='#176b88')
        axes[0].set(xlabel='Mo-plane separation (Angstrom)',ylabel='Interaction energy / area (meV/Angstrom^2)',title='r0: rigid-layer distance samples')
        reg=sorted([r for r in rows if r['setting']=='selected' and r['Mo_plane_distance_angstrom']==6.3],key=lambda r:r['registry_id'])
        axes[1].bar([f"r{r['registry_id']}" for r in reg],[r['gamma_meV_angstrom2'] for r in reg],color=['#176b88','#62a58c','#62a58c'][:len(reg)])
        axes[1].set(xlabel='Registry at D = 6.3 Angstrom',title='Same orientation, fixed geometry')
        for ax in axes:ax.axhline(0,lw=.7,color='#888888');ax.grid(axis='y',alpha=.2)
        fig.suptitle('Pure MoS2 / MoS2 control; each energy uses matched frozen-layer references',fontsize=11)
        fig.text(.5,.01,'Negative values mean lower energy than matched separated layers; they do not establish the force direction or equilibrium.',ha='center',fontsize=8)
        fig.tight_layout(rect=[0,.06,1,.92]);save(fig,'bilayer-controls')

path=ROOT/'reports/week2-xtb-analysis.json'
if path.exists():
    data=json.loads(path.read_text());names=list(data['species']);values=[data['species'][n] for n in names]
    labels={ 'diflunisal':'DIF','caffeine':'CAF','adenine':'ADE','theobromine':'THEO','xanthine':'XAN'}
    x=np.arange(len(names));fig,ax=plt.subplots(figsize=(12,4.5))
    ax.bar(x-.24,[s['started'] for s in values],.24,label='Starts',color='#b5c4d0')
    ax.bar(x,[s['geometry_usable'] for s in values],.24,label='Converged and geometry-usable',color='#176b88')
    ax.bar(x+.24,[s['cluster_count'] for s in values],.24,label='Conservative geometry clusters',color='#62a58c')
    ax.set_xticks(x,['-'.join(labels.get(p,p) for p in n.split('--')) for n in names],rotation=40,ha='right')
    ax.set_ylabel('Number of structures');ax.set_title('GFN2-xTB screening: counts, not verified minima or prediction accuracy')
    ax.legend(fontsize=9);ax.grid(axis='y',alpha=.2);fig.tight_layout();save(fig,'xtb-screening-counts')

(ROOT/'reports/week2-figure-index.json').write_text(json.dumps({'figures':files,'data_source':'saved week2 JSON reports; no simulated results'},indent=2)+'\n')
print(json.dumps({'figures':files}))
