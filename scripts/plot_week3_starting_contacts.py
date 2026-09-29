"""Two orthogonal projections of the actual xTB starts, with measured contact distances."""
import json
import os
from pathlib import Path
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows=json.loads((ROOT/'reports/week2-xtb-dimers.json').read_text())['results']
entries=[('p36','Stacked start selected by xTB',None),
         ('p14','Phenol OH to carbonyl O',(17,25,29)),
         ('p06','Literature COOH to CAF N9',(2,18,32))]
colors={1:'#c1c5ca',6:'#485366',7:'#326ed0',8:'#d94040',9:'#2b9d76'}
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for col,(short,title,contact) in enumerate(entries):
    row=next(r for r in rows if r['candidate_id']=='diflunisal--caffeine-'+short)
    atoms=read(ROOT/row['run_dir']/'xtbopt.xyz');p=atoms.positions-atoms.positions.mean(axis=0)
    _,_,vt=np.linalg.svd(p[atoms.numbers>1]);p=p@vt.T
    extra=''
    if contact:
        donor,h,acceptor=contact;v=p[donor]-p[h];w=p[acceptor]-p[h]
        angle=np.degrees(np.arccos(np.clip(np.dot(v,w)/np.linalg.norm(v)/np.linalg.norm(w),-1,1)))
        extra=f'\nH...A {np.linalg.norm(w):.2f} Angstrom; {angle:.1f} deg'
    for view,dims in enumerate([(0,1),(0,2)]):
        ax=axes[view,col];q=p[:,dims]
        for i,j,_ in row['candidate']['graph']['bonds']:ax.plot(*q[[i,j]].T,c='#8794a3',lw=1.5,zorder=1)
        for z in sorted(set(atoms.numbers)):
            ids=np.flatnonzero(atoms.numbers==z);ax.scatter(*q[ids].T,s=22 if z==1 else 78,c=colors[z],edgecolors='white',lw=.6,zorder=3)
        if contact:
            donor,h,acceptor=contact;ax.plot(*q[[h,acceptor]].T,'--',color='#e69411',lw=2.2,zorder=4)
            if view==0:
                for idx,offset in [(donor,(8,13)),(acceptor,(8,-15))]:
                    symbol='N' if atoms.numbers[idx]==7 else 'O'
                    ax.annotate(f'{symbol}{idx}',q[idx],xytext=offset,textcoords='offset points',fontsize=9,
                                bbox={'facecolor':'white','edgecolor':'none','alpha':.85},zorder=5)
        ax.set_aspect('equal');ax.margins(.14);ax.axis('off')
        if view==0:ax.set_title(f'{short}: {title}'+extra,fontsize=11,pad=16)
        else:ax.set_title('Orthogonal view',fontsize=9,color='#596675')
fig.suptitle('DFT starting contacts: two views of each molecular pair\nAll shown coordinates are xTB outputs, not measured crystal structures',fontsize=13)
out=ROOT/'reports/figures/week3';out.mkdir(parents=True,exist_ok=True)
fig.savefig(out/'starting-contact-review.png',dpi=180)
