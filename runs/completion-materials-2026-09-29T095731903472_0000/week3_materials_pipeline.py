"""Bounded adaptive monolayer/heterobilayer reference; every engine call is supplied by the controller."""
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from scipy.optimize import brentq

from week3_qe import ROOT, job, load, json_file
from research_run import save, sha
from week2_qe_results import load_job as load_week2

REGISTRIES = [(0.,0.), (1/3,2/3), (2/3,1/3)]


def centered(atoms, vacuum):
    atoms=atoms.copy();atoms.center(vacuum=vacuum/2,axis=2);return atoms


def bilayer(bottom, top, rid, distance, vacuum):
    assert np.allclose(bottom.cell[:2],top.cell[:2],atol=1e-9,rtol=0)
    b=bottom.positions.copy();t=top.positions.copy()
    b[:,2]-=b[0,2];t[:,2]-=t[0,2]
    thickness=(np.ptp(b[:,2])+np.ptp(t[:,2]))/2
    cell=bottom.cell.array.copy();cell[2]=[0,0,distance+thickness+vacuum]
    b[:,2]+=cell[2,2]/2-distance/2;t[:,2]+=cell[2,2]/2+distance/2
    r=REGISTRIES[rid];t+=r[0]*cell[0]+r[1]*cell[1]
    atoms=Atoms(bottom.get_chemical_symbols()+top.get_chemical_symbols(),
                positions=np.vstack([b,t]),cell=cell,pbc=True)
    assert abs(atoms.positions[3,2]-atoms.positions[0,2]-distance)<1e-9
    assert abs(cell[2,2]-np.ptp(atoms.positions[:,2])-vacuum)<1e-8
    return atoms


def summarize(r):
    return dict(job_id=r['job_id'],energy_eV=r['energy_eV'],max_force_eV_angstrom=r['max_force_eV_angstrom'],
                lattice_a_angstrom=float(r['atoms'].cell.lengths()[0]),vacuum_gap_angstrom=r['vacuum_gap_angstrom'],
                number_of_electrons=r['nelec'],xml_sha256=r['xml_sha256'])


def execute(run_config, stage):
    report=ROOT/'reports/week3-materials.json'
    result=dict(status='running',scope='C09 frozen-layer interaction reference; no phonon or alloy validation',stages=[])
    save(report,result)

    def run(path):
        run_config(path)
        cfg=json.loads((ROOT/path).read_text());return load(cfg['job_id'])

    stage('semiconductor_lattice_equilibria')
    equilibria={name:run(f'experiments/week3/materials/w3-{name}-lattice-vc.json') for name in ['mos2','mose2']}
    result['equilibria']={name:summarize(r) for name,r in equilibria.items()};save(report,result)

    stage('semiconductor_lattice_and_precision_checks')
    # Local energy brackets at the relaxed z coordinates verify the cell stationary point.
    brackets={};convergence={};settings={'selected':(80,6,20),'middle':(100,9,25),'reference':(120,12,30)}
    for name,eq in equilibria.items():
        a0=eq['atoms'].cell.lengths()[0];rows=[]
        for direction in [-1,1]:
            atoms=eq['atoms'].copy();cell=atoms.cell.array.copy();cell[:2]*=(a0+direction*.005)/a0
            atoms.set_cell(cell,scale_atoms=True)
            r=run(job(f'w3-{name}-lattice-probe-{"minus" if direction<0 else "plus"}',atoms,
                      role='fixed_internal_coordinates_lattice_bracket',delta_a_angstrom=direction*.005))
            delta=r['energy_eV']-eq['energy_eV']
            rows.append(dict(delta_a_angstrom=direction*.005,delta_energy_eV=delta,job_id=r['job_id']))
        assert all(row['delta_energy_eV']>0 for row in rows), 'Cell energy bracket requires review'
        brackets[name]=rows
        computed={label:run(job(f'w3-{name}-precision-{label}',centered(eq['atoms'],v),cut=c,k=k,
                               role='combined_numerical_sensitivity',comparison_setting=label))
                  for label,(c,k,v) in settings.items()}
        anchor=computed['reference'];convergence[name]={}
        for label,r in computed.items():
            de=abs(r['energy_eV']-anchor['energy_eV'])*1000/len(r['atoms'])
            df=float(np.linalg.norm(r['forces']-anchor['forces'],axis=1).max())
            convergence[name][label]=dict(delta_energy_meV_atom=de,max_delta_force_eV_angstrom=df,
                                           passed=de<=1 and df<=.005,job_id=r['job_id'])
    result.update(lattice_brackets=brackets,monolayer_combined_sensitivity=convergence);save(report,result)
    selected=next(label for label in ['selected','middle','reference'] if all(convergence[n][label]['passed'] for n in equilibria))
    cut,k,vac=settings[selected]
    result['setting_selection']=dict(label=selected,ecutwfc_Ry=cut,k=k,vacuum_angstrom=vac,
                                    scope='Relative to 120/12/30 at both relaxed monolayers; heterobilayer difference checked separately')
    a_common=float(np.mean([r['atoms'].cell.lengths()[0] for r in equilibria.values()]))
    result['common_lattice_angstrom']=a_common
    result['strain_fraction']={n:a_common/r['atoms'].cell.lengths()[0]-1 for n,r in equilibria.items()};save(report,result)

    stage('semiconductor_common_lattice_internal_relaxation')
    layers={}
    for name,eq in equilibria.items():
        atoms=eq['atoms'].copy();cell=atoms.cell.array.copy();cell[:2]*=a_common/atoms.cell.lengths()[0]
        atoms.set_cell(cell,scale_atoms=True);atoms=centered(atoms,vac)
        r=run(job(f'w3-{name}-common-lattice-relax',atoms,calculation='relax',cut=cut,k=k,
                  role='strained_monolayer_internal_relaxation',common_lattice_angstrom=a_common))
        layers[name]=r['atoms'];result.setdefault('common_lattice_layers',{})[name]=summarize(r)
    save(report,result)

    stage('semiconductor_heterobilayer_distance_search')
    comparisons=[];cache={}

    def compare(tag,bottom,top,rid,distance,setting):
        c,kk,v=setting;distance=round(float(distance),6)
        key=(tag,rid,distance,c,kk,v)
        if key in cache:return cache[key]
        atoms=bilayer(bottom,top,rid,distance,v);number=len(comparisons)
        values={}
        for part,geometry in [('pair',atoms),('bottom',atoms[:3]),('top',atoms[3:])]:
            cfg=job(f'w3-{tag}-g{number:03d}-{part}',geometry,cut=c,k=kk,role=part,
                    registry_id=rid,registry_fraction=list(REGISTRIES[rid]),Mo_plane_distance_angstrom=distance,
                    bilayer_vacuum_gap_angstrom=v,orientation='R-type; no twist',interaction_group=tag)
            values[part]=run(cfg)
        pair=values['pair']
        # Verify that reference geometries really are the extracted frozen layers.
        for part,indices in [('bottom',slice(0,3)),('top',slice(3,6))]:
            assert np.allclose(values[part]['atoms'].positions,pair['atoms'].positions[indices],rtol=0,atol=1e-7)
            assert np.allclose(values[part]['atoms'].cell,pair['atoms'].cell,rtol=0,atol=1e-7)
        gamma=(pair['energy_eV']-values['bottom']['energy_eV']-values['top']['energy_eV'])*1000/pair['area_angstrom2']
        fz=float(pair['forces'][3:,2].sum())
        assert abs(float(pair['forces'][:,2].sum()))<.002, 'Net z force requires review'
        record=dict(tag=tag,registry_id=rid,Mo_plane_distance_angstrom=distance,setting=[c,kk,v],
                    gamma_meV_angstrom2=gamma,upper_layer_force_z_eV_angstrom=fz,
                    area_angstrom2=pair['area_angstrom2'],jobs={part:r['job_id'] for part,r in values.items()})
        comparisons.append(record);cache[key]=record
        result['comparisons']=comparisons;save(report,result)
        return record

    bottom,top=layers['mos2'],layers['mose2'];chosen=(cut,k,vac)
    distances=[6.0,6.3,6.6,6.9,7.2,7.5,7.8,8.1]
    scan=[];force_bracket=None
    for distance in distances:
        row=compare('hetero',bottom,top,0,distance,chosen);scan.append(row)
        if len(scan)>=2 and scan[-2]['upper_layer_force_z_eV_angstrom']>0 and row['upper_layer_force_z_eV_angstrom']<0:
            force_bracket=[scan[-2]['Mo_plane_distance_angstrom'],distance];break
    assert force_bracket is not None, 'No D0 bracket inside reviewed 6.0-8.1 Angstrom search'
    d0=float(brentq(lambda d:compare('hetero',bottom,top,0,d,chosen)['upper_layer_force_z_eV_angstrom'],
                   *force_bracket,xtol=.002,rtol=1e-8,maxiter=10))
    center=compare('hetero',bottom,top,0,d0,chosen);d0=center['Mo_plane_distance_angstrom']
    assert abs(center['upper_layer_force_z_eV_angstrom'])<=.01, 'D0 force residual requires review'
    result['D0']=dict(Mo_plane_distance_angstrom=d0,force_bracket_angstrom=force_bracket,
                     force_residual_eV_angstrom=center['upper_layer_force_z_eV_angstrom'],
                     scope='Frozen internally relaxed layers at the common in-plane lattice')
    save(report,result)

    stage('semiconductor_heterobilayer_precision_check')
    strict=(120,12,30) if chosen!=(120,12,30) else (140,15,35)
    strict_row=compare('hetero-strict',bottom,top,0,d0,strict)
    sensitivity=abs(strict_row['gamma_meV_angstrom2']-center['gamma_meV_angstrom2'])
    result['heterobilayer_numerical_sensitivity']=dict(delta_gamma_meV_angstrom2=sensitivity,
                                                     threshold_meV_angstrom2=.1,passed=sensitivity<=.1,
                                                     settings=[list(chosen),list(strict)],registry_id=0,distance_angstrom=d0)
    save(report,result)
    assert sensitivity<=.1, 'Representative heterobilayer precision exceeds 0.1 meV/Angstrom^2'

    stage('semiconductor_heterobilayer_nine_point_grid')
    grid=[compare('hetero',bottom,top,rid,d,chosen) for rid in range(3) for d in [d0-.25,d0,d0+.25]]
    assert len(grid)==9
    mid=next(r for r in grid if r['registry_id']==0 and abs(r['Mo_plane_distance_angstrom']-d0)<1e-5)
    sides=[r for r in grid if r['registry_id']==0 and abs(r['Mo_plane_distance_angstrom']-d0)>.1]
    assert all(r['gamma_meV_angstrom2']>mid['gamma_meV_angstrom2'] for r in sides)
    result['nine_point_grid']=grid;save(report,result)

    stage('semiconductor_pure_MoS2_distance_bracket_completion')
    old=json.loads((ROOT/'experiments/week2/bilayer-comparisons.json').read_text())['configurations']
    pure=[]
    for rec in old:
        if rec['setting']!='selected' or rec['registry_id']!=0:continue
        vals={part:load_week2(j) for part,j in rec['jobs'].items()}
        p=vals['pair'];gamma=(p['energy_eV']-vals['bottom']['energy_eV']-vals['top']['energy_eV'])*1000/p['area_angstrom2']
        pure.append(dict(Mo_plane_distance_angstrom=rec['Mo_plane_distance_angstrom'],gamma_meV_angstrom2=gamma,
                         upper_layer_force_z_eV_angstrom=float(p['forces'][3:,2].sum()),jobs=rec['jobs'],source_week=2))
    base=read(ROOT/'structures/optimized/week1/MoS2_fixed_cell_relaxed.cif')
    for d in [6.9,7.2,7.5,7.8]:
        r=compare('pure-mos2',base,base,0,d,(80,6,20));pure.append({**r,'source_week':3})
        energies=[r['gamma_meV_angstrom2'] for r in pure];index=int(np.argmin(energies))
        if 0<index<len(energies)-1 and energies[index]<energies[index-1] and energies[index]<energies[index+1]:break
    index=int(np.argmin([r['gamma_meV_angstrom2'] for r in pure]))
    assert 0<index<len(pure)-1, 'Pure MoS2 sampled minimum remains at boundary'
    result.update(pure_mos2_distance_points=pure,pure_mos2_sample_minimum=pure[index],
                  pure_mos2_scope='Bracketed sampled minimum at original a=3.18 Angstrom; not reused as heterobilayer D0',
                  status='passed',alloy_24_atom_comparison_started=False)
    save(report,result)
    return result
