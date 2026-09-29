"""Hand-solvable controls and invariances, not evidence of physical results."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import completion_drug_contrasts as c


class DrugContrasts(unittest.TestCase):
    def setUp(self):
        # Analytic synthetic energies; no claim that these describe molecules.
        k = c.HARTREE_KCAL
        self.energy = {'diflunisal': -100., 'diflunisal--diflunisal': -200. - 10/k}
        for i, b in enumerate(c.COFORMERS):
            monomer = -10. * (i+1)
            self.energy[b] = monomer
            self.energy[f'{b}--{b}'] = 2*monomer - (2+2*i)/k
            self.energy[f'diflunisal--{b}'] = -100. + monomer - (8-i)/k
        self.context = dict(method='synthetic', environment='synthetic', numerical_protocol_id='fixture', energy_quantity='synthetic_electronic')
        self.rows = [dict(record_id=s, species=s, selection_role='primary', energy_Hartree=e, comparison_context=self.context)
                     for s, e in self.energy.items()]

    def test_analytic_binding_and_self_competition(self):
        scores = c.panel_scores(self.energy)
        self.assertEqual(len(scores), 4)
        for i, row in enumerate(scores):
            self.assertAlmostEqual(row['B']['kcal_mol'], -8+i, places=9)
            self.assertAlmostEqual(row['X']['kcal_mol'], -2+2*i, places=9)

    def test_independent_species_reference_shifts_cancel(self):
        shifted = {}
        offsets = {'diflunisal': 13.12345, **{b: 5.1214*(i+1) for i,b in enumerate(c.COFORMERS)}}
        for species, energy in self.energy.items():
            shifted[species] = energy + sum(offsets[b] for b in species.split('--'))
        for before, after in zip(c.panel_scores(self.energy), c.panel_scores(shifted)):
            for key in ('B','X'):
                self.assertAlmostEqual(before[key]['kcal_mol'], after[key]['kcal_mol'], places=8)

    def test_missing_self_control_is_unknown_not_zero_or_negative_prediction(self):
        del self.energy['caffeine--caffeine']
        row = next(r for r in c.panel_scores(self.energy) if r['coformer']=='caffeine')
        self.assertIsNotNone(row['B']['kcal_mol'])
        self.assertIsNone(row['X']['kcal_mol'])
        self.assertEqual(row['X']['missing_species'], ['caffeine--caffeine'])

    def test_tie_boundaries_and_missing_differences(self):
        result = c.compare_observed_groups(dict(adenine=-2., caffeine=None, theobromine=-1., xanthine=-4.),1.)
        self.assertEqual(result['counts'],dict(consistent=0,reversed=1,tie=1,unavailable=2))
        exact = c.compare_observed_groups(dict(adenine=0.,caffeine=0.,theobromine=0.,xanthine=1.),0.)
        self.assertEqual(exact['counts'],dict(consistent=2,reversed=0,tie=2,unavailable=0))

    def test_primary_scope_does_not_absorb_added_pose(self):
        alternative = copy.deepcopy(next(r for r in self.rows if r['species']=='diflunisal--caffeine'))
        alternative.update(record_id='caf-alt', selection_role='alternative', energy_Hartree=alternative['energy_Hartree']-3/c.HARTREE_KCAL)
        result = c.analyze_uniform_records([*self.rows,alternative], self.context)['scopes']
        self.assertEqual(result['primary_only']['selected_record_ids']['diflunisal--caffeine'],'diflunisal--caffeine')
        self.assertEqual(result['all_selected_poses']['selected_record_ids']['diflunisal--caffeine'],'caf-alt')
        self.assertEqual(len(result['primary_only']['comparisons']['B']),4)

    def test_context_mixing_is_rejected(self):
        for key in c.CONTEXT_KEYS:
            rows=copy.deepcopy(self.rows);rows[-1]['comparison_context']=dict(self.context,**{key:'different'})
            with self.subTest(key=key), self.assertRaises(ValueError):
                c.analyze_uniform_records(rows,self.context)

    def test_nonfinite_and_boolean_energies_rejected(self):
        for value in [float('nan'), float('inf'), -float('inf'), True]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                c.panel_scores({**self.energy,'adenine':value})

    def test_shared_dif_self_reference_shifts_all_X_equally(self):
        before=c.panel_scores(self.energy)
        self.energy['diflunisal--diflunisal'] += .5/c.HARTREE_KCAL
        after=c.panel_scores(self.energy)
        for a,b in zip(before,after):
            self.assertAlmostEqual(b['X']['kcal_mol']-a['X']['kcal_mol'],-.25,places=9)
            self.assertEqual(a['B'],b['B'])
        # Absolute X needs AA, though its common additive term cancels in rankings.
        for width in c.TIE_WIDTHS:
            a=c.compare_observed_groups({r['coformer']:r['X']['kcal_mol'] for r in before},width)
            b=c.compare_observed_groups({r['coformer']:r['X']['kcal_mol'] for r in after},width)
            self.assertEqual(a['counts'],b['counts'])

    def test_real_archived_fixed_geometry_diagnostic_binding(self):
        data=json.loads((ROOT/'reports/completion/precision-diagnostic.json').read_text())
        values={}
        for grid in ('g2','g3','dense'):
            rows={r['species']:r['energy_Hartree'] for r in data['rows'] if r['grid_setting']==grid}
            values[grid]=c.contrast(rows,{'dif2':1.,'dif':-2.})['kcal_mol']
        # Independently recorded prior diagnostic result; this is one geometry,
        # not a completed four-coformer calculation or a minima certificate.
        self.assertAlmostEqual(values['dense']-values['g3'],.004744127722,places=7)


if __name__=='__main__':unittest.main()
