"""Output guards with archived r2SCAN and synthetic SMD text; no engine runs."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import evaluate_completion_orca_compatibility as check


class CompatibilityEvidenceTests(unittest.TestCase):
    def smd(self, final=-9.999, solvent='ETHANOL'):
        # Intentionally artificial numbers, not computed or published energies.
        return f'''****END OF INPUT****
CPCM SOLVATION MODEL
Solvent ... {solvent}
SMD-CDS solvent descriptors
SMD CDS free energy correction energy : 1.882528422 Kcal/mol
Total Energy after SMD CDS correction = -10.000000000 Eh
TOTAL SCF ENERGY
Total Energy : -10.000000000 Eh
CPCM Dielectric : -0.010000000 Eh
SMD CDS (Gcds) : 0.003000000 Eh
Dispersion correction -0.001000000
gCP correction 0.002000000
FINAL SINGLE POINT ENERGY {final:.12f}
'''

    def test_smd_total_already_contains_cds_and_composite_terms_follow(self):
        result = check.smd_energy_evidence(self.smd())
        self.assertAlmostEqual(result['final_energy_Hartree'], -9.999)
        self.assertFalse(result['cds_added_again'])
        self.assertFalse(result['molecular_formation_free_energy'])

    def test_double_cds_addition_and_missing_component_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Composite energy sum'):
            check.smd_energy_evidence(self.smd(final=-9.996))
        with self.assertRaisesRegex(ValueError, 'Missing or nonfinite'):
            check.smd_energy_evidence(self.smd().replace('gCP correction', 'missing correction'))
        with self.assertRaisesRegex(ValueError, 'CDS unit conversion'):
            check.smd_energy_evidence(self.smd().replace('1.882528422', '2.882528422'))

    def test_actual_solvent_header_and_last_step_are_used(self):
        with self.assertRaisesRegex(ValueError, 'not ethanol'):
            check.smd_energy_evidence(self.smd(solvent='WATER'))
        # Earlier valid output must not hide a final inconsistent energy.
        tail = self.smd(final=-9.996).split('****END OF INPUT****', 1)[1]
        with self.assertRaisesRegex(ValueError, 'Composite energy sum'):
            check.smd_energy_evidence(self.smd() + tail)

    def test_archived_r2scan_method_is_read_without_accepting_it_as_wb97(self):
        path = ROOT / 'runs/20260910T092739076636Z-week1-water-opt-freq/output.log'
        before = check.sha(path)
        result = check.method_evidence(path.read_text(), 'r2SCAN-3c')
        self.assertEqual(result['orbital_basis'], 'def2-mTZVPP')
        self.assertTrue(result['detailed_parameter_and_ecp_review_required'])
        with self.assertRaisesRegex(ValueError, 'runtime composite-method'):
            check.method_evidence(path.read_text(), 'wB97X-3c')
        self.assertEqual(check.sha(path), before)

    def test_method_credits_and_echoed_input_are_not_runtime_evidence(self):
        text = '''Author: wB97X-3c, vDZP, DFTD4
Your calculation utilizes the basis: vDZP
The WB97X3C composite method is recognized
****END OF INPUT****
SCF CONVERGED
'''
        with self.assertRaisesRegex(ValueError, 'runtime composite-method'):
            check.method_evidence(text, 'wB97X-3c')


if __name__ == '__main__':
    unittest.main()
