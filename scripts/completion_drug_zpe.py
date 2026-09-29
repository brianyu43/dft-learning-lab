"""Read harmonic zero-point corrections without changing archived results.

The numerical extraction is reusable after a separate minimum/provenance
audit. It does not certify a method, minimum, grid or formation free energy.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import re

import numpy as np
from ase.io import read

from completion_drug_contrasts import HARTREE_KCAL
from prepare_completion_drug_lineage import build as revalidate_lineage
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_frequency_resume import inspect_hessian

# CODATA 2022 hartree-inverse meter relationship, divided by 100.
# Matches the installed scipy.constants value; freeze it for this analysis.
WAVENUMBER_CM1_PER_HARTREE = 219474.63136314
PRINT_COMPARISON_TOLERANCE_HARTREE = 1e-8  # ORCA prints the summary to 8 decimal places.
SOURCES = {
    'frequency_scaling': 'https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html',
    'thermochemistry': 'https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html',
    'conversion': 'https://docs.scipy.org/doc/scipy/reference/constants.html',
}
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def scalar(pattern, text):
    found = re.findall(pattern, text, re.M)
    require(len(found) == 1, 'Missing or repeated scalar: ' + pattern)
    value = float(found[0].replace('D', 'E').replace('d', 'e'))
    require(math.isfinite(value), 'Nonfinite scalar')
    return value


def harmonic_zpe(stored_frequencies_cm1, *, applied_scale_factor):
    """Frequencies from ORCA .hess are already scaled: do not scale twice."""
    freq = np.asarray(stored_frequencies_cm1, dtype=float)
    require(freq.ndim == 1 and freq.size > 0 and np.isfinite(freq).all() and np.all(freq > 0),
            'All internal frequencies must be finite and strictly positive')
    require(not isinstance(applied_scale_factor, bool) and math.isfinite(applied_scale_factor)
            and applied_scale_factor > 0, 'Invalid applied scaling factor')
    value = .5 * math.fsum(freq.tolist()) / WAVENUMBER_CM1_PER_HARTREE
    return dict(zpe_Hartree=value, zpe_kcal_mol=value * HARTREE_KCAL,
                applied_frequency_scale_factor=applied_scale_factor,
                scale_applied_again=False, internal_frequency_count=len(freq),
                minimum_internal_frequency_cm1=float(freq.min()),
                modes_below_50_cm1=int(np.count_nonzero(freq < 50)),
                modes_below_100_cm1=int(np.count_nonzero(freq < 100)))


def extract(log_path, hessian_path, atoms, electronic_energy_Hartree, *, expected_scale_factor=1.):
    """Check raw frequency/thermochemistry agreement for a nonlinear molecule.

    Caller must separately establish the minimum and preserve source hashes.
    Multiple thermochemistry blocks or discarded positive modes require a new
    reviewed extraction, rather than silently selecting a convenient block.
    """
    require(math.isfinite(electronic_energy_Hartree), 'Nonfinite electronic energy')
    log_path, hessian_path = Path(log_path), Path(hessian_path)
    log = log_path.read_text(); lines = hessian_path.read_text().splitlines()
    require('ORCA TERMINATED NORMALLY' in log and 'SCF CONVERGED' in log
            and 'SCF NOT CONVERGED' not in log, 'Raw normal completion is missing')
    require(re.search(r'Program Version\s+6\.1\.1\b', log), 'Unreviewed ORCA version')
    require(log.count('VIBRATIONAL FREQUENCIES') == 1, 'Expected one complete frequency analysis')
    block = log.split('VIBRATIONAL FREQUENCIES', 1)[1].split('NORMAL MODES', 1)[0]
    printed_scale = scalar(r'Scaling factor for frequencies\s*=\s*(' + NUMBER + ')', block)
    require(lines.count('$frequency_scale_factor') == 1, 'Missing or repeated Hessian scale factor')
    stored_scale = float(lines[lines.index('$frequency_scale_factor') + 1])
    require(abs(stored_scale - expected_scale_factor) < 1e-10
            and abs(printed_scale - stored_scale) < 1e-10, 'Frequency scale differs from the specified policy')
    modes = [(int(i), float(f)) for i, f in re.findall(r'^\s*(\d+):\s+([-+]?\d+\.\d+)\s+cm', block, re.M)]
    require([i for i, _ in modes] == list(range(3 * len(atoms)))
            and all(abs(f) <= .001 for _, f in modes[:6]), 'Incomplete or nonstandard nonlinear mode list')
    require(lines.count('$vibrational_frequencies') == 1, 'Missing or repeated Hessian frequency table')
    start = lines.index('$vibrational_frequencies'); n = int(lines[start + 1])
    require(n == 3 * len(atoms)
            and [int(x.split()[0]) for x in lines[start + 2:start + 2 + n]] == list(range(n)),
            'Hessian mode indices are incomplete or duplicated')
    frequencies = inspect_hessian(hessian_path, atoms, [f for _, f in modes[6:]])
    result = harmonic_zpe(frequencies, applied_scale_factor=stored_scale)
    require(log.count('THERMOCHEMISTRY AT') == 1, 'Expected one thermochemistry block')
    thermo = log.split('THERMOCHEMISTRY AT', 1)[1]
    printed_zpe = scalar(r'^Zero point energy\s+\.\.\.\s+(' + NUMBER + r')\s+Eh', thermo)
    printed_nonthermal = scalar(r'^Non-thermal \(ZPE\) correction\s+(' + NUMBER + r')\s+Eh', thermo)
    printed_electronic = scalar(r'^Electronic energy\s+\.\.\.\s+(' + NUMBER + r')\s+Eh', thermo)
    cutoff = scalar(r'^Cut-Off Frequency\s+\.\.\.\s+(' + NUMBER + r')\s+cm', thermo)
    require(result['minimum_internal_frequency_cm1'] > cutoff,
            'A positive mode lies below the thermochemistry cutoff; separate review required')
    require(abs(printed_electronic - electronic_energy_Hartree) <= PRINT_COMPARISON_TOLERANCE_HARTREE,
            'Thermochemistry electronic energy differs from the verified run')
    for value in (printed_zpe, printed_nonthermal):
        require(abs(value - result['zpe_Hartree']) <= PRINT_COMPARISON_TOLERANCE_HARTREE,
                'Frequency sum and printed ZPE disagree')
    start = lines.index('$atoms'); count = int(lines[start + 1])
    masses = [float(x.split()[1]) for x in lines[start + 2:start + 2 + count]]
    require(count == len(atoms) and np.isfinite(masses).all() and all(m > 0 for m in masses), 'Invalid Hessian masses')
    qrrho = re.findall(r'^Quasi RRHO\s+\.\.\.\s+(True|False)', thermo, re.M)
    require(len(qrrho) == 1, 'Missing thermochemistry convention')
    result.update(electronic_energy_Hartree=electronic_energy_Hartree,
        electronic_plus_harmonic_zpe_Hartree=math.fsum([electronic_energy_Hartree, result['zpe_Hartree']]),
        printed_zpe_Hartree=printed_zpe, sum_minus_printed_zpe_Hartree=result['zpe_Hartree'] - printed_zpe,
        printed_energy_comparison_tolerance_Hartree=PRINT_COMPARISON_TOLERANCE_HARTREE,
        hessian_masses_amu=masses, hessian_atomic_numbers=atoms.numbers.tolist(),
        thermochemistry_metadata=dict(temperature_K=scalar(r'^Temperature\s+\.\.\.\s+(' + NUMBER + r')\s+K', thermo),
            pressure_atm=scalar(r'^Pressure\s+\.\.\.\s+(' + NUMBER + r')\s+atm', thermo),
            quasi_rrho=qrrho[0] == 'True', cutoff_frequency_cm1=cutoff),
        energy_quantity='harmonic_zero_point_correction', thermal_or_entropy_terms_added=False,
        frequency_grid_convergence_verified=False, source_sha256={str(p): sha(p) for p in (log_path, hessian_path)})
    return result


def archived():
    lineage_path = ROOT / 'experiments/completion/drug-start-lineage-r1.json'
    lineage = json.loads(lineage_path.read_text())
    require(revalidate_lineage() == lineage, 'Original lineage or minimum evidence changed')
    rows = []
    evidence = dict(lineage['input_hashes'])
    evidence.update({str(lineage_path.relative_to(ROOT)): sha(lineage_path),
                     str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))})
    for row in lineage['rows']:
        minima = [a for a in row['attempts'] if a['archived_minimum_verified_under_own_settings']]
        if not minima:
            continue
        chosen = minima[-1]; meta = json.loads((ROOT / chosen['report']).read_text())
        run = Path(meta['run_dir']); geometry = ROOT / row['source_geometry']
        hessian = run / (Path(meta['input']).stem + '.hess')
        require(sha(hessian) == chosen['hessian_sha256'], 'Changed archived Hessian')
        require(not meta['linear'] and meta['result']['energy_unit'] == 'Hartree', 'Unsupported geometry/energy type')
        result = extract(run / 'output.log', hessian, read(geometry), chosen['energy_Hartree'])
        result['source_sha256'] = {str(Path(p).relative_to(ROOT)): h for p, h in result['source_sha256'].items()}
        evidence.update(result['source_sha256'])
        rows.append(dict(candidate_id=row['candidate_id'], species=row['species'], selection_role=row['selection_role'],
            method=meta['method'], solvent=meta['solvent'], numerical_input_sha256=meta['source_sha256'][meta['input']],
            minimum_report=chosen['report'], geometry=str(geometry.relative_to(ROOT)), geometry_sha256=sha(geometry),
            archived_minimum_verified=True, common_grid_comparison_eligible=False, **result))
    require(len(rows) == lineage['original_minimum_lineages'] == 9, 'Changed archived minimum coverage')
    for rel, expected in evidence.items():
        require(sha(ROOT / rel) == expected, f'Source changed during extraction: {rel}')
    record = dict(status='archived_harmonic_zpe_extracted_not_uniform_panel', rows=rows,
        original_start_count=17, archived_minimum_count=9, source_sha256=evidence,
        wavenumber_cm1_per_Hartree=WAVENUMBER_CM1_PER_HARTREE, Hartree_to_kcal_mol=HARTREE_KCAL,
        sources=SOURCES, source_documentation_checked_on='2026-09-29',
        maximum_absolute_printed_zpe_difference_Hartree=max(abs(r['sum_minus_printed_zpe_Hartree']) for r in rows),
        rankings_computed=False, B_X_scores_computed=False, frequency_grid_convergence_verified=False,
        new_engine_runs=0, scientific_project_complete=False,
        scope='Harmonic ZPE of existing minima under their own archived settings; no common-grid ranking or crystal free energy')
    target = ROOT / 'reports/completion/drug-archived-zpe-r1.json'
    immutable(target, json.dumps(record, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(archived())
