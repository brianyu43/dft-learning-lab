"""Arithmetic for an already provenance-verified, uniform drug energy panel.

This module cannot certify a minimum, a numerical protocol, or a calculation.
The caller must supply energies from one method/environment/protocol at a time.
No engine is launched and no current research result is inferred here.
"""
from __future__ import annotations

import math
from collections.abc import Mapping

HARTREE_KCAL = 627.509474  # Preserve the existing project's conversion convention.
COFORMERS = ("adenine", "caffeine", "theobromine", "xanthine")
OBSERVED = {"adenine": "formed", "caffeine": "formed",
            "theobromine": "not_formed", "xanthine": "not_formed"}
TIE_WIDTHS = (0.0, 0.5, 1.0, 2.0)
CONTEXT_KEYS = ("method", "environment", "numerical_protocol_id", "energy_quantity")


def finite_number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Expected a finite numeric energy, got {value!r}")
    return float(value)


def contrast(energy, terms):
    """Return kcal/mol and exact missing species; never substitute zero."""
    missing = [name for name in terms if name not in energy or energy[name] is None]
    if missing:
        return {"kcal_mol": None, "missing_species": missing}
    value = math.fsum(coefficient * finite_number(energy[name]) for name, coefficient in terms.items())
    return {"kcal_mol": value * HARTREE_KCAL, "missing_species": []}


def panel_scores(energies: Mapping[str, float | None]):
    """The four B/X values, including explicit partial-data states."""
    for value in energies.values():
        if value is not None:
            finite_number(value)
    records = []
    for coformer in COFORMERS:
        cross = f"diflunisal--{coformer}"
        records.append({"coformer": coformer,
            "B": contrast(energies, {cross: 1., "diflunisal": -1., coformer: -1.}),
            "X": contrast(energies, {cross: 1., "diflunisal--diflunisal": -.5, f"{coformer}--{coformer}": -.5})})
    return records


def compare_observed_groups(scores, width):
    """Four shared-molecule contrasts, not four independent experiments.

    Positive delta means the known formed pair has the lower energy score.
    A difference equal to the chosen width is a tie. A separate 1e-8 kcal/mol
    arithmetic guard handles machine rounding at a boundary; it is not the
    scientific numerical uncertainty estimate.
    """
    width = finite_number(width)
    if width < 0:
        raise ValueError("Tie width must be nonnegative")
    for value in scores.values():
        if value is not None:
            finite_number(value)
    guard = 1e-8
    rows = []
    for positive in COFORMERS[:2]:
        for negative in COFORMERS[2:]:
            p, n = scores.get(positive), scores.get(negative)
            delta = None if p is None or n is None else math.fsum([finite_number(n), -finite_number(p)])
            if delta is None:
                outcome = "unavailable"
            elif abs(delta) <= width + guard:
                outcome = "tie"
            else:
                outcome = "consistent" if delta > 0 else "reversed"
            rows.append({"formed_pair_coformer": positive, "not_formed_pair_coformer": negative,
                         "delta_not_formed_minus_formed_kcal_mol": delta, "assessment": outcome})
    return {"tie_width_kcal_mol": width, "boundary_roundoff_guard_kcal_mol": guard,
            "rows": rows, "counts": {k: sum(r['assessment'] == k for r in rows)
                                      for k in ("consistent", "reversed", "tie", "unavailable")},
            "independent_observation_count": None,
            "scope": "four contrasts from four molecular pairs sharing diflunisal; not independent replicates"}


def analyze_uniform_records(records, expected_context):
    """Select by electronic energy only, keeping primary/all scopes separate.

    This is a narrow grouping safeguard, not a replacement for raw output,
    geometry, Hessian and numerical-sensitivity verification. Input records
    must already have been curated from such verified results by the caller.
    """
    if set(expected_context) != set(CONTEXT_KEYS) or any(not expected_context[k] for k in CONTEXT_KEYS):
        raise ValueError("Specify method, environment, numerical protocol and energy quantity")
    allowed_species = {"diflunisal", "diflunisal--diflunisal", *COFORMERS,
                       *(f'diflunisal--{b}' for b in COFORMERS), *(f'{b}--{b}' for b in COFORMERS)}
    seen = set()
    for row in records:
        if row['comparison_context'] != expected_context:
            raise ValueError("Cannot combine different methods, environments, grids or energy quantities")
        if row['record_id'] in seen:
            raise ValueError("Repeated record_id")
        seen.add(row['record_id'])
        if row['species'] not in allowed_species:
            raise ValueError("Species outside the fixed five-monomer/nine-dimer panel")
        if row['selection_role'] not in ('primary', 'alternative'):
            raise ValueError("Selection role must come from the frozen start policy")
        finite_number(row['energy_Hartree'])
    result = {"status": "arithmetic_only_requires_external_provenance_verification",
              "comparison_context": expected_context, "scopes": {}, "class_predictor_fitted": False}
    for scope in ('primary_only', 'all_selected_poses'):
        pool = [r for r in records if scope != 'primary_only' or r['selection_role'] == 'primary']
        by_species = {name: [r for r in pool if r['species'] == name] for name in sorted(allowed_species)}
        best = {name: min(rs, key=lambda r: (r['energy_Hartree'], r['record_id']))
                for name, rs in by_species.items() if rs}
        scores = panel_scores({name: r['energy_Hartree'] for name, r in best.items()})
        comparisons = {score: [compare_observed_groups({r['coformer']: r[score]['kcal_mol'] for r in scores}, width)
                               for width in TIE_WIDTHS] for score in ('B', 'X')}
        result['scopes'][scope] = {"scores": scores, "comparisons": comparisons,
            "selected_record_ids": {name: r['record_id'] for name, r in best.items()},
            "available_record_counts_not_unique_minima": {name: len(rs) for name, rs in by_species.items()},
            "pose_gaps_kcal_mol": {name: [{"record_id": r['record_id'],
                                          "gap_kcal_mol": math.fsum([r['energy_Hartree'], -best[name]['energy_Hartree']]) * HARTREE_KCAL}
                                         for r in sorted(rs, key=lambda r: (r['energy_Hartree'], r['record_id']))]
                                   for name, rs in by_species.items() if rs}}
    return result
