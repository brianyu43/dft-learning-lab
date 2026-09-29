"""Geometry and raw-SCF audit of a fixed-registry 24-atom trajectory.

The reference-derived Mo-X distance graph is a geometric continuity check,
not an electronic bond-order calculation or proof of physical stability.
"""
from __future__ import annotations

import numpy as np

from completion_alloy_geometry import FixMoLayerCentroids, centroids, validate_geometry
from completion_qe_bridge import load_scf, require
from completion_qe_relaxation import settings_fingerprint


class RegistryGeometryAudit:
    def __init__(self, initial):
        validate_geometry(initial)
        self.initial = initial.copy(); self.initial.calc = None; self.initial.set_constraint()
        self.centers = centroids(initial)
        self.mos = np.flatnonzero(initial.numbers == 42)
        self.chalcogens = np.flatnonzero(initial.numbers != 42)
        self.layers = initial.arrays['layer_id'].copy()
        self.faces = np.sign(initial.positions[self.chalcogens, 2]
                             - self.centers[self.layers[self.chalcogens], 2])
        require(np.all(self.faces != 0), 'Reference chalcogen lies on the Mo centroid plane')
        distances = initial.get_all_distances(mic=True)[np.ix_(self.mos, self.chalcogens)]
        expected = np.zeros(distances.shape, dtype=bool)
        for j, atom in enumerate(self.chalcogens):
            same_layer = np.flatnonzero(self.layers[self.mos] == self.layers[atom])
            order = same_layer[np.argsort(distances[same_layer, j])]
            require(distances[order[3], j] - distances[order[2], j] > 1e-8,
                    'Reference nearest-three-Mo shell is ambiguous')
            expected[order[:3], j] = True
        require(np.all(expected.sum(axis=0) == 3) and np.all(expected.sum(axis=1) == 6),
                'Reference lacks the expected Mo-X coordination graph')
        bonded_max, nonbonded_min = float(distances[expected].max()), float(distances[~expected].min())
        require(nonbonded_min - bonded_max > 1e-8, 'No distance gap for a fixed reference graph')
        self.cutoff = .5 * (bonded_max + nonbonded_min)
        self.expected_graph = expected
        all_distances = initial.get_all_distances(mic=True)
        self.minimum_distance = float(all_distances[np.triu_indices(len(initial), 1)].min())
        self.reference_shell = dict(max_neighbor_distance_angstrom=bonded_max,
                                    min_non_neighbor_distance_angstrom=nonbonded_min,
                                    fixed_distance_cutoff_angstrom=self.cutoff)

    def inspect(self, atoms):
        require(np.array_equal(atoms.numbers, self.initial.numbers), 'Atom composition/order changed')
        for key in ('layer_id', 'primitive_index', 'tile_ij'):
            require(key in atoms.arrays and np.array_equal(atoms.arrays[key], self.initial.arrays[key]),
                    'Atom provenance changed: ' + key)
        require(np.array_equal(atoms.pbc, self.initial.pbc)
                and np.allclose(atoms.cell, self.initial.cell, rtol=0, atol=1e-10), 'Fixed cell/periodicity changed')
        require(np.isfinite(atoms.positions).all(), 'Nonfinite relaxed coordinates')
        centers = centroids(atoms)
        drift = float(np.abs(centers - self.centers).max())
        heights = self.faces * (atoms.positions[self.chalcogens, 2]
                               - centers[self.layers[self.chalcogens], 2])
        changed_faces = self.chalcogens[heights <= 1e-8].tolist()
        # Spatial layer separation is evaluated from actual coordinates, not labels alone.
        ranges = [[float(atoms.positions[self.layers == layer, 2].min()),
                   float(atoms.positions[self.layers == layer, 2].max())] for layer in (0, 1)]
        gap = ranges[1][0] - ranges[0][1]
        vacuum = float(atoms.cell[2, 2] - np.ptp(atoms.positions[:, 2]))
        contained = bool(atoms.positions[:, 2].min() > 0 and atoms.positions[:, 2].max() < atoms.cell[2, 2])
        all_distances = atoms.get_all_distances(mic=True)
        minimum_distance = float(all_distances[np.triu_indices(len(atoms), 1)].min())
        distances = all_distances[np.ix_(self.mos, self.chalcogens)]
        graph = distances < self.cutoff

        def edges(mask):
            return [dict(mo_index=int(self.mos[i]), chalcogen_index=int(self.chalcogens[j]),
                         distance_angstrom=float(distances[i, j])) for i, j in np.argwhere(mask)]

        lost, added = edges(self.expected_graph & ~graph), edges(~self.expected_graph & graph)
        reasons = []
        if drift > 1e-8: reasons.append('Mo_centroids_moved')
        if changed_faces: reasons.append('chalcogen_crossed_reference_Mo_face')
        if gap <= 1e-8: reasons.append('spatial_layer_ranges_overlap')
        if vacuum <= 0 or not contained: reasons.append('outside_unwrapped_slab_cell_or_no_vacuum')
        if lost or added: reasons.append('reference_distance_graph_changed')
        if minimum_distance < .5 * self.minimum_distance: reasons.append('severe_atomic_approach_requires_review')
        return dict(status='geometry_continuity_passed' if not reasons else 'requires_review',
                    reasons=reasons, Mo_centroids_angstrom=centers.tolist(),
                    max_centroid_drift_angstrom=drift, changed_face_atom_indices=changed_faces,
                    min_retained_signed_face_height_angstrom=float(heights.min()),
                    layer_z_ranges_angstrom=ranges, interlayer_surface_gap_angstrom=gap,
                    current_vacuum_gap_angstrom=vacuum, coordinates_within_unwrapped_z_cell=contained,
                    minimum_interatomic_distance_angstrom=minimum_distance,
                    severe_approach_threshold_angstrom=.5 * self.minimum_distance,
                    reference_shell=self.reference_shell, lost_distance_graph_edges=lost,
                    added_distance_graph_edges=added,
                    scope='Reference geometry continuity only; no bond-order or dynamical-stability claim')


def audit_qe_trajectory(initial, reference_report, trace, *, method='C09', fmax=.005, maxstep=.1):
    """Re-read every SCF and recompute projected forces before accepting a trace."""
    require(trace and [row['step'] for row in trace] == list(range(len(trace))), 'Empty or discontinuous trajectory')
    require(len(trace) <= 101 and 0 < fmax <= .005 and 0 < maxstep <= .1, 'Unreviewed trajectory limits')
    require(not initial.info.get('fixture_only', False), 'Fixture is not a production trajectory')
    require(all(row.get('scf_report') for row in trace), 'Every production step requires its raw SCF report')
    geometry_audit = RegistryGeometryAudit(initial)
    load_scf(reference_report, expected_method=method, expected_atoms=initial)
    fingerprint = settings_fingerprint(reference_report)
    constraint = FixMoLayerCentroids.from_atoms(initial)
    previous = initial.positions.copy()
    rows = []
    for saved in trace:
        atoms = initial.copy(); atoms.set_constraint(); atoms.calc = None
        positions = np.asarray(saved['positions_angstrom'], dtype=float)
        require(positions.shape == initial.positions.shape and np.isfinite(positions).all(), 'Bad trajectory coordinates')
        atoms.set_positions(positions, apply_constraint=False)
        displacement = float(np.linalg.norm(positions - previous, axis=1).max())
        require(displacement <= (1e-10 if saved['step'] == 0 else maxstep + 1e-8),
                'Initial geometry changed, step exceeded its limit, or coordinates were wrapped')
        result = load_scf(saved['scf_report'], expected_method=method, expected_atoms=atoms)
        require(settings_fingerprint(saved['scf_report']) == fingerprint, 'Trajectory SCF settings changed')
        raw = result['forces'].copy(); projected = raw.copy(); constraint.adjust_forces(atoms, projected)
        values = dict(energy_eV=float(result['energy_eV']),
                      raw_max_force_eV_angstrom=float(np.linalg.norm(raw, axis=1).max()),
                      projected_max_force_eV_angstrom=float(np.linalg.norm(projected, axis=1).max()),
                      reaction_max_force_eV_angstrom=float(np.linalg.norm(raw - projected, axis=1).max()),
                      max_Mo_centroid_drift_angstrom=float(np.abs(centroids(atoms) - centroids(initial)).max()))
        for key, value in values.items():
            require(key in saved and np.isfinite(saved[key]) and abs(saved[key] - value) <= 1e-8,
                    'Saved trace disagrees with raw SCF/geometry: ' + key)
        rows.append(dict(step=saved['step'], **values, maximum_step_displacement_angstrom=displacement,
                         geometry=geometry_audit.inspect(atoms), raw_scf_audit=result['audit']))
        previous = positions.copy()
    geometry_preserved = all(row['geometry']['status'] == 'geometry_continuity_passed' for row in rows)
    force_converged = rows[-1]['projected_max_force_eV_angstrom'] < fmax
    return dict(status='verified_constrained_trajectory' if geometry_preserved and force_converged else 'requires_review',
                geometry_preserved=geometry_preserved, projected_force_converged=bool(force_converged),
                fmax_eV_angstrom=fmax, maxstep_angstrom=maxstep, rows=rows,
                evidence_kind='verified_QE_SCF', separated_final_layers_complete=False,
                scope='Fixed-registry constrained trajectory; final separated-layer energies and sensitivity remain separate')
