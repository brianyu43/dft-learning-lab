"""Geometry-only preparation for the conditional 24-atom interface comparison.

No guessed lattice or distance is promoted to a production DFT input here.
The production caller must supply the verified common-lattice monolayers.
"""
from collections import Counter

import numpy as np
import spglib
from ase.constraints import FixConstraint
from scipy.optimize import linear_sum_assignment

from week3_materials_pipeline import REGISTRIES, bilayer


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _canonical_layer(atoms, chalcogen):
    _require(Counter(atoms.get_chemical_symbols()) == Counter(Mo=1, **{chalcogen: 2}),
             'Expected one Mo and two identical chalcogens per primitive layer')
    _require(np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell).all(),
             'Non-finite geometry')
    _require(atoms.pbc.all() and np.linalg.det(atoms.cell) > 0, 'Expected a right-handed periodic cell')
    _require(np.allclose(atoms.cell[:2, 2], 0, atol=1e-10, rtol=0)
             and np.allclose(atoms.cell[2, :2], 0, atol=1e-10, rtol=0),
             'The layer normal must be the Cartesian z axis')
    mo = int(np.flatnonzero(atoms.numbers == 42)[0])
    others = [i for i in range(3) if i != mo]
    others.sort(key=lambda i: atoms.positions[i, 2], reverse=True)
    _require(atoms.positions[others[0], 2] > atoms.positions[mo, 2]
             > atoms.positions[others[1], 2], 'Mo must lie between its chalcogen planes')
    result = atoms[[mo, *others]]
    result.positions -= result.positions[0]
    return result


def build_family(bottom, top, registry_id, distance, vacuum, top_exchange_tile=(0, 0)):
    """Return C0/CI/CO with identical coordinates; only two atomic species change.

    Site ids are (primitive tile i, j, primitive index), where primitive indices
    are bottom Mo/upper-X/lower-X and top Mo/upper-X/lower-X. Initial exchange
    tiles are (0,0) on both layers. A different top tile is for a separately
    reviewed independent-placement experiment, never selected by energy here.
    """
    _require(registry_id in range(3), 'Unknown registry')
    _require(np.isfinite(distance) and distance > 0 and np.isfinite(vacuum) and vacuum > 0,
             'Distance and vacuum must be positive and finite')
    _require(tuple(top_exchange_tile) in [(0, 0), (0, 1), (1, 0), (1, 1)],
             'Exchange tile lies outside the 2x2 supercell')
    b = _canonical_layer(bottom, 'S')
    t = _canonical_layer(top, 'Se')
    _require(np.allclose(b.cell[:2], t.cell[:2], atol=1e-9, rtol=0), 'Different in-plane cells')
    _require(distance > b.positions[:, 2].max() - t.positions[:, 2].min(), 'Layers overlap in z')
    primitive = bilayer(b, t, registry_id, float(distance), float(vacuum))
    primitive.set_array('layer_id', np.array([0, 0, 0, 1, 1, 1]))
    primitive.set_array('primitive_index', np.arange(6))
    parent = primitive.repeat((2, 2, 1))
    tiles = np.repeat(np.array([(0, 0), (0, 1), (1, 0), (1, 1)]), 6, axis=0)
    parent.set_array('tile_ij', tiles)
    expected = (np.tile(primitive.positions, (4, 1))
                + tiles[:, :1] * primitive.cell[0] + tiles[:, 1:] * primitive.cell[1])
    _require(np.allclose(parent.positions, expected, atol=1e-12, rtol=0), 'Unexpected repeat ordering')

    def site(tile, primitive_index):
        indices = np.flatnonzero((tiles == tile).all(axis=1)
                                 & (parent.arrays['primitive_index'] == primitive_index))
        _require(len(indices) == 1, 'Site is not unique')
        return int(indices[0])

    exchanges = {'C0': [], 'CI': [site((0, 0), 1), site(top_exchange_tile, 5)],
                 'CO': [site((0, 0), 2), site(top_exchange_tile, 4)]}
    family = {}
    for label, indices in exchanges.items():
        atoms = parent.copy()
        if indices:
            atoms.numbers[indices] = atoms.numbers[indices[::-1]]
        atoms.info.update(configuration_id=label, registry_id=int(registry_id),
                          registry_fraction=list(REGISTRIES[registry_id]),
                          primitive_cell_angstrom=primitive.cell.array.tolist(),
                          Mo_plane_distance_angstrom=float(distance),
                          bilayer_vacuum_gap_angstrom=float(vacuum),
                          exchange_indices=indices, top_exchange_tile=list(top_exchange_tile),
                          coordinate_scope='Frozen common-lattice monolayers; species exchange only')
        validate_geometry(atoms)
        _require(np.array_equal(atoms.positions, parent.positions), 'Exchange moved atoms')
        family[label] = atoms
    return family


def mo_groups(atoms):
    _require('layer_id' in atoms.arrays, 'Explicit layer membership is required')
    groups = [np.flatnonzero((atoms.arrays['layer_id'] == layer) & (atoms.numbers == 42))
              for layer in (0, 1)]
    _require(all(len(group) == 4 for group in groups), 'Expected four Mo atoms in each layer')
    return groups


def centroids(atoms):
    return np.array([atoms.positions[indices].mean(axis=0) for indices in mo_groups(atoms)])


def validate_geometry(atoms, tolerance=1e-8):
    """Verify frozen structures, not constrained or unconstrained relaxed ones."""
    _require(len(atoms) == 24 and Counter(atoms.get_chemical_symbols()) == Counter(Mo=8, S=8, Se=8),
             'Total composition is not Mo8S8Se8')
    label = atoms.info['configuration_id']
    _require(label in ('C0', 'CI', 'CO'), 'Unknown configuration')
    required = ([Counter(Mo=4, S=8), Counter(Mo=4, Se=8)] if label == 'C0'
                else [Counter(Mo=4, S=7, Se=1), Counter(Mo=4, S=1, Se=7)])
    for layer, counts in enumerate(required):
        _require(Counter(atoms[atoms.arrays['layer_id'] == layer].get_chemical_symbols()) == counts,
                 'Layer composition changed')
    centers = centroids(atoms)
    shift = centers[1] - centers[0]
    primitive_cell = np.array(atoms.info['primitive_cell_angstrom'])
    r = atoms.info['registry_fraction']
    expected = r[0] * primitive_cell[0] + r[1] * primitive_cell[1]
    expected[2] = atoms.info['Mo_plane_distance_angstrom']
    _require(np.allclose(shift, expected, atol=tolerance, rtol=0), 'Mo registry or distance changed')
    _require(np.allclose(atoms.cell[:2], primitive_cell[:2] * 2, atol=tolerance, rtol=0),
             'Not a 2x2 primitive supercell')
    _require(abs(atoms.cell[2, 2] - np.ptp(atoms.positions[:, 2])
                 - atoms.info['bilayer_vacuum_gap_angstrom']) <= tolerance, 'Incorrect vacuum gap')
    exchanged = [int(i) for i in atoms.info['exchange_indices']]
    actual = np.flatnonzero(((atoms.arrays['layer_id'] == 0) & (atoms.numbers == 34))
                           | ((atoms.arrays['layer_id'] == 1) & (atoms.numbers == 16))).tolist()
    _require(sorted(exchanged) == sorted(actual), 'Exchange annotation disagrees with species')
    if label != 'C0':
        for layer, index in enumerate(exchanged):
            _require(atoms.arrays['layer_id'][index] == layer, 'Wrong exchange layer')
            dz = atoms.positions[index, 2] - centers[layer, 2]
            sign = (1 if layer == 0 else -1) * (1 if label == 'CI' else -1)
            _require(sign * dz > 0, 'Exchange is on the wrong face')
    return dict(configuration=label, total_formula=atoms.get_chemical_formula(),
                layer_formulas=[atoms[atoms.arrays['layer_id'] == i].get_chemical_formula() for i in (0, 1)],
                Mo_centroids_angstrom=centers.tolist(), exchange_indices=exchanged,
                registry_id=int(atoms.info['registry_id']), distance_angstrom=float(shift[2]))


def _match(frac_a, frac_b, types_a, types_b, cell, tolerance):
    delta = frac_a[:, None, :] - frac_b[None, :, :]
    # Exact-equivalence test only: a near-zero periodic displacement necessarily
    # lies near an integer fractional vector, also for this skew hexagonal cell.
    delta -= np.rint(delta)
    distances = np.linalg.norm(delta @ cell, axis=2)
    allowed = (distances <= tolerance) & (types_a[:, None] == types_b[None, :])
    rows, cols = linear_sum_assignment(~allowed)
    if not allowed[rows, cols].all():
        return None
    return cols, float(distances[rows, cols].max())


def decoration_equivalence(first, second, tolerance=1e-5):
    """Check decorations of the same Mo/X parent under its full space group.

    Permutation, periodic wrapping and a different global origin are allowed.
    Unlike swapping each layer's faces separately, every tested operation acts
    on the entire periodic bilayer. Different parent geometries are rejected.
    """
    _require(len(first) == len(second) and np.allclose(first.cell, second.cell, atol=1e-9, rtol=0),
             'Equivalence requires the same cell and atom count')
    _require(first.pbc.all() and second.pbc.all() and 0 < tolerance <= 1e-2,
             'Invalid periodicity or tolerance')
    fa, fb = first.get_scaled_positions(), second.get_scaled_positions()
    pa = np.where(first.numbers == 42, 42, 16)
    pb = np.where(second.numbers == 42, 42, 16)
    anchor = int(np.flatnonzero(pa == 42)[0])
    alignment = None
    for target in np.flatnonzero(pb == pa[anchor]):
        shift = fb[target] - fa[anchor]
        if _match(fa + shift, fb, pa, pb, first.cell.array, tolerance) is not None:
            alignment = shift
            break
    _require(alignment is not None, 'The species-blind parent geometries differ')
    symmetry = spglib.get_symmetry((first.cell.array, fa, pa), symprec=tolerance)
    _require(symmetry is not None, 'Parent symmetry search failed')
    operations = list(zip(symmetry['rotations'], symmetry['translations']))
    for index, (rotation, translation) in enumerate(operations):
        offset = translation + alignment
        match = _match(fa @ rotation.T + offset, fb, first.numbers, second.numbers,
                       first.cell.array, tolerance)
        if match is not None:
            return dict(equivalent=True, parent_operation_count=len(operations),
                        operations_checked=index + 1, tolerance_angstrom=tolerance,
                        permutation=match[0].tolist(), max_residual_angstrom=match[1],
                        rotation_fractional=rotation.tolist(), translation_fractional=offset.tolist())
    return dict(equivalent=False, parent_operation_count=len(operations), operations_checked=len(operations),
                tolerance_angstrom=tolerance)


class FixMoLayerCentroids(FixConstraint):
    """Six linear constraints on the unwrapped Cartesian Mo coordinates.

    All members of each group must have equal mass, hence their centroid is
    also their center of mass. Each projection touches its own group only.
    Chalcogen motions and relative motions among a layer's Mo atoms remain free.
    Use with a fixed cell and never wrap positions between optimization steps.
    """
    def __init__(self, groups, reference_centroids):
        self.groups = [np.array(group, dtype=int) for group in groups]
        self.reference_centroids = np.array(reference_centroids, dtype=float)
        flat = np.concatenate(self.groups)
        _require(len(self.groups) == 2 and all(len(g) == 4 for g in self.groups)
                 and len(np.unique(flat)) == 8 and (flat >= 0).all(), 'Invalid Mo groups')
        _require(self.reference_centroids.shape == (2, 3) and np.isfinite(self.reference_centroids).all(),
                 'Invalid centroid reference')

    @classmethod
    def from_atoms(cls, atoms):
        return cls(mo_groups(atoms), centroids(atoms))

    def _validate(self, atoms):
        for indices in self.groups:
            _require((indices < len(atoms)).all() and (atoms.numbers[indices] == 42).all(),
                     'Constraint atom membership changed')
            masses = atoms.get_masses()[indices]
            _require(np.allclose(masses, masses[0], rtol=0, atol=1e-12), 'Unequal Mo masses')

    def get_removed_dof(self, atoms):
        self._validate(atoms)
        return 6

    def get_indices(self):
        return np.concatenate(self.groups)

    def adjust_positions(self, atoms, newpositions):
        self._validate(atoms)
        for indices, target in zip(self.groups, self.reference_centroids):
            newpositions[indices] += target - newpositions[indices].mean(axis=0)

    def adjust_forces(self, atoms, forces):
        self._validate(atoms)
        for indices in self.groups:
            forces[indices] -= forces[indices].mean(axis=0)

    def adjust_momenta(self, atoms, momenta):
        self.adjust_forces(atoms, momenta)

    def index_shuffle(self, atoms, ind):
        selected = np.arange(len(atoms))[ind]
        _require(len(selected) == len(atoms) and len(np.unique(selected)) == len(atoms),
                 'Constrained structures must not be sliced; extract references before attaching constraints')
        old_to_new = {int(old): new for new, old in enumerate(selected)}
        self.groups = [np.array([old_to_new[int(i)] for i in group]) for group in self.groups]

    def todict(self):
        # Explicit project loader required; not registered in ASE dict2constraint.
        return dict(name=type(self).__name__, kwargs=dict(groups=[g.tolist() for g in self.groups],
                    reference_centroids=self.reference_centroids.tolist()))

    @classmethod
    def fromdict(cls, value):
        _require(value['name'] == cls.__name__, 'Unknown constraint')
        return cls(**value['kwargs'])
