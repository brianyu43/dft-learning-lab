"""Geometry-assigned proton states and whole-complex symmetry comparisons.

These descriptors are not electronic charges, a crystal model, or a polymorph
inventory. Historical energies are not mixed to certify common-grid duplicates.
"""
from itertools import islice, combinations
import io
import json
from pathlib import Path
import re

import networkx as nx
import numpy as np
from ase.io import read
from scipy.optimize import linear_sum_assignment

from completion_qe_bridge import require
from prepare_completion_drug_lineage import build as audit_lineage
from prepare_week3_drug import immutable
from research_run import ROOT, sha, now
from run_week2_xtb_screen import geometry_audit

POLICY = dict(duplicate_all_atom_RMSD_angstrom=.05, duplicate_energy_gap_kcal_mol=.01,
              screening_heavy_RMSD_angstrom=.75, contact_H_acceptor_max_angstrom=2.6,
              contact_DHA_min_angle_degrees=120., hydrogen_owner_ambiguity_margin_angstrom=.1,
              maximum_heavy_graph_mappings=10000, hydrogen_alignment_iterations=20)


def describe(candidate, atoms):
    audit = geometry_audit(candidate, atoms)
    require(audit.get('heavy_connectivity_preserved') and not audit.get('unassigned_hydrogens')
            and not audit.get('separated_geometry', False) and np.isfinite(atoms.positions).all(),
            'Geometry needs chemical/connectivity review')
    z = atoms.numbers; heavy = np.flatnonzero(z > 1); hydrogens = np.flatnonzero(z == 1)
    require(len(heavy) >= 2, 'At least two heavy atoms are required for this drug-state audit')
    owners, distances, ambiguous = {}, [], []
    for h in hydrogens:
        d = np.linalg.norm(atoms.positions[heavy] - atoms.positions[h], axis=1)
        order = np.argsort(d); nearest, second = int(heavy[order[0]]), int(heavy[order[1]])
        require(d[order[0]] > 1e-8, 'Hydrogen overlaps a heavy atom')
        owners[int(h)] = nearest
        gap = float(d[order[1]] - d[order[0]])
        if gap < POLICY['hydrogen_owner_ambiguity_margin_angstrom']: ambiguous.append(int(h))
        distances.append(dict(hydrogen=int(h), nearest_heavy=nearest, nearest_distance_angstrom=float(d[order[0]]),
                              second_heavy=second, second_distance_angstrom=float(d[order[1]]), distance_gap_angstrom=gap))
    counts = candidate.get('fragment_atom_counts', [len(atoms)])
    fragments = np.repeat(np.arange(len(counts)), counts)
    require(len(fragments) == len(atoms), 'Fragment sizes disagree')
    contacts = []
    for h, donor in owners.items():
        if z[donor] not in (7, 8): continue
        for acceptor in heavy:
            if fragments[donor] == fragments[acceptor] or z[acceptor] not in (7, 8): continue
            dh, ah = atoms.positions[donor] - atoms.positions[h], atoms.positions[acceptor] - atoms.positions[h]
            angle = float(np.degrees(np.arccos(np.clip(np.dot(dh, ah) / np.linalg.norm(dh) / np.linalg.norm(ah), -1, 1))))
            if np.linalg.norm(ah) < POLICY['contact_H_acceptor_max_angstrom'] and angle > POLICY['contact_DHA_min_angle_degrees']:
                contacts.append(dict(donor=donor, hydrogen=h, acceptor=int(acceptor),
                                     H_acceptor_distance_angstrom=float(np.linalg.norm(ah)), DHA_angle_degrees=angle))
    return dict(geometry_audit=audit, hydrogen_owners=owners, hydrogen_assignment_distances=distances,
                ambiguous_hydrogen_indices=ambiguous, contacts=contacts,
                scope='Nearest-heavy-atom ownership and geometric N/O-H...N/O contacts; no electronic charge assignment')


def state_graph(candidate, atoms, description):
    graph = nx.Graph()
    counts = candidate.get('fragment_atom_counts', [len(atoms)])
    fragments = np.repeat(np.arange(len(counts)), counts)
    labels = candidate['species'].split('--')
    require(len(labels) == len(counts), 'Fragment species mapping differs')
    for i in np.flatnonzero(atoms.numbers > 1):
        graph.add_node(int(i), element=int(atoms.numbers[i]), fragment_species=labels[fragments[i]],
                       hydrogen_count=sum(owner == i for owner in description['hydrogen_owners'].values()))
    for i, j, order in candidate['graph']['bonds']:
        if atoms.numbers[i] > 1 and atoms.numbers[j] > 1: graph.add_edge(i, j, order=order)
    return graph


def align(p, q):
    pc, qc = p.mean(axis=0), q.mean(axis=0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    translated = (p - pc) @ rotation + qc
    return float(np.sqrt(np.mean(np.sum((translated - q) ** 2, axis=1)))), rotation, pc, qc


def compare(candidate_a, a, candidate_b, b, *, context_a=None, context_b=None, energy_a=None, energy_b=None):
    require(candidate_a['species'] == candidate_b['species'] and len(a) == len(b), 'Different molecular species')
    da, db = describe(candidate_a, a), describe(candidate_b, b)
    require(candidate_a['charge'] == candidate_b['charge'] and candidate_a['multiplicity'] == candidate_b['multiplicity'],
            'Electronic state differs')
    if da['ambiguous_hydrogen_indices'] or db['ambiguous_hydrogen_indices']:
        return dict(status='hydrogen_assignment_requires_review', duplicate_minimum_verified=False)
    ga, gb = state_graph(candidate_a, a, da), state_graph(candidate_b, b, db)
    matcher = nx.algorithms.isomorphism.GraphMatcher(ga, gb,
                node_match=nx.algorithms.isomorphism.categorical_node_match(['element', 'fragment_species', 'hydrogen_count'], [None]*3),
                edge_match=nx.algorithms.isomorphism.categorical_edge_match('order', None))
    cap = POLICY['maximum_heavy_graph_mappings']
    mappings = list(islice(matcher.isomorphisms_iter(), cap + 1))
    complete = len(mappings) <= cap
    if not mappings:
        return dict(status='different_geometry_assigned_state_or_graph', heavy_graph_mapping_count=0,
                    mapping_search_complete=True, duplicate_minimum_verified=False)
    contacts_a = {(c['donor'], c['acceptor']) for c in da['contacts']}
    contacts_b = {(c['donor'], c['acceptor']) for c in db['contacts']}
    heavy = sorted(ga); best = None; minimum_heavy = float('inf')
    for mapping in mappings[:cap]:
        target = [mapping[i] for i in heavy]
        rms, rotation, pc, qc = align(a.positions[heavy], b.positions[target])
        minimum_heavy = min(minimum_heavy, rms)
        contact_match = {(mapping[d], mapping[x]) for d, x in contacts_a} == contacts_b
        previous = None
        for _ in range(POLICY['hydrogen_alignment_iterations']):
            full = dict(mapping)
            transformed = (a.positions - pc) @ rotation + qc
            for owner in heavy:
                ha = [h for h, x in da['hydrogen_owners'].items() if x == owner]
                hb = [h for h, x in db['hydrogen_owners'].items() if x == mapping[owner]]
                require(len(ha) == len(hb), 'Hydrogen-count mapping mismatch')
                if ha:
                    costs = np.sum((transformed[ha, None, :] - b.positions[None, hb, :]) ** 2, axis=2)
                    ii, jj = linear_sum_assignment(costs)
                    full.update({ha[i]: hb[j] for i, j in zip(ii, jj)})
            atom_map = tuple(full[i] for i in range(len(a)))
            value, rotation, pc, qc = align(a.positions, b.positions[list(atom_map)])
            if best is None or (contact_match, -value) > (best['contact_match'], -best['all_atom_RMSD_upper_bound_angstrom']):
                best = dict(all_atom_RMSD_upper_bound_angstrom=value, atom_map_a_to_b=list(atom_map),
                            contact_match=contact_match)
            if atom_map == previous: break
            previous = atom_map
    comparable_energy = bool(context_a is not None and context_a == context_b
                         and energy_a is not None and energy_b is not None
                         and np.isfinite([energy_a, energy_b]).all())
    gap = abs(energy_a - energy_b) * 627.509474 if comparable_energy else None
    close_geometry = best['contact_match'] and best['all_atom_RMSD_upper_bound_angstrom'] <= POLICY['duplicate_all_atom_RMSD_angstrom']
    # This bound can prove separation even without exhaustive H permutations.
    lower_bound = minimum_heavy * np.sqrt(len(heavy) / len(a)) if complete else None
    separated = lower_bound is not None and lower_bound > POLICY['duplicate_all_atom_RMSD_angstrom']
    duplicate = bool(close_geometry and gap is not None and gap <= POLICY['duplicate_energy_gap_kcal_mol'])
    status = ('duplicate_under_declared_tolerances' if duplicate else
              'geometrically_separated_under_declared_tolerance' if separated else
              'geometry_match_energy_context_missing_or_gap_exceeded' if close_geometry else 'requires_pair_review')
    return dict(status=status, **best, minimum_heavy_RMSD_angstrom=minimum_heavy,
                all_atom_RMSD_lower_bound_angstrom=lower_bound, heavy_graph_mapping_count=min(len(mappings), cap),
                mapping_search_complete=complete, energy_comparable=comparable_energy, energy_gap_kcal_mol=gap,
                duplicate_minimum_verified=duplicate, proper_rotation_only=True, fragments_aligned_together=True,
                scope='Threshold-based geometric duplicate check; hydrogen assignment is a candidate alignment, not a proven global all-atom minimum')


def inline_orca_geometry(path):
    text = Path(path).read_text()
    match = re.search(r'^\s*\*\s*xyz\s+0\s+1\s*\n(.*?)^\s*\*\s*$', text, re.M | re.S)
    require(match is not None, 'Expected one neutral singlet inline XYZ input')
    lines = [line.strip() for line in match.group(1).splitlines() if line.strip()]
    return read(io.StringIO(str(len(lines)) + '\ninline ORCA input\n' + '\n'.join(lines) + '\n'), format='xyz')


def historical_report():
    path = ROOT / 'experiments/completion/drug-start-lineage-r1.json'
    frozen = json.loads(path.read_text()); refreshed = audit_lineage()
    require(json.loads(json.dumps(refreshed)) == frozen, 'Frozen lineage audit changed')
    rows, minima = [], []
    for row in frozen['rows']:
        candidate = row['candidate']; cfg = json.loads((ROOT / row['original_config']).read_text())
        inputs = [source for source, destination in cfg.get('assets', {}).items() if destination == 'initial.xyz']
        if inputs:
            initial_path = ROOT / inputs[0]
            initial = read(initial_path)
        else:
            original_meta = json.loads((ROOT / row['attempts'][0]['report']).read_text())
            initial_path = Path(original_meta['run_dir']) / Path(original_meta['input']).name
            initial = inline_orca_geometry(initial_path)
        initial_rel = str(initial_path.relative_to(ROOT))
        require(sha(initial_path) == frozen['input_hashes'][initial_rel], 'Unverified original input geometry')
        initial_state = describe(candidate, initial)
        attempts = [a for a in row['attempts'] if a['archived_minimum_verified_under_own_settings']]
        result = dict(candidate_id=row['candidate_id'], species=row['species'], selection_role=row['selection_role'],
                      initial_geometry_source=initial_rel, initial_geometry_source_sha256=sha(initial_path),
                      initial=initial_state, archived_minimum_verified=False, final=None,
                      net_hydrogen_owner_changes=None)
        if attempts:
            chosen = attempts[-1]; atoms = read(ROOT / chosen['source_xyz']); final = describe(candidate, atoms)
            net = [dict(hydrogen=h, initial_owner=owner, final_owner=final['hydrogen_owners'][h])
                   for h, owner in initial_state['hydrogen_owners'].items() if final['hydrogen_owners'][h] != owner]
            result.update(archived_minimum_verified=True, final=final, net_hydrogen_owner_changes=net,
                          source_geometry=chosen['source_xyz'], source_geometry_sha256=chosen['source_xyz_sha256'],
                          minimum_report=chosen['report'], minimum_frequency_cm1=chosen['lowest_internal_frequency_cm1'])
            minima.append((row, atoms))
        rows.append(result)
    pairs = []
    for (ra, a), (rb, b) in combinations(minima, 2):
        if ra['species'] != rb['species']: continue
        pairs.append(dict(candidate_a=ra['candidate_id'], candidate_b=rb['candidate_id'], species=ra['species'],
                          comparison=compare(ra['candidate'], a, rb['candidate'], b)))
    record = dict(checked_at=now(), status='archived_geometry_review_not_uniform_panel', policy=POLICY,
                  lineage_sha256=sha(path), lineage_source_files_reverified=len(frozen['input_hashes']),
                  original_starts=len(rows), archived_minimum_lineages=len(minima), rows=rows, pair_comparisons=pairs,
                  original_sources_overwritten=False, new_engine_runs=0, common_grid_duplicates_certified=False,
                  cocrystal_polymorph_count=None, B_X_scores_computed=False,
                  script_sha256=sha(Path(__file__)), networkx_version=nx.__version__,
                  scope='Existing gas-phase molecular geometries only; original 17 starts and all failed attempts remain preserved')
    target = ROOT / 'reports/completion/drug-archived-geometry-r1.json'
    immutable(target, json.dumps(record, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(historical_report())
