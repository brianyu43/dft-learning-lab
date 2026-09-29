"""Read the official MACE alloy archive without executing code or model pickle."""
from collections import Counter
import hashlib
import io
import json
import re
import zipfile

import numpy as np

from research_run import ROOT, now, save, sha
from prepare_week3_drug import immutable


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def frames(raw):
    stream = io.StringIO(raw.decode('utf8'))
    result = []
    while True:
        line = stream.readline()
        if not line:
            break
        if not line.strip():
            continue
        n = int(line)
        header = stream.readline()
        properties = re.search(r'Properties=(\S+)', header).group(1)
        assert properties.startswith('species:S:1:pos:R:3'), 'Unexpected coordinate column schema'
        cell = np.fromstring(re.search(r'Lattice="([^"]+)"', header).group(1), sep=' ').reshape(3, 3)
        rows = [stream.readline() for _ in range(n)]
        assert all(row for row in rows)
        tokens = [row.split() for row in rows]
        symbols = [row[0] for row in tokens]
        positions = np.array([[float(x) for x in row[1:4]] for row in tokens])
        assert np.isfinite(positions).all() and np.isfinite(cell).all()
        assert set(symbols) <= {'Mo', 'W', 'S', 'Se'}
        normalized = np.round(np.vstack((cell, positions)), 8)
        normalized[normalized == 0] = 0.0
        key = digest(json.dumps([symbols, normalized.tolist()], separators=(',', ':')).encode())
        metal = np.array([symbol in ('Mo', 'W') for symbol in symbols])
        zspan = float(np.ptp(positions[metal, 2])) if metal.any() else None
        result.append(dict(index=len(result), n_atoms=n, formula_counts=dict(sorted(Counter(symbols).items())),
                           ordered_geometry_sha256=key, raw_frame_sha256=digest((line + header + ''.join(rows)).encode()),
                           metal_z_span_angstrom=zspan, total_z_span_angstrom=float(np.ptp(positions[:, 2])),
                           cell_c_angstrom=float(np.linalg.norm(cell[2])), properties=properties))
    return result


def main():
    base = ROOT / 'data/completion/materials-models'
    archive_path = base / 'quaternary_repo.zip'
    provenance = json.loads((base / 'quaternary-archive-source.json').read_text())
    assert sha(archive_path) == provenance['sha256']
    assert archive_path.stat().st_size == 55899962
    assert hashlib.md5(archive_path.read_bytes()).hexdigest() == '201aa56c63ec99d0e966134a59b4793c'
    prefix = 'quaternary_repo/'
    members_to_read = ['readme.txt', 'training/run_mace_training', 'training/make_training_data.py',
                       'training/MoWSSe_MACE_201_run-123.log', 'binding_energy/get_binding_e_data.py']
    files, datasets, raw_splits = {}, {}, {}
    with zipfile.ZipFile(archive_path) as archive:
        for name in members_to_read:
            raw = archive.read(prefix + name)
            path = base / 'quaternary-inspected-text' / name
            immutable(path, raw.decode('utf8'))
            files[prefix + name] = dict(sha256=digest(raw), local=str(path.relative_to(ROOT)))
        for split, name in [('train', 'MoWSSe_allxyT.xyz'), ('test', 'MoWSSe_allxyT_test.xyz')]:
            raw = archive.read(prefix + 'training/' + name)
            raw_splits[split] = raw
            rows = frames(raw)
            dataset = dict(archive_member=prefix + 'training/' + name, file_sha256=digest(raw),
                           frame_count=len(rows), atom_count_distribution=dict(Counter(row['n_atoms'] for row in rows)),
                           maximum_metal_z_span_angstrom=max(row['metal_z_span_angstrom'] or 0 for row in rows),
                           maximum_all_atom_z_span_angstrom=max(row['total_z_span_angstrom'] for row in rows),
                           frames=rows)
            datasets[split] = dataset
        model_name = prefix + 'training/MoWSSe_MACE_201_run-123.model'
        model_raw = archive.read(model_name)
        model = dict(archive_member=model_name, bytes=len(model_raw), sha256=digest(model_raw),
                     model_deserialized=False, model_inference_performed=False)
    geometry = {split: Counter(row['ordered_geometry_sha256'] for row in data['frames'] if row['n_atoms'] > 1)
                for split, data in datasets.items()}
    raw_frames = {split: Counter(row['raw_frame_sha256'] for row in data['frames'] if row['n_atoms'] > 1)
                  for split, data in datasets.items()}
    overlap = dict(scope='Non-isolated frames; exact same atom ordering, cell and Cartesian positions rounded to 8 decimals',
                   test_is_exact_byte_prefix_of_train=raw_splits['train'].startswith(raw_splits['test']),
                   trailing_train_bytes=len(raw_splits['train']) - len(raw_splits['test']),
                   unique_ordered_geometries={s: len(c) for s, c in geometry.items()},
                   shared_unique_ordered_geometries=len(geometry['train'].keys() & geometry['test'].keys()),
                   shared_geometry_frame_multiplicity=sum((geometry['train'] & geometry['test']).values()),
                   shared_raw_frame_multiplicity=sum((raw_frames['train'] & raw_frames['test']).values()))
    output = dict(status='archive_scope_audited', checked_at=now(), provenance=provenance,
                  script_sha256=sha(ROOT / 'scripts/audit_completion_materials_archive.py'),
                  inspected_text_files=files, datasets=datasets, train_test_overlap=overlap, model=model,
                  recorded_training=dict(mace_version='0.2.0', cutoff_angstrom=5.0,
                                         interaction_blocks=2, elements=['S', 'Se', 'Mo', 'W']),
                  interpretation=dict(
                      generated_geometry='Single MX2 sheet; binding_energy script measures substitution-pair energetics within it',
                      interlayer_training_or_validation_evidence_found=False,
                      suitable_for_primary_CI_CO_interlayer_model_comparison=False,
                      decision='Exclude this released monolayer model from in-domain interlayer performance claims; retain as a possible monolayer-alloy follow-up',
                      dft_reference='Code variables name optB88; exact electronic settings and pseudopotentials not established by this archive audit',
                      limits='Dataset membership and recorded config inspected; model architecture/weights not executed; this is not a general assessment of MACE'),
                  source_url='https://zenodo.org/records/10912291')
    save(ROOT / 'reports/completion/materials-mace-archive-audit.json', output)
    print(json.dumps(dict(status=output['status'], model=model, train_test_overlap=overlap,
                          datasets={s: {k: v for k, v in d.items() if k != 'frames'} for s, d in datasets.items()}), indent=2))


if __name__ == '__main__':
    main()
