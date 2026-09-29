"""Verify the immutable file inventory of the public snapshot, without engines."""
import hashlib
import json
from pathlib import Path
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    manifest_path = root / 'publication-manifest.json'
    if not manifest_path.is_file():
        raise SystemExit('publication-manifest.json is missing')
    manifest = json.loads(manifest_path.read_text())
    errors = []
    total = 0
    for relative, expected in manifest['files'].items():
        path = root / relative
        if not path.resolve().is_relative_to(root.resolve()) or path.is_symlink():
            errors.append(f'Unsafe path: {relative}')
            continue
        if not path.is_file():
            errors.append(f'Missing: {relative}')
            continue
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if len(data) != expected['bytes'] or digest != expected['sha256']:
            errors.append(f'Changed: {relative}')
        total += len(data)
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        raise SystemExit(1)
    print(f"Verified {len(manifest['files'])} files ({total:,} bytes); scientific completion is not implied.")


if __name__ == '__main__':
    main()
