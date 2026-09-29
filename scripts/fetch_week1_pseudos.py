"""Download three versioned official PseudoDojo files; TLS verification stays on."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
target = ROOT / 'data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'
target.mkdir(parents=True, exist_ok=True)
records = []
for element in ('Mo', 'S', 'Se'):
    url = f'https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard/{element}.upf.gz'
    archive = target / f'{element}.upf.gz'
    subprocess.run(['curl', '-fsSL', '--max-time', '60', '-o', str(archive), url], check=True)
    raw = archive.read_bytes()
    data = gzip.decompress(raw)
    if b'<UPF' not in data or b'PP_HEADER' not in data:
        raise ValueError(f'Not UPF: {element}')
    (target / f'{element}.upf').write_bytes(data)
    records.append({'element': element, 'url': url,
                    'compressed_sha256': hashlib.sha256(raw).hexdigest(),
                    'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})
    print(f'{element}: {len(data)} bytes', flush=True)
(target / 'provenance.json').write_text(json.dumps({
    'retrieved_at': datetime.now(timezone.utc).isoformat(),
    'source': 'https://www.pseudo-dojo.org/', 'family': 'nc-sr-04_pbe_standard',
    'records': records}, indent=2) + '\n')
