"""Print a hash-verified archived checkpoint without changing this report."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'checkpoint-archive.json').read_text())
rows = {row['original_path']: row for row in manifest['files']}
if len(sys.argv) != 2 or sys.argv[1] not in rows:
    raise SystemExit('Pass one original_path from checkpoint-archive.json')
row = rows[sys.argv[1]]
packed = (root / row['compressed_path']).read_bytes()
assert hashlib.sha256(packed).hexdigest() == row['compressed_sha256']
data = gzip.decompress(packed)
assert len(data) == row['original_bytes']
assert hashlib.sha256(data).hexdigest() == row['original_sha256']
sys.stdout.buffer.write(data)
