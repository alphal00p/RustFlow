#!/usr/bin/env python3
"""Verify the additive archive map and all unchanged inner manifest entries."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--restore', action='store_true',
                        help='Restore missing logical files byte-exactly; never overwrite existing bytes.')
    args = parser.parse_args()
    record = json.loads((HERE / 'archive-map.json').read_bytes())
    indexed = {r['logical_path']: r for r in record['entries']}
    for row in record['entries']:
        data = (ROOT / row['archive_path']).read_bytes()
        assert sha(data) == row['archive_sha256'] and len(data) == row['archive_bytes']
        raw = gzip.decompress(data)
        assert sha(raw) == row['sha256'] and len(raw) == row['bytes']
        path = ROOT / row['logical_path']
        if path.exists():
            assert path.read_bytes() == raw
        elif args.restore:
            with path.open('xb') as handle:
                handle.write(raw)
    count = 0
    for name, expected in record['historical_manifests'].items():
        manifest = ROOT / name
        assert sha(manifest.read_bytes()) == expected
        rows = json.loads(manifest.read_bytes())['files']
        if isinstance(rows, dict):
            rows = [{'path': key, **value} for key, value in rows.items()]
        for row in rows:
            path = Path(row['path'])
            if not path.is_absolute():
                path = (ROOT if path.parts[0] == 'reports' else manifest.parent) / path
            if path.exists():
                data = path.read_bytes()
            else:
                item = indexed[str(path.relative_to(ROOT))]
                data = gzip.decompress((ROOT / item['archive_path']).read_bytes())
            assert sha(data) == row['sha256'], path
            if 'bytes' in row:
                assert len(data) == row['bytes'], path
            count += 1
    print(json.dumps({'result': 'PASS', 'historical_manifests': len(record['historical_manifests']),
                      'historical_entries': count, 'restorations': len(indexed),
                      'restored_files': args.restore}))


if __name__ == '__main__':
    main()
