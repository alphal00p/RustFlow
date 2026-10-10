#!/usr/bin/env python3
"""Add a lossless storage layer without changing any historical manifest bytes."""
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPORTS = ROOT / 'reports/validation'
ROOT_NAMES = [
    'targeted-native-depth-controls', 'targeted-portfolio-controls',
    'mixed-energy-polynomial-source-experiment',
    'singleton-energy-ward-first-controls', 'singleton-energy-ward-native-controls',
    'factorized-compact-ward-assessment', 'global-euler-source-experiment',
    'occupied-coordinate-control', 'three-loop-single0-flow-check',
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    assert not path.exists(), path
    path.write_text(json.dumps(value, indent=2) + '\n')


def entries(manifest):
    rows = json.loads(manifest.read_bytes())['files']
    if isinstance(rows, dict):
        rows = [{'path': name, **value} for name, value in rows.items()]
    for row in rows:
        path = Path(row['path'])
        if not path.is_absolute():
            path = (ROOT if path.parts[0] == 'reports' else manifest.parent) / path
        yield path, row


def main():
    assert not (HERE / 'archive-map.json').exists(), 'preserve the completed storage record'
    roots = [REPORTS / ('2026-10-10-' + name) for name in ROOT_NAMES]
    manifests = []
    for directory in roots:
        candidates = [p for p in [directory / 'manifest.json', directory / 'artifact-manifest.json'] if p.exists()]
        assert len(candidates) == 1, (directory, candidates)
        manifests.extend(candidates)
    before = {str(m.relative_to(ROOT)): sha(m.read_bytes()) for m in manifests}
    manifest_entries = 0
    for manifest in manifests:
        for path, row in entries(manifest):
            data = path.read_bytes()
            assert sha(data) == row['sha256'], (manifest, path, 'pre-archive hash')
            if 'bytes' in row:
                assert len(data) == row['bytes'], (manifest, path, 'pre-archive size')
            manifest_entries += 1
    # These report owners explicitly released their readers. The allowlist excludes
    # all active roots, all native binaries, and the latest larger-frontier pairs.
    candidates = sorted(p for d in roots for p in d.rglob('*.json')
                        if p.stat().st_size > 256 * 1024
                        and 'manifest' not in p.name and 'archive' not in p.name)
    rows = []
    for path in candidates:
        data = path.read_bytes()
        json.loads(data)  # Every selected file is an intact JSON document.
        archive = path.with_suffix(path.suffix + '.gz')
        assert not archive.exists(), archive
        compressed = gzip.compress(data, compresslevel=9, mtime=0)
        assert gzip.decompress(compressed) == data
        archive.write_bytes(compressed)
        assert gzip.decompress(archive.read_bytes()) == data
        rows.append({'logical_path': str(path.relative_to(ROOT)),
                     'sha256': sha(data), 'bytes': len(data),
                     'archive_path': str(archive.relative_to(ROOT)),
                     'archive_sha256': sha(compressed), 'archive_bytes': len(compressed),
                     'gzip_mtime': 0})
        path.unlink()
    record = {'schema': 'additive-completed-evidence-storage-v1',
              'scope': 'Storage only. Historical source, proof, numerical and report content is byte-exact after decompression. Historical manifests are unchanged.',
              'reader_release': 'Completed owners confirmed no readers. Active longer-singleton, partial-placement and native-program-work-profile roots are excluded. Native binary corpora and latest larger-frontier checkpoint pairs remain raw.',
              'historical_manifests': before, 'entries': rows,
              'raw_bytes': sum(r['bytes'] for r in rows),
              'archive_bytes': sum(r['archive_bytes'] for r in rows)}
    save(HERE / 'archive-map.json', record)
    by_path = {r['logical_path']: r for r in rows}
    verified = 0
    for manifest in manifests:
        assert sha(manifest.read_bytes()) == before[str(manifest.relative_to(ROOT))]
        for path, old in entries(manifest):
            if path.exists():
                data = path.read_bytes()
            else:
                row = by_path[str(path.relative_to(ROOT))]
                compressed = (ROOT / row['archive_path']).read_bytes()
                assert sha(compressed) == row['archive_sha256']
                data = gzip.decompress(compressed)
                assert sha(data) == row['sha256'] and len(data) == row['bytes']
            assert sha(data) == old['sha256'], (path, 'post-archive logical hash')
            if 'bytes' in old:
                assert len(data) == old['bytes']
            verified += 1
    assert verified == manifest_entries
    save(HERE / 'verification.json', {
        'result': 'PASS', 'historical_manifests_unchanged': len(manifests),
        'historical_manifest_entries_verified': verified,
        'archive_restorations_verified': len(rows),
        'raw_bytes': record['raw_bytes'], 'archive_bytes': record['archive_bytes'],
        'archive_map_sha256': sha((HERE / 'archive-map.json').read_bytes()),
    })
    print(json.dumps({'status': 'PASS', 'archives': len(rows), 'manifest_entries': verified,
                      'raw_bytes': record['raw_bytes'], 'archive_bytes': record['archive_bytes']}))


if __name__ == '__main__':
    main()
