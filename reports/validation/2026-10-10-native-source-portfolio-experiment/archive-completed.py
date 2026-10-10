#!/usr/bin/env python3
"""Losslessly archive completed verbose controls, preserving resume inputs."""
import datetime
import gzip
import hashlib
import json
from pathlib import Path

base = Path(__file__).resolve().parent
sha = lambda data: hashlib.sha256(data).hexdigest()
protected = set()
for arm, checkpoint in [('baseline',7),('proposal',6)]:
    directory = base/'active-resume'/f'polynomial62-{arm}'
    binding = json.loads((directory/'run-binding.json').read_text())
    for name, digest in binding['frozen_parent_and_build_artifact_sha256'].items():
        path = Path(name)
        assert path.is_file() and sha(path.read_bytes()) == digest, name
        if path.is_relative_to(base): protected.add(path)
    protected.update(directory/name for name in ['closed.bin','result.json','configuration.json',
                                                  'progress.json','resume-state.json','run-binding.json',
                                                  'resources.json','resources.provenance.json'])
    # Preserve both original and copied committed parent marker pairs.
    for phase in ['active-pilot','active-resume']:
        for ext in ['json','bin']:
            path = base/phase/f'polynomial62-{arm}'/f'round-{checkpoint:03}.{ext}'
            assert path.is_file()
            protected.add(path)
before = {str(p.relative_to(base)):sha(p.read_bytes()) for p in sorted(protected)}
map_path = base/'archive-map.json'
assert not map_path.exists(), 'preserve previous archive map'
selected = []
for path in sorted(base.rglob('*')):
    if not path.is_file() or path in protected or '__pycache__' in path.parts: continue
    if (path.suffix == '.bin' or (path.name.startswith('round-') and path.suffix == '.json')
            or path.name.endswith('-19-points.json')
            or (path.suffix == '.log' and path.stat().st_size > 65536)):
        selected.append(path)
records = []
for path in selected:
    data = path.read_bytes()
    archived = gzip.compress(data, compresslevel=9, mtime=0)
    target = path.with_name(path.name+'.gz')
    assert not target.exists()
    target.write_bytes(archived)
    assert gzip.decompress(target.read_bytes()) == data
    records.append({'original_path':str(path.relative_to(base)), 'archive_path':str(target.relative_to(base)),
                    'original_bytes':len(data),'original_sha256':sha(data),
                    'archive_bytes':len(archived),'archive_sha256':sha(archived),'byte_roundtrip_verified':True})
    path.unlink()
assert all((base/name).is_file() and sha((base/name).read_bytes()) == digest for name,digest in before.items())
map_path.write_text(json.dumps({'format':'gzip;mtime=0;level=9','scope':'Completed local verbose native programs/results only. Preserve raw closed programs, final result JSON, authenticated parent checkpoint pairs and every hash-bound local resume input. Frozen historical references retain original paths; restore before reuse.',
                               'protected_raw_sha256':before,'files':records},indent=2)+'\n')
files=[]
for path in sorted(base.rglob('*')):
    if path.is_file() and path.name != 'manifest.json' and '__pycache__' not in path.parts:
        data=path.read_bytes();files.append({'path':str(path.relative_to(base)), 'bytes':len(data),'sha256':sha(data)})
(base/'manifest.json').write_text(json.dumps({'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'excluded':['manifest.json','**/__pycache__/**'],'files':files},indent=2)+'\n')
assert all(sha((base/x['path']).read_bytes())==x['sha256'] for x in files)
print(json.dumps({'archived_files':len(records),'original_bytes':sum(x['original_bytes']for x in records),
    'compressed_bytes':sum(x['archive_bytes']for x in records),'protected_raw_files':len(before),
    'manifest_files':len(files),'manifest_sha256':sha((base/'manifest.json').read_bytes())},indent=2))
