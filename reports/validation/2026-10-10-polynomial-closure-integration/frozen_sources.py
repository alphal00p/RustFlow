"""Source inventory/verification only; importing this module captures nothing."""
import hashlib,json,subprocess
from pathlib import Path
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PREDECESSOR=ROOT/'reports/validation/2026-10-10-native-zero-domain-projection/guard-point-source-hashes.json'

def digest(path):
    before=path.stat();h=hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda:source.read(1024*1024),b''):h.update(block)
    after=path.stat()
    if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino):
        raise ValueError('source changed while hashing: '+str(path))
    return h.hexdigest()

def inventory(extra=()):
    previous=json.loads(PREDECESSOR.read_text())['files']
    paths={name for name in previous if (ROOT/name).is_file()}
    found=subprocess.check_output(['rg','--files','src','tests','build_support','examples/finite_density','vendor/rustred'],cwd=ROOT,text=True).splitlines()
    for name in found:
        p=Path(name)
        if p.suffix=='.rs' or p.name=='Cargo.toml' or (name.startswith('tests/') and p.suffix=='.py') or (name.startswith('examples/finite_density/') and p.suffix=='.json'):
            paths.add(name)
    for name in extra:
        p=Path(name)
        if p.is_absolute() or '..' in p.parts or not (ROOT/p).is_file():raise ValueError('invalid additional source: '+name)
        paths.add(name)
    return sorted(paths)

def verify_snapshot(path):
    snapshot=json.loads(path.read_text())
    current=inventory(snapshot.get('additional_sources',[]))
    if set(current)!=set(snapshot['files']):raise ValueError('source inventory changed after freeze')
    changed=[name for name in current if digest(ROOT/name)!=snapshot['files'][name]]
    if changed:raise ValueError('source content changed after freeze: '+repr(changed))
    return snapshot
