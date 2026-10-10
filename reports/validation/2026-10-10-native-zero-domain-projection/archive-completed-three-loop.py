"""Losslessly archive the explicitly released completed full-three-loop corpus.

Scope is fixed: this script cannot archive the massless diagnostic programs.
Run only after every reader has released this exact corpus and sibling log.
"""
import gzip,hashlib,json,os,shutil
from pathlib import Path
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
manifest=BASE/'zero-projection-three-loop-archives.json'
if manifest.exists():raise ValueError('preserve previous archive manifest')
paths=sorted(p for p in (BASE/'zero-projection-three-loop-fixed/native-closure').rglob('*') if p.is_file())
paths.append(BASE/'zero-projection-three-loop-fixed-resources.log')
if not paths or any(p.suffix=='.gz' for p in paths):raise ValueError('unexpected archive scope')
report={'schema':1,'scope':'Completed failed full-three-loop native checkpoints and sibling raw log only. Inputs/configuration/failure/resources/bindings remain plain. Massless programs remain raw for diagnosis.','reader_release':'Parent, guarded dependency owner, and boundary owner explicitly confirmed no remaining readers.','compression':'gzip level9, mtime0, empty embedded filename','complete':False,'records':[]}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def save():
 tmp=manifest.with_suffix('.json.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(manifest)
for path in paths:
 before=path.stat();original=sha(path);destination=Path(str(path)+'.gz');temporary=Path(str(destination)+'.tmp')
 if destination.exists() or temporary.exists():raise ValueError('archive destination exists')
 with path.open('rb') as source,temporary.open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=9,mtime=0) as compressed:shutil.copyfileobj(source,compressed,1024*1024)
 check=hashlib.sha256();size=0
 with gzip.open(temporary,'rb') as source:
  for block in iter(lambda:source.read(1024*1024),b''):check.update(block);size+=len(block)
 after=path.stat()
 assert (before.st_ino,before.st_size,before.st_mtime_ns)==(after.st_ino,after.st_size,after.st_mtime_ns)
 assert check.hexdigest()==original and size==before.st_size
 temporary.replace(destination)
 report['records'].append({'original_path':str(path.relative_to(ROOT)),'original_sha256':original,'original_bytes':before.st_size,'archive_path':str(destination.relative_to(ROOT)),'archive_sha256':sha(destination),'archive_bytes':destination.stat().st_size,'decompressed_sha256_verified':True})
 save();path.unlink()
report['complete']=True;report['file_count']=len(paths);report['original_bytes']=sum(x['original_bytes'] for x in report['records']);report['archive_bytes']=sum(x['archive_bytes'] for x in report['records']);report['script_sha256']=sha(Path(__file__));save()
print(json.dumps({k:report[k] for k in ['complete','file_count','original_bytes','archive_bytes']}))
