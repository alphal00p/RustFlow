"""Losslessly archive completed verbose checkpoints; preserve final proof raw."""
import gzip,hashlib,json,os
from pathlib import Path
BASE=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 assert read(BASE/'run-binding.json')['exit_code']==0
 assert read(BASE/'selected-reference-comparison.json')['independent_reference_comparison_count']==2
 output=BASE/'archive-map.json';assert not output.exists()
 native=BASE/'single3/native-closure'
 preserve={'round-018-closed.json','round-018-closed.bin','round-018-provisional.json','round-018-provisional.bin','requested-discovery-0019.json','retired-unresolved.json'}
 selected=[p for p in native.iterdir()if p.is_file()and p.name not in preserve]+[BASE/'resources.log']
 rows=[]
 for p in sorted(selected):
  target=p.with_name(p.name+'.gz');assert not target.exists()
  rawsha=sha(p);size=p.stat().st_size
  with p.open('rb')as source,target.open('wb')as dest:
   with gzip.GzipFile(filename='',mode='wb',fileobj=dest,mtime=0,compresslevel=6)as compressed:
    while block:=source.read(1024*1024):compressed.write(block)
  with gzip.open(target,'rb')as restored:
   restoredsha=hashlib.file_digest(restored,'sha256').hexdigest()
  assert rawsha==restoredsha
  rows.append({'logical_raw_path':str(p.relative_to(BASE)),'raw_sha256':rawsha,'raw_bytes':size,'gzip_path':str(target.relative_to(BASE)),'gzip_sha256':sha(target),'gzip_bytes':target.stat().st_size,'restored_sha256':restoredsha})
  p.unlink()
 report={'status':'all deterministic gzip roundtrips passed','gzip_mtime':0,'scope':'Completed control only. Final provisional/closed native proof and transaction kept raw; no active readers or production inputs modified. Prior logical-path hash bindings remain verifiable by restoration.','files':rows,'raw_bytes':sum(r['raw_bytes']for r in rows),'gzip_bytes':sum(r['gzip_bytes']for r in rows),'preserved_raw':{str((native/p).relative_to(BASE)):sha(native/p)for p in sorted(preserve)}}
 output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'files':len(rows),'raw_bytes':report['raw_bytes'],'gzip_bytes':report['gzip_bytes']}))
if __name__=='__main__':main()
