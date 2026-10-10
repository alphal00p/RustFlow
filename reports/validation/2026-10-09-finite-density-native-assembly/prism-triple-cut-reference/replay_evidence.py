"""Recheck saved exact pole/tail evidence against each archived numerical contour.

No quadrature, AMF, supplied oracle or production source is read or changed.
"""
from pathlib import Path
from fractions import Fraction
import gzip,hashlib,json,subprocess,sys,tempfile,time
BASE=Path(__file__).resolve().parent
PROFILES=['laurent-baseline','laurent-nearby-grid','laurent-fine-baseline','laurent-fine-precision','laurent-fine-quadrature']
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(b):return hashlib.sha256(b).hexdigest()
start=time.monotonic(); records=[]
with tempfile.TemporaryDirectory(prefix='rustflow-prism-evidence-') as tmp:
 tmp=Path(tmp)
 for profile in PROFILES:
  for i in range(1,13):
   name=f'sample-{i:02d}'; directory=BASE/profile; sample=read(directory/(name+'.json')); provenance=read(directory/(name+'-resources.provenance.json'))
   archives=read(directory/(name+'-archives.json')); archive=next(a for a in archives if a['file']==name+'-contours.json.gz'); path=directory/archive['file']; data=gzip.decompress(path.read_bytes()); assert sha(path)==archive['sha256_gzip']; assert digest(data)==archive['sha256_uncompressed']
   contours=json.loads(data); assert Fraction(contours['epsilon'])==Fraction(sample['epsilon'])
   numerical_input=[v for k,v in provenance['argument_files'].items() if Path(k).name==name+'-contours.json']; assert len(numerical_input)==1 and numerical_input[0]['sha256']==digest(data)
   input_path=tmp/'contours.json'; input_path.write_bytes(data); checks={}
   for kind,script in [('regulator','prism_barnes_regulator_analysis.py'),('tails','prism_barnes_tail_check.py')]:
    source=BASE/'generator'/script; output=tmp/(kind+'.json'); saved=directory/(name+'-'+kind+'.json'); assert Fraction(read(saved)['epsilon'])==Fraction(sample['epsilon'])
    completed=subprocess.run([sys.executable,str(source),str(input_path),str(output)],capture_output=True,text=True,timeout=30); assert completed.returncode==0,completed.stdout+completed.stderr; assert read(output)==read(saved)
    checks[kind]={'script':str(source.relative_to(BASE)),'script_sha256':sha(source),'saved_report':str(saved.relative_to(BASE)),'saved_sha256':sha(saved),'replayed_report_equal':True}
   records.append({'profile':profile,'sample':name,'epsilon':sample['epsilon'],'contour_archive':str(path.relative_to(BASE)),'contour_sha256':digest(data),'contour_archive_sha256':sha(path),'numerical_command_input_sha256':numerical_input[0]['sha256'],'checks':checks})
  print('replayed',profile,flush=True)
report={'schema_version':1,'status':'passed','scope':__doc__.strip(),'sample_count':len(records),'exact_reports_replayed':2*len(records),'script_sha256':sha(Path(__file__)),'wall_seconds':time.monotonic()-start,'records':records}
assert len(records)==60
(BASE/'evidence-replay.json').write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps({k:v for k,v in report.items()if k!='records'}))
