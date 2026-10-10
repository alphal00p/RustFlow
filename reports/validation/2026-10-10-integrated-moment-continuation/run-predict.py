from pathlib import Path
import hashlib,json,os,subprocess,sys,time
B=Path(__file__).resolve().parent
def bind(p):p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
assert not (B/'prediction-run.json').exists()
plan=json.loads((B/'plan.json').read_text());inputs={n:bind(B/n) for n in ['plan.json','predict.py','run-predict.py']};command=[sys.executable,str(B/'predict.py')];start=time.monotonic()
record={'command':command,'timeout_seconds':plan['limits']['numerical_prediction_timeout_seconds'],'inputs':inputs,'python':bind(sys.executable),'reference_values_passed':False}
(B/'prediction-launch.json').write_text(json.dumps(record,indent=2)+'\n')
with (B/'prediction.log').open('wb') as log:
 try:result=subprocess.run(command,env=os.environ,stdout=log,stderr=subprocess.STDOUT,timeout=record['timeout_seconds']);code=result.returncode
 except subprocess.TimeoutExpired:code=124
record.update(exit_code=code,wall_seconds=time.monotonic()-start,inputs_unchanged=all(bind(v['path'])==v for v in inputs.values()))
if (B/'predictions.json').exists():record['prediction']=bind(B/'predictions.json')
(B/'prediction-run.json').write_text(json.dumps(record,indent=2)+'\n');assert record['inputs_unchanged'];print(json.dumps(record));raise SystemExit(code)
