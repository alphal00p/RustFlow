"""Evaluate saved formal master identities with frozen native owners only."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
def bind(p):
 p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
def main():
 if len(sys.argv)!=3:raise SystemExit('run-normalizer.py SAVED_COEFFICIENTS OUTPUT_DIRECTORY')
 coefficient=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir()
 build_path=BASE/'build-01/binding.json';build=json.loads(build_path.read_text());assert build['exit_code']==0
 executable=build['executable'];assert bind(executable['path'])==executable
 assert bind(build['source']['path'])==build['source']
 for item in build['dependencies'].values():assert bind(item['path'])==item
 data=json.loads(coefficient.read_text());assert data['schema']=='integrated-gaussian-rational-prefix.v1'
 input_path=ROOT/'examples/finite_density/massless_three_loop_chain.json'
 source_inputs={k:bind(p) for k,p in [('coefficients',coefficient),('input',input_path),('native_build',build_path),('normalizer_source',build['source']['path']),('plan',BASE/'plan.json'),('launcher',Path(__file__))]}
 (out/'coefficients.json').write_bytes(coefficient.read_bytes());(out/'input.json').write_bytes(input_path.read_bytes())
 command=[executable['path'],str(out/'coefficients.json'),str(out/'input.json'),str(out/'normalization.json')]
 record={'scope':'Only saved master channels and source-bound native owners; no reference data read','inputs':source_inputs,'executable':executable,'command':command,'timeout_seconds':60,'producer_known_terms':data['known_prefix_length']}
 (out/'launch.json').write_text(json.dumps(record,indent=2)+'\n')
 start=time.monotonic()
 with (out/'stdout.log').open('wb') as log,(out/'stderr.log').open('wb') as err:
  try:result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=err,timeout=60);code=result.returncode
  except subprocess.TimeoutExpired:code=124
 record.update(exit_code=code,wall_seconds=time.monotonic()-start)
 record['inputs_unchanged']=all(bind(v['path'])==v for v in source_inputs.values())
 record['executable_unchanged']=bind(executable['path'])==executable
 record['dependencies_unchanged']=all(bind(v['path'])==v for v in build['dependencies'].values())
 if (out/'normalization.json').exists():record['output']=bind(out/'normalization.json')
 (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
 assert record['inputs_unchanged'] and record['executable_unchanged'] and record['dependencies_unchanged']
 print(json.dumps({k:record[k] for k in ['exit_code','wall_seconds','producer_known_terms']}));raise SystemExit(code)
if __name__=='__main__':main()
