"""Compare complete joint predictions only after static producer authentication."""
import argparse,json,subprocess,sys
from pathlib import Path
from joint_provenance import BASE,ROOT,read,digest,verify_joint

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['sample','laurent']);args=parser.parse_args()
 output=BASE/(args.mode+'-comparison.json');binding_path=BASE/(args.mode+'-comparison-binding.json')
 assert not output.exists()and not binding_path.exists(),'preserve previous comparisons'
 native,inputs=verify_joint(args.mode)
 plan_path=BASE/'comparison-plan.json';plan=read(plan_path)
 for key in ['frozen_comparator_and_reference_artifacts_sha256','report_local_metadata_validators_sha256']:
  for name,expected in plan[key].items():
   path=ROOT/name;assert digest(path)==expected,'changed frozen comparison input: '+name;inputs.append(path)
 inputs.extend([plan_path,Path(__file__)])
 spec=plan['gates']['three-loop-'+('fixed'if args.mode=='sample'else'laurent')]
 command=[sys.executable,str(ROOT/spec['comparator']),args.mode,'--predictions',str(BASE/'predictions'),'--reference',str(ROOT/spec['reference']),'--output',str(output),'--run-provenance',str(BASE/'resources.provenance.json'),'--source-snapshot',str(BASE/'evidence/singleton-physical-zero-v2-source-hashes.json')]
 if args.mode=='laurent':command+=['--sample-comparison',str(BASE/'sample-comparison.json')]
 def name(p):
  try:return str(p.relative_to(ROOT))
  except ValueError:return str(p)
 hashes={name(p):digest(p)for p in inputs}
 evidence={'status':'comparison_running','scope':plan['scope'],'command':command,'native_provenance':native,'inputs_sha256':hashes,'definitions_and_tolerances_changed':False,'production_modified':False}
 binding_path.write_text(json.dumps(evidence,indent=2)+'\n')
 result=subprocess.run(command,cwd=ROOT)
 for n,h in hashes.items():assert digest(ROOT/n if not Path(n).is_absolute()else Path(n))==h
 evidence.update(exit_code=result.returncode,status='passed'if result.returncode==0 else'failed')
 if output.exists():
  evidence['comparison_sha256']=digest(output)
  if result.returncode==0:assert read(output)['status']=='passed'
 else:assert result.returncode!=0
 binding_path.write_text(json.dumps(evidence,indent=2)+'\n')
 raise SystemExit(result.returncode)
if __name__=='__main__':main()
