from pathlib import Path
from time import perf_counter_ns
import hashlib, importlib.util, json, os, shutil, statistics, sys
import symbolica
from symbolica.community.hep.integration import BoundaryCache

root=Path('/common/dev/amflow'); out=root/'target/notebook-cache-hit-benchmark'
original=root/'target/gg-hg-notebook-acceptance/publication-native30'
controllers={
 'baseline':Path('/common/dev/symbolica-community/loop-integration-publication/examples/hep/gg_hg_support.py'),
 'candidate':Path('/common/dev/symbolica-community/notebook-cache-hit-persistence/examples/hep/gg_hg_support.py'),
}
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
files={name: original/name/'physical-boundaries.bin' for name in ('seeds','transport')}
before={name:sha(path) for name,path in files.items()}
os.sched_setaffinity(0,{44})
report={'schema':1,'scope':'Matched notebook repeated exact cache-hit stage, unchanged native publication wheel; no ODE work or upstream timing.',
 'python':sys.executable,'cpu_affinity':sorted(os.sched_getaffinity(0)),
 'original_bank_sha256':before,'controller_sha256':{name:sha(path) for name,path in controllers.items()},'measurements':[]}
core=Path(symbolica.core.__file__)
report['native_extension']={'path':str(core),'sha256':sha(core)}
sessions={}; modules={}; expected=None
save_calls=[]; original_save=BoundaryCache.save

def traced_save(cache,directory):
 save_calls.append(str(directory))
 return original_save(cache,directory)
BoundaryCache.save=traced_save

def result_evidence(results):
 return {name:(r.identity,r.leading_power,r.verified_digits,r.input_verified_digits,r.provenance,
              r.coordinates,r.starting_coordinates,r.root_sheets,r.coefficients,r.comparison_errors)
         for name,r in results.items()}
try:
 for name,path in controllers.items():
  spec=importlib.util.spec_from_file_location('cache_hit_'+name,path)
  module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); modules[name]=module
  directory=out/name; directory.mkdir(exist_ok=False)
  for bank,source in files.items():
   (directory/bank).mkdir(); shutil.copy2(source,directory/bank/source.name)
  start=perf_counter_ns()
  sessions[name]=module.CalculationSession(path.parent/'data/gg_hg/native-model.json',directory)
  report.setdefault('initialization_ns',{})[name]=perf_counter_ns()-start
  report.setdefault('cache_entries',{})[name]=len(sessions[name].cache)
  print('initialized',name,'entries',len(sessions[name].cache),flush=True)
 baseline_evidence=modules['baseline'].boundary_evidence(sessions['baseline'].cache)
 assert baseline_evidence==modules['candidate'].boundary_evidence(sessions['candidate'].cache)
 for repetition in range(3):
  order=('baseline','candidate') if repetition%2==0 else ('candidate','baseline')
  for name in order:
   start_calls=len(save_calls); start=perf_counter_ns(); result=sessions[name].transport(); elapsed=perf_counter_ns()-start
   assert len(result)==16
   assert all(r.cache_hit and r.steps==0 and r.inserted_points==0 for r in result.values())
   evidence=result_evidence(result)
   if expected is None: expected=evidence
   else: assert evidence==expected
   count=sum(len(row) for r in result.values() for row in r.coefficients)
   assert count==4360
   record={'controller':name,'repetition':repetition,'elapsed_ns':elapsed,'save_calls':len(save_calls)-start_calls,
           'cache_hits':16,'ode_steps':0,'new_points':0,'coefficients_equal':count,'provenance_and_error_evidence_equal':True}
   report['measurements'].append(record)
   (out/'partial.json').write_text(json.dumps(report,indent=2)+'\n')
   print(record,flush=True)
 for name,session in sessions.items():
  restored=BoundaryCache.load(out/name/'transport')
  assert modules[name].boundary_evidence(restored)==baseline_evidence
 report['native_values_and_evidence_preserved_after_reload']=True
 report['original_banks_unchanged']={name:sha(path) for name,path in files.items()}==before
 assert report['original_banks_unchanged']
 report['median_elapsed_ns']={name:statistics.median(r['elapsed_ns'] for r in report['measurements'] if r['controller']==name) for name in controllers}
 report['status']='passed'
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(report['median_elapsed_ns'],flush=True)
finally:
 BoundaryCache.save=original_save
 for session in sessions.values(): session.close()
