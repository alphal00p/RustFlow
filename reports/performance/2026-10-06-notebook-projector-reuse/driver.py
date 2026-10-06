from pathlib import Path
from time import perf_counter_ns
import hashlib, importlib.util, json, os, shutil, statistics
import symbolica
from symbolica.community.hep import integration as n
root=Path('/common/dev/amflow'); raw=root/'target/notebook-amplitude-performance'
source=root/'target/gg-hg-notebook-acceptance/publication-native30'
controllers={
 'baseline':Path('/common/dev/symbolica-community/loop-integration-publication/examples/hep/gg_hg_support.py'),
 'candidate':Path('/common/dev/symbolica-community/notebook-amplitude-performance/examples/hep/gg_hg_support.py'),
}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
affinity=os.sched_getaffinity(0);cpu=44 if 44 in affinity else min(affinity);os.sched_setaffinity(0,{cpu})
original={x:sha(source/x/'physical-boundaries.bin') for x in ('seeds','transport')}
report={'schema':1,'scope':'Warm forced amplitude assembly with unchanged release extension and supplied native completed banks; no fresh boundaries, no ODE work, no upstream speed ratio.','cpu':cpu,'original_bank_sha256':original,'native_extension_sha256':sha(Path(symbolica.core.__file__)),'controller_sha256':{k:sha(p) for k,p in controllers.items()},'measurements':[]}
sessions={};modules={}
def evidence(s,r):
 return ([(k,v.values,v.absolute_errors,v.verified_relative_digits,v.provenance) for k,v in sorted(s.form_factors.items())],r.values,r.absolute_errors,r.arithmetic_changes,r.verified_relative_digits,r.working_bits,r.provenance)
try:
 for name,path in controllers.items():
  spec=importlib.util.spec_from_file_location('amplitude_'+name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);modules[name]=m
  bank=raw/(name+'-bank');bank.mkdir(exist_ok=False)
  for kind in original:
   (bank/kind).mkdir();shutil.copy2(source/kind/'physical-boundaries.bin',bank/kind/'physical-boundaries.bin')
  sessions[name]=m.CalculationSession(path.parent/'data/gg_hg/native-model.json',bank)
 # Keep one shared immutable native amplitude to isolate warm assembly. Both
 # controllers normally retain this kernel after their first assembly too.
 start=perf_counter_ns();kernel=n.HiggsJetAmplitude(sessions['baseline'].model);report['native_amplitude_preparation_ns']=perf_counter_ns()-start
 for s in sessions.values():s.amplitude=kernel
 for point in ('original','nearby_cached'):
  transport_evidence=None
  for name,s in sessions.items():
   result=s.transport(nearby=point=='nearby_cached')
   assert len(result)==16 and all(r.cache_hit and r.steps==0 and r.inserted_points==0 for r in result.values())
   current={k:(v.identity,v.coordinates,v.coefficients,v.comparison_errors,v.root_sheets,v.provenance) for k,v in result.items()}
   if transport_evidence is None:transport_evidence=current
   else:assert current==transport_evidence
  warmup=[evidence(s,s.assemble(recompute=True)) for s in sessions.values()]
  assert warmup[0]==warmup[1]
  for repetition in range(3):
   expected=None
   order=('baseline','candidate') if repetition%2==0 else ('candidate','baseline')
   for name in order:
    s=sessions[name];start=perf_counter_ns();value=s.assemble(recompute=True);elapsed=perf_counter_ns()-start
    actual=evidence(s,value)
    if expected is None:expected=actual
    else:assert actual==expected
    record={'controller':name,'point':point,'repetition':repetition,'elapsed_ns':elapsed,'form_factors':8,'observables':3,'exact_values_errors_precision_and_provenance_equal':True}
    report['measurements'].append(record);(raw/'amplitude-partial.json').write_text(json.dumps(report,indent=2)+'\n');print(record,flush=True)
 report['median_ns']={point:{name:statistics.median(r['elapsed_ns'] for r in report['measurements'] if r['point']==point and r['controller']==name) for name in controllers} for point in ('original','nearby_cached')}
 report['original_banks_unchanged']=original=={x:sha(source/x/'physical-boundaries.bin') for x in original};assert report['original_banks_unchanged']
 report['status']='passed';(raw/'amplitude-report.json').write_text(json.dumps(report,indent=2)+'\n');print(report['median_ns'],flush=True)
finally:
 for s in sessions.values():s.close()
