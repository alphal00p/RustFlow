from pathlib import Path
from time import perf_counter_ns
import cProfile, hashlib, importlib.util, json, os, pstats, shutil
import symbolica
from symbolica.community.hep import integration as n
raw=Path('/common/dev/amflow/target/notebook-amplitude-performance')
source=Path('/common/dev/amflow/target/gg-hg-notebook-acceptance/publication-native30')
controller=Path('/common/dev/symbolica-community/notebook-amplitude-performance/examples/hep/gg_hg_support.py')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
affinity=os.sched_getaffinity(0); cpu=44 if 44 in affinity else min(affinity);os.sched_setaffinity(0,{cpu})
original={x:sha(source/x/'physical-boundaries.bin') for x in ('seeds','transport')}
copy=raw/'profile-bank-completed';copy.mkdir(exist_ok=False)
for name in original:
 (copy/name).mkdir();shutil.copy2(source/name/'physical-boundaries.bin',copy/name/'physical-boundaries.bin')
spec=importlib.util.spec_from_file_location('amplitude_profile_controller',controller);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
records=[]
def constructor(name):
 original=getattr(n,name)
 def traced(*args,**kwargs):
  start=perf_counter_ns();value=original(*args,**kwargs);record={'operation':name,'elapsed_ns':perf_counter_ns()-start};records.append(record);print(record,flush=True);return value
 setattr(n,name,traced)
 return original
constructors={name:constructor(name) for name in ('HiggsJetAmplitude','HiggsJetFormFactorProjector')}
session=m.CalculationSession(controller.parent/'data/gg_hg/native-model.json',copy)
try:
 start=perf_counter_ns(); session.transport();records.append({'operation':'exact_transport','elapsed_ns':perf_counter_ns()-start})
 for i in range(3):
  profile=cProfile.Profile();start=perf_counter_ns();profile.enable();result=session.assemble(recompute=True);profile.disable()
  row={'operation':'assemble','repetition':i,'elapsed_ns':perf_counter_ns()-start};records.append(row);print(row,flush=True)
  profile.dump_stats(str(raw/f'assemble-{i}.prof'))
  with (raw/f'assemble-{i}.txt').open('w') as f:pstats.Stats(profile,stream=f).sort_stats('cumulative').print_stats(35)
 assert original=={x:sha(source/x/'physical-boundaries.bin') for x in original}
 (raw/'profile.json').write_text(json.dumps({'controller_sha256':sha(controller),'native_extension_sha256':sha(Path(symbolica.core.__file__)),'cpu':cpu,'records':records,'original_banks_unchanged':True},indent=2)+'\n')
finally:
 session.close()
 for name,value in constructors.items():setattr(n,name,value)
