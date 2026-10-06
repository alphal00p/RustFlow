import pathlib,subprocess,json,time,datetime,re,hashlib,sys
w=pathlib.Path('/common/dev/amflow-sheet-constrained-endpoints'); t=w/'target/sheet-constraint-validation'
manifest=json.loads((t/'final-source.json').read_text())
def recheck():
 return [p for p,v in manifest['files'].items() if hashlib.sha256((w/p).read_bytes()).hexdigest()!=v['sha256']]
assert not recheck()
targets=['constrained_singular_endpoints','supplied_singular_endpoints','algebraic_endpoints','frobenius','prepared_frobenius','asymptotic','analytic_origin','supplied_flow','physical_transport','transport_cache','prescribed_cache','algebraic_quotients','algebraic_cache','algebraic_kinematic','epsilon_algebraic']
commands=[('final-lib',['cargo','test','--release','--lib','--','--nocapture']),('final-public',['cargo','test','--release']+sum((['--test',x] for x in targets),[])+['--','--nocapture']),('final-fmt',['cargo','fmt','--all','--','--check'])]
records=[]
for label,command in commands:
 started=datetime.datetime.now(datetime.timezone.utc).isoformat(); tic=time.monotonic(); log=t/(label+'.log')
 print('START',label,started,flush=True)
 with log.open('w') as out: result=subprocess.run(command,cwd=w,stdout=out,stderr=subprocess.STDOUT)
 text=log.read_text(); checks={}
 for chunk in re.split(r'(?m)^\s*Running ',text)[1:]:
  name=chunk.split(' (',1)[0].strip(); summaries=re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out; finished in ([^\n]+)',chunk)
  if summaries:
   last=summaries[-1]; checks[name]={'passed':int(last[1]),'failed':int(last[2]),'ignored':int(last[3]),'seconds':last[6],'all_summaries_including_child_processes':summaries}
 record={'label':label,'command':command,'start_utc':started,'elapsed_seconds':time.monotonic()-tic,'exit_code':result.returncode,'checks':checks,'source_recheck_changed':recheck()}
 records.append(record);(t/'final-run-progress.json').write_text(json.dumps(records,indent=2)+'\n');print('END',label,record['exit_code'],record['elapsed_seconds'],flush=True)
 if result.returncode or record['source_recheck_changed']:sys.exit(result.returncode or 100)
(t/'final-run-complete.json').write_text(json.dumps(records,indent=2)+'\n')
