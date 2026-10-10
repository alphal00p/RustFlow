#!/usr/bin/env python3
"""One predeclared theta attempt; does not stat/read/hash holdout before freeze."""
import datetime,hashlib,json,pathlib,resource,subprocess,sys,time
if len(sys.argv)!=5:raise SystemExit('usage: run-frozen-fit.py EXE TRAIN HOLDOUT NEW_OUTPUT_DIR')
exe,train,holdout=map(pathlib.Path,sys.argv[1:4]);out=pathlib.Path(sys.argv[4]);out.mkdir(exist_ok=False)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(name,value):
    with (out/name).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
start=time.monotonic();runs=[]
# Intentionally bind only executable, this launcher and training before candidate selection.
write('input-binding.json',{'started_utc':now(),'primary_family':'theta','selection':'smallest training-only shape (unknowns,theta order,t degree); ambiguity stops','max_total_fit_seconds':300,'executable':{'path':str(exe.resolve()),'sha256':sha(exe)},'launcher_sha256':sha(pathlib.Path(__file__)),'training':{'path':str(train.resolve()),'sha256':sha(train)},'holdout_path_not_opened':str(holdout),'identity_certified':False})
def run(name,args):
    before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic();remaining=300-(t-start)
    if remaining<=0:raise TimeoutError('cumulative fit wall cap before next command')
    timeout=False
    try:p=subprocess.run([str(exe),*map(str,args)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=remaining);output=p.stdout;code=p.returncode
    except subprocess.TimeoutExpired as e:output=(e.stdout or b'');output=output.decode() if isinstance(output,bytes) else output;code=124;timeout=True
    (out/(name+'.log')).write_text(output)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    record={'stage':name,'args':list(map(str,args)),'exit_code':code,'timeout':timeout,'wall_seconds':time.monotonic()-t,'user_seconds':after.ru_utime-before.ru_utime,'system_seconds':after.ru_stime-before.ru_stime,'cumulative_child_peak_rss_kib':after.ru_maxrss}
    runs.append(record);write(name+'-resources.json',record)
    if code:raise RuntimeError(f'{name} exit {code}')
    return output.strip()
status='incomplete';holdout_opened=False
try:
    candidate=out/'candidate.json';run('fit',['fit',train,candidate,'theta'])
    fitted=json.loads(candidate.read_text())
    if fitted['candidate'] is None:status=fitted['audit']['status']
    else:
        digest=run('digest',['digest',candidate])
        if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):raise RuntimeError('invalid candidate BLAKE3')
        write('candidate-freeze.json',{'frozen_utc':now(),'sha256':sha(candidate),'blake3':digest,'candidate_path':str(candidate.resolve()),'training_sha256':sha(train),'heldout_read':False,'selection_from_training_only':True,'identity_certified':False})
        # The immutable candidate freeze exists before any call receives the heldout file.
        holdout_opened=True
        run('validate',['validate',candidate,digest,train,holdout,out/'holdout-result.json'])
        status=json.loads((out/'holdout-result.json').read_text())['status']
        write('heldout-binding-after-validation.json',{'captured_utc':now(),'path':str(holdout.resolve()),'sha256':sha(holdout),'candidate_freeze_sha256':sha(out/'candidate-freeze.json')})
except Exception as exc:
    status='incomplete';write('failure.json',{'error':repr(exc),'holdout_access_started':holdout_opened})
finally:
    write('run-summary.json',{'finished_utc':now(),'status':status,'cumulative_wall_seconds':time.monotonic()-start,'holdout_access_started':holdout_opened,'identity_certified':False,'runs':runs,'no_retry_after_holdout':True,'algebra_caps':'200000 retained monomials / 128MiB canonical payload; not CAS scratch'})
print(json.dumps({'status':status,'report':str(out),'identity_certified':False}))
if status=='incomplete':raise SystemExit(2)
