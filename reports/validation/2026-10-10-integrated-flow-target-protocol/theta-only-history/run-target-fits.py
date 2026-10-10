#!/usr/bin/env python3
"""Training-only per-target theta fits; all candidates frozen before any holdout."""
import datetime,hashlib,json,pathlib,resource,subprocess,sys,time
if len(sys.argv)!=4:raise SystemExit('usage: run-target-fits.py EXE CONFIG_JSON NEW_OUTPUT_DIR')
exe=pathlib.Path(sys.argv[1]);config_path=pathlib.Path(sys.argv[2]);out=pathlib.Path(sys.argv[3]);out.mkdir(exist_ok=False)
MAX_BYTES=128*1024*1024
start=time.monotonic();runs=[];any_holdout=False

def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,value):
    with path.open('x')as f:json.dump(value,f,indent=2);f.write('\n')
def read(p):
    if p.stat().st_size>MAX_BYTES:raise RuntimeError('metadata/input byte cap')
    return json.loads(p.read_text())
def run(stage,args,folder):
    t=time.monotonic();remaining=300-(t-start)
    if remaining<=0:raise TimeoutError('aggregate fit300s exhausted')
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        p=subprocess.run([str(exe),*map(str,args)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=remaining);log=p.stdout;code=p.returncode
    except subprocess.TimeoutExpired as e:
        log=e.stdout or b'';log=log.decode()if isinstance(log,bytes)else log;code=124
    (folder/(stage+'.log')).write_text(log)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    record={'stage':stage,'args':list(map(str,args)),'exit_code':code,'wall_seconds':time.monotonic()-t,'user_seconds':after.ru_utime-before.ru_utime,'system_seconds':after.ru_stime-before.ru_stime,'cumulative_child_peak_rss_kib':after.ru_maxrss}
    write(folder/(stage+'-resources.json'),record);runs.append(record)
    if code:raise RuntimeError(f'{stage} exit{code}')
    return log.strip()
status='incomplete';selected=[];validation=[]
try:
    cfg=read(config_path)
    if cfg.get('schema')!='integrated-per-target-fit-plan.v1':raise RuntimeError('unexpected config schema')
    targets=cfg.get('targets',[])
    if not 1<=len(targets)<=8:raise RuntimeError('target count outside declared1..8 scope')
    ids=[t['id']for t in targets]
    if any(not isinstance(i,int)or i<0 for i in ids)or len(set(ids))!=len(ids):raise RuntimeError('invalid/duplicate target IDs')
    targets=sorted(targets,key=lambda t:t['id'])
    bindings=[]
    for target in targets:
        train=pathlib.Path(target['training']);meta=read(train)
        if meta.get('target_index')!=target['id']:raise RuntimeError('training target identity mismatch')
        bindings.append({'id':target['id'],'training_path':str(train.resolve()),'training_sha256':sha(train),'series_identity':meta.get('series_identity'),'holdout_path_not_opened':target['holdout']})
    write(out/'input-binding.json',{'captured_utc':now(),'config_sha256':sha(config_path),'launcher_sha256':sha(pathlib.Path(__file__)),'executable':{'path':str(exe.resolve()),'sha256':sha(exe)},'targets':bindings,'primary_family':'theta per target; joint formal channels within target','aggregate_fit_seconds':300,'no_physical_identity_certified':True})
    for target in targets:
        folder=out/f"target-{target['id']}";folder.mkdir()
        train=pathlib.Path(target['training']);candidate=folder/'candidate.json'
        run('fit',['fit',train,candidate,'theta'],folder)
        result=read(candidate)
        digest=run('digest',['digest',candidate],folder)
        if len(digest)!=64 or any(c not in '0123456789abcdef'for c in digest):raise RuntimeError('invalid candidate digest')
        selected.append({'id':target['id'],'candidate':str(candidate.resolve()),'candidate_sha256':sha(candidate),'candidate_blake3':digest,'training_sha256':sha(train),'candidate_present':result['candidate']is not None,'training_status':result['audit']['status']})
    # This freeze lists every selected operator/no-candidate outcome. No heldout file has been accessed.
    write(out/'all-candidates-freeze.json',{'frozen_utc':now(),'targets':selected,'heldout_access_started':False,'structure':'separate theta operators; no cross-target common annihilator; any future companion is block diagonal','selection_from_training_only':True,'identity_certified':False})
    for target,selection in zip(targets,selected):
        if not selection['candidate_present']:
            validation.append({'id':target['id'],'status':selection['training_status'],'heldout_access_started':False});continue
        folder=out/f"target-{target['id']}";holdout=pathlib.Path(target['holdout']);any_holdout=True
        meta=read(holdout)
        if meta.get('target_index')!=target['id']:raise RuntimeError('holdout target identity mismatch')
        run('validate',['validate',selection['candidate'],selection['candidate_blake3'],target['training'],holdout,folder/'holdout-result.json'],folder)
        result=read(folder/'holdout-result.json')
        write(folder/'heldout-binding-after-freeze.json',{'captured_utc':now(),'path':str(holdout.resolve()),'sha256':sha(holdout),'all_candidates_freeze_sha256':sha(out/'all-candidates-freeze.json')})
        validation.append({'id':target['id'],'status':result['status'],'heldout_access_started':True})
    status='all_targets_heldout_pass_candidate_only'if all(v['status']=='heldout_pass_candidate_only'for v in validation)else 'bounded_candidate_attempt_incomplete_or_rejected'
except Exception as exc:
    write(out/'failure.json',{'error':repr(exc),'holdout_access_started':any_holdout})
finally:
    write(out/'summary.json',{'finished_utc':now(),'status':status,'selected':selected,'validation':validation,'aggregate_wall_seconds':time.monotonic()-start,'holdout_access_started':any_holdout,'runs':runs,'identity_certified':False,'no_retry_or_family_change_after_holdout':True})
print(json.dumps({'status':status,'out':str(out),'identity_certified':False}))
if status=='incomplete':raise SystemExit(2)
