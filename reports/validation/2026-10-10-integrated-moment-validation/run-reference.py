"""Budgeted reference-only driver. Requires frozen integrated predictions first.

A release JSON binds the actual conversion audit and saved integrated coefficient
prediction artifacts. This driver hashes prediction bytes, but never parses their
values or passes their paths/content to the quadrature child.
"""
import argparse, datetime, fcntl, hashlib, json, os
from pathlib import Path
import subprocess, sys, time

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PYTHON=Path('/nix/store/0bq2kpdgy8kimy82pzpjjw787h1sx6a9-python3-3.14.6-env/bin/python3')
MPMATH=Path('/nix/store/vm12q3zvz6yp0wklj2b76smgdknnc7wz-python3.13-mpmath-1.3.0/lib/python3.13/site-packages')
BUDGET_SECONDS=600
PER_RUN_SECONDS=180

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,data):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)
def bound(record):
 p=Path(record['path']).resolve();assert p.is_file(),p
 assert sha(p)==record['sha256'],p
 return p

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--release',required=True,type=Path)
 p.add_argument('--eta',required=True,choices=['8','2','0.5'])
 p.add_argument('--nodes',required=True,type=int,choices=[16,24,32])
 p.add_argument('--digits',required=True,type=int,choices=[50,80])
 p.add_argument('--tag',required=True)
 a=p.parse_args()
 assert a.tag and all(c.isalnum() or c in '-_' for c in a.tag)
 release=json.loads(a.release.read_text())
 assert release['schema']==1 and release['status']=='frozen producer integrated predictions; reference generation released'
 assert release['dimension']=='6.5' and release['mu']=='1'
 assert release['predictions_are_integrated_coefficients'] is True
 assert release['equation_may_still_be_pending'] is True
 assert a.eta in release['authorized_reference_eta']
 audit_path=bound(release['native_conversion_audit']);audit=json.loads(audit_path.read_text())
 assert audit['status']=='passed' and audit['special_case_N_eta_equals_N'] is True
 metadata=bound(release['producer_metadata'])
 assert sha(metadata)==audit['metadata_sha256']
 actual_input=bound(release['input'])
 assert sha(actual_input)==audit['input_sha256']
 assert actual_input==ROOT/'examples/finite_density/massless_three_loop_chain.json'
 predictions=release['producer_integrated_predictions']
 assert predictions and all(r['kind']=='integrated moment coefficients' for r in predictions)
 for record in predictions:bound(record)
 for record in release['producer_source_bindings']:bound(record)
 assert release['producer_source_bindings']
 amendment=bound(release['prospective_validation_amendment'])
 amendment_data=json.loads(amendment.read_text())
 assert amendment_data['reference_eta']==[8,2,'1/2']
 assert amendment_data['reference_runtime_seconds']=={'one_run':180,'cumulative':600}
 evaluator=BASE/'finite_eta_reference.py'
 run_dir=BASE/'runs'/a.tag
 assert not run_dir.exists(),'Never overwrite a completed or failed attempt'
 ledger_path=BASE/'reference-budget.json'
 lock_path=BASE/'reference-budget.lock'
 with lock_path.open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  ledger=json.loads(ledger_path.read_text()) if ledger_path.exists() else {'schema':1,'maximum_cumulative_seconds':600,'maximum_run_seconds':180,'runs':[]}
  assert ledger['maximum_cumulative_seconds']==BUDGET_SECONDS
  charged=sum(row.get('charged_seconds',row['reservation_seconds']) for row in ledger['runs'])
  remaining=BUDGET_SECONDS-charged
  # Reserve timeout plus one second of termination/accounting overhead. An
  # interrupted driver keeps its reservation charged rather than granting time.
  timeout=min(PER_RUN_SECONDS-1,remaining-1)
  assert timeout>=10,'Insufficient pre-authorized numerical reference budget'
  run_dir.mkdir(parents=True)
  profile={'dimension':'6.5','mu':'1','eta':a.eta,'nodes':a.nodes,'digits':a.digits}
  receipt={'schema':1,'created_utc':now(),'prediction_artifacts_frozen_before_reference':True,
    'profile':profile,'release_sha256':sha(a.release),'reference_evaluator_sha256':sha(evaluator),
    'native_conversion_audit_sha256':sha(audit_path),'producer_metadata_sha256':sha(metadata),
    'prediction_sha256':[r['sha256'] for r in predictions],
    'values_or_coefficients_passed_to_reference_child':False,'timeout_seconds':timeout}
  write(run_dir/'receipt.json',receipt)
  record={'tag':a.tag,'started_utc':now(),'status':'running','reservation_seconds':timeout+1,
          'receipt_sha256':sha(run_dir/'receipt.json'),'profile':profile}
  ledger['runs'].append(record);write(ledger_path,ledger)
  # Serial work is intentional: one ledger lock remains held for the whole run.
  command=[str(PYTHON),str(evaluator),'--receipt',str(run_dir/'receipt.json'),'--output',str(run_dir/'reference.json')]
  env=dict(os.environ);env['PYTHONPATH']=str(MPMATH);env['PYTHONDONTWRITEBYTECODE']='1'
  provenance={'command':command,'cwd':str(ROOT),'python_sha256':sha(PYTHON.resolve()),
   'mpmath_init_sha256':sha(MPMATH/'mpmath/__init__.py'),'driver_sha256':sha(__file__),
   'evaluator_sha256':sha(evaluator),'release':str(a.release.resolve()),'release_sha256':sha(a.release),
   'profile':profile,'scope':'Reference-only calculation after producer artifact freeze; no producer invocation or equation input',
   'producer_bindings':release,'timeout_seconds':timeout,'cumulative_charged_before_seconds':charged}
  write(run_dir/'provenance.json',provenance)
  started=time.monotonic()
  with (run_dir/'stdout.log').open('wb') as out,(run_dir/'stderr.log').open('wb') as err:
   try:
    process=subprocess.run(command,cwd=ROOT,env=env,stdout=out,stderr=err,timeout=timeout)
    code=process.returncode;status='completed' if code==0 else 'failed'
   except subprocess.TimeoutExpired:
    code=124;status='timed out'
   except BaseException as exc:
    record.update(status='driver interrupted; full reservation remains charged',error=repr(exc))
    write(ledger_path,ledger);raise
  elapsed=time.monotonic()-started
  record.update(status=status,exit_code=code,elapsed_seconds=elapsed,charged_seconds=elapsed,finished_utc=now())
  # Source/prediction bindings are checked again; changes fail the gate, never
  # silently redefine its saved reference.
  unchanged=True
  for item in [release['native_conversion_audit'],release['producer_metadata'],release['input'],*predictions,*release['producer_source_bindings'],release['prospective_validation_amendment']]:
   unchanged=unchanged and sha(item['path'])==item['sha256']
  unchanged=unchanged and sha(evaluator)==provenance['evaluator_sha256'] and sha(__file__)==provenance['driver_sha256']
  record['within_individual_budget']=elapsed<=PER_RUN_SECONDS
  record['bound_inputs_unchanged']=unchanged
  if (run_dir/'reference.json').exists():record['reference_sha256']=sha(run_dir/'reference.json')
  write(run_dir/'resources.json',record);write(ledger_path,ledger)
  assert unchanged,'A frozen input changed during reference evaluation'
  assert sum(r.get('charged_seconds',r['reservation_seconds']) for r in ledger['runs'])<=BUDGET_SECONDS
  raise SystemExit(code)

if __name__=='__main__':main()
