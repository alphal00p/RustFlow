from pathlib import Path
import subprocess,sys,json,concurrent.futures
p=Path(__file__).resolve().parent
runner=p.parent/'run-resource-command.py'
def run(item):
    cap,hints=item
    name=f'budget-{cap:02d}-anchor-{str(hints).lower()}'
    command=[sys.executable,str(runner),str(p/f'{name}-resources.json'),'timeout','60s','env',f'PROBE_DOMAINS={cap}','PROBE_DEPTH=3',f'PROBE_PERSONAL_ANCHOR={str(hints).lower()}','/tmp/three_loop_personal_anchor_probe','/tmp/three-loop-double-bc3c6bd.bin','target-rays',str(p/f'{name}.json')]
    result=subprocess.run(command,check=False)
    return {'case':name,'runner_exit':result.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    completed=list(pool.map(run,[(n,h) for n in [1,2,4,8,16,32] for h in [False,True]]))
(p/'runner-results.json').write_text(json.dumps(completed,indent=2)+'\n')
print(json.dumps(completed,indent=2))
