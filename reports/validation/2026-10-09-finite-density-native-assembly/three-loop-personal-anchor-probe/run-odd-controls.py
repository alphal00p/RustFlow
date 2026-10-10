from pathlib import Path
import subprocess,sys,json,concurrent.futures
p=Path(__file__).resolve().parent
runner=p.parent/'run-resource-command.py'
def run(item):
    power,cap,hints,mode=item
    target=[1,power,2,1,0,-1,0,0,0,0,0,0,0,0,0,0]
    name=f'odd-{power}-{mode}-budget-{cap:02d}-anchor-{str(hints).lower()}'
    command=[sys.executable,str(runner),str(p/f'{name}-resources.json'),'timeout','60s','env',f'PROBE_DOMAINS={cap}','PROBE_DEPTH=3',f'PROBE_PERSONAL_ANCHOR={str(hints).lower()}',f'PROBE_TARGET={json.dumps(target)}','/tmp/three_loop_personal_anchor_flexible_probe','/tmp/three-loop-double-bc3c6bd.bin',mode,str(p/f'{name}.json')]
    result=subprocess.run(command,check=False)
    return {'case':name,'runner_exit':result.returncode}
cases=[(power,n,h,'target-rays') for power in [2,3] for n in [1,2,4,8,16,32] for h in [False,True]]+[(power,32,False,'point') for power in [2,3]]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    completed=list(pool.map(run,cases))
(p/'odd-runner-results.json').write_text(json.dumps(completed,indent=2)+'\n')
print(json.dumps(completed,indent=2))
