import json,pathlib,subprocess,sys
r=pathlib.Path(__file__).resolve().parent;data=r/'synthetic-fallback';data.mkdir(exist_ok=False)
exe=sys.argv[1];targets=[]
for target,ratio,beta in [(0,2,'-1/4'),(1,3,'1/2'),(2,0,'0')]:
    common={'schema':'integrated-rational-prefix.v1','dimension':'13/2','q':1,'grading_proof_identity':'synthetic-grading-only','series_identity':f'synthetic-target-{target}','target_index':target}
    for name,start,end in [('train',0,48),('holdout',48,64)]:
        v=dict(common);v['first_index']=start;v['channels']=[{'id':'formal-normalization','beta':beta,'log_power':0,'coefficients':[str(ratio**n)if ratio else f'1/{(n+2)**4}'for n in range(start,end)]}]
        (data/f'{target}-{name}.json').write_text(json.dumps(v,indent=2)+'\n')
    targets.append({'id':target,'training':str(data/f'{target}-train.json'),'holdout':str(data/f'{target}-holdout.json')})
config=data/'config.json';config.write_text(json.dumps({'schema':'integrated-per-target-fit-plan.v1','targets':targets},indent=2)+'\n')
p=subprocess.run([sys.executable,str(r/'run-target-fits.py'),exe,str(config),str(data/'output')],check=True,capture_output=True,text=True)
s=json.loads((data/'output/summary.json').read_text());assert s['status']=='all_targets_heldout_pass_candidate_only'
a=json.loads((data/'output/target-0/theta-candidate.json').read_text())['candidate'];b=json.loads((data/'output/target-1/theta-candidate.json').read_text())['candidate'];assert a['order']==b['order']==1 and a['coefficients']!=b['coefficients']
f=json.loads((data/'output/all-candidates-freeze.json').read_text());assert len(f['targets'])==3 and f['heldout_access_started']==False
assert [x['stage']for x in s['runs']]==['fit-theta','digest','fit-theta','digest','fit-theta','fit-recurrence','digest','validate','validate','validate']
assert s['selected'][2]['family']=='recurrence' and len(s['selected'][2]['attempts'])==2
assert json.loads((data/'output/target-2/recurrence-candidate.json').read_text())['candidate']['degree']==4
(data/'tests.json').write_text(json.dumps({'checks_passed':5,'checks':['per-target first-order operators','different operators remain separate','both candidates frozen before any holdout','all target holdouts pass only as candidates','theta miss uses declared recurrence training-only fallback'],'scope':'synthetic only','stdout':p.stdout},indent=2)+'\n')
print('PASS5 per-target protocol checks; synthetic only')
