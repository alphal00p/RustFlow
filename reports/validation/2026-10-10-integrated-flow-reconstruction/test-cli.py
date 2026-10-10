#!/usr/bin/env python3
import json,pathlib,subprocess,sys
r=pathlib.Path(__file__).resolve().parent;folder=r/(sys.argv[2] if len(sys.argv)>2 else 'synthetic-cli');folder.mkdir(exist_ok=False)
common={'schema':'integrated-rational-prefix.v1','dimension':'13/2','q':1,'grading_proof_identity':'synthetic-exact-grading','series_identity':'synthetic-geometric-no-physics'}
def data(begin,end,bad=False):
 v=dict(common);v['first_index']=begin;seq=[str(2**n+(1 if bad and n==48 else 0)) for n in range(begin,end)];v['channels']=[{'id':'unit-normalization','beta':'-1/4','log_power':0,'coefficients':seq}];return v
for name,v in [('train',data(0,48)),('good',data(48,64)),('bad',data(48,64,True))]:(folder/(name+'.json')).write_text(json.dumps(v,indent=2)+'\n')
results=[]
for case in ['good','bad']:
 cmd=[sys.executable,str(r/'run-frozen-fit.py'),sys.argv[1],str(folder/'train.json'),str(folder/(case+'.json')),str(folder/case)]
 p=subprocess.run(cmd,check=True,capture_output=True,text=True)
 summary=json.loads((folder/case/'run-summary.json').read_text());results.append({'case':case,'summary':summary,'stdout':p.stdout})
 assert summary['status']==('heldout_pass_candidate_only' if case=='good' else 'heldout_rejected_protocol_stops')
 assert not summary['identity_certified']
 freeze=json.loads((folder/case/'candidate-freeze.json').read_text());assert freeze['heldout_read']==False
# Candidate selection cannot depend on heldout content.
assert (folder/'good/candidate.json').read_bytes()==(folder/'bad/candidate.json').read_bytes()
(folder/'result.json').write_text(json.dumps({'passed':3,'checks':['good holdout','false holdout rejected','identical frozen candidates'],'identity_certified':False,'results':results},indent=2)+'\n')
print('PASS 3 CLI checks; synthetic data only')
