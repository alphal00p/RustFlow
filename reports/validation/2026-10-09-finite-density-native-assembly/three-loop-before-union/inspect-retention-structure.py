"""Read-only metadata coverage; this does not replace native rule replay."""
import gzip, hashlib, importlib.util, json, pathlib, struct
ROOT=pathlib.Path('/common/dev/rustflow_fermi/reports/validation/2026-10-09-finite-density-native-assembly')
spec=importlib.util.spec_from_file_location('guarded_structure',ROOT/'E7-origin-native-coverage-probe/inspect-generated-structure.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
base=ROOT/'three-loop-before-union/native-closure/cut-0'
def read(path):
 if path.exists():return path.read_bytes()
 return gzip.decompress(path.with_name(path.name+'.gz').read_bytes())
def parse(path):
 raw=read(path);assert raw[:8]==b'RRPBIN\r\n';n=struct.unpack_from('<I',raw,16)[0];p=20;data=None
 for _ in range(n):
  tag,res,length=struct.unpack_from('<HHQ',raw,p);p+=12
  if tag==4:data=raw[p:p+length]
  p+=length
 assert p==len(raw) and data is not None
 r=helper.Reader(data)
 x={'schema':r.string(),'measure':r.string(),'roles':r.vec(r.byte),'indices':r.vec(r.u),'sources':r.vec(r.source),'zero_domains':r.vec(r.domain),'rules':r.vec(r.rule),'terminals':r.vec(lambda:r.vec(r.i))}
 assert r.p==len(data)
 return x,hashlib.sha256(raw).hexdigest()
frontier=json.loads(read(base/'round-003-provisional.json'))['frontier']
rounds=[];matches=[]
for round in range(4):
 path=base/f'round-{round:03d}-provisional.bin';program,digest=parse(path)
 covered=set();first=[]
 for ordinal,rule in enumerate(program['rules']):
  for index,point in enumerate(frontier):
   if helper.contains(rule['domain'],point) and helper.targetmatch(rule['target'],point):
    covered.add(index)
    if index==0:first.append({'rule':ordinal,'domain':rule['domain'],'target':rule['target'],'sources':rule['sources']})
 matches.append(covered)
 rounds.append({'round':round,'program':str(path.relative_to(ROOT)),'raw_sha256':digest,'rules':len(program['rules']),'sources':len(program['sources']),'terminals':len(program['terminals']),'latest_frontier_structurally_matched':len(covered),'first_latest_frontier_matches':first})
old=set.union(*matches[:-1]);last=matches[-1];lost=old-last
result={'scope':'Read-only structure of saved native programs. A matching domain and target pattern is only a candidate; coefficients, conditions, descent and proofs are not evaluated here. This is not a native replay or new coverage certificate.','latest_frontier_count':len(frontier),'first_latest_frontier':frontier[0],'rounds':rounds,'matched_by_any_previous_program':len(old),'matched_by_latest_program':len(last),'matched_by_previous_but_not_latest':len(lost),'previous_only_example_labels':[frontier[i] for i in sorted(lost)[:12]]}
out=ROOT/'three-loop-before-union/retention-structure.json';out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['rounds','scope']}))
for row in rounds:print('round',row['round'],'rules',row['rules'],'latest-frontier matches',row['latest_frontier_structurally_matched'],'first',len(row['first_latest_frontier_matches']))
