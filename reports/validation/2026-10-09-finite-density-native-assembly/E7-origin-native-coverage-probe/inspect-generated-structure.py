import gzip,json,struct,pathlib,hashlib
ROOT=pathlib.Path('/common/dev/rustflow_fermi/reports/validation/2026-10-09-finite-density-native-assembly')
class Reader:
 def __init__(self,b): self.b=b;self.p=0
 def take(self,n): r=self.b[self.p:self.p+n];assert len(r)==n;self.p+=n;return r
 def u(self):
  n=self.take(1)[0]
  if n<251:return n
  assert n in (251,252,253,254)
  return int.from_bytes(self.take({251:2,252:4,253:8,254:16}[n]),'little')
 def i(self):
  n=self.u();return (n>>1)^-(n&1)
 def byte(self):return self.take(1)[0]
 def opt(self,f):
  t=self.byte();assert t in(0,1);return f() if t else None
 def vec(self,f):return [f() for _ in range(self.u())]
 def string(self):return self.take(self.u()).decode()
 def domain(self):return self.vec(lambda:[self.opt(self.i),self.opt(self.i)])
 def integral(self):return self.vec(lambda:[self.byte(),self.i()])
 def term(self):return {'integral':self.integral(),'coefficient':self.u()}
 def source(self):return {'id':self.string(),'domain':self.domain(),'conditions':self.vec(self.u),'terms':self.vec(self.term)}
 def seed(self):return {'row':self.u(),'integral':self.integral(),'shifts':self.vec(self.i)}
 def rule(self):return {'fixed':self.vec(lambda:self.opt(self.i)),'target':self.integral(),'rhs':self.vec(self.term),'sources':self.vec(self.seed),'domain':self.domain(),'discovery_domain':self.domain(),'conditions':self.vec(self.u),'sector':self.vec(self.byte),'permutation':self.opt(lambda:self.vec(self.u))}
def parse(path):
 raw=gzip.decompress(path.read_bytes());assert raw[:8]==b'RRPBIN\r\n';n=struct.unpack_from('<I',raw,16)[0];p=20
 for _ in range(n):
  tag,res,length=struct.unpack_from('<HHQ',raw,p);p+=12
  if tag==4: data=raw[p:p+length]
  p+=length
 assert p==len(raw)
 r=Reader(data);out={'schema':r.string(),'measure':r.string(),'roles':r.vec(r.byte),'indices':r.vec(r.u),'sources':r.vec(r.source),'zero_domains':r.vec(r.domain),'rules':r.vec(r.rule),'terminals':r.vec(lambda:r.vec(r.i))};assert r.p==len(data)
 return out,hashlib.sha256(raw).hexdigest()
def contains(domain,x):return all((lo is None or lo<=v) and(hi is None or v<=hi) for(lo,hi),v in zip(domain,x))
def targetmatch(integral,x):return all(s or v==p for(s,p),v in zip(integral,x))
if __name__=='__main__':
 result=[]
 for name,target in [('single',[1,1,0,1,1,1,2,-2,0,0,0,0,0,0,0,0]),('double',[1,1,1,2,1,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0])]:
  for round in range(3):
   path=ROOT/f'E7-origin-{name}/native-closure/round-{round:03d}-provisional.bin.gz';r,h=parse(path)
   domain_matches=[i for i,row in enumerate(r['rules']) if contains(row['domain'],target)]
   candidate_matches=[i for i,row in enumerate(r['rules']) if targetmatch(row['target'],target)]
   both=sorted(set(domain_matches)&set(candidate_matches))
   domain_misses={}
   for i in candidate_matches:
    misses=tuple(j for j,(bound,v) in enumerate(zip(r['rules'][i]['domain'],target)) if not contains([bound],[v]));domain_misses[misses]=domain_misses.get(misses,0)+1
   nearest=sorted(candidate_matches,key=lambda i:sum(not contains([bound],[v]) for bound,v in zip(r['rules'][i]['domain'],target)))[:3]
   entry={'case':name,'round':round,'target':target,'program':str(path.relative_to(ROOT)),'raw_sha256':h,'rules':len(r['rules']),'sources':len(r['sources']),'domain_matches':len(domain_matches),'candidate_target_matches':len(candidate_matches),'applicable_structure_matches':both,'source_domain_contains_target':sum(contains(s['domain'],target) for s in r['sources']),'is_zero_domain':any(contains(d,target) for d in r['zero_domains']),'candidate_match_domain_miss_counts':[{'axes':list(axes),'count':n} for axes,n in sorted(domain_misses.items(),key=lambda x:-x[1])[:8]],'nearest_rules':[{'ordinal':i,**{k:r['rules'][i][k] for k in ('target','domain','discovery_domain','sources')}} for i in nearest]}
   result.append(entry)
 print(json.dumps({'scope':'Read-only structural inspection of native serialized generated programs; not an independent replay or discovery certificate. Coefficient algebra is not decoded.','cases':result},indent=2))
