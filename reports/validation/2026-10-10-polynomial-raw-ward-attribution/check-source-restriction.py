"""Exact rational polynomial source-row comparison; no integral reduction."""
import ast,collections,datetime,fractions,hashlib,json,pathlib
F=fractions.Fraction
P=pathlib.Path(__file__).resolve().parent
normal_path=P.parent/'2026-10-10-normalized-boost-reciprocal-experiment/native-source-probe/normalized-last/result.json'
raw_path=P/'native-source-probe/unrecentered-boost-ward/result.json'
normal=json.loads(normal_path.read_text());raw=json.loads(raw_path.read_text())
def add(a,b):
 out=dict(a)
 for m,c in b.items():out[m]=out.get(m,F(0))+c
 return {m:c for m,c in out.items() if c}
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def mul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():
   powers=collections.Counter(dict(m));powers.update(dict(n));key=tuple(sorted(powers.items()));out[key]=out.get(key,F(0))+c*d
 return {m:c for m,c in out.items() if c}
def poly(text,fixed):
 def visit(n):
  if isinstance(n,ast.Constant) and isinstance(n.value,int):return {():F(n.value)} if n.value else {}
  if isinstance(n,ast.Name):return {():F(fixed[n.id])} if n.id in fixed and fixed[n.id] else ({} if n.id in fixed else {((n.id,1),):F(1)})
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return scale(visit(n.operand),-1)
  if isinstance(n,ast.BinOp):
   a,b=visit(n.left),visit(n.right)
   if isinstance(n.op,ast.Add):return add(a,b)
   if isinstance(n.op,ast.Sub):return add(a,scale(b,-1))
   if isinstance(n.op,ast.Mult):return mul(a,b)
   if isinstance(n.op,ast.Div):assert len(b)==1 and () in b and b[()];return scale(a,1/b[()])
  raise AssertionError(('unsupported polynomial AST',ast.dump(n)))
 return visit(ast.parse(text.replace('^','**'),mode='eval').body)
def serial(row):return [{'shift':list(k),'coefficient':[(list(m),str(c))for m,c in sorted(v.items())]}for k,v in sorted(row.items())]
def subset(child,parent):
 for (lo,hi),(pl,ph) in zip(child,parent):
  assert pl is None or lo is not None and lo>=pl
  assert ph is None or hi is not None and hi<=ph
rows=[]
for r in raw['additional_sources']:
 if not r['id'].startswith('polynomial-singleton-raw-normalized-boost/'):continue
 positive='positive-upper' in r['id'];candidates=[n for n in normal['normalized_boost_sources'] if n['domain'][9]==r['domain'][9] and n['domain'][10]==r['domain'][10]];assert len(candidates)==1;n=candidates[0]
 assert r['domain'][0]==[1,1] and all(r['domain'][i]==[0,0]for i in [6,7,8,10,11]);subset(r['domain'],n['domain']);fixed={f'a_{i}':v for i,v in [(0,1),(6,0),(7,0),(8,0),(10,0),(11,0)]}
 if not positive:fixed['a_9']=0
 projected={};dropped=[]
 for t in n['terms']:
  assert all(s for s,_ in t['powers']);shift=tuple(v for _,v in t['powers']);c=poly(t['coefficient'],fixed)
  reason=None
  if not c:reason='exact coefficient zero under guard'
  elif 1+shift[0]<=0:reason='required-cut C_n zero for n<=0'
  elif 1+shift[0]>=1 and shift[10]>=1:reason='existing finite-label lower-contact origin zero'
  if reason:dropped.append({'shift':list(shift),'coefficient':t['coefficient'],'reason':reason});continue
  projected[shift]=add(projected.get(shift,{}),c)
 expected={}
 for t in r['terms']:
  key=tuple(t['shift']);expected[key]=add(expected.get(key,{}),poly(t['coefficient'],fixed))
 projected={k:v for k,v in projected.items()if v};expected={k:v for k,v in expected.items()if v}
 assert projected=={k:scale(v,8)for k,v in expected.items()}
 normal_conditions=[poly(c,fixed)for c in n['conditions']];raw_conditions=[poly(c,fixed)for c in r['conditions']]
 def nonconstant(values):
  for v in values:
   assert v, "vanished source condition"
  return [v for v in values if not(len(v)==1 and () in v)]
 assert nonconstant(normal_conditions)==nonconstant(raw_conditions)
 rows.append({'branch':'positive upper'if positive else'bulk','normalized_id':n['id'],'raw_id':r['id'],'raw_guard':r['domain'],'normalized_guard':n['domain'],'raw_guard_is_subset':True,'constant_nonzero_scale':8,'coefficient_and_shift_equality':True,'raw_row':serial(expected),'restricted_normalized_row':serial(projected),'dropped':dropped,'nonconstant_conditions_match':True,'normalized_declared_conditions_retained':n['conditions'],'raw_declared_conditions_retained':r['conditions'],'extra_conditions_are_verified_nonzero_rational_constants':True})
assert len(rows)==2
meta=lambda f:{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
out={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Exact symbolic rational-polynomial comparison of source coefficients and shifts after explicit domain restriction and existing zero projection; no elimination or period evaluation','normal_input':meta(normal_path),'raw_input':meta(raw_path),'checker':meta(pathlib.Path(__file__)),'rows':rows,'physical_proof_caveat':'The polynomial raw Ward source is independently justified by the shell cutoff theorem. This comparison does not replay an original polynomial source at a forbidden reciprocal seed or widen production admission.'}
(P/'source-restriction-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: two source branches, exact shift/guard/coefficient equality up to fixed scale8')
