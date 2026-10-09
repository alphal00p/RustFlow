"""Compare native HEPKit support polynomials with independent routed arithmetic.

Usage: python compare_native.py NATIVE_REPORT_DIRECTORY [OUTPUT_JSON]
This compares exact coefficients and supports, not numerical integral values.
"""
from pathlib import Path
from fractions import Fraction
import gzip,hashlib,json,sys
base=Path(__file__).resolve().parent
independent=base.parent/'prism-endpoint-proof/rustflow_prism_parametric.json'
reference=json.loads(independent.read_text());native_dir=Path(sys.argv[1]);runs=[]
def polynomial(terms):
 return {tuple(t['exponents']):Fraction(t['coefficient'])for t in terms}
def restrict(terms,positions):
 omitted=set(range(6))-set(positions);out={}
 for term in terms:
  if any(term['powers'][j]for j in omitted):continue
  power=tuple(term['powers'][j]for j in positions)
  out[power]=out.get(power,Fraction(0))+Fraction(term['coefficient'])
 return {p:c for p,c in out.items()if c}
for ref in reference:
 a,b=ref['cuts'];path=native_dir/f'prism-cuts-{a}-{b}-supports.json'
 if not path.exists():path=path.with_suffix('.json.gz')
 payload=gzip.decompress(path.read_bytes())if path.suffix=='.gz'else path.read_bytes()
 data=json.loads(payload);supports=data['supports'];assert len(supports)==64
 for support in supports:
  positions=[ref['parameter_slots'].index(s)for s in support['active_slots']]
  expected_u=restrict(ref['U'],positions);expected_v=restrict(ref['V'],positions)
  assert polynomial(support['u'])==expected_u,('U',ref['cuts'],support['active_slots'])
  if expected_u:
   assert polynomial(support['v'])==expected_v,('V',ref['cuts'],support['active_slots'])
  else:
   assert support['active_virtual_rank']<support['virtual_loops']
  assert support['routed_rows']==ref['routed_rows'],('routing',ref['cuts'])
  assert [s for s,c in support['pure_transfer_slots']]==ref['pure_transfer_slots']
 runs.append({'cuts':ref['cuts'],'supports_compared':64,'native_file':str(path),'native_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'compiletime_module_blake3':data['source_blake3'],'compiletime_fixture_blake3':data['input_fixture_blake3']})
result={'scope':'Exact native HEPKit versus independent rational routed U/V comparison, with every active support and exact route; no period values or production admission','passed':192,'failed':0,'reference_sha256':hashlib.sha256(independent.read_bytes()).hexdigest(),'runs':runs}
out=Path(sys.argv[2])if len(sys.argv)>2 else native_dir/'independent-polynomial-comparison.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'passed':192,'failed':0,'output':str(out)}))
