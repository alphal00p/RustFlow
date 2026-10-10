"""Read-only exact helper-body, history, archive and immutable assessment check."""
from pathlib import Path
import json, hashlib, gzip, subprocess
B=Path(__file__).resolve().parent;R=B.parents[2]
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def function(s,name):
 pos=s.index('fn '+name+'(');start=s.rfind('\n',0,pos)+1;at=s.index('{',pos);depth=1;i=at+1
 while depth:
  if s[i]=='{':depth+=1
  if s[i]=='}':depth-=1
  i+=1
 return s[start:i]
checks=[]
for name in ['scalar_symbols','imaginary_parameter','encode_complex','convert','convert_at_epsilon']:
 original=function((R/'src/family.rs').read_text(),name).replace('pub(crate) fn','fn')
 copy=function((B/'harness/family.rs').read_text(),name).replace('pub(crate) fn','fn')
 # Trait declarations precede the converted bodies: locate implementation first.
 if name in ['convert','convert_at_epsilon']:
  copy=function((B/'harness/family.rs').read_text().split('impl DraftConvert for IntegralFamily {',1)[1],name)
 assert original==copy,name
 checks.append('exact family '+name)
for name in ['validate_distribution','polynomial_terms','laurent_terms']:
 assert function((R/'src/finite_density/boundary.rs').read_text(),name)==function((B/'harness/boundary.rs').read_text(),name),name
 checks.append('exact boundary '+name)
for name in ['rational_expression','rational_parts','canonical_conditions','rational_denominator_conditions']:
 assert function((R/'src/physical_conditions.rs').read_text(),name)==function((B/'harness/physical_conditions.rs').read_text(),name),name
 checks.append('exact conditions '+name)
for module,stop in [('algebra','/// Rational eigenvalues'),('coefficient','#[cfg(test)]'),('cut_regions','#[cfg(test)]')]:
 assert (R/f'src/{module}.rs').read_text().split(stop)[0]==(B/f'harness/{module}.rs').read_text()
 checks.append('exact '+module+' module prefix')
source=B/'source/src/finite_density/boundary/virtual_soft.rs';harness=B/'harness/virtual_soft.rs'
assert harness.read_text().replace('use crate::family::DraftConvert;\n','')==source.read_text()
checks.append('draft/test module identical except extension trait import')
for attempt in ['attempt1','attempt2','attempt3']:
 j=json.loads((B/(attempt+'-build-binding.json')).read_text())
 for rel,sha in j['inputs_sha256'].items():assert h(B/'history'/attempt/rel)==sha
 checks.append('historical exact '+attempt+' inputs')
for row in json.loads((B/'archive-map.json').read_text())['entries']:
 raw=gzip.decompress((B/row['archive']).read_bytes());assert hashlib.sha256(raw).hexdigest()==row['sha256'];assert h(B/row['archive'])==row['archive_sha256']
 checks.append('lossless archive '+row['path'])
A=R/'reports/validation/2026-10-10-partial-placement-boundary-assessment'
for row in json.loads((A/'artifact-manifest.json').read_text())['files']:assert h(A/row['path'])==row['sha256']
checks.append('frozen original assessment remains byte-exact')
subprocess.run(['git','apply','--check',str(B/'revised-virtual-soft.patch')],cwd=R,check=True)
checks.append('unapplied production patch applies cleanly')
print(json.dumps({'checks':checks,'passed':len(checks),'production_edits':False},indent=2))
