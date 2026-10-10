from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib
root=Path('reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference')
paths=[root/'laurent-baseline/sample-01.json',root/'laurent-refined-pilot/sample-01.json'];j=[json.loads(p.read_text())for p in paths]
with localcontext()as c:
 c.prec=100;a,b=[Decimal(v['full_msbar_reference_normalization'][1:].split('+',1)[0])for v in j];delta=abs(a-b);relative=delta/max(abs(a),abs(b))
report={'status':'finite-epsilon quadrature/precision refinement only; no Laurent or native AMF acceptance','target':j[0]['target'],'epsilon':j[0]['epsilon'],'inputs':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'digits':v['digits'],'step':v['step'],'cutoff':v['truncation']}for p,v in zip(paths,j)],'observed_absolute_difference':str(delta),'observed_relative_difference':str(relative),'per_monomial_regulator_residue_checks':[len(v['rho_residue_checks'])for v in j],'uncertainty':'empirical difference, not rigorous enclosure','note':'The pilot has only one epsilon sample; its generic controller interpolation output is not a Laurent fit and is not used.'}
(root/'near-zero-pilot-comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
