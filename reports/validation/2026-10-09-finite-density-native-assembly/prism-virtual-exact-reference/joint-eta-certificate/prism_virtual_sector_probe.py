from itertools import permutations,product
from collections import Counter
from pathlib import Path
import json
names='abcdef';U=['ab','ad','af','bc','cd','cf','be','de','ef','bf','df'];Phi=['abc','acd','acf','bcf','abd','bcd','bde','bdf','adf'];Psi=['abe','ade','aef','bef']
def mono(s):return tuple(s.count(x)for x in names)
def exponents(s,order):
 e=mono(s);return tuple(sum(e[i]for i in order[j:])for j in range(1,6))
def chart(poly,order):
 terms=[exponents(s,order)for s in poly];minimum=tuple(min(t[j]for t in terms)for j in range(5));return minimum,minimum in terms,[[t[j]-minimum[j]for j in range(5)]for t in terms]
rows=[]
for order in permutations(range(6)):
 um,uu,up=chart(U,order);fm,fu,fp=chart(Phi,order);pm,pu,pp=chart(Psi,order)
 rows.append({'order':[names[i]for i in order],'U':{'monomial':um,'unit':uu,'residual':up},'Phi':{'monomial':fm,'unit':fu,'residual':fp},'Psi':{'monomial':pm,'unit':pu,'residual':pp}})
report={'scope':'exact positive-polynomial primary Hepp chart probe for coupled two-virtual prism coefficient family; no values/contour admission','charts':len(rows),'unit_counts':{name:sum(r[name]['unit']for r in rows)for name in ['U','Phi','Psi']},'both_U_Phi_units':sum(r['U']['unit']and r['Phi']['unit']for r in rows),'charts_data':rows}
Path('/tmp/prism-virtual-sector-probe.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items()if k!='charts_data'})
