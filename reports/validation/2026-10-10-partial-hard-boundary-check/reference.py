#!/usr/bin/env python3
"""Independent convergent one-mass Gaussian/Schwinger reference; no predictions."""
from decimal import Decimal as D, localcontext
from pathlib import Path
import hashlib,json,time
from gamma_reference import gamma_positive
BASE=Path(__file__).resolve().parent
assert not (BASE/'reference.json').exists()
rows=[]
for digits in (35,60,85):
    start=time.monotonic()
    with localcontext() as ctx:
        ctx.prec=digits+20
        dimension=D(12)/5
        arguments=[dimension/2-1,2-dimension/2,3-dimension,dimension/2]
        gamma=[gamma_positive(x,digits) for x in arguments]
        result=gamma[0][0]**2*gamma[1][0]*gamma[2][0]/gamma[3][0]
        rows.append({'digits':digits,'gamma_arguments':[str(x) for x in arguments],
            'gamma_values':[str(g[0]) for g in gamma],
            'gamma_log_truncation_bounds':[str(g[1]) for g in gamma],
            'native_minkowski':str(-result),'euclidean_pi_normalized':str(result),
            'wall_seconds':time.monotonic()-start})
with localcontext() as ctx:
    ctx.prec=110
    reference=D(rows[-1]['native_minkowski'])
    changes=[abs(D(r['native_minkowski'])-reference)/abs(reference) for r in rows]
    assert max(changes)<D('1e-33')
    out={'scope':'Independent validation only. Positive convergent Gaussian/Schwinger reference; no native predictions or oracle constants read.',
        'dimension':'12/5','epsilon':'4/5','masses_squared':['1','0','0'],'powers':[1,1,1],
        'formula':'Gamma(D/2-1)^2 Gamma(2-D/2) Gamma(3-D) / Gamma(D/2)',
        'absolute_convergence_domain':'2 < Re(D) < 3',
        'normalization':'Euclidean measures d^Dk/pi^(D/2); native Minkowski sign (-1)^3.',
        'reference_native':str(reference),'reference_absolute_imaginary_part':'0',
        'accepted_comparison_relative_tolerance':'1e-10','profiles':rows,
        'relative_changes_from_selected':[str(x) for x in changes],
        'oracle_records_read':0,'native_prediction_files_read':0,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'gamma_implementation_sha256':hashlib.sha256((BASE/'gamma_reference.py').read_bytes()).hexdigest()}
(BASE/'reference.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'reference_native':out['reference_native'],'profiles':len(rows)}))
