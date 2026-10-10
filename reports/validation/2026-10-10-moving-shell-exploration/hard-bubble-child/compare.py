"""Compare saved ordinary AMF child predictions; no integral evaluation is run."""
import hashlib,json
from pathlib import Path
import mpmath as mp
mp.mp.dps=100
HERE=Path(__file__).resolve().parent
prediction_path=HERE/'result/predictions.json'
reference_path=HERE.parent/'radial-reference.json'
pred=json.loads(prediction_path.read_text()); ref=json.loads(reference_path.read_text())
assert pred['schema']==1 and pred['status']=='predictions saved before comparison'
assert pred['reference_values_read'] is False and pred['bubble_subloops'] is False
assert pred['terminal_policy']=='TadpolesOnly' and pred['D']=='12/5' and pred['epsilon']=='4/5'
assert ref['status']=='passed' and ref['production_reference_values_loaded'] is False
assert len(pred['records'])==3
profiles={(r['digits'],r['guard_digits'],r['series_order']) for r in pred['records']}
assert profiles=={(20,40,60),(30,40,60),(30,40,80)}
assert len(profiles)==len(pred['records'])
# Two Minkowski denominators give a positive Wick phase. The ordinary native
# measure d^Dk/(i*pi^(D/2)) is (4*pi)^(D/2) times d^Dk_E/(2*pi)^D.
factor=(4*mp.pi)**(mp.mpf(6)/5)
raw=ref['leading_hard_coefficient']; assert mp.mpf(raw['im'])==0
expected=factor*mp.mpf(raw['re'])
rows=[]
for r in pred['records']:
 assert len(r['values'])==1
 value=mp.mpc(r['values'][0]['re'],r['values'][0]['im'])
 assert mp.isfinite(value.real) and mp.isfinite(value.imag)
 tolerance=mp.mpf('1e-18') if r['digits']==20 else mp.mpf('1e-28')
 relative=abs(value.real-expected)/abs(expected)
 imaginary=abs(value.imag)
 assert relative<tolerance and imaginary<tolerance,(r,relative,imaginary)
 rows.append({'profile':{k:r[k] for k in ('digits','guard_digits','series_order')},'relative_error':mp.nstr(relative,95),'imaginary_absolute_error':mp.nstr(imaginary,95),'tolerance':str(tolerance)})
finest=next(r for r in pred['records'] if r['digits']==30 and r['series_order']==80)
best=mp.mpc(finest['values'][0]['re'],finest['values'][0]['im']); refinements=[]
for r in pred['records']:
 if r is finest:continue
 value=mp.mpc(r['values'][0]['re'],r['values'][0]['im'])
 tolerance=mp.mpf('1e-18') if r['digits']==20 else mp.mpf('1e-28')
 error=abs(value-best)/abs(best);assert error<tolerance
 refinements.append({'digits':r['digits'],'series_order':r['series_order'],'relative_error':mp.nstr(error,95),'tolerance':str(tolerance)})
files=[prediction_path,reference_path,HERE/'probe.rs',HERE/'build-binding.json',HERE/'resources.json',HERE/'resources.provenance.json',HERE/'result/system.json',Path(__file__)]
report={'schema':1,'status':'passed','scope':'Ordinary native AMF hard-child feasibility only; not the whole moving-shell sunset, an integrated boundary adapter, or three/four-loop acceptance.',
 'reference_conversion':'positive (4*pi)^(D/2) for two virtual propagators, D=12/5',
 'source_reference_leading_coefficient':raw,'normalization_factor':mp.nstr(factor,95),'native_normalized_reference':mp.nstr(expected,95),
 'reference_comparisons':rows,'refinement_comparisons':refinements,
 'maximum_relative_reference_error':mp.nstr(max(mp.mpf(x['relative_error']) for x in rows),95),
 'maximum_relative_refinement_error':mp.nstr(max(mp.mpf(x['relative_error']) for x in refinements),95),
 'source_and_input_hashes':{str(p.relative_to(HERE.parents[3])):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
 'uncertainty_scope':'Independent saved 90-digit quadrature and empirical native precision/order refinements; no rigorous full integration error bound claimed.',
 'production_seed_values_supplied':False,'analytic_bubble_shortcut':False,'three_loop_acceptance':False,'four_loop_acceptance':False}
(HERE/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'passed','reference_comparisons':3,'refinement_comparisons':2,'maximum_relative_reference_error':report['maximum_relative_reference_error']}))
