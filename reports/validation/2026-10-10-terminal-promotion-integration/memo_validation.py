"""Check saved memo configuration/counters without inferring proof or total work."""
CAPS={'max_entries':64,'max_output_terms':4096,'max_conditions':32768,'max_polynomial_terms':262144,'max_owned_bytes':67108864}
PEAKS={'unit_memo_peak_entries':'max_entries','unit_memo_peak_output_terms':'max_output_terms','unit_memo_peak_conditions':'max_conditions','unit_memo_peak_polynomial_terms':'max_polynomial_terms','unit_memo_peak_owned_bytes':'max_owned_bytes'}
COUNTERS=['native_rule_applications','performed_native_rule_applications','memoized_rule_applications','unit_memo_hits','unit_memo_misses',*PEAKS]
def validate_caps(value,enabled):
 if not enabled:
  assert value is None,'disabled memo must have no configured caps'
  return
 assert isinstance(value,dict)
 assert {k:value[k]for k in CAPS}==CAPS
 assert all(type(value[k])is int for k in CAPS)
 assert set(value) in [set(CAPS)|{'metric'},set(CAPS)|{'owned_bytes_metric'}]
 text=value.get('metric',value.get('owned_bytes_metric'))
 assert isinstance(text,str)and 'scratch' in text and 'allocator' in text

def validate_memo(config,closure,enabled,cuts):
 validate_caps(config['unit_reduction_memo'],enabled)
 assert isinstance(closure,list)and [r['cut_slots']for r in closure]==cuts
 rows=[]
 for entry in closure:
  d=entry['diagnostics']
  assert all(type(d[k])is int and d[k]>=0 for k in COUNTERS)
  assert d['native_rule_applications']==d['performed_native_rule_applications']+d['memoized_rule_applications']
  if enabled:
   for peak,cap in PEAKS.items():assert d[peak]<=CAPS[cap]
   assert d['unit_memo_misses']>0
  else:
   assert all(d[k]==0 for k in COUNTERS if k not in ['native_rule_applications','performed_native_rule_applications'])
  rows.append({'cut_slots':entry['cut_slots'],**{k:d[k]for k in COUNTERS}})
 if enabled:assert sum(r['unit_memo_hits']for r in rows)>0,'numeric control did not exercise a completed-call hit'
 return {'enabled':enabled,'configured_caps':CAPS if enabled else None,'cuts':rows,'logical_equals_performed_plus_memoized':True,'scope':'Native unit reduction applications only; full proof replay, terminal-admission checks, encoding, decoding and arithmetic scratch are not measured by these counters. Quotas bound retained payload, not total process memory or wall time.'}
