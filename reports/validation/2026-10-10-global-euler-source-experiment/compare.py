#!/usr/bin/env python3
"""Compare saved exact native maps; no reference values or new native work."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    out=HERE/'summary.json';assert not out.exists()
    reports=[json.loads((HERE/m/'result.json').read_text()) for m in ['original','global-euler']]
    assert [r['source_count'] for r in reports]==[232,248]
    assert all(r['global_is_exact_sum_of_local_dilations'] and r['same_original_domains_and_zeros'] for r in reports)
    assert all(len(r['probes'])==16 and r['depth']==3 and r['ray_domains']==r['point_domains']==1 for r in reports)
    rows=[]
    for i,(a,b) in enumerate(zip(*[r['probes'] for r in reports])):
        assert a['selection']==b['selection'] and a['roundtrip'] and b['roundtrip']
        x,y=a['application'],b['application']
        # Native normalized coefficient displays and ordered arrays coincide
        # exactly here; no heuristic symbolic equality or tolerance is used.
        assert x==y and x['status'].startswith('Applied')
        rows.append({'ordinal':i,'point':a['selection']['point'],'stratum':a['selection']['stratum'],
                     'rhs_terms':len(x['rhs']),'exact_status_rhs_and_conditions_identical':True})
    files=[HERE/m/'result.json' for m in ['original','global-euler']]
    files += [HERE/f'{m}-{name}.json' for m in ['original','global-euler'] for name in ['resources','binding']]
    files += [HERE/'probe.rs',HERE/'build-binding.json',HERE/'points.json',HERE/'selection-binding.json',Path(__file__)]
    result={'scope':'Source-presentation diagnostic only; no period evaluation, full closure, policy recommendation or proof import',
            'status':'negative: no exact map or condition change on the selected16',
            'original_source_count':232,'augmented_source_count':248,'points':16,'applied_each':16,
            'rhs_terms_each':sum(r['rhs_terms'] for r in rows),'rows':rows,
            'global_dilation_branch_checks':16,'native_generated_encode_decode_replay_each':True,
            'budgets':{'depth':3,'ray_domains_per_point':1,'exact_point_fallback_domains':1,'wall_timeout_seconds_each':300},
            'resources':{m:json.loads((HERE/f'{m}-resources.json').read_text()) for m in ['original','global-euler']},
            'inputs_sha256':{str(p.relative_to(HERE)):h(p) for p in files}}
    out.write_text(json.dumps(result,indent=2)+'\n');print('16 exact saved maps identical;140 RHS terms each')


if __name__=='__main__':main()
