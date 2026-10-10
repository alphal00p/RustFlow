#!/usr/bin/env python3
"""Audit exact saved point maps, retain failed attempts, archive and freeze."""
import collections
import gzip
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DEPTH=HERE.parent/'2026-10-10-targeted-native-depth-controls'
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not (HERE/'artifact-manifest.json').exists()
    bindings={}
    def bind(path):
        bindings[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)]=sha(path)
    def check_meta(item):
        path=Path(item['path']);assert sha(path)==item['sha256'];bind(path)
    for path in [HERE/'build-binding.json',HERE/'source-only/build-binding.json']:
        build=read(path);assert build['exit_code']==0 and build['cargo_invoked']is False;bind(path)
        for key in ['source','launcher','parent_source','parent_probe_source','parent_probe_build','parent_portfolio_build','executable']:
            if key in build:check_meta(build[key])
        for item in build['dependencies'].values():check_meta(item)
        for item in build.get('prior_failed_attempts',{}).values():check_meta(item)
    failures=[];cases=[];completions=collections.Counter();all_trials=[];total_terms=0
    for sector in ['single3','double']:
        old=HERE/sector;case=HERE/'source-only'/sector
        failed=read(old/'run-binding.json');res=read(old/'resources.json')
        assert failed['exit_code']==res['exit_code']==101 and failed['inputs_unchanged']
        assert 'unsupported guarded proof schema' in (old/'resources.log').read_text()
        assert not (old/'proposal-depth-3/result.json').exists()
        for item in failed['inputs'].values():check_meta(item)
        for key in ['build_binding','executable']:check_meta(failed[key])
        failures.append({'sector':sector,'exit_code':101,'reason':'unsupported guarded proof schema',
            'wall_seconds':res['wall_seconds'],'scope':'Explicit production/isolated proof-schema mismatch while checking historical selection; no fresh point result was produced.'})
        run=read(case/'run-binding.json');res=read(case/'resources.json')
        assert run['exit_code']==res['exit_code']==0 and run['inputs_unchanged']
        for name,digest in run['inputs'].items():
            path=Path(name);assert sha(path)==digest;bind(path)
        for key in ['build_binding','executable']:check_meta(run[key])
        for name in ['input.bin','input.json','points.json']:
            assert sha(case/name)==sha(old/name)==sha(DEPTH/sector/name);bind(DEPTH/sector/name)
        baseline_path=DEPTH/sector/'depth-3/result.json';baseline=read(baseline_path);bind(baseline_path)
        assert sha(case/'baseline-result.json')==sha(baseline_path)
        proposal_path=case/'proposal-depth-3/result.json';proposal=read(proposal_path)
        for key in ['depth','ray_domains','point_domains','source_count','source_identity','source_schema','selection_observations']:
            assert baseline[key]==proposal[key],key
        assert proposal['depth']==3 and proposal['ray_domains']==proposal['point_domains']==1
        assert proposal['all_rhs_labels_admitted']
        assert len(proposal['probes'])==len(baseline['probes'])==(4 if sector=='single3'else 3)
        points=[]
        for i,(b,p)in enumerate(zip(baseline['probes'],proposal['probes'])):
            assert b['selection']==p['selection']
            assert b['application']==p['application'],'exact ordered RHS/coefficient/condition comparison'
            assert b['ray']==p['ray'] and b['point']==p['point']
            assert p['roundtrip'] and p['point_fallback']
            assert p['ray']['status']=='Unresolved(NoApplicableRule)'
            assert p['application']['status'].startswith('Applied')
            assert len(p['native_portfolio'])==1
            stats=p['native_portfolio'][0];assert stats is not None
            assert stats['selected_arm']==0 and stats['completion']=='baseline-retained'
            assert stats['selected_rhs_terms']==stats['baseline_rhs_terms']==len(p['application']['rhs'])
            assert len(stats['trials'])==2
            for trial in stats['trials']:
                assert trial['attempted_rows']<=2048 and trial['accepted_rows']<=512
                assert trial['attempted_rows']==trial['accepted_rows']+trial['guard_rejected_rows']
                assert trial['exact_trace_rows']<=64 and trial['exact_trace_terms']<=16384
                assert trial['empty_rows']<=trial['accepted_rows']
                if trial['completion']=='attempted-row-budget':assert trial['attempted_rows']==2048
                completions[trial['completion']]+=1;all_trials.append(trial)
            proof=case/'proposal-depth-3'/f'point-{i:02}.bin';bind(proof)
            points.append({'selection':p['selection'],'status':p['application']['status'],
                'rhs_terms':len(p['application']['rhs']),'exact_rhs_and_conditions_match_baseline':True,
                'conditions':p['application']['conditions'],'ray_missed_point_applied':True,
                'native_roundtrip_replay_passed':True,'portfolio':stats,'proof_sha256':sha(proof)})
            total_terms+=len(p['application']['rhs'])
        cases.append({'sector':sector,'source_count':proposal['source_count'],'point_count':len(points),
            'successful_process_exit':0,'wall_seconds':res['wall_seconds'],'peak_child_rss_kib':res['peak_child_rss_kib'],
            'baseline_wall_seconds':read(DEPTH/sector/'depth-3-resources.json')['wall_seconds'],
            'rhs_terms':sum(p['rhs_terms']for p in points),'points':points})
        for folder in (old,case):
            for path in folder.rglob('*'):
                if path.is_file():bind(path)
    assert total_terms==45 and len(all_trials)==14
    assert dict(completions)=={'attempted-row-budget':11,'not-strictly-shorter':2,'rejected-new-condition-locus':1}
    summary={'status':'selected_points_replayed_no_portfolio_improvement','failed_historical_schema_attempts':failures,
        'source_only_cases':cases,'selected_points':7,'exact_application_matches':7,'baseline_rhs_terms':45,'proposal_rhs_terms':45,
        'ray_misses_recovered_by_points':7,'portfolio_baseline_retained':7,'trial_completion_counts':dict(completions),
        'trial_limits':{'attempted_rows':2048,'accepted_rows':512,'exact_trace_rows':64,'exact_trace_terms':16384,'depth':3},
        'source_import_audit':{'full_ordered_source_rows':True,'source_domains_and_nonzero_conditions':True,'all_zero_domains':True,
            'original_roles_and_index_variables':True,'original_selection_verified_by_bound_baseline':True,
            'historical_rules_used_in_source_only_run':False,'foreign_proof_schema_bypassed':False},
        'scope':'Seven saved labels only (six historical frontier leaves and one already-covered control). Source-only import rebuilds the complete native source context, then freshly searches/replays. All final maps and parameter conditions exactly match baseline. Trial budgets and one new-condition-locus rejection remain explicit. No general coverage, closure, period, production-change or speed claim.',
        'production_modified':False,'full_closure_proved':False,'numerical_predictions':0,
        'bindings_sha256':bindings}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Input/proof binaries are restored byte-exactly through the map. Small logs,
    # including both explicit failure logs, and all result JSON remain readable.
    entries=[]
    for path in sorted(HERE.rglob('*.bin')):
        raw=path.read_bytes();target=path.with_name(path.name+'.gz');assert not target.exists()
        target.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0));assert gzip.decompress(target.read_bytes())==raw
        entries.append({'original_path':str(path.relative_to(HERE)),'archive_path':str(target.relative_to(HERE)),
            'original_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':sha(target),
            'original_bytes':len(raw),'archive_bytes':target.stat().st_size})
        path.unlink()
    archives={'scope':'Completed copied inputs and native point proofs only; all failures, scripts, bindings, selected points and results remain readable. No original external baseline is changed.',
        'gzip_mtime':0,'restore_checks_passed':True,'entries':entries,'file_count':len(entries),
        'original_bytes':sum(x['original_bytes']for x in entries),'archive_bytes':sum(x['archive_bytes']for x in entries)}
    (HERE/'archive-map.json').write_text(json.dumps(archives,indent=2)+'\n')
    mapping={str(HERE/row['original_path']):row for row in entries}
    for name,digest in bindings.items():
        path=Path(name)
        if not path.is_absolute():path=ROOT/path
        raw=path.read_bytes()if path.exists()else gzip.decompress((HERE/mapping[str(path)]['archive_path']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==digest
    files={str(p.relative_to(HERE)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(HERE.rglob('*'))if p.is_file()and'__pycache__'not in p.parts}
    manifest={'status':'frozen','scope':summary['scope'],'file_count':len(files),'files':files,'archive_and_binding_verification':'passed'}
    (HERE/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'frozen','files':len(files),'archives':len(entries),'raw_bytes':archives['original_bytes'],
        'gzip_bytes':archives['archive_bytes'],'manifest_sha256':sha(HERE/'artifact-manifest.json')}))

if __name__=='__main__':main()
