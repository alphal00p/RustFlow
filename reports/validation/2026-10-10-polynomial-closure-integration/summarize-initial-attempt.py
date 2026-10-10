#!/usr/bin/env python3
"""Preserve the first integration attempt with its failed fixture admission tests."""
import hashlib
import json
import re
from pathlib import Path

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PREFIX='polynomial-closure'

def main():
    summary_path=BASE/(PREFIX+'-attempt-summary.json')
    manifest_path=BASE/(PREFIX+'-attempt-manifest.json')
    assert not summary_path.exists()and not manifest_path.exists()
    read=lambda path:json.loads(path.read_text())
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    build_path=BASE/(PREFIX+'-build-provenance.json')
    snapshot_path=BASE/(PREFIX+'-source-hashes.json')
    build=read(build_path)
    assert build['source_snapshot_sha256']==sha(snapshot_path)
    rows=[]
    gates=['lib','flow-boundary','massless-sources','source-fingerprint','python','capacity','polynomial-sources']
    for gate in gates:
        stem=PREFIX+'-'+gate
        resources=read(BASE/(stem+'-resources.json'))
        binding=read(BASE/(stem+'-binding.json'))
        log=(BASE/(stem+'-resources.log')).read_text()
        matches=re.findall(r'test result: (?:ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out; finished in ([0-9.]+)s',log)
        assert len(matches)==1
        passed,failed,ignored,measured,filtered=map(int,matches[0][:5])
        assert binding['source_snapshot_sha256']==sha(snapshot_path)
        assert binding['build_provenance_sha256']==sha(build_path)
        assert resources['exit_code']==binding['exit_code']==(101 if gate=='lib'else 0)
        assert failed==(2 if gate=='lib'else 0)
        rows.append({'gate':gate,'passed':passed,'failed':failed,'ignored':ignored,
            'filtered':filtered,'exit_code':resources['exit_code'],
            'harness_seconds':float(matches[0][5]),'wall_seconds':resources['wall_seconds']})
    correction_path=BASE/'fixture-correction-v2/correction-binding.json'
    correction=read(correction_path)
    snapshot=read(snapshot_path)
    for name,info in correction['files'].items():
        assert snapshot['files'][name]==info['before_sha256']
        assert sha(BASE/'fixture-correction-v2/failed-source'/name)==info['before_sha256']
    assert sha(BASE/'fixture-correction-v2/fixture-only.patch')==correction['patch_sha256']
    source_report_path=BASE/(PREFIX+'-actual62-source-equivalence.json')
    source_report=read(source_report_path)
    assert source_report['source_count']==62 and source_report['all_original_legacy_rows_preserved']
    assert source_report['source_context_roundtrip']and not source_report['historical_rules_imported']
    native_path=BASE/(PREFIX+'-native-gate-reuse.json')
    native=read(native_path)
    assert native['reused_passed_tests']==84 and native['new_native_test_invocations']==0
    summary={'schema':1,'status':'failed_fixture_admission_tests',
        'scope':'First optional polynomial source policy integration attempt. Compilation and exact factory/source equivalence passed; two new test fixtures failed the existing graph-admission checks before source generation. Numerical runs were not started.',
        'build':{'path':str(build_path.relative_to(ROOT)),'sha256':sha(build_path),'compile_identities':build['compile_identities']},
        'source_snapshot':{'path':str(snapshot_path.relative_to(ROOT)),'sha256':sha(snapshot_path),'files':len(snapshot['files'])},
        'new_passed_tests':sum(row['passed']for row in rows),'new_failed_tests':sum(row['failed']for row in rows),
        'reused_native_passed_tests':84,'new_native_invocations':0,'gates':rows,
        'failed_tests':[
            'finite_density::massless_endpoint::raw_ward_draft_tests::raw_ward_maps_actual_energy_and_checks_every_base_guard',
            'finite_density::measure::polynomial_closure_draft_tests::exact_energy_certificate_distinguishes_rescaled_compact_from_mixed_or_virtual_energy'],
        'failure':'Unsupported: species 0 does not form directed nonbranching fermion cycles at vertex 0',
        'exact_source_gate':{'source_count':62,'original_legacy_rows_preserved':True,'source_context_roundtrip':True,
            'historical_rules_imported':False,'report_sha256':sha(source_report_path)},
        'next_attempt':{'directory':'reports/validation/2026-10-10-polynomial-closure-integration-v2',
            'source_policy_remains':'polynomial-closure-v1','change_scope':'cfg(test) fixture routing/charge assignments and associated assertions only; production source policy unchanged',
            'correction_binding_sha256':sha(correction_path),'correction_patch_sha256':correction['patch_sha256'],
            'validation_status':'separate build and gates pending at this checkpoint'},
        'native_full_amplitude_runs':0,'new_three_loop_numerical_acceptance':False,'new_four_loop_numerical_acceptance':False,
        'summary_script_sha256':sha(Path(__file__))}
    assert summary['new_passed_tests']==103 and summary['new_failed_tests']==2
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    files={str(path.relative_to(BASE)):{'sha256':sha(path),'bytes':path.stat().st_size}
           for path in sorted(BASE.rglob('*'))if path.is_file()and'__pycache__'not in path.parts and path!=manifest_path}
    manifest={'schema':1,'scope':'Frozen complete first-attempt artifacts, including failed-fixture sources and correction evidence. No passing release or new numerical acceptance claim.',
              'files':files,'file_count':len(files),'bytes':sum(item['bytes']for item in files.values()),
              'summary_sha256':sha(summary_path)}
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    assert all(sha(BASE/name)==info['sha256']for name,info in files.items())
    print(json.dumps({'status':summary['status'],'passed':103,'failed':2,'reused_native':84,'manifest_files':len(files)}))

if __name__=='__main__':main()
