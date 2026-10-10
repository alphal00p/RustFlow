#!/usr/bin/env python3
"""Compare only complete successful saved predictions, with frozen reference hashes.

Existing comparator implementations and tolerances are reused unchanged. This
wrapper adds current build/source/policy binding, complete profile checks and a
matching fixed-dimension prerequisite for Laurent construction provenance.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
from history_coverage import verify_history_coverage
from certificate_validation import load_certificates, validate_sample

PREFIX='singleton-physical-zero-v2'
PROFILES=[[18,60,8],[28,60,8],[28,80,8],[28,80,12]]

def read(path):return json.loads(path.read_text())

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gate',choices=['two-loop-fixed','two-loop-laurent','three-loop-fixed','three-loop-laurent','massive-regression'])
    args=parser.parse_args()
    stem=PREFIX+'-'+args.gate
    directory=BASE/stem
    output=BASE/(stem+'-comparison.json')
    evidence_path=BASE/(stem+'-comparison-binding.json')
    assert not output.exists()and not evidence_path.exists(),'preserve existing comparisons'
    snapshot_path=BASE/(PREFIX+'-source-hashes.json')
    snapshot=verify_snapshot(snapshot_path)
    build_path=BASE/(PREFIX+'-build-provenance.json')
    build=read(build_path)
    assert build['source_snapshot_sha256']==digest(snapshot_path)
    binding_path=BASE/(stem+'-binding.json')
    resources_path=BASE/(stem+'-resources.json')
    provenance_path=BASE/(stem+'-resources.provenance.json')
    binding,resources,provenance=map(read,[binding_path,resources_path,provenance_path])
    assert binding['exit_code']==resources['exit_code']==0,'native run must finish successfully before comparison'
    assert binding['post_run_source_and_executable_unchanged'] is True
    assert binding['source_snapshot_sha256']==digest(snapshot_path)
    assert binding['build_provenance_sha256']==digest(build_path)
    assert binding['command']==resources['command']==provenance['command']
    assert not(directory/'failure.json').exists(),'failure artifact prevents numerical acceptance'
    captured=[v for k,v in provenance['argument_files'].items()if Path(k).name.startswith('finite_density_runtime_flow-')]
    assert len(captured)==1
    names=[name for name in build['files']if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert len(names)==1
    executable=names[0]
    assert captured[0]['sha256']==build['files'][executable]['sha256']==digest(ROOT/executable)
    config_path=directory/'configuration.json'
    config=read(config_path)
    massive=args.gate=='massive-regression'
    three=args.gate.startswith('three')
    laurent=args.gate.endswith('laurent')
    expected_history_limit=None if massive else 0
    assert binding['settings'].get('RUSTFLOW_WEIGHTED_MAX_HISTORY_CANDIDATE_MAPS')==(None if massive else '0')
    assert binding['requested_ray_point_schedule_enabled']is three
    if three:
        expected_schedule={'RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT':'true',
            'RUSTFLOW_WEIGHTED_RAY_DOMAINS':'1','RUSTFLOW_WEIGHTED_POINT_DOMAINS':'1',
            'RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'0','RUSTFLOW_WEIGHTED_GUARD_PASSES':'0',
            'RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'8192'}
        assert all(binding['settings'].get(key)==value for key,value in expected_schedule.items())
        assert not binding['experimental_programs_or_rules_imported']
    else:
        assert binding['settings'].get('RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT','false')=='false'
    expected_input='examples/finite_density/'+('massive_two_loop_sunset'if massive else'massless_three_loop_chain'if three else'massless_two_loop_sunset')+'.json'
    input_path=ROOT/expected_input
    assert config['input_path']==expected_input
    assert digest(input_path)==snapshot['files'][expected_input]
    assert read(directory/'input.json')==read(input_path)
    assert config['profiles']==([[28,80,12]]if massive else PROFILES)
    assert config['guard_digits']==40 and config['independent_reference_comparisons']==0
    assert config['source_options']=={'policy':'polynomial-closure-v1',
        'positive_compact_energy_powers':False,'free_virtual_zero_sectors':not massive}
    assert config['epsilon']==('4/5'if massive else'-5/4')
    assert config['initial_epsilon_grid']==1000
    profiles=[(28,80,12)]if massive else[tuple(p)for p in PROFILES]
    if laurent:profiles=[(*p,1000)for p in profiles]+[(28,80,12,2000)]
    prediction_paths=[directory/('prediction-'+'-'.join(map(str,p[:3]))+(f'-grid-{p[3]}'if laurent else'')+'.json')for p in profiles]
    assert set(directory.glob('prediction-*.json'))==set(prediction_paths),'complete exact profile set required'
    cuts=[[],[0],[3],[0,3]]if three else[[],[0],[1],[0,1]]
    certificate_path=directory/'physical-zero-certificates.json'
    certificates=read(certificate_path)
    if massive:assert certificates==[], 'massive gate must not receive a massless physical-zero certificate'
    else:certificates,_=load_certificates(directory,read(input_path),cuts)
    for path,profile in zip(prediction_paths,profiles):
        row=read(path)
        assert row['full_amplitude']is True and row['independent_reference_comparisons']==0
        assert row['normalization']=='unscaled Euclidean amplitude'
        assert [row[key]for key in ['digits','series_order','occupied_start_scale']]==list(profile[:3])
        assert row['guard_digits']==40
        if laurent:
            assert row['epsilon_grid_denominator']==profile[3]
            assert len(row['expansions'])==2
        else:
            assert row['epsilon']==config['epsilon']
            assert [item['cut_slots']for item in row['contributions']]==cuts
            assert [item['cut_slots']for item in row['occupied_reports']]==cuts[1:]
            if massive:
                assert row['physical_zero_certificates']==[]
                for item in row['occupied_reports']:
                    evaluation=item.get('report',item)
                    assert evaluation['construction']=='weighted_amf'
                    assert evaluation['source_options']==config['source_options']
            else:validate_sample(row,certificates,read(input_path),cuts,config['source_options'])
    tool_record_path=BASE/'native-digest-build.json'
    tool_record=read(tool_record_path)
    assert tool_record['exit_code']==0
    assert tool_record['source_snapshot_sha256']==digest(snapshot_path)
    assert tool_record['build_provenance_sha256']==digest(build_path)
    digest_executable=Path(tool_record['executable']['path'])
    assert digest(digest_executable)==tool_record['executable']['sha256']
    def native_blake3(path):
        before=path.stat()
        fields=subprocess.check_output([str(digest_executable),str(path)],text=True).strip().split()
        after=path.stat()
        assert len(fields)==2 and re.fullmatch('[0-9a-f]{64}',fields[0])
        assert int(fields[1])==before.st_size==after.st_size
        assert (before.st_mtime_ns,before.st_ino)==(after.st_mtime_ns,after.st_ino)
        return fields[0]
    checked=[tool_record_path,BASE/'native-digest.rs',BASE/'build-native-digest.py',BASE/'history_coverage.py',
             snapshot_path,build_path,binding_path,resources_path,provenance_path,
             config_path,input_path,directory/'input.json',*prediction_paths,Path(__file__),certificate_path,BASE/'certificate_validation.py']
    metadata=[]
    if not massive:
        for cut in cuts[1:-1]:
            assert not list((directory/'native-closure'/('cut-'+'-'.join(map(str,cut)))).glob('round-*-closed.*')), 'a physical-zero certificate must not pretend to close a native program'
    for cut in (cuts[1:] if massive else [cuts[-1]]):
        paths=sorted((directory/'native-closure'/('cut-'+'-'.join(map(str,cut)))).glob('round-*-closed.json'))
        assert paths,'every occupied cut needs saved successful native closure evidence'
        for path in paths:
            row=read(path)
            assert row['status']=='closed'
            assert row['discovery_schedule']==('requested-ray-point-v1'if three else'legacy')
            assert row['requested_ray_point']==({'max_domains_per_ray':1,'max_domains_per_point':1}if three else None)
            proof_path=path.with_suffix('.bin')
            assert proof_path.is_file(),'closed metadata must retain its native proof binary'
            assert proof_path.stat().st_size==row['program_bytes']
            assert native_blake3(proof_path)==row['program_blake3']
            transaction_binding=None
            coverage_binding=None
            assert row['max_history_candidate_maps']==expected_history_limit
            assert row['active_target_closure'] is (not massive)
            if row['active_target_closure']:
                provisional_path=path.with_name(path.name.replace('-closed.json','-provisional.json'))
                provisional=read(provisional_path)
                provisional_proof=provisional_path.with_suffix('.bin')
                retired_path=path.parent/'retired-unresolved.json'
                coverage_binding=verify_history_coverage(row,provisional,expected_history_limit,read(retired_path))
                assert provisional_proof.stat().st_size==provisional['program_bytes']
                assert native_blake3(provisional_proof)==provisional['program_blake3']
                coverage_binding.update(provisional_metadata_sha256=digest(provisional_path),
                    provisional_program_sha256=digest(provisional_proof),retired_unresolved_sha256=digest(retired_path))
                checked.extend([provisional_path,provisional_proof,retired_path])
            if three:
                link=row['active_state']['requested_discovery']
                assert re.fullmatch(r'requested-discovery-[0-9]+[.]json',link['path'])
                transaction_path=path.parent/link['path']
                assert native_blake3(transaction_path)==link['blake3']
                transaction=read(transaction_path)
                assert transaction['schema']==1 and transaction['status']=='complete'
                assert transaction['schedule']=='requested-ray-point-v1'
                assert transaction['source_measure_id']==row['measure_id']
                assert transaction['physical_arity']==row['physical_arity']
                assert transaction['storage_capacity']==row['storage_capacity']
                assert transaction['policy']==row['requested_ray_point']
                assert transaction['empty_terminals']is True
                assert transaction['completed_search_implies_coverage']is False
                assert provisional['round']==row['round']
                # Exact shared preclosure metadata/state already verified above;
                # only the final optional coverage object is newly present.
                assert provisional_proof.stat().st_size==provisional['program_bytes']==transaction['native_program']['bytes']
                assert native_blake3(provisional_proof)==provisional['program_blake3']==transaction['native_program']['blake3']
                transaction_binding={'path':str(transaction_path.relative_to(ROOT)),'sha256':digest(transaction_path),
                    'blake3':link['blake3'],'provisional_metadata_sha256':digest(provisional_path),
                    'provisional_program_sha256':digest(provisional_proof)}
                checked.extend([transaction_path,provisional_path,provisional_proof])
            fields=row['measure_id'].split(':',3)
            assert fields[:3]==['rustflow-weighted-sources-v1',
                build['compile_identities']['RUSTRED_SOURCE_DIGEST'],build['compile_identities']['DEPENDENCY_SOURCE_DIGEST']]
            metadata.append({'cut_slots':cut,'path':str(path.relative_to(ROOT)),'sha256':digest(path),
                'proof_path':str(proof_path.relative_to(ROOT)),'proof_sha256':digest(proof_path),
                'program_blake3_verified':row['program_blake3'],'program_bytes_verified':row['program_bytes'],
                'requested_discovery':transaction_binding,'historical_candidate_collection':coverage_binding})
            checked.extend([path,proof_path])
    native={'executable_sha256':captured[0]['sha256'],'source_snapshot_sha256':digest(snapshot_path),
        'compile_identities':build['compile_identities'],
        'configuration':{'source_options':config['source_options'],'guard_digits':config['guard_digits'],
            'native_closure':re.sub(r'checkpoints: Some\((?:\\"|\").*?(?:\\"|\")\)','checkpoints: None',config['native_closure'])},
        'closed_metadata':metadata,
        'physical_zero_certificates':certificates,
        'construction_scope':'massive native AMF throughout'if massive else'certified singleton physical zeros plus actual native multicut AMF',
        'native_digest_tool':{'record_sha256':digest(tool_record_path),**tool_record['executable']}}
    if laurent:
        sample=BASE/(PREFIX+'-'+('three'if three else'two')+'-loop-fixed-comparison.json')
        sample_binding=sample.with_name(sample.name.replace('-comparison.json','-comparison-binding.json'))
        assert read(sample)['status']=='passed'
        old=read(sample_binding)
        assert old['status']=='passed'
        for key in ['executable_sha256','source_snapshot_sha256','compile_identities','configuration','physical_zero_certificates']:
            assert old['native_provenance'][key]==native[key],'Laurent/fixed construction provenance differs: '+key
        for name,sha in old['inputs_sha256'].items():assert digest(ROOT/name)==sha
        checked.extend([sample,sample_binding])
    # No independent values are accessed until all native completion, profile,
    # construction and source checks above have passed.
    plan_path=BASE/'comparison-plan.json'
    plan=read(plan_path)
    for name,sha in plan['frozen_comparator_and_reference_artifacts_sha256'].items():
        assert digest(ROOT/name)==sha,'frozen comparison/reference artifact changed: '+name
        checked.append(ROOT/name)
    for name,sha in plan['report_local_metadata_validators_sha256'].items():
        assert digest(ROOT/name)==sha,'report-local construction validator changed: '+name
        checked.append(ROOT/name)
    checked.append(plan_path)
    spec=plan['gates'][args.gate]
    command=[sys.executable,str(ROOT/spec['comparator'])]
    if massive:command+=['--predictions',str(directory),'--output',str(output)]
    else:
        command+=['laurent'if laurent else'sample','--predictions',str(directory),
                  '--reference',str(ROOT/spec['reference']),'--output',str(output)]
        if three:
            command+=['--run-provenance',str(provenance_path),'--source-snapshot',str(snapshot_path)]
            if laurent:command+=['--sample-comparison',str(sample)]
    hashes={str(path.relative_to(ROOT)):digest(path)for path in checked}
    evidence={'status':'comparison_running','scope':spec['scope'],'command':command,
        'native_provenance':native,'inputs_sha256':hashes,
        'definitions_and_tolerances_changed':False,'production_modified':False}
    evidence_path.write_text(json.dumps(evidence,indent=2)+'\n')
    result=subprocess.run(command,cwd=ROOT)
    assert all(digest(ROOT/name)==sha for name,sha in hashes.items())
    assert digest(digest_executable)==tool_record['executable']['sha256']
    verify_snapshot(snapshot_path)
    evidence['exit_code']=result.returncode
    evidence['status']='passed'if result.returncode==0 else'failed'
    if output.exists():
        evidence['comparison_sha256']=digest(output)
        if result.returncode==0:assert read(output)['status']=='passed'
    else:assert result.returncode!=0
    evidence_path.write_text(json.dumps(evidence,indent=2)+'\n')
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
