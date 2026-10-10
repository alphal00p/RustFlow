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

PREFIX='polynomial-closure-v2'
PROFILES=[[18,60,8],[28,60,8],[28,80,8],[28,80,12]]

def read(path):return json.loads(path.read_text())

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gate',choices=['three-loop-fixed','three-loop-laurent'])
    args=parser.parse_args()
    stem=PREFIX+'-'+args.gate+'-guards3-per-residual1'
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
    assert binding['settings']['RUSTFLOW_WEIGHTED_GUARD_PASSES']=='3'
    assert binding['settings']['RUSTFLOW_WEIGHTED_DOMAINS']=='8192'
    assert binding['settings']['RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL']=='1'
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
            for item in row['occupied_reports']:
                evaluation=item.get('report',item)
                assert evaluation['construction']=='weighted_amf'
                assert evaluation['source_options']==config['source_options']
                if not massive:assert evaluation['massless_endpoint']is not None
    checked=[snapshot_path,build_path,binding_path,resources_path,provenance_path,
             config_path,input_path,directory/'input.json',*prediction_paths,Path(__file__)]
    metadata=[]
    for cut in cuts[1:]:
        paths=sorted((directory/'native-closure'/('cut-'+'-'.join(map(str,cut)))).glob('round-*-closed.json'))
        assert paths,'every occupied cut needs saved successful native closure evidence'
        for path in paths:
            row=read(path)
            assert row['status']=='closed'
            proof_path=path.with_suffix('.bin')
            assert proof_path.is_file(),'closed metadata must retain its native proof binary'
            fields=row['measure_id'].split(':',3)
            assert fields[:3]==['rustflow-weighted-sources-v1',
                build['compile_identities']['RUSTRED_SOURCE_DIGEST'],build['compile_identities']['DEPENDENCY_SOURCE_DIGEST']]
            metadata.append({'cut_slots':cut,'path':str(path.relative_to(ROOT)),'sha256':digest(path),
                'proof_path':str(proof_path.relative_to(ROOT)),'proof_sha256':digest(proof_path)})
            checked.extend([path,proof_path])
    native={'executable_sha256':captured[0]['sha256'],'source_snapshot_sha256':digest(snapshot_path),
        'compile_identities':build['compile_identities'],
        'configuration':{'source_options':config['source_options'],'guard_digits':config['guard_digits'],
            'native_closure':re.sub(r'checkpoints: Some\((?:\\"|\").*?(?:\\"|\")\)','checkpoints: None',config['native_closure'])},
        'closed_metadata':metadata}
    if laurent:
        sample=BASE/(PREFIX+'-'+('three'if three else'two')+'-loop-fixed-guards3-per-residual1-comparison.json')
        sample_binding=sample.with_name(sample.name.replace('-comparison.json','-comparison-binding.json'))
        assert read(sample)['status']=='passed'
        old=read(sample_binding)
        assert old['status']=='passed'
        for key in ['executable_sha256','source_snapshot_sha256','compile_identities','configuration']:
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
        'definitions_and_tolerances_changed':False,'production_modified':False,
        'retry_configuration':{'guard_passes':3,'shared_native_domains':8192,'domains_per_residual':1},
        'wrapper_template_sha256':'2257a2b25c3d96cddabcdef8bc4f156b5c5e8216f5f28595f3bba6f061d4c12f'}
    evidence_path.write_text(json.dumps(evidence,indent=2)+'\n')
    result=subprocess.run(command,cwd=ROOT)
    assert all(digest(ROOT/name)==sha for name,sha in hashes.items())
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
