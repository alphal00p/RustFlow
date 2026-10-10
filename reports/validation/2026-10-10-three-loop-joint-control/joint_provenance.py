"""Authenticate the captured standalone producer and complete saved phases.

Uses frozen source-archive/evidence identities, never mutable production sources.
Native proof digests remain mandatory. No reference value is read here.
"""
import hashlib,json,re,subprocess,tarfile
from pathlib import Path
from certificate_validation import load_certificates,validate_sample
from history_coverage import verify_history_coverage
BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
CUTS=[[],[0],[3],[0,3]];PROFILES=[(18,60,8),(28,60,8),(28,80,8),(28,80,12)]
read=lambda p:json.loads(p.read_text())
def digest(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def verify_joint(mode):
    assert mode in ['sample','laurent']
    cap_path=BASE/'capsule-binding.json';cap=read(cap_path);run_path=BASE/'run-binding.json';run=read(run_path)
    resources_path=BASE/'resources.json';resources=read(resources_path);provenance_path=BASE/'resources.provenance.json';provenance=read(provenance_path)
    assert run['exit_code']==resources['exit_code']==0 and run['post_run_static_capsule_unchanged']is True
    assert run['capsule_binding_sha256']==digest(cap_path)
    assert run['launcher_sha256']==digest(BASE/'run.py')
    assert run['numerical_prerequisites_sha256']==digest(BASE/'numerical-prerequisites.json')
    assert cap['prepare_script_sha256']==digest(BASE/'build-capsule.py')
    assert run['command']==resources['command']==provenance['command']
    assert not run['imported_rules_or_checkpoints']and not run['reference_values_loaded_by_producer']
    assert not run['prerequisite_values_passed_to_child']
    paths=[cap_path,run_path,resources_path,provenance_path,BASE/'harness-binding.json',BASE/'build.json',BASE/'numerical-prerequisites.json',BASE/'input.json',Path(__file__),BASE/'certificate_validation.py',BASE/'history_coverage.py',BASE/'run.py',BASE/'build-capsule.py',BASE/'capture-prerequisites.py']
    for row in cap['static_files']:
        path=Path(row['path']);assert digest(path)==row['sha256'];paths.append(path)
    prerequisites=read(BASE/'numerical-prerequisites.json')
    assert prerequisites['all_three_comparisons_passed']is True
    assert prerequisites['capture_script_sha256']==digest(BASE/'capture-prerequisites.py')
    assert prerequisites['source_snapshot_sha256']==cap['source_snapshot_sha256']
    assert prerequisites['build_provenance_sha256']==cap['build_provenance_sha256']
    for row in cap['evidence_copies']+prerequisites['captured_files']:
        path=BASE/row['copy_path'];assert digest(path)==row['sha256'];paths.append(path)
    evidence={row['source_path']:BASE/row['copy_path']for row in cap['evidence_copies']}
    snapshot_path=next(p for name,p in evidence.items()if name.endswith('singleton-physical-zero-v2-source-hashes.json'))
    build_path=next(p for name,p in evidence.items()if name.endswith('singleton-physical-zero-v2-build-provenance.json'))
    gates_path=next(p for name,p in evidence.items()if name.endswith('singleton-physical-zero-v2-root-gates.json'))
    snapshot=read(snapshot_path);build=read(build_path);gates=read(gates_path)
    assert digest(snapshot_path)==cap['source_snapshot_sha256']==run['source_snapshot_sha256']
    assert digest(build_path)==cap['build_provenance_sha256']==run['build_provenance_sha256']
    assert digest(gates_path)==cap['root_gates_sha256']==run['root_gates_sha256']
    assert gates['status']=='passed'and gates['build_provenance']['sha256']==digest(build_path)
    assert gates['compile_identities']==cap['compile_identities']==build['compile_identities']
    archive=BASE/cap['source_archive']['path'];assert digest(archive)==cap['source_archive']['sha256'];paths.append(archive)
    with tarfile.open(archive,'r:gz')as tar:
        members=[m for m in tar.getmembers()if m.isfile()]
        assert {m.name for m in members}==set(snapshot['files'])
        for m in members:assert hashlib.file_digest(tar.extractfile(m),'sha256').hexdigest()==snapshot['files'][m.name]
        original=tar.extractfile('tests/finite_density_runtime_flow.rs').read()
        definition=json.loads(tar.extractfile('examples/finite_density/massless_three_loop_chain.json').read())
    harness=read(BASE/'harness-binding.json');standalone=read(BASE/'build.json')
    assert digest(BASE/'harness-binding.json')==cap['harness_binding_sha256']
    assert digest(BASE/'build.json')==cap['standalone_build_sha256']==run['standalone_harness_build_sha256']
    assert standalone['exit_code']==0 and standalone['cargo_invoked']is False and standalone['shared_target_modified']is False
    assert standalone['production_library_build_provenance_sha256']==digest(build_path)
    assert hashlib.sha256(original).hexdigest()==harness['original_sha256']==snapshot['files']['tests/finite_density_runtime_flow.rs']
    generated=BASE/'joint-runtime.rs';assert digest(generated)==harness['generated_sha256']==standalone['source_sha256'];paths.append(generated)
    assert generated.read_bytes().startswith(original)
    original_text=original.decode();generated_text=generated.read_text()
    for name,key in [('runtime_graph_full_amplitude','fixed_phase_body_sha256'),('runtime_graph_full_laurent','laurent_phase_body_sha256')]:
        marker='fn '+name+'() {';assert original_text.count(marker)==1
        body=original_text.split(marker,1)[1].split('\n#[test]',1)[0].rstrip();assert body.endswith('}')
        body=body[:-1];prefix='\n    let run = Run::new();\n    let flow = run.full();';assert body.startswith(prefix)
        phase=body[len(prefix):];assert hashlib.sha256(phase.encode()).hexdigest()==harness[key]
        assert generated_text[len(original_text):].count(phase)==1
    joined=generated_text[len(original_text):]
    assert joined.count('let run = Run::new();')==joined.count('let flow = run.full();')==1
    assert digest(BASE/'generate-harness.py')==harness['generator_sha256'];paths.append(BASE/'generate-harness.py')
    for name in ['symbolica_amflow','rustred']:
        expected=[v['sha256']for k,v in build['files'].items()if Path(k).name.startswith('lib'+name+'-')and k.endswith('.rlib')]
        assert expected==[standalone['libraries'][name]['sha256']]
    runtime=Path(cap['capsule_path'])/'runtime'
    assert digest(runtime)==run['executable_sha256']
    captures=[v for k,v in provenance['argument_files'].items()if Path(k)==runtime]
    assert len(captures)==1 and captures[0]['sha256']==digest(runtime)
    directory=BASE/'predictions';assert not(directory/'failure.json').exists()
    assert read(directory/'input.json')==definition==read(BASE/'input.json')
    assert digest(BASE/'input.json')==run['input_sha256']
    config_path=directory/'configuration.json';config=read(config_path);paths.extend([config_path,directory/'input.json'])
    assert config['input_path']==str(Path(cap['capsule_path'])/'input.json')
    assert config['profiles']==[list(p)for p in PROFILES]and config['epsilon']=='-5/4'
    assert config['guard_digits']==40 and config['initial_epsilon_grid']==1000 and config['independent_reference_comparisons']==0
    assert config['source_options']=={'policy':'polynomial-closure-v1','positive_compact_energy_powers':False,'free_virtual_zero_sectors':True}
    expected={'RUSTFLOW_WEIGHTED_ACTIVE_TARGET_CLOSURE':'true','RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT':'true','RUSTFLOW_WEIGHTED_RAY_DOMAINS':'1','RUSTFLOW_WEIGHTED_POINT_DOMAINS':'1','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'0','RUSTFLOW_WEIGHTED_GUARD_PASSES':'0','RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'32768','RUSTFLOW_WEIGHTED_ROUNDS':'32','RUSTFLOW_WEIGHTED_FRONTIER':'8192','RUSTFLOW_WEIGHTED_REQUESTED':'16384','RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'8192','RUSTFLOW_WEIGHTED_MAX_HISTORY_CANDIDATE_MAPS':'0','RUSTFLOW_WEIGHTED_SOURCE_POLICY':'polynomial-closure-v1','RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED':'false','RUSTFLOW_WEIGHTED_REQUESTED_RAYS':'true','RUSTFLOW_WEIGHTED_FREE_VIRTUAL_ZEROS':'true','RUSTFLOW_WEIGHTED_ZERO_FACES':'true','RUSTFLOW_DENSITY_FLOW_EPSILON':'-5/4','RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8,28:60:8,28:80:8,28:80:12','RUSTFLOW_DENSITY_FLOW_GRID':'1000'}
    assert all(run['settings'].get(k)==v for k,v in expected.items())
    for filename,fields in [('joint-harness.json',{'prepared_flow_count':1,'phase_order':['fixed_dimension','laurent'],'original_runtime_phase_bodies':True,'persisted_rules_loaded':False,'reference_values_loaded':False}),('fixed-phase-completed.json',{'predictions_saved':True,'reference_comparison_performed':False}),('joint-phases-completed.json',{'both_phases_completed':True,'reference_comparison_performed':False})]:
        p=directory/filename;r=read(p);assert all(r[k]==v for k,v in fields.items());paths.append(p)
    certificates,certpath=load_certificates(directory,definition,CUTS);paths.append(certpath)
    predictions=[]
    for digits,order,start in PROFILES:
        p=directory/f'prediction-{digits}-{order}-{start}.json';r=read(p)
        assert r['full_amplitude']is True and r['independent_reference_comparisons']==0
        assert r['normalization']=='unscaled Euclidean amplitude'and r['epsilon']=='-5/4'
        assert [r[k]for k in ['digits','series_order','occupied_start_scale','guard_digits']]==[digits,order,start,40]
        assert [x['cut_slots']for x in r['contributions']]==CUTS
        validate_sample(r,certificates,definition,CUTS,config['source_options']);predictions.append(p)
    for digits,order,start,grid in [(*p,1000)for p in PROFILES]+[(28,80,12,2000)]:
        p=directory/f'prediction-{digits}-{order}-{start}-grid-{grid}.json';r=read(p)
        assert r['full_amplitude']is True and r['independent_reference_comparisons']==0
        assert [r[k]for k in ['digits','series_order','occupied_start_scale','epsilon_grid_denominator','guard_digits']]==[digits,order,start,grid,40]
        assert len(r['expansions'])==2
        assert all(type(e['verified_digits'])is int and e['verified_digits']>=digits and set(e['coefficients'])=={'-3','-2','-1','0'}for e in r['expansions'])
        predictions.append(p)
    assert set(directory.glob('prediction-*.json'))==set(predictions);paths.extend(predictions)
    tool_record=read(BASE/'comparison-digest-binding.json');tool=Path(tool_record['executable']['path']);assert digest(tool)==tool_record['executable']['sha256'];paths.extend([BASE/'comparison-digest-binding.json',tool])
    source_record=ROOT/tool_record['source_record'];assert digest(source_record)==tool_record['source_record_sha256'];paths.append(source_record)
    def blake3(path):
        fields=subprocess.check_output([str(tool),str(path)],text=True).split();assert int(fields[1])==path.stat().st_size;return fields[0]
    for cut in [[0],[3]]:assert not list((directory/'native-closure'/('cut-'+str(cut[0]))).glob('round-*-closed.*'))
    native_dir=directory/'native-closure/cut-0-3';closed_paths=list(native_dir.glob('round-*-closed.json'));assert len(closed_paths)==1
    closed_path=closed_paths[0];closed=read(closed_path);proof=closed_path.with_suffix('.bin')
    assert closed['status']=='closed'and closed['active_target_closure']is True and closed['max_history_candidate_maps']==0
    assert closed['requested_ray_point']=={'max_domains_per_ray':1,'max_domains_per_point':1}and closed['discovery_schedule']=='requested-ray-point-v1'
    assert proof.stat().st_size==closed['program_bytes']and blake3(proof)==closed['program_blake3']
    provisional_path=closed_path.with_name(closed_path.name.replace('-closed','-provisional'));provisional=read(provisional_path);provisional_proof=provisional_path.with_suffix('.bin')
    retired=native_dir/'retired-unresolved.json';coverage=verify_history_coverage(closed,provisional,0,read(retired))
    assert provisional_proof.stat().st_size==provisional['program_bytes']and blake3(provisional_proof)==provisional['program_blake3']
    link=closed['active_state']['requested_discovery'];assert re.fullmatch(r'requested-discovery-[0-9]+[.]json',link['path'])
    transaction_path=native_dir/link['path'];assert blake3(transaction_path)==link['blake3'];transaction=read(transaction_path)
    assert transaction['schema']==1 and transaction['status']=='complete'and transaction['schedule']=='requested-ray-point-v1'
    assert transaction['source_measure_id']==closed['measure_id']and transaction['physical_arity']==closed['physical_arity']and transaction['storage_capacity']==closed['storage_capacity']
    assert transaction['policy']==closed['requested_ray_point']and transaction['empty_terminals']is True and transaction['completed_search_implies_coverage']is False
    assert transaction['native_program']['bytes']==provisional['program_bytes']and transaction['native_program']['blake3']==provisional['program_blake3']
    assert closed['measure_id'].split(':',3)[:3]==['rustflow-weighted-sources-v1',cap['compile_identities']['RUSTRED_SOURCE_DIGEST'],cap['compile_identities']['DEPENDENCY_SOURCE_DIGEST']]
    paths.extend([closed_path,proof,provisional_path,provisional_proof,retired,transaction_path])
    output={'configuration':{'source_options':config['source_options'],'guard_digits':config['guard_digits'],'native_closure':re.sub(r'checkpoints: Some\((?:\\"|\").*?(?:\\"|\")\)','checkpoints: None',config['native_closure'])},'executable_sha256':digest(runtime),'dependency_compile_ids':[cap['compile_identities']['RUSTRED_SOURCE_DIGEST'],cap['compile_identities']['DEPENDENCY_SOURCE_DIGEST']],'source_snapshot_sha256':digest(snapshot_path),'build_provenance_sha256':digest(build_path),'standalone_harness_build_sha256':digest(BASE/'build.json'),'physical_zero_certificates':certificates,'closed_metadata':[{'path':str(closed_path.relative_to(ROOT)),'sha256':digest(closed_path),'proof_sha256':digest(proof),'program_blake3':closed['program_blake3'],'requested_discovery_sha256':digest(transaction_path),'history_coverage':coverage}], 'source_binding_scope':'restored static source archive and coherent copied evidence; no later live ROOT source reads','construction_limit':'One fresh public PreparedDensityFlow; separate singleton physical zeros plus actual native double-cut AMF. Laurent shares the same prepared owner; no individual sector Laurent coefficient report is inferred.'}
    if mode=='laurent':
        fixed=BASE/'sample-comparison.json';assert read(fixed)['status']=='passed'
        previous=read(fixed)['native_provenance']
        for k in ['configuration','executable_sha256','dependency_compile_ids','source_snapshot_sha256','standalone_harness_build_sha256','physical_zero_certificates']:assert previous[k]==output[k]
        paths.append(fixed)
    return output,list(dict.fromkeys(paths))
