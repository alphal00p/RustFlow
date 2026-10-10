#!/usr/bin/env python3
"""Compare complete saved selected-sector partial-flow predictions only.

No full-amplitude or Laurent acceptance follows from this selected [0,3] cut.
The report validator binds the native owners; it is not a second theorem prover.
"""
import json,re,subprocess,sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from frozen_sources import BASE,ROOT,digest,verify_snapshot
from history_coverage import verify_history_coverage
from memo_validation import CAPS,COUNTERS,PEAKS,validate_caps
from compare_three_loop_reference import reference_adapter,TARGETS,value,difference

PREFIX='partial-flow'
STEM=PREFIX+'-three-loop-double-partial'
PROFILES=[[18,60,8],[28,60,8],[28,80,8],[28,80,12]]
CUTS=[0,3]
SHIFTED=[1,2]
SOURCE_OPTIONS={'policy':'polynomial-closure-v1','positive_compact_energy_powers':False,'free_virtual_zero_sectors':True}
def read(path):return json.loads(path.read_text())
def bound(path):return {'path':str(path.relative_to(ROOT)),'sha256':digest(path)}
def hexhash(s):return isinstance(s,str) and re.fullmatch('[0-9a-f]{64}',s) is not None

def validate_partial(evaluation,closure):
    assert evaluation['construction']=='weighted_amf'
    assert evaluation['source_options']==SOURCE_OPTIONS
    assert evaluation['shifted_slots']==SHIFTED
    assert evaluation['massless_endpoint'] is None,'partial proof must remain distinct from all-shift evidence'
    assert evaluation['basis_size']==len(closure['basis'])>0
    assert evaluation['native_storage_capacity']==16
    arity=evaluation['physical_arity']
    assert type(arity)is int and 0<arity<=16
    assert all(len(label)==arity for label in closure['basis'])
    assert evaluation['nonzero_conditions']==closure['nonzero_conditions']
    p=evaluation['partial_continuation'];assert isinstance(p,dict)
    proof=p['proof'];audit=p['physical_consumers'];sources=p['universal_source_class']
    assert proof['version']=='partial-continuation-h1-k2-v1'
    assert proof['shifted_slots']==SHIFTED and proof['source_options']==SOURCE_OPTIONS
    assert isinstance(proof['input_identity'],str) and proof['input_identity']
    assert isinstance(proof['family_signature'],list) and proof['family_signature']
    assert all(isinstance(s,str)and s for s in proof['family_signature'])
    for key in ('source_identity','origin_identity','endpoint_identity'):assert isinstance(proof[key],str)and proof[key]
    assert evaluation['contour_admission']==proof['source_identity']
    for key in ('native_closure','transport_performed','raw_ward_authority','infinity_boundary_authority','internal_cas_work_bound'):assert proof[key]is False
    assert proof['cooperative_cancellation']is True
    requirements=proof['endpoint_requirements']
    for key in ('audit_every_target_basis_candidate_and_derivative_label','require_universal_source_image_class_certificate_for_native_intermediates','retain_raw_target_conditions_before_cancellation_or_zero_classification','check_original_conditions_at_each_epsilon_sample','audit_combined_rational_reconstruction_with_all_eta_poles','require_actual_regular_singular_frobenius_representation_and_sufficient_series','retain_native_joint_target_divergence_checks','right_half_plane_eta_only'):assert requirements[key]is True
    assert requirements['origin_identity']==proof['origin_identity']
    assert audit['continuation_identity']==proof['source_identity']
    assert hexhash(audit['reduced_system_signature'])
    assert len(audit['combined_target_signature'])==2
    for key in ('native_weighted_reconstruction_required','actual_symbolic_frobenius_required','combined_target_projection_required','universal_source_image_certificate_required','boundary_region_certificates_required'):assert audit[key]is True
    assert set(proof['raw_nonzero_conditions'])<=set(audit['nonzero_conditions'])<=set(closure['nonzero_conditions'])
    endpoint=audit['endpoint']
    assert endpoint['endpoint_identity']==proof['endpoint_identity']
    assert endpoint['native_mode_authority']is False
    assert endpoint['labels']
    labels=set()
    for item in endpoint['labels']:
        label=item['label'];assert len(label)==arity and all(type(n)is int for n in label)
        labels.add(tuple(label))
        assert Fraction(item['sigma'])>0
        assert item['strict_margins'] and all(Fraction(m['value'])>0 for m in item['strict_margins'])
    assert set(map(tuple,closure['basis']))<=labels
    assert sources['physical_arity']==arity and sources['storage_capacity']==16
    assert len(sources['admitted_bounds'])==16
    for key in ('binding_blake3','source_only_program_blake3'):assert hexhash(sources[key])
    for key in ('individual_application_trace','one_dimension_for_unbounded_domain','physical_identity_validity_certified_by_this_checker','internal_cas_work_or_heap_bound'):assert sources[key]is False
    assert sources['counts']['sources']>0
    return p

def main():
    directory=BASE/STEM;output=BASE/(STEM+'-comparison.json');evidence_path=BASE/(STEM+'-comparison-binding.json')
    assert not output.exists()and not evidence_path.exists(),'preserve prior comparisons'
    snapshot_path=BASE/(PREFIX+'-source-hashes.json');snapshot=verify_snapshot(snapshot_path)
    build_path=BASE/(PREFIX+'-build-provenance.json');build=read(build_path)
    binding_path=BASE/(STEM+'-binding.json');resources_path=BASE/(STEM+'-resources.json');provenance_path=BASE/(STEM+'-resources.provenance.json')
    binding,resources,provenance=map(read,(binding_path,resources_path,provenance_path))
    assert binding['exit_code']==resources['exit_code']==0
    assert binding['post_run_source_and_executable_unchanged']is True
    assert binding['source_snapshot_sha256']==build['source_snapshot_sha256']==digest(snapshot_path)
    assert binding['build_provenance_sha256']==digest(build_path)
    assert binding['command']==resources['command']==provenance['command']
    assert binding['command'][1]=='1800' and 'runtime_graph_occupied_flow'in binding['command']
    assert binding['full_amplitude']is False and binding['cut_slots']==CUTS and binding['shifted_slots']==SHIFTED
    assert binding['required_fixed_profiles']==PROFILES
    assert binding['experimental_programs_or_rules_imported']is False and binding['predictions_before_reference_access']is True
    assert not(directory/'failure.json').exists()
    executable=binding['executable'];assert Path(executable).name.startswith('finite_density_runtime_flow-')
    assert digest(ROOT/executable)==binding['executable_sha256']==build['files'][executable]['sha256']
    captured=[v for k,v in provenance['argument_files'].items()if Path(k).name.startswith('finite_density_runtime_flow-')]
    assert len(captured)==1 and captured[0]['sha256']==binding['executable_sha256']
    checked=[snapshot_path,build_path,binding_path,resources_path,provenance_path,Path(__file__),BASE/'history_coverage.py',BASE/'memo_validation.py',BASE/'compare_three_loop_reference.py']
    for prerequisite in binding['prerequisite_comparisons'].values():
        path=ROOT/prerequisite['path'];assert digest(path)==prerequisite['sha256'] and read(path)['status']=='passed';checked.append(path)
    assert set(binding['prerequisite_comparisons'])=={'two-loop-fixed','two-loop-laurent','massive-regression'}
    settings=binding['settings']
    expected_settings={'RUSTFLOW_WEIGHTED_INPUT':'examples/finite_density/massless_three_loop_chain.json','RUSTFLOW_WEIGHTED_SOURCE_POLICY':'polynomial-closure-v1','RUSTFLOW_WEIGHTED_SHIFTED_SLOTS':'[1,2]','RUSTFLOW_WEIGHTED_CUTS':'0,3','RUSTFLOW_WEIGHTED_UNIT_MEMO':'true','RUSTFLOW_WEIGHTED_REQUESTED_RAY_POINT':'true','RUSTFLOW_WEIGHTED_RAY_DOMAINS':'1','RUSTFLOW_WEIGHTED_POINT_DOMAINS':'1','RUSTFLOW_WEIGHTED_DOMAINS_PER_RESIDUAL':'0','RUSTFLOW_WEIGHTED_GUARD_PASSES':'0','RUSTFLOW_WEIGHTED_DEPTH':'3','RUSTFLOW_WEIGHTED_DOMAINS':'32768','RUSTFLOW_WEIGHTED_ROUNDS':'32','RUSTFLOW_WEIGHTED_FRONTIER':'8192','RUSTFLOW_WEIGHTED_REQUESTED':'16384','RUSTFLOW_WEIGHTED_REUSED_RULES':'65536','RUSTFLOW_WEIGHTED_DIRECT_ZERO_ATTEMPTS':'8192','RUSTFLOW_WEIGHTED_MAX_HISTORY_CANDIDATE_MAPS':'0','RUSTFLOW_DENSITY_FLOW_EPSILON':'-5/4','RUSTFLOW_DENSITY_FLOW_PROFILES':'18:60:8,28:60:8,28:80:8,28:80:12'}
    assert all(settings[k]==v for k,v in expected_settings.items())
    config_path=directory/'configuration.json';config=read(config_path)
    input_path=ROOT/expected_settings['RUSTFLOW_WEIGHTED_INPUT'];saved_input=directory/'input.json'
    assert config['input_path']==expected_settings['RUSTFLOW_WEIGHTED_INPUT']
    assert digest(input_path)==snapshot['files'][config['input_path']]==binding['input_sha256']
    assert read(saved_input)==read(input_path)
    assert config['profiles']==PROFILES and config['guard_digits']==40 and config['epsilon']=='-5/4'
    assert config['source_options']==SOURCE_OPTIONS and config['mass_mode']=='Propagators([1, 2])'
    assert config['independent_reference_comparisons']==0
    sector_path=directory/'sector.json';sector=read(sector_path)
    assert sector=={'cut_slots':CUTS,'storage_capacity':16,'source_options':SOURCE_OPTIONS,'full_amplitude':False}
    closure_path=directory/'closure.json';closure=read(closure_path)
    validate_caps(config['unit_reduction_memo'],True)
    diagnostics=closure['diagnostics']
    assert all(type(diagnostics[k])is int and diagnostics[k]>=0 for k in COUNTERS)
    assert diagnostics['native_rule_applications']==diagnostics['performed_native_rule_applications']+diagnostics['memoized_rule_applications']
    assert all(diagnostics[peak]<=CAPS[cap]for peak,cap in PEAKS.items())
    paths=[directory/('prediction-'+'-'.join(map(str,p))+'.json')for p in PROFILES]
    assert set(directory.glob('prediction-*.json'))==set(paths)
    predictions=[];partial=None
    for path,profile in zip(paths,PROFILES):
        row=read(path)
        assert row['full_amplitude']is False and row['cut_slots']==CUTS
        assert row['epsilon']=='-5/4'and row['normalization']=='unscaled Euclidean cut contribution'
        assert row['independent_reference_comparisons']==0 and row['guard_digits']==40
        assert [row[k]for k in ('digits','series_order','occupied_start_scale')]==profile
        p=validate_partial(row['report'],closure)
        assert partial is None or partial==p,'physical proof binding changed between profiles'
        partial=p;assert len(row['report']['values'])==2;predictions.append(row)
    checked.extend([config_path,input_path,saved_input,sector_path,closure_path,*paths])
    tool_path=BASE/'native-digest-build.json';tool=read(tool_path);digest_exe=Path(tool['executable']['path'])
    assert tool['exit_code']==0 and tool['source_snapshot_sha256']==digest(snapshot_path) and tool['build_provenance_sha256']==digest(build_path)
    assert digest(digest_exe)==tool['executable']['sha256']
    checked.extend([tool_path,BASE/'native-digest.rs',BASE/'build-native-digest.py'])
    def blake3(path):
        before=path.stat();parts=subprocess.check_output([str(digest_exe),str(path)],text=True).strip().split();after=path.stat()
        assert len(parts)==2 and hexhash(parts[0]) and int(parts[1])==before.st_size==after.st_size
        assert (before.st_ino,before.st_mtime_ns)==(after.st_ino,after.st_mtime_ns)
        return parts[0]
    # Run.occupied writes directly to its native-closure directory, without the
    # full-assembly cut subdirectory. Check only this explicitly selected sector.
    closed_paths=list((directory/'native-closure').glob('round-*-closed.json'))
    assert len(closed_paths)==1
    closed_path=closed_paths[0];closed=read(closed_path);proof_path=closed_path.with_suffix('.bin')
    assert closed['status']=='closed'and closed['active_target_closure']is True
    assert closed['discovery_schedule']=='requested-ray-point-v1'
    assert closed['requested_ray_point']=={'max_domains_per_ray':1,'max_domains_per_point':1}
    validate_caps(closed['unit_reduction_memo'],True)
    assert closed['max_history_candidate_maps']==0 and closed['storage_capacity']==16
    assert closed['physical_arity']==predictions[0]['report']['physical_arity']
    assert proof_path.stat().st_size==closed['program_bytes'] and blake3(proof_path)==closed['program_blake3']
    provisional_path=closed_path.with_name(closed_path.name.replace('-closed.json','-provisional.json'))
    provisional=read(provisional_path);provisional_proof=provisional_path.with_suffix('.bin');retired_path=closed_path.parent/'retired-unresolved.json'
    history=verify_history_coverage(closed,provisional,0,read(retired_path))
    link=closed['active_state']['requested_discovery'];assert re.fullmatch(r'requested-discovery-[0-9]+[.]json',link['path'])
    transaction_path=closed_path.parent/link['path'];transaction=read(transaction_path)
    assert blake3(transaction_path)==link['blake3']
    assert transaction['schema']==1 and transaction['status']=='complete'and transaction['schedule']=='requested-ray-point-v1'
    assert transaction['source_measure_id']==closed['measure_id'] and transaction['physical_arity']==closed['physical_arity'] and transaction['storage_capacity']==16
    assert transaction['policy']==closed['requested_ray_point'] and transaction['empty_terminals']is True and transaction['completed_search_implies_coverage']is False
    assert provisional_proof.stat().st_size==provisional['program_bytes']==transaction['native_program']['bytes']
    assert blake3(provisional_proof)==provisional['program_blake3']==transaction['native_program']['blake3']
    fields=closed['measure_id'].split(':',3)
    assert fields[:3]==['rustflow-weighted-sources-v1',build['compile_identities']['RUSTRED_SOURCE_DIGEST'],build['compile_identities']['DEPENDENCY_SOURCE_DIGEST']]
    checked.extend([closed_path,proof_path,provisional_path,provisional_proof,retired_path,transaction_path])
    # All producer completion/proof/configuration checks precede reading values
    # from the frozen independently generated reference.
    plan_path=BASE/'comparison-plan.json';plan=read(plan_path)
    for name,sha in plan['frozen_comparator_and_reference_artifacts_sha256'].items():assert digest(ROOT/name)==sha
    reference_path=ROOT/plan['gates']['three-loop-fixed']['reference']
    reference,reference_inputs=reference_adapter(reference_path)
    assert reference['source']['input']==config['input_path'] and reference['source']['input_sha256']==digest(input_path)
    matches=[s for s in reference['samples']if Fraction(s['dimension'])==Fraction(13,2) and Fraction(s['chemical_potential'])==1]
    assert len(matches)==1 and len(matches[0]['double_cut'])==2
    expected=list(map(value,matches[0]['double_cut']))
    checked.extend([plan_path,*reference_inputs])
    hashes={str(p.relative_to(ROOT)):digest(p)for p in checked}
    comparisons=[];refinements=[];previous=None
    for row,path,profile in zip(predictions,paths,PROFILES):
        actual=list(map(value,row['report']['values']))
        for i,target in enumerate(TARGETS):
            comparisons.append({'target':target,'cut_slots':CUTS,'configuration':profile,'prediction_path':str(path.relative_to(ROOT)),**difference(actual[i],expected[i]),'imaginary_zero_passed':abs(actual[i][1])<=Decimal('1e-25')})
            if previous is not None:refinements.append({'target':target,'to_configuration':profile,**difference(actual[i],previous[i])})
        previous=actual
    passed=all(r['passed']and r['imaginary_zero_passed']for r in comparisons)and all(r['passed']for r in refinements)
    verify_snapshot(snapshot_path)
    assert all(digest(ROOT/p)==sha for p,sha in hashes.items())
    result={'schema':1,'status':'passed'if passed else'failed','scope':'Selected three-loop double cut [0,3], partial shifted slots [1,2], both original targets, four fixed-D profiles. No full-amplitude three-loop, Laurent or four-loop acceptance.','full_amplitude':False,'epsilon':'-5/4','cut_slots':CUTS,'shifted_slots':SHIFTED,'configurations':PROFILES,'reference_comparison_count':len(comparisons),'independent_refinement_comparison_count':len(refinements),'relative_tolerance':'1e-12','small_magnitude_threshold':'1e-20','small_absolute_tolerance':'1e-25','imaginary_zero_absolute_tolerance':'1e-25','reference_uncertainty':reference['uncertainty'],'reference_error_is_empirical':True,'comparisons':comparisons,'refinements':refinements,'inputs_sha256':hashes,'native_proof':{'closed':bound(closed_path),'program':bound(proof_path),'program_blake3':closed['program_blake3'],'requested_discovery':bound(transaction_path),'historical_candidate_collection':history},'partial_continuation':partial,'memo_counters':{k:diagnostics[k]for k in COUNTERS}}
    output.write_text(json.dumps(result,indent=2)+'\n')
    evidence_path.write_text(json.dumps({'status':result['status'],'inputs_sha256':hashes,'comparison':bound(output),'source_snapshot':bound(snapshot_path),'build_provenance':bound(build_path),'script_sha256':digest(Path(__file__))},indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ('status','reference_comparison_count','independent_refinement_comparison_count')}))
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
