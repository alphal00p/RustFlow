#!/usr/bin/env python3
"""Generate/check acceptance definitions and exact graph evidence; never evaluate integrals.

Only whitelisted input fields leave the oracle fixture. Reference expressions are
not parsed here: their native Symbolica parsing belongs to Rust validation tests.
Default --check compares deterministic artifacts; --write regenerates them.
"""
import argparse
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ORACLE = ROOT / 'fixtures/finite_density/oracle_results_supplied.json'
ROUTINGS = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
            [1,0,0,-1],[0,1,0,-1],[0,0,1,-1],[1,-1,0,0],[1,0,-1,0],[1,-1,-1,0]]
GRAPHS = [
    ('chain_of_three_parallel_pairs', 'I37', [[0,3],[3,2],[2,1],[1,0],[3,0],[2,3],[1,2]]),
    ('five_vertex_eight_edge', 'I91', [[0,4],[4,3],[3,2],[2,0],[1,0],[3,1],[2,3],[4,1]]),
    ('triangular_prism', 'I115', [[0,5],[5,4],[1,0],[4,1],[2,3],[4,2],[3,1],[5,2],[3,0]]),
]


def rank(matrix):
    rows = [[Fraction(v) for v in row] for row in matrix]
    if not rows:
        return 0
    pivot = 0
    for column in range(len(rows[0])):
        found = next((r for r in range(pivot, len(rows)) if rows[r][column]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        q = rows[pivot][column]
        rows[pivot] = [v/q for v in rows[pivot]]
        for r in range(pivot+1, len(rows)):
            q = rows[r][column]
            rows[r] = [v-q*w for v,w in zip(rows[r],rows[pivot])]
        pivot += 1
    return pivot


def det(matrix):
    if not matrix:
        return 1
    return sum((-1)**j * matrix[0][j] * det([row[:j]+row[j+1:] for row in matrix[1:]])
               for j in range(len(matrix)))


def connected(vertices, edges, removed_vertices=(), removed_edges=()):
    nodes = set(range(vertices)) - set(removed_vertices)
    if not nodes:
        return True
    seen = {min(nodes)}
    while True:
        following = set(seen)
        for i,(u,v) in enumerate(edges):
            if i in removed_edges or u not in nodes or v not in nodes:
                continue
            if u in seen or v in seen:
                following.update([u,v])
        if following == seen:
            return seen == nodes
        seen = following


def certificate(edges, routing, signatures):
    e = len(edges)
    v = max(max(edge) for edge in edges)+1
    loops = len(routing[0])
    incidence = [[int(head==i)-int(tail==i) for tail,head in edges] for i in range(v)]
    assert connected(v, edges)
    assert e-v+1 == loops
    assert rank(incidence)==v-1 and rank(routing)==loops
    assert all(sum(incidence[a][j]*routing[j][b] for j in range(e))==0
               for a in range(v) for b in range(loops))
    charges = [sum(a*b for a,b in zip(row,signatures)) for row in routing]
    assert all(abs(c)<=1 for c in charges)
    for row in incidence:
        at_vertex = [a*b for a,b in zip(row,charges)]
        assert sum(at_vertex)==0 and sum(x>0 for x in at_vertex)<=1
    assert all(connected(v,edges,removed_vertices=[i]) for i in range(v))
    assert all(connected(v,edges,removed_edges=[i]) for i in range(e))
    # A split into independent denominator blocks would partition routing rows
    # with additive ranks. Exhaust all partitions, fixing edge0 in the first.
    for bits in range(1, 1 << e):
        if not bits&1 or bits==(1<<e)-1:
            continue
        left=[routing[i] for i in range(e) if bits>>i&1]
        right=[routing[i] for i in range(e) if not bits>>i&1]
        assert rank(left)+rank(right)>loops
    adjacency=[[sum({a,b}=={i,j} for a,b in edges) if i!=j else 0
                for j in range(v)] for i in range(v)]
    canonical=min(tuple(adjacency[p[i]][p[j]] for i in range(v) for j in range(i+1,v))
                  for p in itertools.permutations(range(v)))
    minors=[(list(indices),det([routing[i] for i in indices]))
            for indices in itertools.combinations(range(e),loops)]
    nonzero=[(s,d) for s,d in minors if d]
    assert {abs(d) for _,d in nonzero}=={1}
    # Cut complements are explicitly tested as graphs; rank alone is not used.
    charged=[j for j,c in enumerate(charges) if c]
    cuts=[]
    for k in range(1,min(loops,len(charged))+1):
        for selected in itertools.combinations(charged,k):
            if not connected(v,edges,removed_edges=selected):
                continue
            full,d=next((s,d) for s,d in nonzero if set(selected)<=set(s))
            cuts.append({'slots':list(selected),'completion_to_loop_basis':full,
                         'routing_determinant':d,'absolute_jacobian':1})
    return {
        'incidence_convention':'-1 at tail, +1 at head',
        'vertices':v,'physical_edges':e,'loops':loops,'incidence_matrix':incidence,
        'routing_matrix':routing,'incidence_times_routing_is_zero':True,
        'routing_rank':loops,'incidence_rank':v-1,
        'edge_charges':charges,'conserved_fermion_cycles':True,
        'bridges':[],'articulation_vertices':[],
        'independent_denominator_partition_exists':False,
        'partition_test':'Every nontrivial row partition has rank(A)+rank(B)>L; exact rationals.',
        'degree_multiset':sorted(sum(abs(x) for x in row) for row in incidence),
        'canonical_undirected_adjacency_upper_triangle':list(canonical),
        'canonicalization':'Exhaustive vertex permutations preserving edge multiplicity.',
        'full_rank_routing_minor_absolute_values':[1],
        'connected_admissible_occupied_cutsets':cuts,
        'native_validation_status':'Exact certificate independently checked here; native admission is a separate required gate.',
    }


def artifacts():
    raw=ORACLE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == 'd49acc63e40b9fdccf7fb3eddea6e71f211aba0067021e917925274952ad66d1', 'Supplied oracle must remain byte-for-byte unchanged'
    oracle=json.loads(raw)
    assert oracle['schema_version']==3
    records=oracle['four_loop_oracles']+oracle['lower_loop_oracles']
    assert len(oracle['four_loop_oracles'])==95 and len(oracle['lower_loop_oracles'])==13
    assert len(oracle['constants'])==9
    records_by_id={r['id']:r for r in records}
    # Deliberate allowlist: do not copy conventions.constants/coefficient prose,
    # coefficients, constants, uncertainties, det tags, or whole source records.
    definitions={
        'schema_version':1,'purpose':'Input definitions only; not numerical validation results.',
        'basis_maps':oracle['bases'],
        'conventions':{key:oracle['conventions'][key] for key in [
            'expression_syntax','parameters','momenta','dot_product','measure','measure_factor','normalization']},
        'targets':[{'id':r['id'],'basis':r['basis'],'loop_signatures':r['loop_signatures'],
                    'finite_density_target':r['finite_density_target']} for r in records],
    }
    outputs={'examples/finite_density/oracle_definitions.json':definitions}
    entries=[]
    for name,record_id,edges in GRAPHS:
        record=records_by_id[record_id]
        target=record['finite_density_target']
        e=len(edges)
        signatures=record['loop_signatures']
        graph=certificate(edges,ROUTINGS[:e],signatures)
        raised_slot=next(j for j,charge in enumerate(graph['edge_charges']) if charge)
        raised_powers=[1]*e
        raised_powers[raised_slot]=2
        native={
            'name':name,'loops':4,'vertices':graph['vertices'],
            'edges':[{'vertices':ends,'routing':[str(x) for x in ROUTINGS[j]],
                      'mass_squared':'0','charges':[graph['edge_charges'][j]]}
                     for j,ends in enumerate(edges)],
            'loop_charges':[[x] for x in signatures],'chemical_potentials':['1'],
            'numerator_convention':'shifted_euclidean','laurent_orders':[-4,0],'digits':12,
            'targets':[{'powers':target['denominator_powers'][:e],
                        'numerator':target['numerator_polynomial']},
                       {'powers':raised_powers,'numerator':'g1_2^2+g1_3*g2_4'}],
        }
        file=f'examples/finite_density/{name}.json'
        outputs[file]=native
        entries.append({
            'id':name,'definition':file,'physical_slots_zero_based':list(range(e)),
            'oracle_completion_slots_zero_based':list(range(e,10)),
            'completion_slots_are_physical_edges':False,
            'medium_completion_coordinates':['u1','u2','u3','u4'],
            'graph_certificate':graph,
            'targets':[
                {'input_target':0,'oracle_record':record_id,'kind':'polynomial numerator with raised neutral line' if record_id=='I37' else 'polynomial numerator',
                 'reference':{'status':'available','file':'fixtures/finite_density/oracle_results_supplied.json',
                              'record':record_id,'last_supplied_order':0,
                              'method':'Independent user-supplied reference; retain supplied absolute uncertainties and det labels.',
                              'absolute_uncertainties':record.get('uncertainties', {}),
                              'comparison_precision':{'I37':'Exact symbolic coefficients at supplied orders; 10-digit numerical comparison target.', 'I91':'Finite coefficient has supplied absolute uncertainty 1/1000; independent prediction still requires 10-digit stability.', 'I115':'Finite coefficient has supplied absolute uncertainty 1/10^10; independent prediction still requires 10-digit stability.'}[record_id]}},
                {'input_target':1,'kind':'raised occupied line and polynomial numerator',
                 'raised_physical_slot_zero_based':raised_slot,
                 'reference':{'status':'required_not_generated',
                              'method':'Independent complete thermal-contour/cutting evaluation with independent line masses and distributional mass differentiation; independently implemented sector resolution and arbitrary-precision spatial quadrature after subtraction.',
                              'target_absolute_uncertainty':'1/10^12',
                              'attainable_precision':'Not established. Reference generation and error certification are blocking acceptance gates; no 12-digit availability is claimed.'}},
            ],
            'evaluation_status':'not_evaluated','coefficients_compared':[],
        })
        if record_id == 'I37':
            independent = {
                'status':'available_analytic_reference_not_native_acceptance',
                'file':'reports/validation/2026-10-09-finite-density-native-assembly/independent-e7-reference/reference.json',
                'derivation':'docs/finite-density-e7-reference.md',
                'method':'Independent analytically regulated neutral tensor bubbles and compact beta integrals; raised charged line differentiated before massless limit, retaining nonzero upper surface. Gamma-duplication cross-check and exact Laurent expansion.',
                'last_supplied_order':0,
                'orders_below_minus_two':'zero by the derived meromorphic expression',
                'comparison_precision':'Exact analytic coefficients; native 50/80-digit precision and independent epsilon-extrapolant checks passed. No native AMF predictions compared.',
            }
            entries[-1]['targets'][0]['additional_independent_reference'] = {**independent, 'target_index':0}
            entries[-1]['targets'][1]['reference'] = {**independent, 'target_index':1}
        else:
            entries[-1]['targets'][1]['reference'].update({
                'derivation_and_obstacles':'docs/finite-density-missing-references.md',
                'definition_only_diagnostics':'reports/validation/2026-10-09-finite-density-native-assembly/missing-reference-diagnostics.json',
                'next_gate':('Common thermal matching for shared interior pole/discontinuity loci across all 29 cut sectors, then local UV/IR subtraction; independently assigned denominator PVs are not established.' if record_id=='I91' else 'Complete arithmetic, Barnes-quadrature and epsilon-grid refinements of the independent raised-target Laurent reference; the derived eight-sector finite-epsilon sum is already checked, but no full Laurent acceptance is claimed.'),
            })
    # A small massive nonfactorized graph requires complete 0/1/2-cut assembly.
    small_edges=[[0,1],[1,0],[1,0]]
    small_routing=[[1,0],[0,1],[1,-1]]
    small=certificate(small_edges,small_routing,[1,1])
    small_input={'name':'massive_two_loop_sunset','loops':2,'vertices':2,
                 'edges':[{'vertices':edge,'routing':[str(x) for x in small_routing[i]],
                           'mass_squared':['1/4','1/4','1'][i],'charges':[small['edge_charges'][i]]}
                          for i,edge in enumerate(small_edges)],
                 'loop_charges':[[1],[1]],'chemical_potentials':['1'],
                 'numerator_convention':'shifted_euclidean','laurent_orders':[-2,0],'digits':12,
                 'targets':[{'powers':[1,1,1],'numerator':'1'},
                            {'powers':[2,1,1],'numerator':'g1_2+u1*u2'}]}
    outputs['examples/finite_density/massive_two_loop_sunset.json']=small_input
    multi_edges=[[0,1],[1,0],[2,3],[3,2],[1,2],[3,0]]
    multi_routing=[[1,0,0],[1,0,-1],[0,1,0],[0,1,-1],[0,0,1],[0,0,1]]
    multi_certificates=[certificate(multi_edges,multi_routing,signatures)
                        for signatures in [[1,0,0],[0,1,0]]]
    multi_input={'name':'two_independent_chemical_cycles','loops':3,'vertices':4,
                 'edges':[{'vertices':edge,'routing':[str(x) for x in multi_routing[i]],
                           'mass_squared':['1/4','1/4','4/9','4/9','1','2'][i],
                           'charges':[c['edge_charges'][i] for c in multi_certificates]}
                          for i,edge in enumerate(multi_edges)],
                 'loop_charges':[[1,0],[0,1],[0,0]],'chemical_potentials':['1','3/2'],
                 'numerator_convention':'shifted_euclidean','laurent_orders':[-3,0],'digits':12,
                 'targets':[{'powers':[1,1,1,1,1,1],'numerator':'g1_2+u1*u2'},
                            {'powers':[2,1,2,1,1,1],'numerator':'g1_2+u1*u2'}]}
    outputs['examples/finite_density/two_independent_chemical_cycles.json']=multi_input
    three_input=json.loads((ROOT/'examples/finite_density/massless_three_loop_chain.json').read_text())
    three_certificate=certificate([e['vertices'] for e in three_input['edges']],
                                  [[int(x) for x in e['routing']] for e in three_input['edges']], [1,0,0])
    manifest={
        'schema_version':1,'status':'acceptance_plan_not_numerical_acceptance',
        'source_baseline':'b3a4843e8327835d1ec3ada5a6f32f1841bab2c2',
        'oracle':{'file':'fixtures/finite_density/oracle_results_supplied.json',
                  'sha256':hashlib.sha256(raw).hexdigest(),'schema_version':3,
                  'four_loop_records':95,'lower_loop_records':13,'constants':9,
                  'mandatory_four_loop_subset':['I37','I91','I115'],
                  'mandatory_lower_loop_subset':['L2_1'],
                  'actually_compared':[],
                  'input_only_export':'examples/finite_density/oracle_definitions.json'},
        'precision':{'reference_normalization':'(4*pi)^(-2*L)*(Lambda_bar/2)^(2*L*eps)',
                     'keep_normalization_unexpanded':True,
                     'nonzero_relative_stability':'abs(c_refined-c_baseline)/max(abs(c_refined),abs(c_baseline)) <= 1/10^10',
                     'zero_absolute_stability':'abs(c) <= 1/10^12 in the stated reference normalization, including an error bound <= 1/10^12',
                     'independent_refinements':['working precision','boundary series order','boundary start scale','epsilon sample grid'],
                     'suggested_decimal_working_precision':[40,60,80],
                     'reference_agreement':'abs(prediction-reference) <= supplied_absolute_uncertainty + prediction_absolute_error; a det tag adds no confidence model.',
                     'unknown_higher_orders':'never compare above the greatest supplied order; lower unlisted orders are zero'},
        'native_storage_capacity':{'status':'implemented_and_validated_for_admitted_massive_flows',
                                   'runtime_capacities':[7,9,12,16,20,24,32],
                                   'selection':'Preserve exact physical arities 7 and 9; otherwise choose the smallest compiled capacity that fits. A separately compiled const-generic Rust capacity can exceed the runtime list.',
                                   'physical_owner':'Only real factors enter numerator conversion and source generation; storage-only Ordinary axes have fixed-zero guards and zero shifts. Physical exports verify and trim every tail.',
                                   'cache_identity':'zero-tail-storage-v1:physical=p:capacity=N only for p<N; exact p=N identity preserved',
                                   'default_closure':{'max_rounds':12,'max_depth':3,'max_domains':8192,'search_frontier_sectors':True,'guard_refinement':{'max_passes':3,'max_added_domains':256,'max_interval_width':2}},
                                   'validation':{'physical_arities':[7,9],'compared_capacity':12,'basis_sizes':[6,11],'target_value_comparisons':4,'relative_tolerance':'1e-15','physical_basis_matrix_target_weights_and_conditions':'exactly_identical','report':'reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-physical-comparison.json','resources':'reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-physical-resources.json','wall_seconds':245.70878702914342,'peak_rss_kib':55332,'fast_regression_tests_passed':33,'fast_regression_report':'reports/validation/2026-10-09-finite-density-native-assembly/native-capacity-fast-gates.json'},
                                   'scope':'Storage capacity is a backend resource bound; it does not admit massless flowing graphs, guarantee closure at a chosen budget, or establish a four-loop numerical reference comparison.',
                                   'documentation':'docs/finite-density-native-capacity.md'},
        'four_loop_families':entries,
        'structural_distinction':{'verified':True,'evidence':'Explicit connected, bridge-free, articulation-free incidence certificates with (V,E)=(4,7),(5,8),(6,9), different degree multisets and canonical multigraph adjacency words. Exhaustive exact row-rank partitions rule out denominator-block factorization. Physical edge count is preserved under graph relabeling and invertible loop rerouting; completion slots are excluded.'},
        'massless_structural_admission':{'status':'implemented_singleton_germs_and_independent_rank_one_virtual_blocks',
                                       'scope':'All physical masses zero; positive occupied chemical magnitudes; polynomial completion factors; all uncut physical factors shifted. Singleton sectors use a uniformly UV-continued external-invariant germ at admitted virtual rank; two-cut sectors use independent rank-one virtual blocks and pure spacelike transfers. Exact generated-label origin and endpoint conditions remain mandatory.',
                                       'E7':'All cuts structurally admitted; bounded native closure remains unresolved (single frontier147 and double frontier88); no native four-loop predictions.',
                                       'remaining':'General coupled virtual blocks, E8 common thermal pinches and partial deformation placements remain unadmitted.',
                                       'current_numerical_regression':'reports/validation/2026-10-09-finite-density-native-assembly/massless-blocks-numerical-validation.json',
                                       'current_provenance':'reports/validation/2026-10-09-finite-density-native-assembly/massless-blocks-numerical-provenance.json',
                                       'current_independent_reference_comparisons':80,'current_independent_refinements':54,'historical_same_profile_regressions':10,
                                       'four_loop_native_predictions':0,'supplied_oracle_numerical_records_compared':0},
        'massless_small_graph':{'definition':'examples/finite_density/massless_two_loop_sunset.json',
                               'reference_status':'independent_Gamma_Beta_and_convergent_D7_quadrature_generated',
                               'reference':'reports/validation/2026-10-09-finite-density-native-assembly/independent-massless-reference/independent-massless-reference.json',
                               'native_status':'complete_amplitude_passed_at_D13_over2_D15_over4_and_Laurent_minus2_through_zero',
                               'report':'reports/validation/2026-10-09-finite-density-native-assembly/massless-native-validation.json',
                               'provenance':'reports/validation/2026-10-09-finite-density-native-assembly/massless-native-provenance.json',
                               'sample_dimensions':['13/2','15/4'],'sample_profiles_per_dimension':4,'laurent_profiles':5,
                               'reference_comparisons':110,'independent_refinement_comparisons':84,
                               'occupied_basis_sizes':[2,2,3],
                               'relative_tolerance':'1e-12','small_magnitude_threshold':'1e-20','small_absolute_tolerance':'1e-25','imaginary_zero_absolute_tolerance':'1e-25',
                               'largest_laurent_relative_reference_discrepancy':'2.041136525965978e-44',
                               'largest_laurent_absolute_expected_zero_residual':'8.661497055645898e-56',
                               'native_D7':'additional_indicial_resonance_before_any_prediction; current vacuum-zero pruning and singleton flows succeed, remaining failure is double-cut infinity recurrence; independent convergent_D7_reference_valid',
                               'admitted_scope':'This numerically validated two-loop example lies inside the current sealed singleton-germ and independent-rank-one-block classes; polynomial completions, all uncut physical factors shifted, bound common thermal origin and finite-label endpoint proof.',
                               'binary_provenance_limitation':'Original native executable digest was not captured before relink; source277 snapshot and run/build evidence retained. Final-order executable digest captured and all10 matched-profile values are exact decimal-string equal.',
                               'four_loop_native_predictions':0,'supplied_oracle_numerical_records_compared':0},
        'three_loop_intermediate':{
            'definition':'examples/finite_density/massless_three_loop_chain.json',
            'input_sha256':hashlib.sha256((ROOT/'examples/finite_density/massless_three_loop_chain.json').read_bytes()).hexdigest(),
            'graph_certificate':three_certificate,
            'purpose':'User-requested full three-loop validation checkpoint before further four-loop numerical scale-up; does not replace any four-loop mandatory target.',
            'targets':three_input['targets'], 'required_cut_subsets':[[],[0],[3],[0,3]],
            'reference_status':'independent_analytic_reference_generated_and_checked',
            'reference':'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-independent-reference/reference.json',
            'reference_validation':'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-independent-reference/validation.json',
            'derivation':'docs/finite-density-three-loop-reference.md',
            'reference_method':'Independent Gaussian first moment and compact Beta integrals, with fixed original mixed medium numerator, independent occupied mass jet and separate nonzero moving upper surface. No AMF predictions or supplied oracle answers read.',
            'reference_dimensions':['31/3','13/2','15/4'], 'laurent_orders':[-3,0],
            'reference_precision_digits':[50,80], 'arithmetic_refinement_checks':14,
            'exact_symbolic_reference_identity_checks':5,
            'exact_reference_validation':'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-independent-reference/exact-identities-validation.json',
            'largest_relative_arithmetic_change':'4.5151223184199595e-54',
            'direct_original_kernel_quadrature_orders':[24,40,64,96],
            'largest_final_direct_quadrature_relative_discrepancy':'8.132611430182034e-15',
            'reference_uncertainty':'Exact analytic formula with empirical arithmetic/quadrature validation; no rigorous interval enclosure asserted.',
            'native_status':'full_amplitude_not_yet_validated; bounded rule-union controls remained unclosed and three subsequent residual-focused controls, including the matched requested-index-priority control, timed out at 240 seconds before predictions',
            'native_controls_report':'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-residual32-controls.json',
            'historical_rule_union_controls':'reports/validation/2026-10-09-finite-density-native-assembly/three-loop-rule-union-controls.json',
            'native_prediction_comparisons':0,
            'comparison_script':'tools/finite_density/compare_three_loop_reference.py',
            'required_native_profiles':{'fixed_dimension':4,'laurent':5},
            'relative_tolerance':'1e-12','small_magnitude_threshold':'1e-20','small_absolute_tolerance':'1e-25',
            'production_requirement':'All occupied sectors use genuine native weighted AMF; the reference Gaussian/Beta expressions are validation-only.'},
        'small_graph':{'definition':'examples/finite_density/massive_two_loop_sunset.json',
                       'graph_certificate':small,'required_cut_components':[0,1,2],
                       'vacuum_contribution':'Nonzero for positive masses; must evaluate and compare, not omit as scaleless.',
                       'reference_status':'available_at_D_12_over_5_and_Laurent_minus2_through_zero',
                       'reference_method':'Independent Schwinger vacuum sectors, reference-only Feynman parameters and compact radial/angular quadrature with exact endpoint maps; original numerator and independent mass derivatives retain moving-support surfaces. Laurent extension uses six-sector analytic UV subtraction and exact rational Lagrange coefficients for epsilon² I.',
                       'reference_artifact':'reports/validation/2026-10-09-finite-density-native-assembly/independent-reference/independent-reference.json',
                       'laurent_reference_artifact':'reports/validation/2026-10-09-finite-density-native-assembly/independent-laurent-reference/independent-laurent-reference.json',
                       'laurent_reference_orders':[-2,0],
                       'reference_derivation':'docs/finite-density-reference.md',
                       'attainable_comparison_precision':'All 12 reference components/totals pass independent node and precision refinements at D=12/5. All 108 Laurent component-coefficient refinement comparisons and six analytic UV-residue checks pass; largest relative quadrature change 3.48e-17, epsilon-grid change 1.0042e-22 and precision change 8.98e-53. These empirical changes are not rigorous interval bounds. Native fixed-dimension and complete Laurent comparisons pass and are recorded separately.',
                       'native_comparison':{'status':'full_amplitude_and_all_sectors_passed_at_D_12_over_5','report':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-frontier-reference-comparison.json','unique_sector_target_values':10,'configurations':4,'reference_comparisons':40,'refinement_comparisons':30,'relative_tolerance':'1e-12','small_magnitude_threshold':'1e-20','small_absolute_tolerance':'1e-25','largest_relative_reference_discrepancy':'9.44882265757e-36','largest_relative_refinement_difference':'9.83297e-50','full_amplitude_compared':True,'laurent_coefficients_compared':0,'accuracy_caution':'Observed agreement; independent reference refinement estimates are not rigorous interval bounds.'},
                       'native_laurent_comparison':{'status':'passed_complete_massive_amplitude',
                                                    'report':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-laurent-reference-comparison.json',
                                                    'orders':[-2,0],'targets':2,'configurations':5,
                                                    'reference_comparisons':30,'independent_refinement_comparisons':24,
                                                    'configurations_digits_order_start_grid':[[18,60,8,1000],[28,60,8,1000],[28,80,8,1000],[28,80,12,1000],[28,80,12,2000]],
                                                    'guard_digits':40,'search_frontier_sectors':True,
                                                    'relative_tolerance':'1e-12','small_magnitude_threshold':'1e-18','small_absolute_tolerance':'1e-20','imaginary_zero_absolute_tolerance':'1e-20',
                                                    'largest_relative_reference_discrepancy':'1.429870667895067e-18',
                                                    'largest_relative_refinement_difference':'9.249528747539574e-45',
                                                    'amf_oracle_comparisons':0,
                                                    'accuracy_caution':'Observed agreement against empirically refined independent quadrature references; not a rigorous interval error bound. Four-loop numerical acceptance is separate.'},
                       'native_closure':{'status':'passed_for_requested_targets_and_derivatives',
                                         'cut_basis_sizes':[{'cut_slots':[0],'native_arity':7,'basis_size':6},
                                                            {'cut_slots':[1],'native_arity':7,'basis_size':6},
                                                            {'cut_slots':[0,1],'native_arity':9,'basis_size':11}],
                                         'guard_refinement':{'max_passes':3,'max_added_domains':256,'max_interval_width':2},
                                         'search_frontier_sectors':True,
                                         'proof':'Exact native guarded original-source replay and final requested-target/derivative closure audit; no basis minimality or complete integer-domain coverage claim.',
                                         'evidence':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-frontier-sectors'},
                       'complete_sample_attempt':{'status':'passed_full_amplitude_and_all_cut_sectors',
                                                  'epsilon':'4/5','guard_digits':40,'search_frontier_sectors':True,
                                                  'configurations_digits_order_start':[[18,60,8],[28,60,8],[28,80,8],[28,80,12]],
                                                  'process_wall_seconds':153.2576642697677,'peak_rss_kib':52372,
                                                  'resources':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-frontier-sectors-resources.json',
                                                  'accepted_prediction_snapshots':4,'written_prediction_snapshots':4,'full_amplitude_target_comparisons':8,'amf_oracle_comparisons':0},
                       'controlled_basis_comparison':{'status':'passed_at_one_matching_profile','old_basis_sizes':[7,6,64],'new_basis_sizes':[6,6,11],
                                                       'digits':18,'guard_digits':40,'series_order':60,'occupied_start_scale':8,'epsilon':'4/5',
                                                       'sector_and_total_comparisons':10,'largest_relative_difference':'2.29744838567709e-46',
                                                       'report':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-basis-comparison.json',
                                                       'scope':'Saved native outputs only; not an exact master transformation or a new independent-reference generation.'},
                       'historical_complete_sample_failure':{'status':'failed_at_two_cut_physical_endpoint_with_20_guard_digits',
                                                  'epsilon':'4/5','digits':18,'guard_digits':20,'series_order':60,'occupied_start_scale':8,
                                                  'reported_error':'Numerical("uncancelled physical endpoint divergence")',
                                                  'error_interpretation':'Numerical endpoint reconstruction failed; this is not a proof of divergence of the original positive-mass integral.',
                                                  'test_harness_seconds':121.18,'process_wall_seconds':121.23372492892668,'peak_rss_kib':184372,
                                                  'resources':'reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-guard-refinement-resources.json',
                                                  'accepted_predictions':0,'written_predictions':0,'full_amplitude_comparisons':0,'amf_oracle_comparisons':0},
                       'capabilities':['raised charged line','medium-vector numerator','nonzero upper Fermi surface','full massive assembly'],
                       'evaluation_status':'full_fixed_dimension_and_laurent_amplitude_passed'},
        'remaining_mandatory_cases':[
            {'case':'support_below_at_above_threshold','reference_status':'exact analytic compact one-loop moments and distributions are derivable; implement and test'},
            {'case':'independent_chemical_potentials','reference_status':'required_not_generated','definition':'examples/finite_density/two_independent_chemical_cycles.json','chemical_parameter_points':[['1','3/2'],['5/4','7/4']],'graph_certificates_by_species':multi_certificates,'reference_method':'Independent regulated complete-contour residues and arbitrary-precision subtracted spatial integration; target 12 absolute digits, attainable accuracy not established.'},
            {'case':'singular_target_weights_endpoint','reference_status':'required_exact_integration_test','requirement':'A finite-density target change of basis exposing cancelling singular rational weights, plus divergence and insufficient-depth failures through evaluation.'},
            {'case':'massless_raised_lines','reference_status':'partially_available','requirement':'I37 gives a massless raised-neutral target only; additional raised-occupied targets require independent dimensionally regulated references.'},
        ],
        'required_reports':['dependency/source identities','reduction coverage and guards','boundary and contour provenance','saved prediction before reference comparison','independent refinement results','runtime and peak RSS with measurement scope'],
    }
    outputs['docs/finite-density-acceptance.json']=manifest
    return outputs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    for relative,value in artifacts().items():
        expected=json.dumps(value,indent=2,ensure_ascii=False)+'\n'
        path=ROOT/relative
        if args.write:
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(expected)
        elif not path.exists() or path.read_text()!=expected:
            raise SystemExit(f'Stale generated artifact: {relative}; run with --write')
    print('Validated 3 distinct nonfactorized four-loop graph certificates, one two-loop graph, the three-loop intermediate graph, two independent charge cycles, and 108 definition-only targets. No numerical predictions or comparisons were performed.')


if __name__=='__main__':
    main()
