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
        native={
            'name':name,'loops':4,'vertices':graph['vertices'],
            'edges':[{'vertices':ends,'routing':[str(x) for x in ROUTINGS[j]],
                      'mass_squared':'0','charges':[graph['edge_charges'][j]]}
                     for j,ends in enumerate(edges)],
            'loop_charges':[[x] for x in signatures],'chemical_potentials':['1'],
            'numerator_convention':'shifted_euclidean','laurent_orders':[-4,0],'digits':12,
            'targets':[{'powers':target['denominator_powers'][:e],
                        'numerator':target['numerator_polynomial']},
                       {'powers':[2]+[1]*(e-1),'numerator':'g1_2^2+g1_3*g2_4'}],
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
                 'reference':{'status':'required_not_generated',
                              'method':'Independent complete thermal-contour/cutting evaluation with independent line masses and distributional mass differentiation; independently implemented sector resolution and arbitrary-precision spatial quadrature after subtraction.',
                              'target_absolute_uncertainty':'1/10^12',
                              'attainable_precision':'Not established. Reference generation and error certification are blocking acceptance gates; no 12-digit availability is claimed.'}},
            ],
            'evaluation_status':'not_evaluated','coefficients_compared':[],
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
        'four_loop_families':entries,
        'structural_distinction':{'verified':True,'evidence':'Explicit connected, bridge-free, articulation-free incidence certificates with (V,E)=(4,7),(5,8),(6,9), different degree multisets and canonical multigraph adjacency words. Exhaustive exact row-rank partitions rule out denominator-block factorization. Physical edge count is preserved under graph relabeling and invertible loop rerouting; completion slots are excluded.'},
        'small_graph':{'definition':'examples/finite_density/massive_two_loop_sunset.json',
                       'graph_certificate':small,'required_cut_components':[0,1,2],
                       'vacuum_contribution':'Nonzero for positive masses; must evaluate and compare, not omit as scaleless.',
                       'reference_status':'required_not_generated',
                       'reference_method':'Independent massive Schwinger-parameter vacuum integral plus regulated complete energy-contour residues and dimensionally subtracted spatial integration; masses varied independently before raising line0.',
                       'attainable_comparison_precision':'Not established; target 12 absolute digits after regulator subtraction.',
                       'capabilities':['raised charged line','medium-vector numerator','nonzero upper Fermi surface','full massive assembly'],
                       'evaluation_status':'not_evaluated'},
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
    print('Validated 3 distinct nonfactorized four-loop graph certificates, one two-loop graph, two independent charge cycles, and 108 definition-only targets. No numerical predictions or comparisons were performed.')


if __name__=='__main__':
    main()
