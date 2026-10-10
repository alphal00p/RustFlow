#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{
    BinaryIoLimits, CoefficientId, DecodedCoefficientTable, SectionTag, inspect_program,
};
use rustred::solver::guarded::{
    GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};
use rustred::solver::{Integral, Power, SearchOptions, Term};
use std::sync::Arc;
type Domain = Vec<(Option<i64>, Option<i64>)>;
type Label = Vec<(bool, i16)>;
#[derive(Decode)]
struct Tr {
    integral: Label,
    coefficient: usize,
}
#[derive(Decode)]
struct Sr {
    id: String,
    domain: Domain,
    conditions: Vec<usize>,
    terms: Vec<Tr>,
}
#[derive(Decode)]
struct Sd {
    row: usize,
    integral: Label,
    shifts: Vec<i16>,
}
#[derive(Decode)]
struct Rr {
    fixed: Vec<Option<i16>>,
    target: Label,
    rhs: Vec<Tr>,
    sources: Vec<Sd>,
    domain: Domain,
    discovery_domain: Domain,
    conditions: Vec<usize>,
    sector: Vec<bool>,
    permutation: Option<Vec<usize>>,
}
#[derive(Decode)]
struct Rec {
    schema: String,
    measure: String,
    roles: Vec<u8>,
    indices: Vec<usize>,
    sources: Vec<Sr>,
    zero_domains: Vec<Domain>,
    rules: Vec<Rr>,
    terminals: Vec<Vec<i64>>,
}
fn domain<const N: usize>(v: &Domain) -> IndexDomain<N> {
    IndexDomain::new(
        v.iter()
            .map(|&(lo, hi)| IndexBounds::new(lo, hi).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
    .unwrap()
}
fn label<const N: usize>(v: &Label) -> Integral<N> {
    Integral::new(
        v.iter()
            .map(|&(s, p)| Power::new(s, p).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
}
fn bounds<const N: usize>(v: &IndexDomain<N>) -> Domain {
    v.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}

use std::collections::BTreeMap;
use symbolica::prelude::*;
fn experiment_stats<const N: usize>(
    rule: &rustred::solver::guarded::GuardedRule<N>,
) -> serde_json::Value {
    #[cfg(source_portfolio)]
    {
        rule.candidate().stats.guarded_portfolio.map(|s| serde_json::json!({
        "baseline_rhs_terms":s.baseline_rhs_terms,"baseline_rows":s.baseline_rows,
        "baseline_seeds":s.baseline_seeds,"baseline_independent_rows":s.baseline_independent_rows,
        "baseline_exact_trace_rows":s.baseline_exact_trace_rows,"baseline_direct_hit":s.baseline_direct_hit,
        "baseline_elapsed":s.baseline_elapsed.as_secs_f64(),"baseline_exact_materialization":s.baseline_exact_materialization.as_secs_f64(),
        "selected_rhs_terms":s.selected_rhs_terms,"selected_arm":s.selected_arm,
        "preseal_calls":s.preseal_calls,"preseal_elapsed":s.preseal_elapsed.as_secs_f64(),
        "policy_elapsed":s.policy_elapsed.as_secs_f64(),"total_elapsed":s.total_elapsed.as_secs_f64(),"completion":s.completion,
        "trials":s.trials.iter().map(|t|serde_json::json!({
            "selector":t.selector,"selected_source_rows":t.selected_source_rows,
            "attempted_rows":t.attempted_rows,"accepted_rows":t.accepted_rows,"seeds":t.seeds,
            "independent_rows":t.independent_rows,"exact_materialization":t.exact_materialization.as_secs_f64(),
            "guard_rejected_rows":t.guard_rejected_rows,"empty_rows":t.empty_rows,
            "exact_trace_rows":t.exact_trace_rows,"exact_trace_terms":t.exact_trace_terms,
            "rhs_terms":t.rhs_terms,"nonzero_conditions":t.nonzero_conditions,
            "completion":t.completion,"search_elapsed":t.search_elapsed.as_secs_f64(),
        })).collect::<Vec<_>>()
    })).unwrap_or(serde_json::Value::Null)
    }
    #[cfg(not(source_portfolio))]
    {
        let _ = rule;
        serde_json::Value::Null
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let bytes = std::fs::read(&args[1]).unwrap();
    let limits = BinaryIoLimits::default();
    let envelope = inspect_program(&bytes, limits).unwrap();
    let (record, used): (Rec, usize) = bincode::decode_from_slice(
        envelope.section(SectionTag::PROGRAM).unwrap(),
        bincode::config::standard(),
    )
    .unwrap();
    assert_eq!(used, envelope.section(SectionTag::PROGRAM).unwrap().len());
    assert!(
        matches!(
            record.schema.as_str(),
            "rustred.guarded-source-program.v1" | "rustred.guarded-source-program.v2"
        ),
        "only explicitly known native source schemas may be imported"
    );
    let table = DecodedCoefficientTable::import_generated_normalized(
        envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),
        envelope.section(SectionTag::COEFFICIENTS).unwrap(),
        limits,
    )
    .unwrap();
    let polynomial = |id| {
        let c = table
            .coefficient(CoefficientId::try_from_index(id).unwrap())
            .unwrap()
            .clone();
        assert!(c.denominator.is_one());
        c.numerator
    };
    let roles: [IndexRole; 16] = record
        .roles
        .iter()
        .map(|r| match r {
            0 => IndexRole::Ordinary,
            1 => IndexRole::RequiredCut,
            2 => IndexRole::Occupation,
            _ => panic!(),
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let indices: [usize; 16] = record.indices.clone().try_into().unwrap();
    let mode = &args[4];
    assert!(["original", "extended-energy-Ward"].contains(&mode.as_str()));
    let order: Vec<usize> = (0..record.sources.len()).collect();
    let rows: Vec<GuardedSource<16>> = order
        .iter()
        .map(|&i| {
            let s = &record.sources[i];
            GuardedSource::new(
                s.id.clone(),
                s.terms
                    .iter()
                    .map(|t| Term {
                        integral: label(&t.integral),
                        coefficient: polynomial(t.coefficient),
                    })
                    .collect(),
                domain(&s.domain),
            )
            .with_nonzero_conditions(s.conditions.iter().map(|&id| polynomial(id)).collect())
        })
        .collect();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            record.measure.clone(),
            roles,
            indices,
            rows,
        )
        .unwrap()
        .with_zero_domains(record.zero_domains.iter().map(domain).collect())
        .unwrap(),
    );
    let depth:u32=args[5].parse().unwrap();assert!([3,4,5].contains(&depth));
    let selected:serde_json::Value=serde_json::from_slice(&std::fs::read(&args[2]).unwrap()).unwrap();let selected=selected.as_array().unwrap();assert!((1..=16).contains(&selected.len()));
    let physical=if roles[12]==IndexRole::Occupation{13}else{11};
    let mut bounds0=*IndexDomain::for_roles(&roles).bounds();for b in &mut bounds0[5..9]{*b=IndexBounds::new(None,Some(0)).unwrap()}for b in &mut bounds0[physical..]{*b=IndexBounds::fixed(0)}let admitted=IndexDomain::new(bounds0).unwrap();
    use rustred::solver::guarded::{GuardedApplicationStatus as Status,GuardedApplicationFailure as Failure};
    let historical=GuardedProgram::decode_generated(&bytes,sources.clone(),limits).unwrap();assert!(historical.terminals().is_empty());
    let mut observations=Vec::new();
    for point in selected{
        let target:[i64;16]=serde_json::from_value(point["point"].clone()).unwrap();let applied=historical.apply(&target).unwrap();
        if point["kind"]=="actual-frontier"{assert_eq!(applied.status,Status::Unresolved(Failure::NoApplicableRule))}else{assert!(matches!(applied.status,Status::Applied{..}))}
        observations.push(serde_json::json!({"point":target,"kind":point["kind"],"status":format!("{:?}",applied.status)}));
    }
    drop(historical);
    assert_eq!(sources.sources().len(),62);
    use symbolica_amflow::finite_density::{DensityInput,guarded::{GuardedContext,GuardedMeasureIdentity,GuardedIdentity,GuardedIdentityTerm}};
    let variables=sources.native_sources().rows()[0][0].coefficient.variables().clone();
    let symbols=variables.iter().map(|v|{let PolyVariable::Symbol(s)=v else{panic!("symbolic variable required")};*s}).collect::<Vec<_>>();
    let index_symbols:[Symbol;16]=record.indices.iter().map(|&i|symbols[i]).collect::<Vec<_>>().try_into().unwrap();
    let parameters=symbols.iter().enumerate().filter(|(i,_)|!record.indices.contains(i)).map(|(_,s)|*s).collect::<Vec<_>>();
    let eta=symbol!("rustflow_occupied::eta");let eps=symbol!("rustflow_occupied::epsilon");assert!(parameters.contains(&eta)&&parameters.contains(&eps));
    let input:DensityInput=serde_json::from_slice(&std::fs::read(&args[6]).unwrap()).unwrap();let prepared=input.prepare().unwrap();assert!(record.measure.contains(prepared.identity()));
    let family=prepared.occupied_cut(&[0],4096).unwrap().at_physical_masses();assert_eq!(family.loops(),3);assert_eq!(family.shells().len(),1);
    let shell=&family.shells()[0];assert_eq!(shell.loop_index,0);assert_eq!(shell.physical_slot,0);assert_eq!(shell.upper_slot,9);assert_eq!(shell.lower_slot,10);assert!(shell.mass_squared.is_zero());assert_eq!(shell.chemical_potential,Rational::one());
    let measure=family.deformed_measure::<16>(eta,&[1,2,3,4]).unwrap();assert_eq!(measure.roles(),&roles);assert_eq!(family.physical_slots(),5);assert_eq!(family.input_slots(),9);
    assert!(record.measure.contains("singleton-UV-meromorphic-null-germ"));
    let coordinates=family.coordinates();let energy=coordinates[6].clone();let medium=&coordinates[6..9];
    let mut pure=Vec::new();let mut medium_dependent=Vec::new();
    for slot in family.physical_slots()..family.input_slots(){
        let factor=&measure.factors()[slot];let depends=medium.iter().any(|c|{let AtomView::Var(v)=c.as_view()else{panic!()};!factor.derivative(v.get_symbol()).expand().is_zero()});
        if depends{medium_dependent.push(slot);let ratio=(factor/&energy).together().cancel();if let Ok(c)=Rational::try_from(ratio.as_view()){if !c.is_zero(){assert!((factor-Atom::num(c.clone())*&energy).expand().is_zero());pure.push((slot,c));}}}
    }
    assert_eq!(pure,vec![(6,Rational::one())]);assert_eq!(medium_dependent,vec![6,7,8]);
    let mut new_rows=Vec::new();
    for positive in [false,true]{
        let mut b=*admitted.bounds();b[shell.physical_slot]=IndexBounds::fixed(1);b[shell.lower_slot]=IndexBounds::fixed(0);b[shell.upper_slot]=if positive{IndexBounds::new(Some(1),None).unwrap()}else{IndexBounds::fixed(0)};
        for &slot in &medium_dependent{if !pure.iter().any(|(s,_)|*s==slot){b[slot]=IndexBounds::fixed(0)}}
        let diag=pure.iter().fold(Atom::num(2)-Atom::num(2)*Atom::var(eps),|a,(slot,_)|a-Atom::var(index_symbols[*slot]));
        let mut shift=[0i16;16];shift[shell.upper_slot]=1;shift[pure[0].0]=-1;
        let coefficient=if positive{Atom::var(index_symbols[shell.upper_slot])}else{-Atom::one()};
        new_rows.push(GuardedIdentity{id:format!("extended-polynomial-singleton-energy-Ward-v1/upper-{}",if positive{"positive"}else{"zero"}),domain:IndexDomain::new(b).unwrap(),terms:vec![GuardedIdentityTerm{shift:[0;16],coefficient:diag},GuardedIdentityTerm{shift,coefficient}],nonzero_conditions:vec![]});
    }
    // The new rows coincide with the original energy-free Ward pair at a6=0;
    // only a guarded C1 polynomial-energy degree is added. No Cn>1 assertion.
    for row in &new_rows{
        let suffix=if row.id.ends_with("positive"){ "positive" }else{"zero"};
        let old=sources.sources().iter().zip(sources.native_sources().rows()).find(|(m,_)|m.id==format!("raw-polynomial-singleton-Ward-v1/upper-{suffix}")).unwrap();
        let replacement=BTreeMap::from([(Atom::var(index_symbols[6]),Atom::zero())]);
        let new=BTreeMap::from_iter(row.terms.iter().map(|t|(t.shift,symbolica_amflow::family::substitute(&t.coefficient,&replacement).expand().together().cancel())));
        let old_terms=BTreeMap::from_iter(old.1.iter().map(|t|(t.integral.powers().map(|p|p.value()),t.coefficient.to_expression().expand().together().cancel())));
        assert_eq!(new,old_terms);assert!(old.0.nonzero_conditions.is_empty());
        let mut restricted=*row.domain.bounds();restricted[6]=IndexBounds::fixed(0);assert_eq!(IndexDomain::new(restricted).unwrap(),old.0.domain);
    }
    let audit=new_rows.iter().map(|r|serde_json::json!({"id":r.id,"domain":bounds(&r.domain),"conditions":r.nonzero_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),"terms":r.terms.iter().map(|t|serde_json::json!({"shift":t.shift,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
    let geometry=serde_json::json!({"input_identity":prepared.identity(),"coordinates":coordinates.iter().map(ToString::to_string).collect::<Vec<_>>(),"factors":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"routing":family.inverse_routing().iter().map(|r|r.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),"occupied_energy":energy.to_string(),"pure_energy_slots":pure.iter().map(|(s,c)|serde_json::json!({"slot":s,"coefficient":c.to_string()})).collect::<Vec<_>>(),"medium_dependent_slots":medium_dependent,"shell_mass_squared":"0","chemical_potential":"1","retained_origin":"sealed singleton-UV-meromorphic-null-germ, pointwise finite polynomial index high-D cutoff removal"});
    let mut actual_rows=sources.sources().iter().zip(sources.native_sources().rows()).map(|(m,r)|GuardedSource::new(m.id.clone(),r.clone(),m.domain.clone()).with_nonzero_conditions(m.nonzero_conditions.clone())).collect::<Vec<_>>();
    if mode=="extended-energy-Ward"{
        let ctx=GuardedContext::new_with_physical_arity(GuardedMeasureIdentity{measure:format!("{};extended-polynomial-singleton-energy-Ward-v1",record.measure),support:"existing C1 massless singleton high-D origin; polynomial energy degree only".into(),orientation:"same future occupied loop0, mu1; no mixed or virtual energy insertions".into(),normalization:"homogeneous native sources, no period values".into(),branch:"same regulator order; no inverse-energy or zero-domain widening".into(),deformation:"same uncut1,2,3,4 D-eta".into()},roles,index_symbols,parameters,new_rows,11).unwrap();
        assert_eq!(ctx.sources().native_sources().coefficient_variables(),variables.as_ref().as_slice());
        actual_rows.extend(ctx.sources().sources().iter().zip(ctx.sources().native_sources().rows()).map(|(m,r)|GuardedSource::new(m.id.clone(),r.clone(),m.domain.clone()).with_nonzero_conditions(m.nonzero_conditions.clone())));
    }
    let sources=Arc::new(GuardedSourceSystem::new(format!("{};source-control={mode};extended-polynomial-singleton-energy-Ward-v1",record.measure),roles,record.indices.clone().try_into().unwrap(),actual_rows).unwrap().with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());

    let out=std::path::PathBuf::from(&args[3]);std::fs::create_dir_all(&out).unwrap();
    let describe=|applied:&rustred::solver::guarded::GuardedApplication<16>|serde_json::json!({"status":format!("{:?}",applied.status),"rhs":applied.terms.iter().map(|(i,c)|serde_json::json!({"indices":i.to_vec(),"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":applied.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>()});
    let fallback=|s:&Status|matches!(s,Status::Unresolved(Failure::NoApplicableRule|Failure::ConditionVanished{..}));
    let mut probes=Vec::new();
    for(ordinal,selection)in selected.iter().enumerate(){
        let target:[i64;16]=serde_json::from_value(selection["point"].clone()).unwrap();assert!(admitted.contains(&target));
        let search=|domain|sources.solve_domains(vec![domain],SearchOptions{max_depth:Some(depth),sample_seed:0,..Default::default()},1).unwrap();
        let ray=IndexDomain::new(std::array::from_fn(|i|if target[i]==0{IndexBounds::fixed(0)}else if target[i]>0{IndexBounds::new(Some(target[i]),None).unwrap()}else{IndexBounds::new(None,Some(target[i])).unwrap()})).unwrap().intersection(&admitted).unwrap();
        let ray=search(ray);let mut gaps=ray.unresolved.iter().map(|g|serde_json::json!({"phase":"ray","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
        let mut program=GuardedProgram::new(sources.clone(),ray.rules,[]).unwrap();let ray_result=program.apply(&target).unwrap();let needs_point=fallback(&ray_result.status);assert!(needs_point||matches!(ray_result.status,Status::Applied{..}|Status::Zero));
        let mut point_result=serde_json::Value::Null;
        if needs_point{
            let exact=search(IndexDomain::new(target.map(IndexBounds::fixed)).unwrap());gaps.extend(exact.unresolved.iter().map(|g|serde_json::json!({"phase":"point","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})));
            let exact=GuardedProgram::new(sources.clone(),exact.rules,[]).unwrap();point_result=describe(&exact.apply(&target).unwrap());program=program.union_replayed(exact,[],65536).unwrap();
        }
        let encoded=program.encode_native(limits).unwrap();std::fs::write(out.join(format!("point-{ordinal:02}.bin")),&encoded).unwrap();let replay=GuardedProgram::decode_generated(&encoded,sources.clone(),limits).unwrap();let result=program.apply(&target).unwrap();let check=replay.apply(&target).unwrap();assert_eq!(result.status,check.status);assert_eq!(result.terms,check.terms);assert_eq!(result.nonzero_conditions,check.nonzero_conditions);assert!(result.terms.keys().all(|p|admitted.contains(p)));
        probes.push(serde_json::json!({"selection":selection,"ray":describe(&ray_result),"point_fallback":needs_point,"point":point_result,"application":describe(&result),"roundtrip":true,"rules":program.rules().len(),"gaps":gaps,"proof_sources":program.rules().iter().map(|r|r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>()})).collect::<Vec<_>>()).collect::<Vec<_>>() }));
    }
    let report=serde_json::json!({"scope":"Bounded source control on exact historical NoApplicableRule leaves, 12 eligible pure-energy C1 rows and four out-of-domain controls. Historical program is replayed ONLY to verify selection status, then discarded before fresh discovery; no old rule union. Source definitions unchanged; no periods/closure inference.","mode":mode,"geometry":geometry,"added_sources":if mode=="original"{vec![]}else{audit},"same_zero_domains":true,"energy_free_old_row_equivalence":true,"old_rules_imported":false,"depth":depth,"ray_domains":1,"point_domains":1,"source_count":sources.sources().len(),"source_identity":record.measure,"source_schema":record.schema,"selection_observations":observations,"all_rhs_labels_admitted":true,"probes":probes});
    std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&report).unwrap()).unwrap();println!("depth={depth} points={} source={} selection and fresh replay passed",selected.len(),record.sources.len());
}
