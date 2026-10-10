//! Exploratory native source/closure experiment; no production admission change.
use std::{collections::BTreeSet, sync::Arc};
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    guarded::{GuardedContext, GuardedMeasureIdentity, IndexBounds, IndexDomain, IndexRole},
    measure::WeightedMeasure,
};
use rustred::solver::{SearchOptions, guarded::{GuardedSource, GuardedSourceSystem, GuardedProgram, GuardedApplicationFailure}};

fn main() {
    let out=std::path::PathBuf::from(std::env::args().nth(1).unwrap());
    std::fs::create_dir_all(&out).unwrap();
    let g=Atom::var(symbol!("moving_toy::g"));
    let u=Atom::var(symbol!("moving_toy::u"));
    let sa=symbol!("moving_toy::A");let sb=symbol!("moving_toy::B");let sd=symbol!("moving_toy::D");
    let a=Atom::var(sa);let b=Atom::var(sb);let d=Atom::var(sd);
    let roles=[IndexRole::RequiredCut,IndexRole::Ordinary,IndexRole::Occupation,IndexRole::Occupation];
    let indices=[symbol!("moving_toy::n"),symbol!("moving_toy::p"),symbol!("moving_toy::h"),symbol!("moving_toy::l")];
    // A=m²+eta, B=mu²-m². Eta fixes B. Physical m² varies A and B oppositely.
    let measure=WeightedMeasure::new(vec![g.clone(),u.clone()],
        [&g-&a,u.clone(),&b-&u*&u+&g,u.clone()],roles).unwrap();
    let mut ids=Vec::new();
    ids.extend(measure.ibp("energy-medium",&indices,&[Atom::num(2)*&u*&u,u.clone()],&Atom::num(1),16).unwrap());
    ids.extend(measure.multiplication_sources().unwrap());
    for (name,dir,div) in [
        ("spatial",vec![Atom::num(2)*(&g-&u*&u),Atom::new()],&d-Atom::num(1)),
        ("medium",vec![Atom::num(2)*&u,Atom::num(1)],Atom::new()),
        ("radial",vec![Atom::num(2)*&g,u.clone()],d.clone()),
    ] { ids.extend(measure.ibp(name,&indices,&dir,&div,16).unwrap()); }
    let admitted=IndexDomain::new([IndexBounds::unbounded(),IndexBounds::new(None,Some(0)).unwrap(),IndexBounds::new(Some(0),None).unwrap(),IndexBounds::new(Some(0),None).unwrap()]).unwrap();
    ids=ids.into_iter().filter_map(|mut id|{id.domain=id.domain.intersection(&admitted)?;Some(id)}).collect();
    let context=GuardedContext::new(GuardedMeasureIdentity{
        measure:"one-loop future C_n(g-A) u^(-p) H_h(B-u²+g) H_l(u); A=m²+eta; B=mu²-m²".into(),
        support:"A>0,B>0; real spatial ball; polynomial energy insertions p<=0; lower contacts vanish since u=0,g=A>0 has no real spatial point".into(),
        orientation:"future u>0, Minkowski g=u²-r²".into(),
        normalization:"unscaled d^Dq delta(g-A); d=D-1 spatial dimensions".into(),
        branch:"positive A,B then meromorphic dimension; compact radial total derivatives with all upper surfaces retained".into(),
        deformation:"eta changes only A; physical m² derivative is partial_A-partial_B; independent parameters during source generation".into(),
    },roles,indices,vec![sa,sb,sd],ids).unwrap();
    // Lossless public reconstruction permits the explicit positive-A support
    // zero to join the same source context; no rule or evaluated period added.
    let original=context.sources();
    let rows=original.sources().iter().zip(original.native_sources().rows()).map(|(info,row)|GuardedSource::new(info.id.clone(),row.clone(),info.domain.clone()).with_nonzero_conditions(info.nonzero_conditions.clone())).collect();
    let lower_zero=IndexDomain::new([IndexBounds::new(Some(1),None).unwrap(),IndexBounds::new(None,Some(0)).unwrap(),IndexBounds::new(Some(0),None).unwrap(),IndexBounds::new(Some(1),None).unwrap()]).unwrap();
    let sources=Arc::new(GuardedSourceSystem::new(original.measure_id(),roles,*original.native_sources().index_variables(),rows).unwrap().with_zero_domains(vec![lower_zero]).unwrap());
    let root=[1,0,0,0];let mut frontier=BTreeSet::from([root]);
    let mut retained:Option<GuardedProgram<4>>=None;let mut history=Vec::new();let mut closed=None;
    for round in 0..16 {
        let requested=frontier.iter().copied().chain(frontier.iter().map(|i|{let mut j=*i;j[0]+=1;j})).collect::<BTreeSet<_>>().into_iter().collect::<Vec<_>>();
        let domains=requested.iter().map(|i|IndexDomain::new(i.map(IndexBounds::fixed)).unwrap()).collect();
        let discovered=sources.solve_domains_with_priority_points(domains,&requested,SearchOptions{max_depth:Some(3),sample_seed:0,..Default::default()},2048).unwrap();
        let gaps=discovered.unresolved.iter().map(|g|format!("{:?}: {}",g.reason,g.detail)).collect::<Vec<_>>();
        let fresh=GuardedProgram::new(sources.clone(),discovered.rules,[]).unwrap();
        let program=if let Some(previous)=retained{previous.union_replayed(fresh,[],16384).unwrap()}else{fresh};
        let before=frontier.clone();
        // Native frontier recurrences discovered in this round retire old labels
        // before new derivative obligations are formed. No manual terminals.
        let mut reduced_frontier=BTreeSet::new();
        for i in &before {
            for term in program.reduce(*i,Default::default()).unwrap().unresolved {
                assert_eq!(term.reason,GuardedApplicationFailure::NoApplicableRule);
                if !term.coefficient.is_zero(){reduced_frontier.insert(term.integral);}
            }
        }
        frontier.clear();
        let active=std::iter::once(root).chain(reduced_frontier.iter().map(|i|{let mut j=*i;j[0]+=1;j})).collect::<Vec<_>>();
        for i in &active {
            let reduced=program.reduce(*i,Default::default()).unwrap();
            for term in reduced.unresolved {
                assert_eq!(term.reason,GuardedApplicationFailure::NoApplicableRule);
                if !term.coefficient.is_zero(){assert!(admitted.contains(&term.integral));frontier.insert(term.integral);}
            }
        }
        history.push(serde_json::json!({"round":round,"frontier":frontier,"rules":program.rules().len(),"discovery_gaps":gaps}));
        std::fs::write(out.join("progress.json"),serde_json::to_vec_pretty(&history).unwrap()).unwrap();
        if frontier==before {
            let audited=program.with_terminals_replayed(frontier.iter().copied(),16384).unwrap();
            let mut matrix=Vec::new();
            for i in &frontier {let mut derivative=*i;derivative[0]+=1;let reduced=audited.reduce(derivative,Default::default()).unwrap();assert!(reduced.unresolved.is_empty());matrix.push(serde_json::json!({"basis":i,"cut_derivative_prefactor":i[0],"terms":reduced.terms.iter().map(|(j,c)|serde_json::json!({"indices":j,"coefficient":c.to_expression().to_string()})).collect::<Vec<_>>(),"conditions":reduced.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>() }));}
            assert!(audited.reduce(root,Default::default()).unwrap().unresolved.is_empty());
            std::fs::write(out.join("closed.bin"),audited.encode_native(Default::default()).unwrap()).unwrap();
            closed=Some(serde_json::json!({"basis":frontier,"derivatives":matrix,"rule_count":audited.rules().len()}));break;
        }
        retained=Some(program);
        if frontier.len()>128{break;}
    }
    let result=serde_json::json!({"scope":"exploratory moving-cut native weighted-source closure; automatic provisional frontier searched alongside its derivatives, exact final derivative audit; no production FixedShellDeformation change, no manual rule or period","source_count":sources.sources().len(),"measure_id":sources.measure_id(),"history":history,"closed":closed,"physical_mass_derivative":"partial_m²=partial_A-partial_B; for bulk I, partial_B I=S, so partial_m² I=partial_eta I-S","fixed_shell_comparison":"For this one-loop occupied tadpole there is no uncut denominator, so existing fixed-shell deformation is eta independent. The moving cut introduces a nontrivial compact integral flow; this toy cannot establish an advantage on multiloop virtual reduction."});
    std::fs::write(out.join("result.json"),serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    println!("closed={} source_count={} history_rounds={}",result["closed"].is_object(),sources.sources().len(),result["history"].as_array().unwrap().len());
}
