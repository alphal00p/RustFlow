//! Exercise an existing proof-backed ordinary vacuum reducer above pow8 limits.
use symbolica::prelude::*;
use symbolica_amflow::*;
use symbolica_amflow::reduction::ReductionBackend;
use rustred::family::IntegralKey;
use rustred::foundry::artifact::derive_one_loop_unit_mass_tadpole;
use rustred::reduction::{Reducer,ReductionLimits};

fn main() -> Result<()> {
    let artifact=derive_one_loop_unit_mass_tadpole().map_err(|e|Error::Unsupported(e.to_string()))?;
    let mut reducer=Reducer::with_limits(&artifact,ReductionLimits {
        max_rule_applications:1024,max_cached_integrals:200000,
        max_cached_coefficient_terms:200000,max_cached_coefficient_bytes:128*1024*1024,
        max_pending_frames:1024,..Default::default()
    }).map_err(|e|Error::Unsupported(e.to_string()))?;
    let dimension=artifact.coefficient_context().parameter("d").unwrap().to_expression();
    let AtomView::Var(dsymbol)=dimension.as_view() else{panic!("plain dimension")};
    let epsilon=Rational::from((-5,4));
    let family=IntegralFamily{name:"ordinary_wide_vacuum_check".into(),loops:vec!["k".into()],external:vec![],external_gram:vec![],propagators:vec![Propagator{constant:Atom::num(-1),scalar_products:vec![Atom::one()]}],physical_propagators:1,epsilon:symbol!("ordinary_wide_check::epsilon"),dimension:4};
    let backend=RustRedBackend{bubble_subloops:false,..Default::default()};
    let mut records=vec![];
    for n in 1..=129 {
        let target=IntegralKey::try_new([i64::from(n)]).unwrap();
        let reduced=reducer.reduce_unit_mass(&target).map_err(|e|Error::Unsupported(e.to_string()))?;
        assert_eq!(reduced.terms().len(),1);
        let(master,c)=reduced.terms().first_key_value().unwrap();assert_eq!(master.powers(),[1]);
        let coefficient=c.to_expression().replace(Atom::var(dsymbol.get_symbol())).with(Atom::num((13,2))).together().cancel();
        let rational=Rational::try_from(coefficient.as_view()).unwrap();
        let mut record=serde_json::json!({"power":n,"master":[1],"rational_coefficient":rational.to_string()});
        if n<=32 {
            let key=Integral(vec![n as i16]);
            let old=backend.reduce_at_epsilon(&family,std::slice::from_ref(&key),&epsilon,&RunContext::default())?;
            let row=old.rules.get(&key).cloned().unwrap_or_else(||std::collections::BTreeMap::from([(key,Atom::one())]));
            assert_eq!(row.len(),1);assert_eq!(row.first_key_value().unwrap().0.0,vec![1]);
            assert!((row.first_key_value().unwrap().1-&coefficient).together().cancel().is_zero());
            record["compact_backend_exact_match"]=true.into();
        }
        if [64,80,128,129].contains(&n) {
            let mut checks=vec![];
            for digits in [50,80] {
                let p=Precision::decimal(digits)?;
                let base=gaussian::terminal(&family,&Integral(vec![1]),&epsilon,p)?;
                let direct=gaussian::terminal(&family,&Integral(vec![n as i16]),&epsilon,p)?;
                let prediction=p.mul(&p.rational(&rational),&base);
                let error=p.norm(&p.div(&p.sub(&prediction,&direct),&direct));
                let passed=error<p.parse("1e-35","0")?.re;assert!(passed);
                checks.push(serde_json::json!({"digits":digits,"relative_error":error.to_string(),"passed":passed}));
            }
            record["native_gaussian_checks"]=checks.into();
        }
        records.push(record);
    }
    let stats=reducer.statistics();
    let result=serde_json::json!({"scope":"Existing generic sealed ordinary-vacuum consumer; no power representation or finite-density production change","algorithm":artifact.algorithm_id(),"family_fingerprint":artifact.family_fingerprint(),"source_rows":artifact.validation().source_rows(),"source_rows_replayed":artifact.validation().replayed_source_rows(),"guarded_rules":artifact.validation().guarded_rules(),"powers_checked":[1,129],"exact_compact_backend_comparisons":32,"independent_gaussian_checks":8,"cache_terms":stats.cached_coefficient_terms(),"cache_bytes":stats.cached_coefficient_bytes(),"rule_applications":stats.rule_applications(),"records":records,"passed":true});
    std::fs::write(std::env::args().nth(1).unwrap(),serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    println!("wide powers1..129 PASS; compact matches32; Gaussian checks8");Ok(())
}
