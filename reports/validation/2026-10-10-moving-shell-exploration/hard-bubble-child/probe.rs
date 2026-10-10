//! Exploratory ordinary AMF child; no endpoint/reference data are read.
use serde_json::json;
use std::{path::PathBuf,sync::Arc,time::Instant};
use symbolica_amflow::{symbolica,*,recursive::{RecursiveBoundary,RecursiveTerminalPolicy}};
use symbolica::prelude::*;
fn main() {
    let output=PathBuf::from(std::env::args().nth(1).expect("OUTPUT_DIRECTORY"));
    std::fs::create_dir_all(&output).unwrap();
    let family=IntegralFamily {
        name:"moving_shell_leading_hard_bubble".into(),loops:vec!["k".into()],external:vec!["p".into()],
        external_gram:vec![vec![Atom::num(1)]],
        propagators:vec![
            Propagator::quadratic(&[1],&[0],Atom::num(1),&[vec![Atom::num(1)]]).unwrap(),
            Propagator::quadratic(&[1],&[-1],Atom::num(1),&[vec![Atom::num(1)]]).unwrap(),
        ],physical_propagators:2,epsilon:symbol!("hard_bubble_eps"),dimension:4,
    };
    let options=FlowOptions {digits:20,guard_digits:40,series_order:60,mass_mode:MassMode::All,sampled_reduction:false,..Default::default()};
    let backend=RustRedBackend {bubble_subloops:false,max_depth:3,checkpoints:Some(output.join("native-reduction")),..Default::default()};
    let context=RunContext {progress:Some(Arc::new(|event|println!("{event:?}"))),..Default::default()};
    let start=Instant::now();
    let prepared=match PreparedFlow::new(&family,&[Integral(vec![1,1])],&KinematicPoint::default(),&backend,&options,&context) {
        Ok(prepared)=>prepared,
        Err(error)=>{
            std::fs::write(output.join("failure.json"),serde_json::to_string_pretty(&json!({"stage":"ordinary symbolic preparation","error":error.to_string(),"reference_values_read":false})).unwrap()).unwrap();
            panic!("ordinary child preparation failed: {error}");
        }
    };
    let preparation_seconds=start.elapsed().as_secs_f64();
    let system=json!({"family":"denominators k^2-1 and (k-p)^2-1; p^2=1","dimension":"4-2 epsilon","epsilon":"4/5","target":[1,1],"basis":prepared.reduced.basis.iter().map(|x|x.0.clone()).collect::<Vec<_>>(),
        "matrix":prepared.system.matrix.iter().map(|r|r.iter().map(|a|a.to_string()).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "nonzero_conditions":prepared.reduced.nonzero_conditions.iter().map(|a|a.to_string()).collect::<Vec<_>>(),
        "targets_exact":format!("{:?}",prepared.reduced.targets),"bubble_subloops":false,"terminal_policy":"TadpolesOnly","preparation_seconds":preparation_seconds});
    std::fs::write(output.join("system.json"),serde_json::to_string_pretty(&system).unwrap()+"\n").unwrap();
    let mut records=Vec::new();
    for (digits,order) in [(20,60),(30,60),(30,80)] {
        let profile=FlowOptions {digits,series_order:order,..options.clone()};
        let provider=RecursiveBoundary::new(&backend,&profile,&context).with_terminal_policy(RecursiveTerminalPolicy::TadpolesOnly);
        let now=Instant::now();
        let values=match prepared.evaluate(&Rational::from((4,5)),&profile,&provider,&context) {
            Ok(values)=>values,
            Err(error)=>{
                std::fs::write(output.join("failure.json"),serde_json::to_string_pretty(&json!({"stage":"ordinary AMF endpoint","digits":digits,"order":order,"error":error.to_string(),"reference_values_read":false})).unwrap()).unwrap();
                panic!("ordinary child evaluation failed: {error}");
            }
        };
        let row=json!({"digits":digits,"guard_digits":40,"series_order":order,"wall_seconds":now.elapsed().as_secs_f64(),"values":values.iter().map(|x|json!({"re":x.re.to_string(),"im":x.im.to_string()})).collect::<Vec<_>>()});
        std::fs::write(output.join(format!("prediction-{digits}-{order}.json")),serde_json::to_string_pretty(&row).unwrap()+"\n").unwrap();
        records.push(row);
    }
    let report=json!({"schema":1,"status":"predictions saved before comparison","scope":"ordinary one-loop AMF evaluation of the leading hard child only; not a moving-shell full amplitude or generic boundary-owner implementation",
        "normalization":"d^D k/(i pi^(D/2)); two Minkowski denominators k^2-1+i0 and (k-p)^2-1+i0; p^2=1", "D":"12/5","epsilon":"4/5","reference_values_read":false,"bubble_subloops":false,"terminal_policy":"TadpolesOnly","preparation_seconds":preparation_seconds,"records":records});
    std::fs::write(output.join("predictions.json"),serde_json::to_string_pretty(&report).unwrap()+"\n").unwrap();
    println!("saved three ordinary AMF hard-child profiles; no reference values read");
}
