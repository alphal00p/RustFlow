//! Validation-only transport of the independently derived radial two-state ODE.
//! Native guarded-source closure is audited separately; no reference values read.
use serde_json::json;
use symbolica_amflow::{
    BoundaryData, ComplexFloat as C, DifferentialSystem, FlowOptions, Precision, RunContext,
    symbolica,
};
use symbolica::prelude::*;

fn atom(s: &str) -> Atom {
    Atom::parse(s, "moving_shell_radial", Default::default()).unwrap()
}

fn scalar(p: Precision, expression: &str) -> C {
    p.eval(&atom(expression), &Default::default()).unwrap()
}

fn boundary(p: Precision, d_num: i64, d_den: i64, start: i64, order: usize) -> Vec<C> {
    let a = scalar(p, &format!("{start}+1/4"));
    let b = scalar(p, "3/4");
    let ratio = p.div(&b, &a);
    let pref_i = scalar(p, &format!("(3/4)^({d_num}/(2*{d_den}))/(2*({start}+1/4)^(1/2))"));
    let pref_s = scalar(p, &format!("(3/4)^(({d_num}/{d_den}-2)/2)/(4*({start}+1/4)^(1/2))"));
    let mut term = p.i(1);
    let mut sum_i = p.zero();
    let mut sum_s = p.zero();
    for n in 0..order {
        sum_i = p.add(&sum_i, &p.scale(&term, d_den, d_num + 2*n as i64*d_den));
        sum_s = p.add(&sum_s, &term);
        term = p.mul(&term, &ratio);
        term = p.scale(&term, -(2*n as i64+1), 2*(n as i64+1));
    }
    vec![p.mul(&pref_i, &sum_i), p.mul(&pref_s, &sum_s)]
}

fn main() {
    let output = std::env::args().nth(1).expect("OUTPUT");
    let eta = symbol!("moving_shell_radial::eta");
    let mut results = Vec::new();
    for (d_num, d_den) in [(7,5), (3,1), (9,2)] {
        let system = DifferentialSystem {
            variable: eta,
            matrix: vec![
                vec![atom(&format!("({d_num}/{d_den}-1)/(2*(eta+1/4))")), atom("-(3/4)/(eta+1/4)")],
                vec![Atom::new(), atom("-1/(2*(eta+1))")],
            ],
        };
        for (digits, taylor_order, start, boundary_order) in [
            (60,48,8,32), (60,64,8,48), (80,64,8,48), (80,64,16,48), (80,64,16,64),
        ] {
            let p = Precision::decimal(digits).unwrap();
            let compiled = system.compile(p, &Default::default()).unwrap();
            for path in ["positive-real", "upper-half-plane"] {
                let initial = boundary(p, d_num, d_den, start, boundary_order);
                let waypoints = if path == "positive-real" {
                    vec![p.zero()]
                } else {
                    vec![p.complex(4,1), p.parse("1", "0.333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333").unwrap(), p.zero()]
                };
                let result = compiled.transport(
                    &BoundaryData { point: p.i(start), values: initial.clone() }, &waypoints,
                    &FlowOptions { digits: 30, guard_digits: digits-30, series_order: taylor_order, ..Default::default() },
                    &RunContext::default(),
                ).unwrap();
                let i = &result.values[0];
                let s = &result.values[1];
                let i_prime = p.sub(&p.scale(i, 2*(d_num-d_den), d_den), &p.scale(s,3,1));
                let physical_scalar = p.sub(&i_prime,s);
                let medium = p.add(&p.scale(i,3,2), &p.sub(&p.scale(&i_prime,1,4),&p.scale(s,5,4)));
                let render = |x: &C| json!({"re":x.re.to_string(),"im":x.im.to_string()});
                results.push(json!({"spatial_dimension":format!("{d_num}/{d_den}"),
                    "working_digits":digits,"transport_order":taylor_order,"start_eta":start,"boundary_order":boundary_order,
                    "path":path,"boundary":initial.iter().map(render).collect::<Vec<_>>(),
                    "I":render(i),"S":render(s),"eta_derivative":render(&i_prime),
                    "physical_scalar_mass_derivative":render(&physical_scalar),
                    "physical_medium_mass_derivative":render(&medium),
                    "transport_steps":result.diagnostics.steps}));
            }
        }
    }
    let result = json!({"schema":1,"status":"predictions saved before comparison",
        "scope":"exploratory two-state radial ODE through existing RustFlow Taylor transport; guarded proof separately audited",
        "mass_squared":"1/4","chemical_potential":"1","radius_squared":"3/4",
        "normalization":"angular volume omitted consistently",
        "boundary_provenance":"uniform convergent binomial expansion integrated over the fixed ball, no finite-eta or endpoint values supplied",
        "matrix":[["(d-1)/(2*(eta+a))","-B/(eta+a)"],["0","-1/(2*(eta+a+B))"]],
        "reference_values_read":false,"records":results});
    std::fs::write(output, serde_json::to_string_pretty(&result).unwrap()+"\n").unwrap();
    println!("saved 30 exploratory transport predictions, including independent orders/precision/start/path refinements");
}
