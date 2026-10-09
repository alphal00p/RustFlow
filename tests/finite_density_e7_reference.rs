//! Validation-only analytic reference, independent of production integral owners.
//! No supplied oracle or native AMF prediction is read. See the full derivation
//! and analytic-regulator qualifications in docs/finite-density-e7-reference.md.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat as C, Precision};

fn gamma(p: Precision, x: &C) -> C {
    p.gamma_real(&x.re).unwrap()
}
fn pi(p: Precision) -> C {
    p.eval(&Atom::var(Symbol::PI), &ahash::HashMap::default())
        .unwrap()
}
fn euler_gamma(p: Precision) -> C {
    p.eval(
        &Atom::var(symbolica::transcendental::euler_gamma()),
        &ahash::HashMap::default(),
    )
    .unwrap()
}
fn rational(p: Precision, n: i64, d: i64) -> C {
    p.rational(&Rational::from((n, d)))
}
fn relative(p: Precision, a: &C, b: &C) -> C {
    p.div(
        &C::new(p.norm(&p.sub(a, b)), p.real(0)),
        &C::new(p.norm(b), p.real(0)),
    )
}

/// Unsimplified neutral tensor factors and independently differentiated compact
/// beta integrals, followed by the original four-loop measure conversion.
fn beta_reference(p: Precision, epsilon: &Rational) -> [C; 2] {
    let eps = p.rational(epsilon);
    let dim = p.sub(&p.i(4), &p.scale(&eps, 2, 1));
    let spatial = p.sub(&dim, &p.i(1));
    let half = p.scale(&dim, 1, 2);
    let pi = pi(p);
    let four_pi = p.scale(&pi, 4, 1);
    let bubble = p.div(
        &p.mul(
            &gamma(p, &p.sub(&p.i(2), &half)),
            &p.powi(&gamma(p, &p.sub(&half, &p.i(1))), 2),
        ),
        &p.mul(&p.pow(&four_pi, &half), &gamma(p, &p.sub(&dim, &p.i(2)))),
    );
    let ad = p.div(
        &p.scale(&p.pow(&pi, &p.scale(&spatial, 1, 2)), 2, 1),
        &p.mul(
            &gamma(p, &p.scale(&spatial, 1, 2)),
            &p.pow(&p.scale(&pi, 2, 1), &spatial),
        ),
    );
    let angle = |c: &C| {
        p.div(
            &p.mul(
                &p.pow(&p.i(2), &p.sub(&p.sub(&spatial, &p.i(2)), c)),
                &p.mul(
                    &gamma(p, &p.scale(&spatial, 1, 2)),
                    &gamma(p, &p.sub(&p.scale(&p.sub(&spatial, &p.i(1)), 1, 2), c)),
                ),
            ),
            &p.mul(
                &p.pow(&pi, &rational(p, 1, 2)),
                &gamma(p, &p.sub(&p.sub(&spatial, &p.i(1)), c)),
            ),
        )
    };
    let c1 = p.sub(&p.i(4), &dim);
    let n1 = p.sub(&p.sub(&spatial, &p.i(1)), &c1);
    let simple = p.div(
        &p.mul(
            &p.powi(&ad, 2),
            &p.mul(&p.pow(&p.i(2), &p.neg(&p.add(&c1, &p.i(2)))), &angle(&c1)),
        ),
        &p.powi(&n1, 2),
    );
    let c2 = p.sub(&p.i(3), &dim);
    let n2 = p.sub(&p.sub(&spatial, &p.i(1)), &c2);
    let denominator = p.mul(&n2, &p.sub(&n2, &p.i(2)));
    let common = p.mul(
        &p.powi(&ad, 2),
        &p.pow(&p.i(2), &p.neg(&p.add(&c2, &p.i(3)))),
    );
    let bulk = p.add(
        &p.div(&angle(&c2), &denominator),
        &p.mul(
            &p.mul(&c2, &angle(&p.add(&c2, &p.i(1)))),
            &p.sub(
                &p.div(&p.i(1), &denominator),
                &p.div(&p.i(1), &p.powi(&p.sub(&n2, &p.i(1)), 2)),
            ),
        ),
    );
    let surface = p.div(&angle(&c2), &n2);
    let raised = p.mul(&common, &p.add(&bulk, &surface));
    let prefactor = p.div(&p.powi(&bubble, 2), &p.scale(&p.sub(&dim, &p.i(1)), 16, 1));
    let first = p.mul(&prefactor, &p.mul(&dim, &simple));
    let second = p.mul(
        &prefactor,
        &p.add(
            &p.mul(
                &p.scale(&p.sub(&p.scale(&dim, 2, 1), &p.i(3)), 2, 1),
                &simple,
            ),
            &p.mul(&p.sub(&p.scale(&dim, 3, 1), &p.i(2)), &raised),
        ),
    );
    let normalization = p.mul(
        &p.powi(&four_pi, 8),
        &p.exp(&p.mul(&p.scale(&eps, 4, 1), &p.sub(&euler_gamma(p), &p.log(&pi)))),
    );
    [
        p.mul(&normalization, &first),
        p.mul(&normalization, &second),
    ]
}

fn simplified_reference(p: Precision, epsilon: &Rational) -> [C; 2] {
    let e = p.rational(epsilon);
    let dim = p.sub(&p.i(4), &p.scale(&e, 2, 1));
    let affine = |n| p.sub(&p.i(1), &p.scale(&e, n, 1));
    let g = p.div(
        &p.mul(
            &p.exp(&p.mul(&p.scale(&e, 4, 1), &euler_gamma(p))),
            &p.mul(
                &p.powi(&gamma(p, &p.add(&p.i(1), &e)), 2),
                &p.mul(&p.powi(&gamma(p, &affine(1)), 5), &gamma(p, &affine(3))),
            ),
        ),
        &p.mul(&p.powi(&gamma(p, &affine(2)), 3), &gamma(p, &affine(4))),
    );
    let first = p.div(
        &p.mul(&p.sub(&p.i(2), &e), &g),
        &p.mul(
            &p.powi(&e, 2),
            &p.scale(
                &p.mul(
                    &p.sub(&p.i(3), &p.scale(&e, 2, 1)),
                    &p.mul(&p.powi(&affine(2), 5), &affine(4)),
                ),
                2,
                1,
            ),
        ),
    );
    let ratio = p.div(
        &p.add(
            &p.sub(
                &p.scale(&p.powi(&dim, 3), 18, 1),
                &p.scale(&p.powi(&dim, 2), 109, 1),
            ),
            &p.sub(&p.scale(&dim, 191, 1), &p.i(72)),
        ),
        &p.mul(&dim, &p.sub(&p.scale(&dim, 2, 1), &p.i(5))),
    );
    [first.clone(), p.mul(&ratio, &first)]
}

#[test]
fn e7_reference_binds_the_complete_original_definition() {
    let input: serde_json::Value = serde_json::from_str(include_str!(
        "../examples/finite_density/chain_of_three_parallel_pairs.json"
    ))
    .unwrap();
    assert_eq!(input["loops"], 4);
    assert_eq!(input["vertices"], 4);
    assert_eq!(input["edges"].as_array().unwrap().len(), 7);
    assert_eq!(input["numerator_convention"], "shifted_euclidean");
    assert_eq!(input["chemical_potentials"], serde_json::json!(["1"]));
    assert_eq!(
        input["loop_charges"],
        serde_json::json!([[1], [0], [0], [0]])
    );
    let rows = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, -1],
        [0, 1, 0, -1],
        [0, 0, 1, -1],
    ];
    let incidence = [[0, 3], [3, 2], [2, 1], [1, 0], [3, 0], [2, 3], [1, 2]];
    for (i, edge) in input["edges"].as_array().unwrap().iter().enumerate() {
        assert_eq!(edge["mass_squared"], "0");
        assert_eq!(
            edge["routing"],
            serde_json::json!(rows[i].map(|n| n.to_string()))
        );
        assert_eq!(edge["vertices"], serde_json::json!(incidence[i]));
        assert_eq!(
            edge["charges"],
            serde_json::json!([if [0, 4].contains(&i) { 1 } else { 0 }])
        );
    }
    assert_eq!(
        input["targets"][0],
        serde_json::json!({"powers":[1,1,1,2,1,1,1],"numerator":"g1_2^2"})
    );
    assert_eq!(
        input["targets"][1],
        serde_json::json!({"powers":[2,1,1,1,1,1,1],"numerator":"g1_2^2+g1_3*g2_4"})
    );
}

#[test]
#[ignore = "explicit validation-only native gamma reference generation"]
fn generate_independent_e7_analytic_reference() {
    let started = std::time::Instant::now();
    let mut evaluations = Vec::new();
    for epsilon in [
        Rational::from((1, 5)),
        Rational::from((1, 10)),
        Rational::from((1, 20)),
    ] {
        let low = Precision::decimal(50).unwrap();
        let high = Precision::decimal(80).unwrap();
        let baseline = beta_reference(low, &epsilon);
        let refined = beta_reference(high, &epsilon);
        let simplified = simplified_reference(high, &epsilon);
        for target in 0..2 {
            let precision_change = relative(high, &baseline[target], &refined[target]);
            let algebra_change = relative(high, &simplified[target], &refined[target]);
            assert!(precision_change.re < high.tolerance(40));
            assert!(algebra_change.re < high.tolerance(65));
            evaluations.push(serde_json::json!({
                "epsilon":epsilon.to_string(),"target_index":target,
                "baseline_50_digits":baseline[target].to_string(),
                "refined_80_digits":refined[target].to_string(),
                "simplified_80_digits":simplified[target].to_string(),
                "precision_relative_change":precision_change.re.as_raw().to_string(),
                "independent_expression_relative_change":algebra_change.re.as_raw().to_string(),
            }));
        }
    }
    let p = Precision::decimal(100).unwrap();
    let exact = [
        [
            rational(p, 1, 3),
            rational(p, 85, 18),
            p.sub(&rational(p, 1066, 27), &p.scale(&p.powi(&pi(p), 2), 1, 3)),
        ],
        [
            rational(p, 25, 9),
            rational(p, 617, 18),
            p.sub(&rational(p, 20887, 81), &p.scale(&p.powi(&pi(p), 2), 25, 9)),
        ],
    ];
    let mut expansion_checks = Vec::new();
    for power in [20_u32, 30] {
        let epsilon = Rational::one() / Rational::from(Integer::from(10).pow(power as u64));
        let e = p.rational(&epsilon);
        let values = beta_reference(p, &epsilon);
        for target in 0..2 {
            let finite = p.div(
                &p.sub(
                    &p.sub(&p.mul(&p.powi(&e, 2), &values[target]), &exact[target][0]),
                    &p.mul(&e, &exact[target][1]),
                ),
                &p.powi(&e, 2),
            );
            let change = relative(p, &finite, &exact[target][2]);
            assert!(change.re < p.tolerance(12));
            expansion_checks.push(serde_json::json!({"epsilon":epsilon.to_string(),"target_index":target,"finite_coefficient_extrapolant":finite.to_string(),"exact_expression_value":exact[target][2].to_string(),"relative_difference":change.re.as_raw().to_string()}));
        }
    }
    let coefficients = exact
        .iter()
        .map(|values| {
            [-2, -1, 0]
                .iter()
                .zip(values)
                .map(|(power, value)| (power.to_string(), value.to_string()))
                .collect::<BTreeMap<_, _>>()
        })
        .collect::<Vec<_>>();
    let output = serde_json::json!({
        "schema":1,"status":"reference_checks_passed","oracle_records_read":0,"feature_predictions_read":0,
        "method":"independent regulated neutral tensor bubbles plus raised compact beta integrals, retaining the upper surface; gamma-duplication cross-check and explicit Laurent expansion",
        "normalization":"MSbar original loop measure divided by unexpanded (4*pi)^(-8)*(Lambda_bar/2)^(8*epsilon)",
        "definition":serde_json::from_str::<serde_json::Value>(include_str!("../examples/finite_density/chain_of_three_parallel_pairs.json")).unwrap(),
        "analytic_coefficients":[{"-2":"1/3","-1":"85/18","0":"1066/27-pi^2/3"},{"-2":"25/9","-1":"617/18","0":"20887/81-25*pi^2/9"}],
        "orders_below_minus_two":"zero by the meromorphic Gamma expression; higher positive orders are not supplied",
        "reference_error_model":"exact analytic coefficient expressions; finite-precision decimal evaluations checked separately",
        "native_100_digit_coefficients":coefficients,"evaluations":evaluations,"expansion_checks":expansion_checks,
        "source_blake3":blake3::hash(include_bytes!("finite_density_e7_reference.rs")).to_hex().to_string(),
        "cargo_lock_blake3":blake3::hash(include_bytes!("../Cargo.lock")).to_hex().to_string(),
        "wall_seconds":started.elapsed().as_secs_f64(),
    });
    let path = std::env::var_os("RUSTFLOW_DENSITY_E7_REFERENCE_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-density-e7-reference.json"));
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    std::fs::write(&path, serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    println!(
        "Independent E7 analytic reference saved to {}",
        path.display()
    );
}
