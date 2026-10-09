//! Independent validation-only massless sunset reference. No finite-density
//! production owner, supplied oracle or native prediction is imported.
//! See docs/finite-density-massless-reference.md for the common prescription.
use std::collections::BTreeMap;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat as C, Precision};

fn rational(p: Precision, n: i64, d: i64) -> C {
    p.rational(&Rational::from((n, d)))
}
fn gamma(p: Precision, x: &C) -> C {
    p.gamma_real(&x.re).unwrap()
}
fn pi(p: Precision) -> C {
    p.eval(&Atom::var(Symbol::PI), &ahash::HashMap::default())
        .unwrap()
}
fn beta(p: Precision, a: &C, b: &C) -> C {
    p.div(&p.mul(&gamma(p, a), &gamma(p, b)), &gamma(p, &p.add(a, b)))
}
fn ad(p: Precision, dim: &C) -> C {
    let d = p.sub(dim, &p.i(1));
    p.div(
        &p.pow(&p.i(2), &p.sub(&p.i(1), &d)),
        &p.mul(
            &p.pow(&pi(p), &p.scale(&d, 1, 2)),
            &gamma(p, &p.scale(&d, 1, 2)),
        ),
    )
}
const LABELS: [&str; 6] = [
    "scalar",
    "raised_bulk_measure",
    "raised_bulk_numerator",
    "raised_bulk_denominator",
    "raised_upper_surface",
    "raised_total",
];
fn total(p: Precision, pieces: &mut [C; 6]) {
    pieces[5] = (1..=4).fold(p.zero(), |s, i| p.add(&s, &pieces[i]));
}
fn relative(p: Precision, actual: &C, expected: &C) -> C {
    p.div(
        &C::new(p.norm(&p.sub(actual, expected)), p.real(0)),
        &C::new(p.norm(expected), p.real(0)),
    )
}
fn assert_relative(p: Precision, actual: &C, expected: &C, digits: u32) {
    assert!(p.norm(expected) > p.real(0));
    let error = relative(p, actual, expected);
    assert!(
        p.norm(&error) < p.tolerance(digits),
        "relative error {error}; actual {actual}; expected {expected}"
    );
}
fn encoded(values: &[C; 6]) -> BTreeMap<String, String> {
    LABELS
        .iter()
        .zip(values)
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

/// Analytically integrated radial factors and independent angular Beta factors.
/// Every raised bulk term and the moving upper surface remain separate.
fn beta_reference(p: Precision, dim: &C, mu: &C) -> [C; 6] {
    let n = p.sub(dim, &p.i(2));
    let alpha = p.scale(&n, 1, 2);
    let k = |b: i64| {
        p.div(
            &beta(p, &p.add(&alpha, &p.i(b)), &alpha),
            &beta(p, &alpha, &alpha),
        )
    };
    let km1 = k(-1);
    let km2 = k(-2);
    let den = p.mul(&n, &p.sub(&n, &p.i(2)));
    let common = p.scale(
        &p.mul(
            &p.powi(&ad(p, dim), 2),
            &p.pow(mu, &p.sub(&p.scale(dim, 2, 1), &p.i(6))),
        ),
        1,
        64,
    );
    let mut out = std::array::from_fn(|_| p.zero());
    out[0] = p.mul(
        &common,
        &p.div(&p.scale(&km1, 4, 1), &p.powi(&p.sub(&n, &p.i(1)), 2)),
    );
    out[1] = p.mul(
        &common,
        &p.div(&p.scale(&p.add(&km1, &p.i(2)), -2, 1), &den),
    );
    out[2] = p.mul(&common, &p.div(&p.scale(&km1, 4, 1), &den));
    out[3] = p.mul(
        &common,
        &p.mul(
            &p.sub(
                &p.div(&p.i(1), &p.powi(&p.sub(&n, &p.i(1)), 2)),
                &p.div(&p.i(1), &den),
            ),
            &p.add(&km2, &p.scale(&km1, 2, 1)),
        ),
    );
    out[4] = p.mul(&common, &p.div(&p.scale(&p.add(&km1, &p.i(2)), -2, 1), &n));
    total(p, &mut out);
    out
}
fn simplified_reference(p: Precision, dim: &C, mu: &C) -> [C; 2] {
    let common = p.mul(
        &p.powi(&ad(p, dim), 2),
        &p.pow(mu, &p.sub(&p.scale(dim, 2, 1), &p.i(6))),
    );
    let den = p.mul(&p.sub(dim, &p.i(3)), &p.sub(dim, &p.i(4)));
    let poly = p.sub(
        &p.add(
            &p.sub(
                &p.scale(&p.powi(dim, 3), 2, 1),
                &p.scale(&p.powi(dim, 2), 25, 1),
            ),
            &p.scale(dim, 98, 1),
        ),
        &p.i(119),
    );
    [
        p.div(&common, &p.scale(&den, 8, 1)),
        p.neg(&p.div(
            &p.mul(&common, &poly),
            &p.scale(
                &p.mul(&den, &p.mul(&p.sub(dim, &p.i(2)), &p.sub(dim, &p.i(6)))),
                16,
                1,
            ),
        )),
    ]
}
fn laurent(p: Precision, mu: &C) -> [[C; 3]; 2] {
    let euler = p
        .eval(
            &Atom::var(symbolica::transcendental::euler_gamma()),
            &ahash::HashMap::default(),
        )
        .unwrap();
    let l = p.sub(
        &p.add(
            &p.sub(&p.i(4), &p.scale(&euler, 2, 1)),
            &p.scale(&p.log(&pi(p)), 2, 1),
        ),
        &p.scale(&p.log(mu), 4, 1),
    );
    let pole = p.neg(&p.div(&p.powi(mu, 2), &p.scale(&p.powi(&pi(p), 4), 64, 1)));
    let raised_pole = p.scale(&pole, 1, 8);
    [
        [p.zero(), pole.clone(), p.mul(&pole, &p.add(&l, &p.i(2)))],
        [
            p.zero(),
            raised_pole.clone(),
            p.mul(&raised_pole, &p.add(&l, &p.i(14))),
        ],
    ]
}
fn legendre(p: Precision, n: usize, x: &C) -> (C, C) {
    let mut previous = p.i(1);
    let mut current = x.clone();
    for k in 2..=n {
        let next = p.scale(
            &p.sub(
                &p.scale(&p.mul(x, &current), (2 * k - 1) as i64, 1),
                &p.scale(&previous, (k - 1) as i64, 1),
            ),
            1,
            k as i64,
        );
        previous = current;
        current = next;
    }
    let derivative = p.scale(
        &p.div(
            &p.sub(&p.mul(x, &current), &previous),
            &p.sub(&p.mul(x, x), &p.i(1)),
        ),
        n as i64,
        1,
    );
    (current, derivative)
}

/// Gauss-Legendre nodes on [0,1]. Binary64 supplies Newton seeds only; nodes,
/// weights, residuals and moment checks use fixed native working precision.
fn quadrature(order: usize, digits: u32, p: Precision) -> Vec<(C, C)> {
    assert!(order >= 2 && digits >= 30);
    let mut nodes = Vec::with_capacity(order);
    for i in 0..order.div_ceil(2) {
        let seed = (std::f64::consts::PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
        let mut x = p.parse(&format!("{seed:.17}"), "0").unwrap();
        let mut converged = false;
        for _ in 0..100 {
            let (value, derivative) = legendre(p, order, &x);
            let delta = p.div(&value, &derivative);
            x = p.sub(&x, &delta);
            if p.norm(&delta) < p.tolerance(digits - 5) {
                converged = true;
                break;
            }
        }
        assert!(converged, "Legendre Newton iteration failed");
        let (residual, derivative) = legendre(p, order, &x);
        assert!(p.norm(&residual) < p.tolerance(digits - 8));
        let weight = p.div(
            &p.i(1),
            &p.mul(
                &p.sub(&p.i(1), &p.mul(&x, &x)),
                &p.mul(&derivative, &derivative),
            ),
        );
        nodes.push((p.scale(&p.sub(&p.i(1), &x), 1, 2), weight.clone()));
        if i != order - 1 - i {
            nodes.push((p.scale(&p.add(&p.i(1), &x), 1, 2), weight));
        }
    }
    nodes.sort_by(|a, b| a.0.re.partial_cmp(&b.0.re).unwrap());
    assert_eq!(nodes.len(), order);
    // Includes normalization and both even and odd polynomial moments.
    for degree in 0..=12.min(2 * order - 1) {
        let actual = nodes.iter().fold(p.zero(), |sum, (x, w)| {
            p.add(&sum, &p.mul(w, &p.powi(x, degree as i64)))
        });
        assert!(p.close(&actual, &rational(p, 1, degree as i64 + 1), digits - 10));
    }
    nodes
}

/// Direct original differentiated integrand at D=7, before radial/angular
/// analytic integration. t=u^2/(u^2+(1-u)^2) removes the square-root endpoints.
fn direct_d7(digits: u32, radial_order: usize, angular_order: usize, mu: &Rational) -> [C; 6] {
    let p = Precision::decimal(digits).unwrap();
    let mu = p.rational(mu);
    let radial = quadrature(radial_order, digits, p);
    let angular = quadrature(angular_order, digits, p);
    let ad = ad(p, &p.i(7));
    let bn = beta(p, &rational(p, 5, 2), &rational(p, 5, 2));
    let mut out = std::array::from_fn(|_| p.zero());
    for (u, uw) in angular {
        let v = p.sub(&p.i(1), &u);
        let q = p.add(&p.powi(&u, 2), &p.powi(&v, 2));
        let t = p.div(&p.powi(&u, 2), &q);
        let z = p.sub(&p.i(1), &p.scale(&t, 2, 1));
        let aw = p.div(
            &p.scale(&p.mul(&uw, &p.mul(&p.powi(&u, 4), &p.powi(&v, 4))), 2, 1),
            &p.mul(&p.powi(&q, 5), &bn),
        );
        for (x2, w2) in &radial {
            let r2 = p.mul(&mu, x2);
            let rw2 = p.scale(&p.mul(&ad, &p.mul(&mu, &p.mul(w2, &p.powi(&r2, 4)))), 1, 2);
            // Original moving upper boundary, before angular or radial integration.
            let hs = p.scale(&p.mul(&mu, &p.mul(&r2, &t)), 4, 1);
            let ns = p.mul(&p.mul(&mu, &r2), &p.sub(&z, &p.i(2)));
            let sw = p.scale(&p.mul(&ad, &p.powi(&mu, 3)), 1, 4);
            out[4] = p.add(
                &out[4],
                &p.mul(&aw, &p.mul(&sw, &p.mul(&rw2, &p.div(&ns, &hs)))),
            );
            for (x1, w1) in &radial {
                let r1 = p.mul(&mu, x1);
                let rw1 = p.scale(&p.mul(&ad, &p.mul(&mu, &p.mul(w1, &p.powi(&r1, 4)))), 1, 2);
                let h = p.scale(&p.mul(&r1, &p.mul(&r2, &t)), 4, 1);
                let numerator = p.mul(&p.mul(&r1, &r2), &p.sub(&z, &p.i(2)));
                let ha = p.sub(&p.div(&r2, &r1), &p.i(1));
                let weight = p.mul(&aw, &p.mul(&rw1, &rw2));
                let kernels = [
                    p.div(&p.i(1), &h),
                    p.div(&numerator, &p.scale(&p.mul(&p.powi(&r1, 2), &h), 2, 1)),
                    p.div(&r2, &p.mul(&r1, &h)),
                    p.div(&p.mul(&numerator, &ha), &p.powi(&h, 2)),
                ];
                for i in 0..4 {
                    out[i] = p.add(&out[i], &p.mul(&weight, &kernels[i]));
                }
            }
        }
    }
    total(p, &mut out);
    out
}

#[test]
fn massless_reference_binds_original_input() {
    let input: serde_json::Value = serde_json::from_str(include_str!(
        "../examples/finite_density/massless_two_loop_sunset.json"
    ))
    .unwrap();
    assert_eq!(input["loops"], 2);
    assert_eq!(input["vertices"], 2);
    assert_eq!(input["numerator_convention"], "shifted_euclidean");
    assert_eq!(input["chemical_potentials"], serde_json::json!(["1"]));
    assert_eq!(input["loop_charges"], serde_json::json!([[1], [1]]));
    let edges = input["edges"].as_array().unwrap();
    assert_eq!(edges.len(), 3);
    for (i, (routing, vertices, charge)) in [
        (["1", "0"], [0, 1], 1),
        (["0", "1"], [1, 0], 1),
        (["1", "-1"], [1, 0], 0),
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(edges[i]["routing"], serde_json::json!(routing));
        assert_eq!(edges[i]["vertices"], serde_json::json!(vertices));
        assert_eq!(edges[i]["charges"], serde_json::json!([charge]));
        assert_eq!(edges[i]["mass_squared"], "0");
    }
    assert_eq!(
        input["targets"],
        serde_json::json!([
        {"powers":[1,1,1],"numerator":"1"},
        {"powers":[2,1,1],"numerator":"g1_2+u1*u2"}])
    );
    assert_eq!(input["laurent_orders"], serde_json::json!([-2, 0]));
}

#[test]
fn independent_beta_moments_match_closed_expression_and_scaling() {
    let p = Precision::decimal(80).unwrap();
    for dim in [
        Rational::from(7),
        Rational::from((13, 2)),
        Rational::from((15, 4)),
    ] {
        for mu in [Rational::one(), Rational::from((3, 2))] {
            let d = p.rational(&dim);
            let m = p.rational(&mu);
            let beta = beta_reference(p, &d, &m);
            let simple = simplified_reference(p, &d, &m);
            assert_relative(p, &beta[0], &simple[0], 65);
            assert_relative(p, &beta[5], &simple[1], 65);
            let unit = simplified_reference(p, &d, &p.i(1));
            let scale = p.pow(&m, &p.sub(&p.scale(&d, 2, 1), &p.i(6)));
            for i in 0..2 {
                assert_relative(p, &simple[i], &p.mul(&unit[i], &scale), 65);
            }
        }
    }
    let at7 = simplified_reference(p, &p.i(7), &p.i(1));
    assert_relative(
        p,
        &at7[0],
        &p.div(&p.i(1), &p.scale(&p.powi(&pi(p), 6), 393216, 1)),
        65,
    );
    assert_relative(
        p,
        &at7[1],
        &p.div(&p.i(-7), &p.scale(&p.powi(&pi(p), 6), 983040, 1)),
        65,
    );
}

#[test]
fn independent_laurent_poles_and_constants_match_small_epsilon_limits() {
    let p = Precision::decimal(100).unwrap();
    for mu in [Rational::one(), Rational::from((3, 2))] {
        let m = p.rational(&mu);
        let coefficients = laurent(p, &m);
        for power in [20, 30] {
            let eps = p.powi(&p.i(10), -power);
            let values = simplified_reference(p, &p.sub(&p.i(4), &p.scale(&eps, 2, 1)), &m);
            for i in 0..2 {
                assert_eq!(p.norm(&coefficients[i][0]), p.real(0));
                assert_relative(p, &p.mul(&eps, &values[i]), &coefficients[i][1], 17);
                let finite = p.sub(&values[i], &p.div(&coefficients[i][1], &eps));
                assert_relative(p, &finite, &coefficients[i][2], 16);
            }
        }
    }
}

#[test]
fn convergent_direct_quadrature_keeps_the_raised_upper_surface() {
    let p = Precision::decimal(60).unwrap();
    let actual = direct_d7(60, 6, 32, &Rational::one());
    let expected = beta_reference(p, &p.i(7), &p.i(1));
    for i in 0..6 {
        assert_relative(p, &actual[i], &expected[i], 17);
    }
    assert!(actual[4].re < p.real(0));
    // At D=7 the bulk sum cancels: omitting the upper surface loses the answer.
    let bulk = (1..4).fold(p.zero(), |s, i| p.add(&s, &actual[i]));
    assert!(p.norm(&bulk) < p.tolerance(25));
}

fn sha256(path: &str) -> String {
    let out = std::process::Command::new("sha256sum")
        .arg(path)
        .output()
        .unwrap();
    assert!(out.status.success());
    String::from_utf8(out.stdout)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .to_owned()
}

#[test]
#[ignore = "explicit independent reference generation with precision and quadrature refinements"]
fn generate_independent_massless_sunset_reference() {
    let start = Instant::now();
    let p = Precision::decimal(90).unwrap();
    let profiles = [(50, 6, 32), (80, 6, 32), (80, 8, 48)];
    let mut quadratures = Vec::new();
    let mut previous: Option<[C; 6]> = None;
    let expected = beta_reference(p, &p.i(7), &p.i(1));
    for (digits, radial_order, angular_order) in profiles {
        let values = direct_d7(digits, radial_order, angular_order, &Rational::one());
        let errors: Vec<_> = (0..6)
            .map(|i| relative(p, &values[i], &expected[i]).to_string())
            .collect();
        for i in 0..6 {
            assert_relative(p, &values[i], &expected[i], 17);
        }
        let refinements = previous.as_ref().map(|old| {
            (0..6)
                .map(|i| relative(p, &values[i], &old[i]).to_string())
                .collect::<Vec<_>>()
        });
        if let Some(old) = &previous {
            for i in 0..6 {
                assert_relative(p, &values[i], &old[i], 17);
            }
        }
        quadratures.push(serde_json::json!({"digits":digits,"radial_order":radial_order,"angular_order":angular_order,"dimension":"7","chemical_potential":"1","values":encoded(&values),"relative_reference_differences":errors,"relative_previous_differences":refinements}));
        previous = Some(values);
    }
    let mut samples = Vec::new();
    for dim in [
        Rational::from(7),
        Rational::from((13, 2)),
        Rational::from((15, 4)),
    ] {
        for mu in [Rational::one(), Rational::from((3, 2))] {
            let values = beta_reference(p, &p.rational(&dim), &p.rational(&mu));
            let lowp = Precision::decimal(60).unwrap();
            let low = beta_reference(lowp, &lowp.rational(&dim), &lowp.rational(&mu));
            for i in 0..6 {
                assert_relative(p, &low[i], &values[i], 50);
            }
            samples.push(serde_json::json!({"dimension":dim.to_string(),"chemical_potential":mu.to_string(),"values":encoded(&values),"vacuum":["0","0"],"single_cut_0":["0","0"],"single_cut_1":["0","0"],"double_cut":[values[0].to_string(),values[5].to_string()],"total":[values[0].to_string(),values[5].to_string()],"precision_refinement_digits":[60,90]}));
        }
    }
    let coefficients = laurent(p, &p.i(1));
    let coefficient_strings: Vec<Vec<_>> = coefficients
        .iter()
        .map(|xs| xs.iter().map(ToString::to_string).collect())
        .collect();
    let input = "examples/finite_density/massless_two_loop_sunset.json";
    let source = "tests/finite_density_massless_reference.rs";
    let report = serde_json::json!({"schema_version":1,"status":"independent_reference_generated","normalization":"raw_euclidean_dDp_over_2pi_power_D_per_loop","scope":"massless two-loop sunset, original scalar and raised polynomial numerator, common regulated thermal dimensional continuation","method":"convergent D>6 radial/angular Beta integrals with original mass derivative and moving upper surface; meromorphic continuation after combination","prescription":"vacuum scaleless; each one-cut amplitude is zero only under the common one-virtual-loop thermal cone and finite-jet endpoint prescription","oracle_records_read":0,"native_predictions_read":0,"source":{"input":input,"input_sha256":sha256(input),"reference":source,"reference_sha256":sha256(source),"cargo_lock_sha256":sha256("Cargo.lock")},"samples":samples,"laurent":{"dimension":"4-2*epsilon","chemical_potential":"1","orders":[-2,-1,0],"target_order":["scalar","raised_original_numerator"],"total":coefficient_strings},"direct_convergent_quadrature":quadratures,"uncertainty":"empirical precision and quadrature refinement, checked against independent exact Beta formulas; no rigorous numerical error interval","acceptance":"17 relative digits for direct D7 integrals; 50 digits for independent working precision change; numerical native acceptance is separate","wall_seconds":start.elapsed().as_secs_f64()});
    let output = std::env::var("RUSTFLOW_DENSITY_MASSLESS_REFERENCE_REPORT")
        .unwrap_or_else(|_| "/tmp/rustflow-massless-sunset-reference.json".into());
    if let Some(parent) = std::path::Path::new(&output).parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    std::fs::write(
        &output,
        format!("{}\n", serde_json::to_string_pretty(&report).unwrap()),
    )
    .unwrap();
    println!("independent reference saved to {output}");
}
