//! Independent, validation-only quadrature; no AMF owner or oracle answers.
//! The derivation and endpoint maps are recorded in finite-density-reference.md.
use std::collections::BTreeMap;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat as C, Precision};

type Contributions = BTreeMap<String, C>;

fn rational(p: Precision, n: i64, d: i64) -> C {
    p.rational(&Rational::from((n, d)))
}
fn power(p: Precision, x: &C, n: i64, d: i64) -> C {
    p.pow(x, &rational(p, n, d))
}
fn gamma(p: Precision, n: i64, d: i64) -> C {
    p.gamma_real(&rational(p, n, d).re).unwrap()
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

#[derive(Clone)]
struct RadialNode {
    r: C,
    e: C,
    raw: C,
    occupied: C,
}

fn radial_nodes(p: Precision, nodes: &[(C, C)], mass: &C, ad: &C) -> (C, Vec<RadialNode>) {
    let gap = p.sub(&p.i(1), mass);
    assert!(gap.re > p.real(0));
    let radius = power(p, &gap, 1, 2);
    let prefactor = p.scale(&p.mul(ad, &power(p, &radius, 7, 5)), 5, 1);
    let values = nodes
        .iter()
        .map(|(s, w)| {
            let r = p.mul(&radius, &p.powi(s, 5));
            let e = power(p, &p.add(mass, &p.mul(&r, &r)), 1, 2);
            let raw = p.mul(&prefactor, &p.mul(w, &p.powi(s, 6)));
            let occupied = p.div(&raw, &p.scale(&e, 2, 1));
            RadialNode {
                r,
                e,
                raw,
                occupied,
            }
        })
        .collect();
    (radius, values)
}

fn reference(order: usize, digits: u32) -> Contributions {
    let p = Precision::decimal(digits).unwrap();
    let nodes = quadrature(order, digits, p);
    let pi = p
        .eval(&Atom::var(Symbol::PI), &ahash::HashMap::default())
        .unwrap();
    let four_pi = p.scale(&pi, 4, 1);
    let a = rational(p, 1, 4);
    let b = rational(p, 1, 4);
    let c = p.i(1);
    let masses = [&a, &b, &c];
    let vacuum_scale = p.div(&gamma(p, 3, 5), &power(p, &four_pi, 12, 5));
    let bubble_scale = p.div(&gamma(p, 4, 5), &power(p, &four_pi, 6, 5));
    let ad = p.div(
        &p.scale(&power(p, &pi, 7, 10), 2, 1),
        &p.mul(&gamma(p, 7, 10), &power(p, &p.scale(&pi, 2, 1), 7, 5)),
    );
    let angle_normalization = p.div(
        &gamma(p, 7, 10),
        &p.mul(&power(p, &pi, 1, 2), &gamma(p, 1, 5)),
    );
    let (rf1, first) = radial_nodes(p, &nodes, &a, &ad);
    let (_, second) = radial_nodes(p, &nodes, &b, &ad);
    let mut angles = Vec::with_capacity(2 * order);
    for (s, w) in &nodes {
        let fifth = p.powi(s, 5);
        let z = p.sub(&p.i(1), &fifth);
        let weight = p.scale(
            &p.mul(
                &angle_normalization,
                &p.mul(w, &power(p, &p.sub(&p.i(2), &fifth), -4, 5)),
            ),
            5,
            1,
        );
        angles.push((z.clone(), weight.clone()));
        angles.push((p.neg(&z), weight));
    }
    let angle_sum = angles.iter().fold(p.zero(), |v, (_, w)| p.add(&v, w));
    assert!(
        p.close(&angle_sum, &p.i(1), 15),
        "angular normalization failed: {angle_sum}"
    );

    let mut vacuum_scalar = p.zero();
    let mut vacuum_raised = p.zero();
    for [j, k, l] in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for (s, ws) in &nodes {
            let t = p.powi(s, 5);
            let prefactor = p.scale(&p.mul(ws, &p.powi(s, 3)), 5, 1);
            for (w, ww) in &nodes {
                let h = p.add(&p.add(&p.i(1), w), &p.mul(&t, w));
                let mass = p.add(
                    masses[j],
                    &p.mul(&t, &p.add(masses[k], &p.mul(w, masses[l]))),
                );
                let kernel = p.mul(
                    &p.mul(&prefactor, ww),
                    &p.mul(&power(p, &h, -6, 5), &power(p, &mass, -3, 5)),
                );
                let ratio_numerator = if j == 1 {
                    p.mul(&t, w)
                } else if k == 1 {
                    w.clone()
                } else {
                    p.i(1)
                };
                vacuum_scalar = p.add(&vacuum_scalar, &kernel);
                vacuum_raised = p.add(
                    &vacuum_raised,
                    &p.mul(&kernel, &p.div(&ratio_numerator, &h)),
                );
            }
        }
    }
    vacuum_scalar = p.mul(&vacuum_scale, &vacuum_scalar);
    vacuum_raised = p.scale(&p.mul(&vacuum_scale, &vacuum_raised), 17, 10);

    let mut bubbles = [p.zero(), p.zero()];
    let mut vectors = [p.zero(), p.zero()];
    let mut derivative1 = p.zero();
    let mut minus_derivative2 = p.zero();
    for (x, w) in &nodes {
        let one_minus = p.sub(&p.i(1), x);
        let product = p.mul(x, &one_minus);
        let f = [
            p.sub(
                &p.add(&p.mul(x, &b), &p.mul(&one_minus, &c)),
                &p.mul(&product, &a),
            ),
            p.sub(
                &p.add(&p.mul(x, &a), &p.mul(&one_minus, &c)),
                &p.mul(&product, &b),
            ),
        ];
        for i in 0..2 {
            let kernel = p.mul(w, &power(p, &f[i], -4, 5));
            bubbles[i] = p.add(&bubbles[i], &kernel);
            vectors[i] = p.add(&vectors[i], &p.mul(&one_minus, &kernel));
        }
        derivative1 = p.add(
            &derivative1,
            &p.mul(
                w,
                &p.mul(&product, &p.mul(&one_minus, &power(p, &f[0], -9, 5))),
            ),
        );
        minus_derivative2 = p.add(
            &minus_derivative2,
            &p.mul(w, &p.mul(&product, &power(p, &f[1], -9, 5))),
        );
    }
    for i in 0..2 {
        bubbles[i] = p.mul(&bubble_scale, &bubbles[i]);
        vectors[i] = p.mul(&bubble_scale, &vectors[i]);
    }
    derivative1 = p.scale(&p.mul(&bubble_scale, &derivative1), 4, 5);
    minus_derivative2 = p.scale(&p.mul(&bubble_scale, &minus_derivative2), 4, 5);
    let mut compact = [p.zero(), p.zero()];
    let mut tensor = [p.zero(), p.zero()];
    let mut tensor_derivative_bulk = p.zero();
    for (i, radial) in [&first, &second].into_iter().enumerate() {
        for node in radial {
            compact[i] = p.add(&compact[i], &node.occupied);
            tensor[i] = p.add(
                &tensor[i],
                &p.mul(&node.occupied, &p.add(masses[i], &p.mul(&node.e, &node.e))),
            );
            if i == 0 {
                tensor_derivative_bulk = p.add(
                    &tensor_derivative_bulk,
                    &p.mul(
                        &node.raw,
                        &p.div(
                            &p.sub(&p.scale(&p.mul(&node.e, &node.e), 3, 1), &a),
                            &p.scale(&p.powi(&node.e, 3), 4, 1),
                        ),
                    ),
                );
            }
        }
    }
    let surface_factor = p.scale(&p.mul(&ad, &power(p, &rf1, -3, 5)), 1, 4);
    let tensor_derivative_surface = p.neg(&p.mul(&surface_factor, &p.add(&a, &p.i(1))));

    let mut full_scalar = p.zero();
    let mut full_raised_bulk = p.zero();
    for n1 in &first {
        for n2 in &second {
            let radial_weight = p.mul(&n1.occupied, &n2.occupied);
            let ee = p.mul(&n1.e, &n2.e);
            let rr = p.mul(&n1.r, &n2.r);
            let ha = p.sub(&p.div(&n2.e, &n1.e), &p.i(1));
            for (z, angular_weight) in &angles {
                let rz = p.mul(&rr, z);
                let h = p.add(&p.sub(&p.sub(&c, &a), &b), &p.scale(&p.sub(&ee, &rz), 2, 1));
                assert!(h.re > p.real(0));
                let n = p.sub(&rz, &p.scale(&ee, 2, 1));
                let weight = p.mul(&radial_weight, angular_weight);
                full_scalar = p.add(&full_scalar, &p.div(&weight, &h));
                // Exact combination of -d/da of the shell Jacobian and N/h.
                // n+2E1E2=rr*z avoids subtracting two equal energy terms.
                let bulk = p.add(
                    &p.div(&rz, &p.scale(&p.mul(&p.mul(&n1.e, &n1.e), &h), 2, 1)),
                    &p.div(&p.mul(&n, &ha), &p.mul(&h, &h)),
                );
                full_raised_bulk = p.add(&full_raised_bulk, &p.mul(&weight, &bulk));
            }
        }
    }
    let mut full_raised_surface = p.zero();
    for n2 in &second {
        for (z, angular_weight) in &angles {
            let rz = p.mul(&p.mul(&rf1, &n2.r), z);
            let h = p.add(
                &p.sub(&p.sub(&c, &a), &b),
                &p.scale(&p.sub(&n2.e, &rz), 2, 1),
            );
            let n = p.sub(&rz, &p.scale(&n2.e, 2, 1));
            full_raised_surface = p.add(
                &full_raised_surface,
                &p.mul(&p.mul(&n2.occupied, angular_weight), &p.div(&n, &h)),
            );
        }
    }
    full_raised_surface = p.mul(&surface_factor, &full_raised_surface);
    let mut out: Contributions = BTreeMap::from([
        ("scalar/vacuum".into(), vacuum_scalar),
        (
            "scalar/cut_0".into(),
            p.neg(&p.mul(&compact[0], &bubbles[0])),
        ),
        (
            "scalar/cut_1".into(),
            p.neg(&p.mul(&compact[1], &bubbles[1])),
        ),
        ("scalar/cut_01".into(), full_scalar),
        ("raised_numerator/vacuum".into(), vacuum_raised),
        (
            "raised_numerator/cut_0_bulk".into(),
            p.neg(&p.add(
                &p.mul(&derivative1, &tensor[0]),
                &p.mul(&vectors[0], &tensor_derivative_bulk),
            )),
        ),
        (
            "raised_numerator/cut_0_surface".into(),
            p.neg(&p.mul(&vectors[0], &tensor_derivative_surface)),
        ),
        (
            "raised_numerator/cut_1".into(),
            p.mul(&minus_derivative2, &tensor[1]),
        ),
        ("raised_numerator/cut_01_bulk".into(), full_raised_bulk),
        (
            "raised_numerator/cut_01_surface".into(),
            full_raised_surface,
        ),
    ]);
    for target in ["scalar", "raised_numerator"] {
        let total = out
            .iter()
            .filter(|(key, _)| key.starts_with(target))
            .fold(p.zero(), |sum, (_, value)| p.add(&sum, value));
        out.insert(format!("{target}/total"), total);
    }
    for (name, value) in &out {
        assert!(
            p.finite(value) && value.im.is_zero(),
            "nonreal or nonfinite {name}: {value}"
        );
    }
    assert!(p.norm(&out["scalar/vacuum"]) > p.tolerance(10));
    assert!(p.norm(&out["raised_numerator/cut_0_surface"]) > p.tolerance(10));
    assert!(p.norm(&out["raised_numerator/cut_01_surface"]) > p.tolerance(10));
    out
}

#[test]
fn validation_quadrature_reproduces_polynomial_moments() {
    for digits in [40, 60] {
        quadrature(16, digits, Precision::decimal(digits).unwrap());
    }
}

#[test]
#[ignore = "explicit independent massive sunset reference generation with node/precision refinement"]
fn generate_independent_massive_sunset_reference() {
    let orders = std::env::var("RUSTFLOW_DENSITY_REFERENCE_ORDERS")
        .unwrap_or_else(|_| "32,48,64".into())
        .split(',')
        .map(|v| v.parse::<usize>().unwrap())
        .collect::<Vec<_>>();
    assert!(orders.len() >= 3 && orders.windows(2).all(|v| v[0] < v[1]));
    let directory = std::env::var_os("RUSTFLOW_DENSITY_REFERENCE_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-density-independent-reference"));
    std::fs::create_dir_all(&directory).unwrap();
    let input: serde_json::Value = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    assert_eq!(input["targets"][0]["powers"], serde_json::json!([1, 1, 1]));
    assert_eq!(input["targets"][1]["powers"], serde_json::json!([2, 1, 1]));
    assert_eq!(input["targets"][1]["numerator"], "g1_2+u1*u2");
    assert_eq!(input["edges"][0]["mass_squared"], "1/4");
    assert_eq!(input["edges"][1]["mass_squared"], "1/4");
    assert_eq!(input["edges"][2]["mass_squared"], "1");
    assert_eq!(input["chemical_potentials"], serde_json::json!(["1"]));
    assert_eq!(input["numerator_convention"], "shifted_euclidean");
    assert_eq!(input["loops"], 2);
    assert_eq!(input["vertices"], 2);
    assert_eq!(input["targets"][0]["numerator"], "1");
    assert_eq!(input["loop_charges"], serde_json::json!([[1], [1]]));
    for (slot, (vertices, routing, charge)) in [
        ([0, 1], ["1", "0"], 1),
        ([1, 0], ["0", "1"], 1),
        ([1, 0], ["1", "-1"], 0),
    ]
    .into_iter()
    .enumerate()
    {
        assert_eq!(input["edges"][slot]["vertices"], serde_json::json!(vertices));
        assert_eq!(input["edges"][slot]["routing"], serde_json::json!(routing));
        assert_eq!(input["edges"][slot]["charges"], serde_json::json!([charge]));
    }
    assert_eq!(input["edges"].as_array().unwrap().len(), 3);
    assert_eq!(input["targets"].as_array().unwrap().len(), 2);
    let mut snapshots = Vec::new();
    let mut values = Vec::new();
    for (order, digits) in orders
        .iter()
        .copied()
        .map(|n| (n, 40))
        .chain(std::iter::once((*orders.last().unwrap(), 60)))
    {
        let started = Instant::now();
        let contributions = reference(order, digits);
        let wall_seconds = started.elapsed().as_secs_f64();
        eprintln!("reference order={order}, digits={digits}, wall={wall_seconds:.3}s");
        let snapshot = serde_json::json!({"order":order,"digits":digits,"wall_seconds":wall_seconds,
            "contributions":contributions.iter().map(|(k,v)|(k.clone(),v.re.to_string())).collect::<BTreeMap<_,_>>()});
        std::fs::write(
            directory.join(format!("quadrature-{order}-{digits}.json")),
            serde_json::to_vec_pretty(&snapshot).unwrap(),
        )
        .unwrap();
        snapshots.push(snapshot);
        values.push(contributions);
    }
    let p = Precision::decimal(80).unwrap();
    let baseline = &values[values.len() - 3];
    let refined = &values[values.len() - 2];
    let high_precision = &values[values.len() - 1];
    let mut comparisons = BTreeMap::new();
    let mut passed = true;
    for (name, value) in refined {
        let node_delta = p.norm(&p.sub(value, &baseline[name]));
        let precision_delta = p.norm(&p.sub(value, &high_precision[name]));
        let scale = p.norm(&high_precision[name]);
        let absolute_ok = node_delta <= p.tolerance(12) && precision_delta <= p.tolerance(12);
        let relative_ok = scale == p.real(0)
            || (node_delta <= scale.clone() * p.tolerance(10)
                && precision_delta <= scale.clone() * p.tolerance(10));
        passed &= absolute_ok && relative_ok;
        comparisons.insert(name, serde_json::json!({"node_absolute_change":node_delta.to_string(),
            "precision_absolute_change":precision_delta.to_string(),"passed":absolute_ok && relative_ok}));
    }
    let report = serde_json::json!({"schema":1,"status":if passed {"refinement_passed"}else{"refinement_failed"},
        "definition":input,"dimension":"12/5","epsilon":"4/5","normalization":"unscaled Euclidean amplitude",
        "method":"independent Schwinger/Feynman parameter and compact radial/angular quadrature; six vacuum sectors; exact fifth-power endpoint maps; native MP Gauss-Legendre",
        "production_integral_owners_used":[],"native_arithmetic_owner":"symbolica_amflow::Precision",
        "oracle_records_read":0,"feature_predictions_read":0,"feature_comparisons":0,
        "source_blake3":blake3::hash(include_bytes!("finite_density_reference.rs")).to_hex().to_string(),
        "cargo_lock_blake3":blake3::hash(include_bytes!("../Cargo.lock")).to_hex().to_string(),
        "error_estimate":"independent node and precision changes, not a rigorous interval bound",
        "snapshots":snapshots,"comparisons":comparisons});
    std::fs::write(
        directory.join("independent-reference.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    assert!(
        passed,
        "independent reference failed node/precision refinement; see {}",
        directory.display()
    );
}
