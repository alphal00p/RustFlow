//! Shared PH1-to-PH6 fixture transport for the pinned planar benchmarks.
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, RootSeed, SquareRoot};
use symbolica_amflow::{
    ComplexFloat, Error, FlowOptions, Precision, Prescription, Progress, Result, RunContext,
    diffexp, kinematics,
};
fn atom(s: &str) -> Result<Atom> {
    Atom::parse(s, "fivepoint", Default::default()).map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(s: &str) -> Symbol {
    if let AtomView::Var(v) = atom(s).unwrap().as_view() {
        v.get_symbol()
    } else {
        panic!("symbol")
    }
}
fn scalar(p: Precision, v: &Value) -> Result<ComplexFloat> {
    p.parse(
        v["real"].as_str().unwrap(),
        v["imaginary"].as_str().unwrap(),
    )
}
fn maximum_error(p: Precision, values: &[Vec<ComplexFloat>], reference: &Value) -> Result<Float> {
    let mut maximum = p.real(0);
    for (k, row) in values.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let e = p.norm(&p.sub(v, &scalar(p, &reference[k][j])?));
            if e > maximum {
                maximum = e;
            }
        }
    }
    Ok(maximum)
}
pub fn evaluate(
    fixture: &Value,
    root_names: &[&str],
    digits: u32,
    order: usize,
) -> Result<Vec<Vec<ComplexFloat>>> {
    let width_divisor = 5;
    let dimension = fixture["dimension"].as_u64().unwrap() as usize;
    let geometry = &fixture["contour_geometry"];
    let p = Precision::decimal(digits)?;
    let start_clock = Instant::now();
    let parameter = symbol("x");
    let epsilon = symbol("eps");
    let variables = fixture["coordinate_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| symbol(n.as_str().unwrap()))
        .collect::<Vec<_>>();
    let mut coordinates = BTreeMap::new();
    for name in fixture["coordinate_names"].as_array().unwrap() {
        let name = name.as_str().unwrap();
        let a = atom(fixture["source_point"][name].as_str().unwrap())?;
        let b = atom(fixture["destination_point"][name].as_str().unwrap())?;
        coordinates.insert(symbol(name), &a + (b - &a) * Atom::var(parameter));
    }
    let path = kinematics::KinematicPath {
        parameter,
        coordinates,
    };
    let used = fixture["sparse_rational_coefficients"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["indices"][0].as_u64().unwrap() as usize)
        .collect::<BTreeSet<_>>();
    let mut letters = Vec::new();
    let mut matrices = Vec::new();
    for index in used {
        letters.push(atom(fixture["letters"][index - 1].as_str().unwrap())?);
        let mut m = vec![vec![Atom::new(); dimension]; dimension];
        for entry in fixture["sparse_rational_coefficients"].as_array().unwrap() {
            let indices = entry["indices"].as_array().unwrap();
            if indices[0].as_u64().unwrap() as usize == index {
                m[indices[1].as_u64().unwrap() as usize - 1]
                    [indices[2].as_u64().unwrap() as usize - 1] =
                    atom(entry["coefficient"].as_str().unwrap())?;
            }
        }
        matrices.push(m);
    }
    let mut roots = Vec::new();
    for &name in root_names {
        let definition = fixture["root_definitions"][name].as_str().unwrap();
        let polynomial = definition
            .strip_prefix("sqrt(")
            .unwrap()
            .strip_suffix(')')
            .unwrap();
        roots.push(SquareRoot {
            symbol: symbol(name),
            radicand: atom(polynomial)?,
        });
    }
    eprintln!("construct canonical physical system, digits={digits},order={order}");
    let symbolic_clock = Instant::now();
    let system =
        AlgebraicKinematicSystem::canonical_dlog(epsilon, &variables, &letters, &matrices, roots)?;
    let pulled = system.pullback(&path, 4)?;
    let symbolic_seconds = symbolic_clock.elapsed().as_secs_f64();
    eprintln!("symbolic pullback complete in {symbolic_seconds:.3}s");
    let compile_clock = Instant::now();
    let compiled = pulled.compile(p)?;
    let compile_seconds = compile_clock.elapsed().as_secs_f64();
    eprintln!(
        "compiled in {compile_seconds:.3}s with {} singularity candidates",
        compiled.singularities().len()
    );
    let mut prescriptions = Vec::new();
    for record in geometry["records"].as_array().unwrap() {
        for zero in record["interior_simple_zeros"].as_array().unwrap() {
            let raw = zero["x"].as_str().unwrap().split('`').next().unwrap();
            prescriptions.push((
                p.parse(raw, "0")?,
                zero["required_parameter_imaginary_sign_for_polynomial_plus_i0"]
                    .as_i64()
                    .unwrap(),
                record["name"].as_str().unwrap(),
            ));
        }
    }
    let merge_tolerance = p.tolerance(25);
    let mut real_poles = compiled
        .singularities()
        .iter()
        .filter(|z| {
            z.im < merge_tolerance
                && z.im > -merge_tolerance.clone()
                && z.re > p.real(0)
                && z.re < p.real(1)
        })
        .cloned()
        .collect::<Vec<_>>();
    real_poles.sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());
    real_poles.dedup_by(|a, b| p.norm(&p.sub(a, b)) < merge_tolerance);
    let mut waypoints = Vec::new();
    let mut matched_prescriptions = BTreeSet::new();
    for root in real_poles {
        let pole = ComplexFloat::new(root.re.clone(), p.real(0));
        let mut clearance = p.real(1) / 20;
        for candidate in compiled.singularities() {
            let d = p.norm(&p.sub(&pole, candidate));
            if d > merge_tolerance && d < clearance {
                clearance = d;
            }
        }
        for endpoint in [p.zero(), p.i(1)] {
            let d = p.norm(&p.sub(&pole, &endpoint));
            if d < clearance {
                clearance = d;
            }
        }
        let radius = clearance / width_divisor;
        let matched = prescriptions
            .iter()
            .enumerate()
            .filter(|(_, (at, _, _))| p.norm(&p.sub(&pole, at)) < merge_tolerance)
            .collect::<Vec<_>>();
        let side = matched.first().map(|(_, (_, s, _))| *s).unwrap_or(1);
        for (index, (_, required_side, name)) in matched {
            assert_eq!(side, *required_side, "conflicting fixture cut {name}");
            matched_prescriptions.insert(index);
        }
        let left = ComplexFloat::new(pole.re.clone() - radius.clone(), p.real(0));
        let right = ComplexFloat::new(pole.re.clone() + radius.clone(), p.real(0));
        let high_left = ComplexFloat::new(left.re.clone(), radius.clone() * side);
        let high_right = ComplexFloat::new(right.re.clone(), radius.clone() * side);
        waypoints.extend([left, high_left, high_right, right]);
    }
    assert_eq!(
        matched_prescriptions.len(),
        prescriptions.len(),
        "every prescribed fixture zero must have a contour detour"
    );
    waypoints.push(p.i(1));
    let coefficients = fixture["boundary_by_epsilon"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| scalar(p, v))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let boundary = diffexp::EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients,
    };
    let seeds = root_names
        .iter()
        .map(|&name| (symbol(name), RootSeed::I0(Prescription::PlusI0)))
        .collect::<BTreeMap<_, _>>();
    let options = FlowOptions {
        digits: 20,
        guard_digits: digits - 20,
        series_order: order,
        max_steps: 5000,
        ..Default::default()
    };
    let context = RunContext {
        progress: Some(Arc::new(|event| {
            if let Progress::Step { index } = event
                && index % 25 == 0
            {
                eprintln!("step {index}");
            }
        })),
        ..Default::default()
    };
    eprintln!("transport {} explicit waypoints", waypoints.len());
    let transport_clock = Instant::now();
    let answer = compiled.transport(&boundary, &waypoints, &seeds, &options, &context, false)?;
    let transport_seconds = transport_clock.elapsed().as_secs_f64();
    let source_error = maximum_error(
        p,
        &answer.solution.coefficients,
        &fixture["reference_by_epsilon"],
    )?;
    let original_error = maximum_error(
        p,
        &answer.solution.coefficients,
        &fixture["original_by_epsilon"],
    )?;
    let pass = source_error < p.tolerance(20) && original_error < p.tolerance(20);
    eprintln!(
        "PH1 to PH6: symbolic {symbolic_seconds:.3}s, compilation {compile_seconds:.3}s, transport {transport_seconds:.3}s, total {:.3}s, {} accepted/{} rejected",
        start_clock.elapsed().as_secs_f64(),
        answer.solution.diagnostics.steps,
        answer.solution.diagnostics.rejected_steps
    );
    if !pass {
        return Err(Error::Accuracy(format!(
            "native endpoint disagreement: source={source_error}, original={original_error}"
        )));
    }
    let half = p.rational(&Rational::from((1, 2)));
    for (root, value) in &answer.branches.roots {
        let name = root_names
            .iter()
            .copied()
            .find(|name| symbol(name) == *root)
            .unwrap();
        let record = geometry["records"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["name"].as_str() == Some(name))
            .unwrap();
        let radicand = p.eval(
            &atom(record["end_value"].as_str().unwrap())?,
            &Default::default(),
        )?;
        let expected = p.pow(&radicand, &half);
        assert!(
            p.norm(&p.sub(value, &expected)) < p.tolerance(30),
            "physical root sheet {}",
            Atom::var(*root)
        );
    }
    Ok(answer.solution.coefficients)
}
