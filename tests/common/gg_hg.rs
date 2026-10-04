//! Offline loader for independently serialized mixed QCD–EW Higgs-plus-jet data.
use serde_json::Value;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};
static TEMPORARY_SEQUENCE: AtomicU64 = AtomicU64::new(0);
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{CanonicalAlgebraicSystem, SquareRoot};
use symbolica_amflow::physical_transport::RustFlow;
use symbolica_amflow::transport_cache::{
    BoundaryAccuracy, CachedBoundary, CachedPoint, EpsilonRange, PointKind, RootGerm, RootSheet,
    RustFlowCache, ScaledDistance,
};
use symbolica_amflow::{
    ComplexFloat, Error, FlowOptions, Precision, Prescription, Result, RunContext,
};

pub fn fixture() -> Value {
    serde_json::from_str(include_str!("../../fixtures/gg-hg/full-systems.json")).unwrap()
}
fn atom(s: &str) -> Result<Atom> {
    Atom::parse(s, "gg_hg_offline", Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(s: &str) -> Symbol {
    if let AtomView::Var(v) = atom(s).unwrap().as_view() {
        v.get_symbol()
    } else {
        panic!("fixture identifier must be a symbol")
    }
}
fn numbers(p: Precision, x: &Value) -> Result<Vec<Vec<ComplexFloat>>> {
    x.as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|z| {
                    p.parse(
                        z["real"].as_str().unwrap(),
                        z["imaginary"].as_str().unwrap(),
                    )
                })
                .collect()
        })
        .collect()
}
fn coordinates(f: &Value, name: &str) -> Result<BTreeMap<Symbol, Atom>> {
    ["s", "t", "b"]
        .into_iter()
        .zip(f[name].as_array().unwrap())
        .map(|(s, a)| Ok((symbol(s), atom(a.as_str().unwrap())?)))
        .collect()
}
fn system(f: &Value) -> Result<CanonicalAlgebraicSystem> {
    let n = f["dimension"].as_u64().unwrap() as usize;
    let letters = f["letters"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| atom(a.as_str().unwrap()))
        .collect::<Result<Vec<_>>>()?;
    let mut matrices = vec![vec![vec![Atom::new(); n]; n]; letters.len()];
    for row in f["matrix_entries"].as_array().unwrap() {
        let letter = row["letter"].as_u64().unwrap() as usize - 1;
        for e in row["entries"].as_array().unwrap() {
            matrices[letter][e[0].as_u64().unwrap() as usize - 1]
                [e[1].as_u64().unwrap() as usize - 1] = atom(e[2].as_str().unwrap())?;
        }
    }
    let roots = f["roots"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| {
            Ok(SquareRoot {
                symbol: symbol(r["name"].as_str().unwrap()),
                radicand: atom(r["radicand"].as_str().unwrap())?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    CanonicalAlgebraicSystem::new(
        symbol("eps"),
        &[symbol("s"), symbol("t"), symbol("b")],
        &letters,
        &matrices,
        roots,
    )
}

/// Evaluate a supplied canonical family with its recorded physical root germ.
/// Source precision is capped by the inherited integration error, not its mantissa.
pub fn evaluate_case(
    fixture: &Value,
    index: usize,
    guard: u32,
    order: usize,
) -> Result<CachedBoundary> {
    let case = &fixture["cases"][index];
    let f = fixture["systems"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["family"] == case["family"])
        .unwrap();
    let n = f["dimension"].as_u64().unwrap() as usize;
    let name = case["label"].as_str().unwrap().replace('-', "_");
    let basis = (1..=n)
        .map(|j| atom(&format!("{name}_I{j}")))
        .collect::<Result<Vec<_>>>()?;
    let flow = RustFlow::new_canonical(
        system(f)?,
        &basis,
        &Atom::one(),
        Prescription::PlusI0,
        "Supplied mixed QCD-EW gg-to-Hg grid, exact regular chart, explicit original root germ",
    )?;
    let germ = RootGerm {
        sheets: f["roots"]
            .as_array()
            .unwrap()
            .iter()
            .zip(case["source_root_germs"].as_array().unwrap())
            .map(|(r, g)| {
                (
                    symbol(r["name"].as_str().unwrap()),
                    match g.as_str().unwrap() {
                        "principal" => RootSheet::Principal,
                        "opposite_principal" => RootSheet::Opposite,
                        _ => panic!("unknown root germ"),
                    },
                )
            })
            .collect(),
    };
    let p = Precision::decimal(140)?;
    let range = EpsilonRange::new(0, 4)?;
    let source = CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(coordinates(case, "start")?).with_root_germ(germ.clone())?,
        kind: PointKind::Physical,
        range,
        coefficients: numbers(p, &case["source_values"])?,
        accuracy: BoundaryAccuracy::supplied(
            24,
            p.bits,
            case["source_errors"]
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|x| p.parse(x.as_str().unwrap(), "0").map(|z| z.re))
                        .collect()
                })
                .collect::<Result<Vec<_>>>()?,
            "Inherited original delta rounded upward plus decimal serialization reserve; cap24",
        )?,
    };
    let mut bank = RustFlowCache::default();
    bank.insert(source)?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: guard,
        series_order: order,
        ..Default::default()
    };
    let destination = coordinates(case, "destination")?;
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = flow.evaluate_to(
        &mut bank,
        &destination,
        Some(&germ),
        range,
        &options,
        &RunContext::default(),
        &policy,
    )?;
    assert!(result.boundary.accuracy.verified_digits() >= 20);
    assert!(result.boundary.accuracy.verified_digits() <= 24);
    let reference = numbers(p, &case["reference_values"])?;
    for (a, b) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(reference.iter().flatten())
    {
        assert!(
            p.norm(&p.sub(a, b)) < p.tolerance(24),
            "{} original mismatch",
            case["label"]
        );
    }
    let path = std::env::temp_dir().join(format!(
        "rustflow-gg-hg-{}-{}-{name}-{guard}-{order}.bin",
        std::process::id(),
        TEMPORARY_SEQUENCE.fetch_add(1, Ordering::Relaxed)
    ));
    bank.save(&path)?;
    let mut restored = RustFlowCache::load(&path)?;
    std::fs::remove_dir_all(path)?;
    let hit = flow.evaluate_to(
        &mut restored,
        &destination,
        Some(&germ),
        range,
        &options,
        &RunContext::default(),
        &policy,
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.inserted_points, 0);
    assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
    Ok(result.boundary)
}

fn scalar(p: Precision, z: &Value) -> Result<ComplexFloat> {
    p.parse(
        z["real"].as_str().unwrap(),
        z["imaginary"].as_str().unwrap(),
    )
}
fn coefficient(p: Precision, text: &str) -> Result<ComplexFloat> {
    let value = atom(text)?
        .replace(atom("imaginary_unit")?)
        .with(Atom::num(Complex::new(
            Rational::from(0),
            Rational::from(1),
        )));
    p.eval(&value, &Default::default())
}
fn project<'a>(
    p: Precision,
    weights: &Value,
    coefficients: &mut BTreeMap<String, ComplexFloat>,
    lookup: impl Fn(&Value) -> &'a CachedBoundary,
) -> Result<(ComplexFloat, Float, u32)> {
    let mut value = p.zero();
    let mut error = p.real(0);
    let mut magnitude = p.real(0);
    let mut cap = u32::MAX;
    for weight in weights.as_array().unwrap() {
        let boundary = lookup(weight);
        let epsilon = weight["epsilon"].as_u64().unwrap() as usize;
        let index = weight["index"].as_u64().unwrap() as usize - 1;
        let text = weight["coefficient"].as_str().unwrap();
        if !coefficients.contains_key(text) {
            coefficients.insert(text.to_owned(), coefficient(p, text)?);
        }
        let c = &coefficients[text];
        if weight.get("coefficient_value").is_some() {
            assert!(p.close(c, &scalar(p, &weight["coefficient_value"])?, 90));
        }
        let term = p.mul(c, &boundary.coefficients[epsilon][index]);
        value = p.add(&value, &term);
        error += p.norm(c) * &boundary.accuracy.comparison_errors()[epsilon][index];
        magnitude += p.norm(&term);
        cap = cap.min(boundary.accuracy.verified_digits());
    }
    // This exceeds rounding at the160-decimal-digit projection precision.
    error += p.tolerance(130) * magnitude;
    Ok((value, error, cap))
}

/// Exact FF and supplied diagnostic projections, with coefficient-weighted
/// transported uncertainty. The source cap does not imply significant digits
/// for small values; relative FF accuracy is checked separately.
pub fn check_projections(fixture: &Value, boundaries: &[CachedBoundary]) -> Result<()> {
    let projections: Value =
        serde_json::from_str(include_str!("../../fixtures/gg-hg/projections.json")).unwrap();
    let cases = fixture["cases"].as_array().unwrap();
    assert_eq!(cases.len(), boundaries.len());
    let data = cases
        .iter()
        .zip(boundaries)
        .map(|(case, boundary)| (case["label"].as_str().unwrap(), boundary))
        .collect::<BTreeMap<_, _>>();
    let p = Precision::decimal(160)?;
    let mut coefficients = BTreeMap::new();
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let loop_factor = p.neg(&p.div(&p.i(1), &p.powi(&p.scale(&pi, 4, 1), 4)));
    let mut ff_count = 0;
    for mass in projections["form_factors"].as_array().unwrap() {
        let name = mass["mass"].as_str().unwrap();
        let normalization = p.mul(
            &loop_factor,
            &p.powi(&coefficient(p, mass["point"][2].as_str().unwrap())?, 2),
        );
        for row in mass["form_factors"].as_array().unwrap() {
            let (value, error, cap) = project(p, &row["weights"], &mut coefficients, |w| {
                data[format!(
                    "{name}-{}-{}",
                    w["family"].as_str().unwrap(),
                    w["permutation"].as_u64().unwrap()
                )
                .as_str()]
            })?;
            let reference = scalar(p, &row["original_exact_i0_limit"])?;
            assert!(p.norm(&p.sub(&value, &reference)) <= error.clone() + p.tolerance(70));
            let normalized = p.mul(&normalization, &value);
            let normalized_error = p.norm(&normalization) * error;
            let digits = (0..=cap)
                .rev()
                .find(|&digits| normalized_error <= p.tolerance(digits) * p.norm(&normalized))
                .unwrap_or(0);
            // The current24-digit source contract supports19–20 relative FF
            // digits; the uniform20-significant-digit gap is documented.
            assert!(
                digits >= 19,
                "{name} form factor{} relative accuracy",
                row["form_factor"]
            );
            eprintln!(
                "{name} F{}: {digits} relative digits; |error|={normalized_error}",
                row["form_factor"]
            );
            ff_count += 1;
        }
    }
    assert_eq!(ff_count, 8);
    let mut row_count = 0;
    for case in projections["precanonical"].as_array().unwrap() {
        let boundary = data[case["label"].as_str().unwrap()];
        for row in case["rows"].as_array().unwrap() {
            let (value, error, _) = project(p, &row["weights"], &mut coefficients, |_| boundary)?;
            let reference = scalar(p, &row["original_value"])?;
            assert!(p.norm(&p.sub(&value, &reference)) <= error + p.tolerance(70));
            row_count += 1;
        }
    }
    assert_eq!(row_count, 4376);
    Ok(())
}
