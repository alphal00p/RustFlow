//! Published input families. These constructors contain no numerical answers.
use crate::{Integral, IntegralFamily, Propagator, Result};
use symbolica::prelude::*;

/// The four rank-three two-loop targets in arXiv:2201.11669, section 4.
pub fn paper_two_loop() -> Result<(IntegralFamily, Vec<Integral>)> {
    let s = Atom::num(30);
    let t = Atom::num((-10, 3));
    let m = Atom::num(1);
    let a = s.clone() / Atom::num(2);
    let b = (&t - &m) / Atom::num(2);
    let c = (-&s - &t + &m) / Atom::num(2);
    let gram = vec![
        vec![Atom::new(), a.clone(), b.clone()],
        vec![a, Atom::new(), c.clone()],
        vec![b, c, m.clone()],
    ];
    let routing = [
        ([1, 0], [0, 0, 0], 0),
        ([1, 0], [1, 0, 0], 0),
        ([1, 0], [1, 1, 0], 0),
        ([0, 1], [0, 0, 0], 0),
        ([0, 1], [0, 0, 1], 1),
        ([0, 1], [-1, -1, 0], 0),
        ([1, 1], [0, 0, 0], 0),
        ([1, 0], [0, 0, -1], 0),
        ([0, 1], [1, 0, 0], 0),
    ];
    let propagators = routing
        .iter()
        .map(|(l, e, m)| Propagator::quadratic(l, e, Atom::num(*m), &gram))
        .collect::<Result<_>>()?;
    let family = IntegralFamily {
        name: "paper_tt".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec!["p1".into(), "p2".into(), "p3".into()],
        external_gram: gram,
        propagators,
        physical_propagators: 7,
        epsilon: symbol!("symbolica_amflow::eps"),
        dimension: 4,
    };
    let targets = (0..=3)
        .map(|i| Integral(vec![1, 1, 1, 1, 1, 1, 1, -3 + i, -i]))
        .collect();
    Ok((family, targets))
}

#[derive(serde::Deserialize)]
struct ReferenceComponent {
    re: String,
    im: String,
    re_precision: Option<String>,
    im_precision: Option<String>,
}
#[derive(serde::Deserialize)]
struct ReferenceTarget {
    indices: Vec<i16>,
    coefficients: std::collections::BTreeMap<i32, ReferenceComponent>,
}
#[derive(serde::Deserialize)]
struct ReferenceData {
    upstream_commit: String,
    targets: Vec<ReferenceTarget>,
}

/// Compare against the pinned upstream data at its recorded component precision.
/// Requires independent 20-digit stability of the calculated coefficients too.
/// Fixture mantissa length is never used as an accuracy estimate.
pub fn validate_paper_two_loop(
    targets: &[Integral],
    values: &[crate::LaurentExpansion],
) -> Result<()> {
    let reference: ReferenceData =
        serde_json::from_str(include_str!("../fixtures/amflow-2.0/two_loop.json"))
            .map_err(|e| crate::Error::InvalidInput(e.to_string()))?;
    if reference.upstream_commit != "26005517a288086c4cb4d1b26d829691bc088485"
        || targets.len() != 4
        || values.len() != 4
    {
        return Err(crate::Error::InvalidInput(
            "two-loop reference provenance/dimensions".into(),
        ));
    }
    let p = crate::Precision::decimal(80)?;
    for (i, ((target, value), expected)) in targets
        .iter()
        .zip(values)
        .zip(reference.targets)
        .enumerate()
    {
        if target.0 != expected.indices {
            return Err(crate::Error::InvalidInput(
                "two-loop target ordering".into(),
            ));
        }
        if value.verified_digits.unwrap_or(0) < 20 {
            return Err(crate::Error::Accuracy(format!(
                "target {i} lacks independent 20-digit stability"
            )));
        }
        for (power, coefficient) in expected.coefficients {
            let actual = value
                .coefficients
                .get(&power)
                .ok_or_else(|| crate::Error::Accuracy(format!("missing epsilon power {power}")))?;
            let change = value.comparison_errors.get(&power).ok_or_else(|| {
                crate::Error::Accuracy("missing independent coefficient comparison".into())
            })?;
            let magnitude = p.norm(actual);
            let scale = if magnitude <= p.tolerance(20) {
                p.real(1)
            } else {
                magnitude
            };
            if !change.is_finite() || *change < p.real(0) || *change > p.tolerance(20) * scale {
                return Err(crate::Error::Accuracy(format!(
                    "target {i}, epsilon^{power} lacks relative 20-digit stability"
                )));
            }
            for (name, actual, text, precision) in [
                ("real", &actual.re, coefficient.re, coefficient.re_precision),
                (
                    "imaginary",
                    &actual.im,
                    coefficient.im,
                    coefficient.im_precision,
                ),
            ] {
                let expected = p.parse(&text, "0")?;
                let actual = crate::ComplexFloat::new(actual.clone(), p.real(0));
                let tolerance = if let Some(precision) = precision {
                    // A five-unit allowance at the recorded precision covers
                    // the reference's last reported digit and rounding.
                    p.norm(&expected) * p.pow(&p.i(10), &p.neg(&p.parse(&precision, "0")?)).re * 5
                } else {
                    p.tolerance(20)
                };
                if !p.finite(&actual) || p.norm(&p.sub(&actual, &expected)) > tolerance {
                    return Err(crate::Error::Accuracy(format!(
                        "reference mismatch: target {i}, epsilon^{power}, {name} component"
                    )));
                }
            }
        }
    }
    Ok(())
}
