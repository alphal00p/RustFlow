//! Complete 15-master unequal banana, seeded by the analytic equal-mass limit.
#[path = "common/banana.rs"]
mod banana;
use banana::{atom, evaluate_equal};
use symbolica::prelude::*;
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::*;

fn transport_line(
    p: Precision,
    input: EpsilonBoundary,
    line: &serde_json::Value,
    options: &FlowOptions,
) -> Result<EpsilonBoundary> {
    let matrices = line["matrix_by_epsilon"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            m.as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|a| atom(a.as_str().unwrap()))
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let system = EpsilonSystem {
        variable: symbol!("banana_proto::x"),
        matrices,
    };
    let compiled = system.compile(p, &Default::default())?;
    let start = p.eval(&atom(line["start"].as_str().unwrap())?, &Default::default())?;
    let end = p.eval(&atom(line["end"].as_str().unwrap())?, &Default::default())?;
    let route = compiled.plan_path(&start, &end, 1)?;
    // Each exact coordinate map has its own line parameter. Coefficients at
    // the shared physical point are carried into the next parameterization.
    let answer = compiled.transport(
        &EpsilonBoundary {
            point: start,
            ..input
        },
        &route,
        options,
        &RunContext::default(),
        false,
    )?;
    Ok(EpsilonBoundary {
        point: answer.point,
        leading: answer.leading,
        coefficients: answer.coefficients,
    })
}

fn evaluate(digits: u32, order: usize, shortcut: bool) -> Result<EpsilonBoundary> {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../fixtures/diffexp/banana-unequal.json")).unwrap();
    let p = Precision::decimal(digits)?;
    let endpoint = if shortcut {
        Rational::from(50)
    } else {
        Rational::from((1, 2))
    };
    let (_, equal) = evaluate_equal(digits, order, endpoint)?;
    // Exact equal-mass symmetry: six double-raised integrals, four single-
    // raised integrals, the scalar banana, and four factorized tadpoles.
    let mapping = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 3, 3, 3, 3];
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: equal
            .coefficients
            .iter()
            .map(|row| mapping.iter().map(|&i| row[i].clone()).collect())
            .collect(),
    };
    let options = FlowOptions {
        digits: 20,
        guard_digits: digits - 20,
        series_order: order,
        ..Default::default()
    };
    let line = if shortcut {
        "mass_deformation_50"
    } else {
        "mass_deformation_half"
    };
    let boundary = transport_line(p, boundary, &fixture["lines"][line], &options)?;
    if shortcut {
        Ok(boundary)
    } else {
        transport_line(p, boundary, &fixture["lines"]["physical_psq"], &options)
    }
}

#[test]
fn unequal_banana_full_route_refinement_and_independent_contour() -> Result<()> {
    let low = evaluate(60, 80, false)?;
    let high = evaluate(80, 112, false)?;
    let shortcut = evaluate(80, 112, true)?;
    let p = Precision::decimal(90)?;
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../fixtures/diffexp/banana-unequal.json")).unwrap();
    for (k, row) in high.coefficients.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert!(
                p.norm(&p.sub(value, &low.coefficients[k][j])) < p.tolerance(40),
                "precision/order refinement eps{k} master{j}"
            );
            assert!(
                p.norm(&p.sub(value, &shortcut.coefficients[k][j])) < p.tolerance(40),
                "independent contour eps{k} master{j}"
            );
            let original = &fixture["reference"]["original_endpoint_by_epsilon"][k][j];
            let expected = p.parse(
                original["real"].as_str().unwrap(),
                original["imaginary"].as_str().unwrap(),
            )?;
            assert!(
                p.norm(&p.sub(value, &expected)) < p.tolerance(20),
                "original DiffExp eps{k} master{j}: {value} != {expected}"
            );
        }
    }
    // Four factorized sectors are analytic products of three tadpoles:
    // exp(3*EulerGamma*eps)*Gamma(1+eps)^3 * Q^(-eps).
    let pi = ComplexFloat::new(p.log(&p.i(-1)).im, p.real(0));
    let zeta3 = ComplexFloat::new(
        Float::with_val(p.bits, p.real(3).as_raw().clone().zeta()),
        p.real(0),
    );
    let gamma = [
        p.i(1),
        p.zero(),
        p.scale(&p.powi(&pi, 2), 1, 4),
        p.neg(&zeta3),
        p.scale(&p.powi(&pi, 4), 19, 480),
    ];
    for (j, q) in [(2, 1), (8, 3), (3, 1), (4, 1)].into_iter().enumerate() {
        let logarithm = p.log(&p.rational(&Rational::from(q)));
        let mut exponential = vec![p.i(1)];
        for k in 1..5 {
            exponential.push(p.scale(&p.mul(&exponential[k - 1], &p.neg(&logarithm)), 1, k as i64));
        }
        for (k, row) in high.coefficients.iter().enumerate() {
            let expected = (0..=k).fold(p.zero(), |sum, l| {
                p.add(&sum, &p.mul(&gamma[l], &exponential[k - l]))
            });
            assert!(
                p.norm(&p.sub(&row[j + 11], &expected)) < p.tolerance(50),
                "analytic tadpole product eps{k} master{}",
                j + 11
            );
        }
    }
    // The printed paper table is an additional independent check of B11.
    for (k, row) in high.coefficients.iter().enumerate() {
        let reference = &fixture["reference"]["paper_B11"][k];
        let expected = p.parse(
            reference["real"].as_str().unwrap(),
            reference["imaginary"].as_str().unwrap(),
        )?;
        assert!(p.norm(&p.sub(&row[10], &expected)) < p.tolerance(20));
    }
    Ok(())
}
