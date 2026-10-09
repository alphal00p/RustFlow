use ahash::HashMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::boundary::OccupiedBoundaryLimits;
use symbolica_amflow::finite_density::flow_boundary::OccupiedFlowBoundary;
use symbolica_amflow::finite_density::geometry::OccupiedCutFamily;
use symbolica_amflow::frobenius::{FrobeniusBasis, FrobeniusColumn};
use symbolica_amflow::*;

fn family(cuts: &[usize]) -> OccupiedCutFamily {
    let input: DensityInput = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/finite_density/massive_two_loop_sunset.json"
    )))
    .unwrap();
    input
        .prepare()
        .unwrap()
        .occupied_cut(cuts, 16)
        .unwrap()
        .at_physical_masses()
}

fn shift(family: &OccupiedCutFamily) -> Vec<bool> {
    (0..family.factors().len())
        .map(|slot| {
            slot < family.physical_slots()
                && !family
                    .shells()
                    .iter()
                    .any(|shell| shell.physical_slot == slot)
        })
        .collect()
}

fn reference(expression: &str, p: Precision) -> ComplexFloat {
    p.eval(&parse!(expression), &HashMap::default()).unwrap()
}

#[test]
fn occupied_matching_plans_and_integrates_subleading_coefficients() {
    let family = family(&[0, 1]);
    let original = family.targets()[0].keys().next().unwrap().clone();
    let mut pinched = original.clone();
    pinched.0[2] = 0;
    let basis = [pinched, original];
    let p = Precision::decimal(60).unwrap();
    // This deliberately simple test-only basis isolates the matching owner:
    // [J,I] = c0*[1,-z] + c1*[0,z^2] + O(z^3).
    // Its second constant cannot be fixed from the leading region terms.
    // The actual boundary coefficients still come from the physical sunset
    // integrand, graph routing, compact moments and spatial tensor projection.
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![
            FrobeniusColumn {
                exponent: Atom::zero(),
                coefficients: vec![
                    vec![vec![p.i(1), p.zero()]],
                    vec![vec![p.zero(), p.i(-1)]],
                    vec![vec![p.zero(), p.zero()]],
                ],
            },
            FrobeniusColumn {
                exponent: Atom::num(2),
                coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
            },
        ],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let eps = symbol!("occupied_flow_boundary_eps");
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        2,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let report = boundary
        .constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    assert!(p.close(&report.constants[0], &reference("1/(4*pi)", p), 40));
    assert!(p.close(&report.constants[1], &reference("13/(32*pi)", p), 40));
    assert_eq!(report.provenance.regions.len(), 1);
    assert_eq!(
        report.provenance.regions[0].exponents,
        vec![Atom::zero(), Atom::one()]
    );
    assert_eq!(
        report.provenance.regions[0].half_orders,
        vec![Some(0), Some(2)]
    );
    assert_eq!(report.provenance.integrated_coefficients, 4);
    assert_eq!(report.provenance.scaleless_noncompact_products, 0);

    let shallow = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        1,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    assert!(matches!(
        shallow.constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::IncompleteReduction(_))
    ));
    let mut short_series = solutions.clone();
    short_series.columns[0].coefficients.pop();
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &short_series,
        ),
        Err(Error::IncompleteReduction(_))
    ));
}

#[test]
fn mixed_boundary_keeps_symbolic_epsilon_power_and_deformation_roles() {
    let family = family(&[0]);
    let basis = vec![family.targets()[0].keys().next().unwrap().clone()];
    let p = Precision::decimal(60).unwrap();
    let eps = symbol!("occupied_mixed_flow_boundary_eps");
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![FrobeniusColumn {
            exponent: Atom::var(eps),
            coefficients: vec![vec![vec![p.i(1)]]],
        }],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        2,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let shifted = shift(&family);
    let report = boundary
        .constants(
            &family,
            &basis,
            &shifted,
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    assert!(p.close(&report.constants[0], &p.scale(&p.i(1), 1, 2), 40));
    assert!(report.provenance.regions.iter().any(|region| {
        region.exponents[0] == Atom::var(eps) && region.half_orders[0] == Some(0)
    }));
    assert!(
        report
            .provenance
            .regions
            .iter()
            .all(|region| !region.hard[0])
    );
    let mut moving_shell = shifted.clone();
    moving_shell[0] = true;
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &moving_shell,
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::Unsupported(_))
    ));
    let mut moving_occupation = shifted;
    moving_occupation[family.shells()[0].upper_slot] = true;
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &moving_occupation,
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::Unsupported(_))
    ));
}
