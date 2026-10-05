use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{diffexp::EpsilonBoundary, kinematics::KinematicSystem, *};

fn system(matrix: Vec<Vec<Atom>>) -> KinematicSystem {
    KinematicSystem {
        epsilon: symbol!("shear::eps"),
        derivatives: BTreeMap::from([(symbol!("shear::s"), matrix)]),
    }
}
fn parse(text: &str) -> Atom {
    Atom::parse(text, "shear", Default::default()).unwrap()
}
fn zero(n: usize) -> Vec<Vec<Atom>> {
    vec![vec![Atom::new(); n]; n]
}

#[test]
fn exact_common_shearing_preserves_both_physical_derivatives_and_gaussian_coefficients() {
    let mut matrix = zero(3);
    matrix[0][1] = parse("1i/eps^2");
    matrix[1][2] = parse("1/(eps*(1-eps))");
    let mut source = system(matrix);
    let mut second = zero(3);
    second[0][2] = parse("3i/eps^2");
    source.derivatives.insert(symbol!("shear::t"), second);
    let (regular, shear) = EpsilonShearing::regularize(&source, &RunContext::default()).unwrap();
    assert_eq!(shear.weights(), &[-3, -1, 0]);
    for (&coordinate, old) in &source.derivatives {
        let expected = DifferentialSystem {
            variable: coordinate,
            matrix: old.clone(),
        }
        .change_basis(&shear.transformation())
        .unwrap();
        assert_eq!(expected.matrix, regular.derivatives[&coordinate]);
        diffexp::EpsilonSystem::from_differential_system(&expected, source.epsilon, 5).unwrap();
    }
    // Opposite physical partials may cancel on a special path, but each still
    // constrains the common gauge before that path is chosen.
    let mut opposite = zero(2);
    opposite[0][1] = parse("-1/eps");
    let mut a = zero(2);
    a[0][1] = parse("1/eps");
    let mut source = system(a);
    source.derivatives.insert(symbol!("shear::t"), opposite);
    assert_eq!(
        EpsilonShearing::regularize(&source, &RunContext::default())
            .unwrap()
            .1
            .weights(),
        &[-1, 0]
    );
}

#[test]
fn negative_cycles_are_typed_and_zero_sum_cycles_are_admitted() {
    let source = system(vec![
        vec![Atom::new(), parse("1/eps")],
        vec![parse("eps"), Atom::new()],
    ]);
    assert_eq!(
        EpsilonShearing::regularize(&source, &RunContext::default())
            .unwrap()
            .1
            .weights(),
        &[-1, 0]
    );
    for matrix in [
        vec![vec![parse("1/eps")]],
        vec![
            vec![Atom::new(), parse("1/eps")],
            vec![parse("1/eps"), Atom::new()],
        ],
        vec![
            vec![parse("1/eps"), parse("1/eps")],
            vec![parse("-1/eps"), parse("-1/eps")],
        ],
    ] {
        assert!(
            matches!(EpsilonShearing::regularize(&system(matrix),&RunContext::default()),Err(Error::Unsupported(message)) if message.contains("negative constraint cycle"))
        );
    }
}

#[test]
fn original_denominator_holes_and_generic_epsilon_domains_survive() {
    let mut a = zero(2);
    a[0][1] = parse("(s^2-1)/(eps*(s-1)) + 1/(s-eps)");
    let source = system(a);
    let (_, shear) = EpsilonShearing::regularize(&source, &RunContext::default()).unwrap();
    let guards = shear.nonzero_conditions();
    assert!(guards.iter().any(|g| (g - parse("s-1")).expand().is_zero()));
    assert!(guards.iter().any(|g| (g - parse("s")).expand().is_zero()));
    assert!(
        guards
            .iter()
            .any(|g| (g - parse("s-eps")).expand().is_zero())
    );
    let labels = [symbol!("shear::I").call(1), symbol!("shear::I").call(2)];
    let scaled = shear.scaled_basis_labels(&labels).unwrap();
    assert_eq!(scaled[0], parse("eps") * &labels[0]);
    assert_eq!(scaled[1], labels[1]);
    let old = transport_cache::BoundaryIdentity::new(
        &source,
        &labels,
        &Atom::one(),
        Prescription::PlusI0,
        "shearing test",
    )
    .unwrap();
    let (regular, _) = EpsilonShearing::regularize(&source, &RunContext::default()).unwrap();
    let new = transport_cache::BoundaryIdentity::with_conditions(
        &regular,
        &scaled,
        &Atom::one(),
        Prescription::PlusI0,
        "shearing test",
        guards,
    )
    .unwrap();
    assert_ne!(old.key(), new.key());
}

#[test]
fn mapping_preserves_errors_and_requires_extra_positive_orders() {
    let a = system(vec![
        vec![Atom::new(), parse("1/eps")],
        vec![parse("eps"), Atom::new()],
    ]);
    let (_, shear) = EpsilonShearing::regularize(&a, &RunContext::default()).unwrap();
    assert_eq!(shear.required_sheared_range(-1, 0).unwrap(), (-1, 1));
    let p = Precision::decimal(60).unwrap();
    let original = EpsilonBoundary {
        point: p.i(0),
        leading: -1,
        coefficients: vec![
            vec![p.i(3), p.i(5)],
            vec![p.i(7), p.i(11)],
            vec![p.i(13), p.i(17)],
        ],
    };
    let errors = vec![
        vec![p.real(1) / 100, p.real(2) / 100],
        vec![p.real(3) / 100, p.real(4) / 100],
        vec![p.real(5) / 100, p.real(6) / 100],
    ];
    let (shifted, shifted_errors) = shear.to_sheared(&original, &errors, 1).unwrap();
    assert_eq!(shifted.coefficients[0], vec![p.zero(), p.i(5)]);
    assert_eq!(shifted_errors[0][0], p.real(0));
    assert_eq!(shifted.coefficients[1][0], original.coefficients[0][0]);
    assert_eq!(shifted_errors[1][0], errors[0][0]);
    let (restored, restored_errors) = shear.to_original(&shifted, &shifted_errors, -1, 0).unwrap();
    assert_eq!(restored.coefficients, original.coefficients[..2]);
    assert_eq!(restored_errors, errors[..2]);
    assert!(matches!(
        shear.to_original(&shifted, &shifted_errors, -1, 1),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.to_sheared(&original, &errors, 2),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.required_sheared_range(i32::MAX, i32::MAX),
        Err(Error::Limit(_))
    ));
    assert!(matches!(
        shear.required_sheared_range(-1023, 0),
        Err(Error::Limit(_))
    ));
}

#[test]
fn shifted_transport_restores_pole_and_finite_terms_analytically() {
    let source = system(vec![
        vec![Atom::new(), parse("1/eps")],
        vec![parse("eps"), Atom::new()],
    ]);
    let (regular, shear) = EpsilonShearing::regularize(&source, &RunContext::default()).unwrap();
    let p = Precision::decimal(80).unwrap();
    // I(0)=[0,1+eps]; I_1(x)=(eps^-1+1)*sinh(x).
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: -1,
        coefficients: vec![
            vec![p.zero(), p.zero()],
            vec![p.zero(), p.i(1)],
            vec![p.zero(), p.i(1)],
        ],
    };
    let errors = vec![vec![p.real(0); 2]; 3];
    let (initial, _) = shear.to_sheared(&boundary, &errors, 1).unwrap();
    let differential = DifferentialSystem {
        variable: symbol!("shear::s"),
        matrix: regular.derivatives[&symbol!("shear::s")].clone(),
    };
    let expanded =
        diffexp::EpsilonSystem::from_differential_system(&differential, source.epsilon, 2).unwrap();
    let options = FlowOptions {
        digits: 20,
        guard_digits: 60,
        series_order: 80,
        ..Default::default()
    };
    let compiled = expanded.compile(p, &Default::default()).unwrap();
    let result = compiled
        .transport(&initial, &[p.i(1)], &options, &RunContext::default(), false)
        .unwrap();
    let at_end = EpsilonBoundary {
        point: result.point,
        leading: result.leading,
        coefficients: result.coefficients,
    };
    let (restored, _) = shear.to_original(&at_end, &errors, -1, 0).unwrap();
    let ex = p.exp(&p.i(1));
    let em = p.exp(&p.i(-1));
    let sinh = p.scale(&p.sub(&ex, &em), 1, 2);
    let cosh = p.scale(&p.add(&ex, &em), 1, 2);
    assert!(p.close(&restored.coefficients[0][0], &sinh, 40));
    assert!(p.close(&restored.coefficients[1][0], &sinh, 40));
    assert!(p.close(&restored.coefficients[1][1], &cosh, 40));
}

#[test]
fn zero_identity_fractional_powers_and_cancelled_context() {
    let (regular, shear) =
        EpsilonShearing::regularize(&system(zero(3)), &RunContext::default()).unwrap();
    assert_eq!(shear.weights(), &[0, 0, 0]);
    assert_eq!(regular.derivatives.values().next().unwrap(), &zero(3));
    assert!(matches!(
        EpsilonShearing::regularize(
            &system(vec![vec![parse("eps^(1/2)")]]),
            &RunContext::default()
        ),
        Err(Error::Unsupported(_))
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        EpsilonShearing::regularize(&system(zero(1)), &context),
        Err(Error::Cancelled)
    ));
}

#[test]
fn covariance_cancellation_and_large_span_are_exact() {
    let mut a = zero(3);
    a[0][1] = parse("1i/eps^2");
    a[1][2] = parse("1/eps");
    let (regular, shear) =
        EpsilonShearing::regularize(&system(a.clone()), &RunContext::default()).unwrap();
    let permutation = [2, 0, 1];
    let permuted = permutation
        .iter()
        .map(|&i| permutation.iter().map(|&j| a[i][j].clone()).collect())
        .collect();
    let (other, permuted_shear) =
        EpsilonShearing::regularize(&system(permuted), &RunContext::default()).unwrap();
    assert_eq!(
        permuted_shear.weights(),
        permutation.map(|i| shear.weights()[i])
    );
    let coordinate = symbol!("shear::s");
    for (i, &old_i) in permutation.iter().enumerate() {
        for (j, &old_j) in permutation.iter().enumerate() {
            assert_eq!(
                other.derivatives[&coordinate][i][j],
                regular.derivatives[&coordinate][old_i][old_j]
            );
        }
    }
    let cancelled = system(vec![vec![parse("(1i*(eps+s)-(1i*eps+1i*s))/eps^3")]]);
    let (regular, shear) = EpsilonShearing::regularize(&cancelled, &RunContext::default()).unwrap();
    assert_eq!(shear.weights(), &[0]);
    assert!(regular.derivatives[&coordinate][0][0].is_zero());
    assert!(matches!(
        EpsilonShearing::regularize(
            &system(vec![
                vec![Atom::new(), parse("1/eps^10001")],
                vec![Atom::new(), Atom::new()]
            ]),
            &RunContext::default()
        ),
        Err(Error::Limit(_))
    ));
}

#[test]
fn malformed_boundary_and_error_coverage_cannot_be_promoted() {
    let (_, shear) = EpsilonShearing::regularize(&system(zero(1)), &RunContext::default()).unwrap();
    let p = Precision::decimal(40).unwrap();
    let b = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(3)]],
    };
    assert!(matches!(
        shear.to_sheared(&b, &[], 0),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.to_sheared(&b, &[vec![p.real(-1)]], 0),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.to_sheared(&b, &[vec![p.real(0), p.real(0)]], 0),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.scaled_basis_labels(&[]),
        Err(Error::InvalidInput(_))
    ));
    let empty = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![],
    };
    assert!(matches!(
        shear.to_sheared(&empty, &[], 0),
        Err(Error::Limit(_))
    ));
    let infinity = Float::with_val(
        p.bits,
        symbolica::domains::float::Float::parse("inf", Some(p.bits))
            .unwrap()
            .as_raw(),
    );
    let mut invalid = b.clone();
    invalid.point.re = infinity.clone();
    assert!(matches!(
        shear.to_sheared(&invalid, &[vec![p.real(0)]], 0),
        Err(Error::InvalidInput(_))
    ));
    let mut invalid = b.clone();
    invalid.coefficients[0][0].im = infinity.clone();
    assert!(matches!(
        shear.to_sheared(&invalid, &[vec![p.real(0)]], 0),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        shear.to_sheared(&b, &[vec![infinity]], 0),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn shared_fuchsian_constraint_core_preserves_physical_derivative_term() {
    let variable = symbol!("shear::s");
    let source = DifferentialSystem {
        variable,
        matrix: vec![
            vec![Atom::new(), parse("1/s^3")],
            vec![parse("s"), Atom::new()],
        ],
    };
    let (regular, weights) = source.diagonal_fuchsian_form().unwrap();
    assert_eq!(weights, vec![0, 2]);
    let expected = vec![
        vec![Atom::new(), parse("1/s")],
        vec![parse("1/s"), parse("-2/s")],
    ];
    for (row, expected) in regular.matrix.iter().zip(expected) {
        for (entry, wanted) in row.iter().zip(expected) {
            assert!((entry - wanted).expand().is_zero());
        }
    }
}
