use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{reduction::Reduction, *};

fn atom(expression: &str) -> Atom {
    Atom::parse(expression, "scoped_reductions", Default::default()).unwrap()
}
fn epsilon() -> Symbol {
    symbol!("scoped_reductions::eps")
}
fn tadpole(mass: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "supplied_tadpole".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -mass,
            scalar_products: vec![Atom::one()],
        }],
        physical_propagators: 1,
        epsilon: epsilon(),
        dimension: 4,
    }
}
fn rule(coefficient: Atom) -> Reduction {
    Reduction {
        rules: BTreeMap::from([(
            Integral(vec![2]),
            BTreeMap::from([(Integral(vec![1]), coefficient)]),
        )]),
        residuals: vec![Integral(vec![1])],
        nonzero_conditions: vec![],
    }
}

#[test]
fn supplied_physical_table_closes_actual_mass_derivative() {
    let mass = symbol!("scoped_reductions::mass_squared");
    let family = tadpole(Atom::var(mass));
    let coefficient = (Atom::one() - Atom::var(epsilon())) / Atom::var(mass);
    let backend =
        ScopedTableBackend::from_table("physical", &family, rule(coefficient.clone())).unwrap();
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[mass],
        &backend,
        &FlowOptions::default(),
        "positive mass",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![1])]);
    assert!(
        (&prepared.flow().system().derivatives[&mass][0][0] - coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        prepared
            .nonzero_conditions()
            .iter()
            .any(|a| a == &Atom::var(mass))
    );
}

#[test]
fn algebraically_equal_auxiliary_coefficients_share_scope() {
    let family = tadpole(Atom::num(2));
    let eta = symbol!("symbolica_amflow::eta");
    let (deformed, _) = family.deform(eta, &MassMode::All).unwrap();
    // HEPKit conversion normalizes rational coefficients. Native auxiliary
    // deformation produces the same polynomial as an unfactored sum.
    let mut normalized = deformed.clone();
    normalized.propagators[0].constant = normalized.propagators[0].constant.together().cancel();
    assert_ne!(
        normalized.propagators[0].constant,
        deformed.propagators[0].constant
    );
    let coefficient = (Atom::one() - Atom::var(epsilon())) / (Atom::num(2) + Atom::var(eta));
    let mut backend =
        ScopedTableBackend::from_table("normalized", &normalized, rule(coefficient.clone()))
            .unwrap();
    let reduced = backend
        .reduce(&deformed, &[Integral(vec![2])], &RunContext::default())
        .unwrap();
    assert!(
        (&reduced.expand(&Integral(vec![2])).unwrap()[&Integral(vec![1])] - &coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert_eq!(
        backend.identity(),
        ScopedTableBackend::from_table("normalized", &deformed, rule(coefficient.clone()))
            .unwrap()
            .identity()
    );
    assert!(matches!(
        backend.insert(&deformed, rule(coefficient)),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn normalized_family_scopes_retain_uncancelled_domain_restrictions() {
    let guarded = tadpole(atom("(m^2-1)/(m-1)"));
    let unrestricted = tadpole(atom("m+1"));
    assert!(
        (&guarded.propagators[0].constant - &unrestricted.propagators[0].constant)
            .together()
            .cancel()
            .is_zero()
    );
    let coefficient = atom("(1-eps)/(m+1)");
    let backend =
        ScopedTableBackend::from_table("original_domain", &guarded, rule(coefficient.clone()))
            .unwrap();
    let reduced = backend
        .reduce(&guarded, &[Integral(vec![2])], &RunContext::default())
        .unwrap();
    assert!(reduced.nonzero_conditions.iter().any(|condition| {
        condition
            .replace(atom("m"))
            .with(Atom::one())
            .together()
            .cancel()
            .is_zero()
    }));
    assert!(matches!(
        backend.reduce(&unrestricted, &[Integral(vec![2])], &RunContext::default()),
        Err(Error::IncompleteReduction(_))
    ));
    assert_ne!(
        backend.identity(),
        ScopedTableBackend::from_table("original_domain", &unrestricted, rule(coefficient))
            .unwrap()
            .identity()
    );
}

#[test]
fn supplied_deformed_table_runs_automatic_amf_boundary_and_transport() {
    let family = tadpole(Atom::num(2));
    let options = FlowOptions {
        digits: 24,
        guard_digits: 35,
        series_order: 72,
        ..Default::default()
    };
    let context = RunContext::default();
    let eta = symbol!("symbolica_amflow::eta");
    let (deformed, _) = family.deform(eta, &options.mass_mode).unwrap();
    let backend = ScopedTableBackend::from_table(
        "deformed",
        &deformed,
        rule((Atom::one() - Atom::var(epsilon())) / (Atom::num(2) + Atom::var(eta))),
    )
    .unwrap();
    let prepared = PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    assert_eq!(prepared.reduced.basis, [Integral(vec![1])]);
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let sample = Rational::from((1, 97));
    let value = prepared
        .evaluate(&sample, &options, &boundary, &context)
        .unwrap()
        .remove(0);
    let p = Precision::decimal(59).unwrap();
    let expected = vacuum::tadpole(
        1,
        &p.i(2),
        &(Rational::from(4) - &sample * &Rational::from(2)),
        p,
    )
    .unwrap();
    assert!(p.close(&value, &expected, 24));
}

#[test]
fn undeformed_and_missing_recursive_family_tables_are_not_reused() {
    let family = tadpole(Atom::num(2));
    let backend =
        ScopedTableBackend::from_table("original_only", &family, rule(atom("(1-eps)/2"))).unwrap();
    let options = FlowOptions::default();
    let context = RunContext::default();
    let error = match PreparedFlow::new(
        &family,
        &[Integral(vec![1])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    ) {
        Ok(_) => panic!("undeformed table was used for auxiliary flow"),
        Err(error) => error,
    };
    assert!(
        matches!(error, Error::IncompleteReduction(ref message) if message.contains("auxiliary deformations"))
    );

    let gram = vec![vec![Atom::num(-1)]];
    let child = IntegralFamily {
        name: "recursive_bubble_child".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::num(2), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1], Atom::num(3), &gram).unwrap(),
        ],
        external_gram: gram,
        physical_propagators: 2,
        epsilon: epsilon(),
        dimension: 4,
    };
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let error = boundary
        .evaluate(
            &child,
            &Integral(vec![1, 1]),
            &Rational::from((1, 97)),
            Precision::decimal(50).unwrap(),
        )
        .unwrap_err();
    assert!(
        matches!(error, Error::IncompleteReduction(ref message) if message.contains("recursive boundary families"))
    );
}

#[test]
fn exact_scope_preserves_parameters_dimension_routing_slots_and_roles() {
    let gram = vec![vec![atom("s")]];
    let family = IntegralFamily {
        name: "scope_collision".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        propagators: vec![
            Propagator::quadratic(&[1], &[0], atom("m"), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1], atom("n"), &gram).unwrap(),
        ],
        external_gram: gram,
        physical_propagators: 2,
        epsilon: epsilon(),
        dimension: 4,
    };
    let residual = Integral(vec![1, 0]);
    let reduction = Reduction {
        residuals: vec![residual.clone()],
        ..Default::default()
    };
    let backend = ScopedTableBackend::from_table("scope", &family, reduction.clone()).unwrap();
    let mut variants = Vec::new();
    let mut changed = family.clone();
    changed.propagators[0].constant = -atom("other_namespace::m");
    variants.push(changed);
    let mut changed = family.clone();
    changed.propagators[0].constant = -atom("m+1");
    variants.push(changed);
    let mut changed = family.clone();
    changed.dimension = 6;
    variants.push(changed);
    let mut changed = family.clone();
    changed.epsilon = symbol!("scoped_reductions::different_eps");
    variants.push(changed);
    let mut changed = family.clone();
    changed.propagators[1].scalar_products[1] = Atom::num(-2);
    variants.push(changed);
    let mut changed = family.clone();
    changed.propagators.swap(0, 1);
    variants.push(changed);
    let mut changed = family.clone();
    changed.physical_propagators = 1;
    variants.push(changed);
    let mut changed = family.clone();
    changed.external_gram[0][0] = atom("s+1");
    variants.push(changed);
    for changed in variants {
        assert!(matches!(
            backend.reduce(
                &changed,
                std::slice::from_ref(&residual),
                &RunContext::default()
            ),
            Err(Error::IncompleteReduction(_))
        ));
        let other = ScopedTableBackend::from_table("scope", &changed, reduction.clone()).unwrap();
        assert_ne!(backend.identity(), other.identity());
    }
}

#[test]
fn admission_preserves_cancelled_denominators_and_rejects_sampled_poles() {
    let family = tadpole(atom("m"));
    let backend = ScopedTableBackend::from_table(
        "guarded",
        &family,
        rule(atom("(m^2-1)/(m-1)+(eps^2-eps)/(eps-1)")),
    )
    .unwrap();
    let reduction = backend
        .reduce(&family, &[Integral(vec![2])], &RunContext::default())
        .unwrap();
    for variable in [atom("m"), Atom::var(epsilon())] {
        assert!(reduction.nonzero_conditions.iter().any(|condition| {
            condition
                .replace(variable.clone())
                .with(Atom::one())
                .together()
                .cancel()
                .is_zero()
        }));
    }
    assert!(matches!(
        backend.reduce_at_epsilon(
            &family,
            &[Integral(vec![2])],
            &Rational::from(1),
            &RunContext::default()
        ),
        Err(Error::Reduction(_))
    ));
    let sampled = backend
        .reduce_at_epsilon(
            &family,
            &[Integral(vec![2])],
            &Rational::from((1, 97)),
            &RunContext::default(),
        )
        .unwrap();
    assert!(sampled.nonzero_conditions.iter().any(|a| {
        a.get_all_symbols(false)
            .contains(&symbol!("scoped_reductions::m"))
    }));
    assert!(ScopedTableBackend::from_table("foreign", &family, rule(atom("D/(2*m)"))).is_err());
    let invalid = Reduction {
        nonzero_conditions: vec![Atom::Zero],
        ..rule(Atom::one())
    };
    assert!(ScopedTableBackend::from_table("zero_guard", &family, invalid).is_err());
}

#[test]
fn table_graph_validation_and_collection_admission_are_transactional() {
    let first = tadpole(Atom::num(2));
    let second = tadpole(Atom::num(3));
    let mut backend = ScopedTableBackend::new("collection");
    backend.insert(&first, rule(atom("(1-eps)/2"))).unwrap();
    let before = backend.identity();
    let missing = Reduction {
        residuals: vec![],
        ..rule(Atom::one())
    };
    assert!(matches!(
        backend.insert(&second, missing),
        Err(Error::IncompleteReduction(_))
    ));
    assert_eq!(backend.identity(), before);
    let cyclic = Reduction {
        rules: BTreeMap::from([
            (
                Integral(vec![1]),
                BTreeMap::from([(Integral(vec![2]), Atom::one())]),
            ),
            (
                Integral(vec![2]),
                BTreeMap::from([(Integral(vec![1]), Atom::one())]),
            ),
        ]),
        ..Default::default()
    };
    assert!(backend.insert(&second, cyclic).is_err());
    assert!(backend.insert(&first, rule(Atom::one())).is_err());
    assert_eq!(backend.identity(), before);
    backend.insert(&second, rule(atom("(1-eps)/3"))).unwrap();
    assert_eq!(backend.len(), 2);
    assert!(!backend.is_empty());
    let mut reverse = ScopedTableBackend::new("collection");
    reverse.insert(&second, rule(atom("(1-eps)/3"))).unwrap();
    reverse.insert(&first, rule(atom("(1-eps)/2"))).unwrap();
    assert_eq!(backend.identity(), reverse.identity());
}

#[test]
fn supplied_tables_do_not_require_scalar_product_completion_or_native_dispatch() {
    let mut family = tadpole(Atom::num(2));
    family.external = vec!["p".into()];
    family.external_gram = vec![vec![Atom::num(-1)]];
    family.propagators[0].scalar_products.push(Atom::Zero);
    // One denominator for two scalar products: valid table metadata, although
    // automatic derivative construction will need a completed family later.
    let backend =
        ScopedTableBackend::from_table("incomplete_basis", &family, rule(atom("(1-eps)/2")))
            .unwrap();
    assert!(
        backend
            .reduce(&family, &[Integral(vec![2])], &RunContext::default())
            .is_ok()
    );

    let width = rustred::compiled_runtime_arities()
        .iter()
        .copied()
        .max()
        .unwrap()
        + 1;
    family.external = (1..width).map(|i| format!("p{i}")).collect();
    family.external_gram = (0..width - 1)
        .map(|i| {
            (0..width - 1)
                .map(|j| Atom::num(i64::from(i == j)))
                .collect()
        })
        .collect();
    family.propagators = (0..width)
        .map(|i| Propagator {
            constant: Atom::num(-1),
            scalar_products: (0..width).map(|j| Atom::num(i64::from(i == j))).collect(),
        })
        .collect();
    family.physical_propagators = width;
    let target = Integral(vec![1; width]);
    let backend = ScopedTableBackend::from_table(
        "above_dispatch",
        &family,
        Reduction {
            residuals: vec![target.clone()],
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        backend
            .reduce(
                &family,
                std::slice::from_ref(&target),
                &RunContext::default()
            )
            .unwrap()
            .expand(&target)
            .unwrap(),
        BTreeMap::from([(target, Atom::one())])
    );
}
