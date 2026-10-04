//! The host follows RustRed's compiled registry for a physical five-loop family.
use symbolica::prelude::*;
use symbolica_amflow::{
    Integral, IntegralFamily, Propagator, Result, RunContext,
    reduction::{ReductionBackend, RustRedBackend},
};

#[test]
fn fifteen_slot_five_loop_tadpole_product_obeys_exact_euler_ibp() -> Result<()> {
    // A custom RustRed build may intentionally omit this arity. The same host
    // must then report that capability limit, rather than hard-code a new list.
    let epsilon = symbol!("native_arity::epsilon");
    let mass_squared = Atom::var(symbol!("native_arity::mass_squared"));
    let coordinates = (0..5)
        .flat_map(|i| (i..5).map(move |j| (i, j)))
        .collect::<Vec<_>>();
    let order = coordinates
        .iter()
        .enumerate()
        .filter(|(_, (i, j))| i == j)
        .chain(coordinates.iter().enumerate().filter(|(_, (i, j))| i != j));
    let propagators = order
        .map(|(index, (i, j))| Propagator {
            constant: if i == j {
                -mass_squared.clone()
            } else {
                Atom::new()
            },
            scalar_products: (0..15)
                .map(|column| Atom::num(i64::from(column == index)))
                .collect(),
        })
        .collect();
    let family = IntegralFamily {
        name: "five_massive_tadpoles_with_complete_isp_basis".into(),
        loops: (0..5).map(|i| format!("k{i}")).collect(),
        external: Vec::new(),
        external_gram: Vec::new(),
        propagators,
        physical_propagators: 5,
        epsilon,
        dimension: 4,
    };
    let mut powers = vec![0; 15];
    powers[..5].fill(1);
    let master = Integral(powers.clone());
    powers[0] = 2;
    let raised = Integral(powers);
    let expected = (Atom::one() - Atom::var(epsilon)) / mass_squared;
    for factorized in [true, false] {
        let backend = RustRedBackend {
            factorized,
            bubble_subloops: false,
            parametric_rules: false,
            symmetry_rules: false,
            max_depth: 1,
            max_targets: 1024,
            ..Default::default()
        };
        let result = backend.reduce(
            &family,
            std::slice::from_ref(&raised),
            &RunContext::default(),
        );
        if !rustred::compiled_runtime_arities().contains(&15) {
            assert!(matches!(
                result,
                Err(symbolica_amflow::Error::Unsupported(_))
            ));
            continue;
        }
        let reduction = result?;
        let rule = reduction.expand(&raised)?;
        assert_eq!(rule.len(), 1);
        assert!((&rule[&master] - &expected).together().cancel().is_zero());
        assert!(reduction.residuals.contains(&master));
        let sampled = backend.reduce_at_epsilon(
            &family,
            std::slice::from_ref(&raised),
            &Rational::from((1, 13)),
            &RunContext::default(),
        )?;
        let sample_rule = sampled.expand(&raised)?;
        let expected_sample = expected
            .replace(Atom::var(epsilon))
            .with(Atom::num((1, 13)));
        assert_eq!(sample_rule.len(), 1);
        assert!(
            (&sample_rule[&master] - expected_sample)
                .together()
                .cancel()
                .is_zero()
        );
    }
    Ok(())
}

#[test]
fn uncompiled_arity_is_a_typed_capability_limit_in_both_backends() -> Result<()> {
    let count = (1..)
        .find(|n| !rustred::compiled_runtime_arities().contains(n))
        .unwrap();
    let family = IntegralFamily {
        name: "uncompiled_complete_family".into(),
        loops: vec!["k".into()],
        external: (1..count).map(|i| format!("p{i}")).collect(),
        external_gram: (1..count)
            .map(|i| (1..count).map(|j| Atom::num(i64::from(i == j))).collect())
            .collect(),
        propagators: (0..count)
            .map(|i| Propagator {
                constant: Atom::num(if i == 0 { -1 } else { 0 }),
                scalar_products: (0..count).map(|j| Atom::num(i64::from(i == j))).collect(),
            })
            .collect(),
        physical_propagators: 1,
        epsilon: symbol!("unsupported_arity::epsilon"),
        dimension: 4,
    };
    let mut powers = vec![0; count];
    powers[0] = 2;
    for factorized in [true, false] {
        let backend = RustRedBackend {
            factorized,
            ..Default::default()
        };
        assert!(matches!(
            backend.reduce(&family, &[Integral(powers.clone())], &RunContext::default()),
            Err(symbolica_amflow::Error::Unsupported(_))
        ));
    }
    Ok(())
}
