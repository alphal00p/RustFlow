use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{reduction::Reduction, *};

fn tadpole(mass: Atom, scale: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "physical_family_identity_test".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -&scale * mass,
            scalar_products: vec![scale],
        }],
        physical_propagators: 1,
        epsilon: symbol!("physical_family_test::epsilon"),
        dimension: 4,
    }
}
fn table(
    rules: BTreeMap<Integral, BTreeMap<Integral, Atom>>,
    conditions: Vec<Atom>,
) -> TableBackend {
    TableBackend {
        name: "physical_family_exact_tadpole".into(),
        reduction: Reduction {
            rules,
            residuals: vec![Integral(vec![1])],
            nonzero_conditions: conditions,
        },
    }
}
#[test]
fn independent_partial_derivatives_close_even_when_their_sum_cancels() {
    let (s, t) = symbol!("physical_family_test::s", "physical_family_test::t");
    let mass = Atom::var(s) - Atom::var(t);
    let family = tadpole(mass.clone(), Atom::one());
    let coefficient = (Atom::one() - Atom::var(family.epsilon)) / &mass;
    let backend = table(
        BTreeMap::from([(
            Integral(vec![2]),
            BTreeMap::from([(Integral(vec![1]), coefficient.clone())]),
        )]),
        vec![mass],
    );
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[s, t],
        &backend,
        &FlowOptions::default(),
        "one mass sheet",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![1])]);
    let derivatives = &prepared.flow().system().derivatives;
    assert!(
        (&derivatives[&s][0][0] - &coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        (&derivatives[&t][0][0] + &coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        (&derivatives[&s][0][0] + &derivatives[&t][0][0])
            .together()
            .cancel()
            .is_zero()
    );
}
#[test]
fn hidden_normalization_parameters_are_rejected_before_reduction() {
    struct ForbiddenBackend;
    impl ReductionBackend for ForbiddenBackend {
        fn identity(&self) -> String {
            "must not run".into()
        }
        fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
            panic!("hidden family parameter must fail before reduction")
        }
    }
    let (s, u) = symbol!("physical_family_test::s", "physical_family_test::hidden_u");
    let family = tadpole(Atom::one(), Atom::var(u));
    let result = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[s],
        &ForbiddenBackend,
        &FlowOptions::default(),
        "one sheet",
        &RunContext::default(),
    );
    assert!(matches!(result, Err(Error::InvalidInput(_))));
}
#[test]
fn identical_connections_with_distinct_family_normalizations_have_distinct_identities() {
    let variable = symbol!("physical_family_test::s");
    let backend = table(BTreeMap::new(), vec![]);
    let prepare = |scale| {
        PreparedPhysicalFamily::new(
            &tadpole(Atom::one(), Atom::num(scale)),
            &[Integral(vec![1])],
            &[variable],
            &backend,
            &FlowOptions::default(),
            "one sheet",
            &RunContext::default(),
        )
        .unwrap()
    };
    let first = prepare(1);
    let second = prepare(2);
    assert_eq!(
        first.flow().system().derivatives,
        second.flow().system().derivatives
    );
    assert_ne!(
        first.flow().identity().key(),
        second.flow().identity().key()
    );
}
