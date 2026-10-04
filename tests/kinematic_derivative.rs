use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::kinematic_derivative::KinematicDerivative;
use symbolica_amflow::reduction::{LinearCombination, Reduction};
use symbolica_amflow::*;

fn equal(actual: &Atom, expected: &Atom) {
    assert!(
        (actual - expected).together().cancel().is_zero(),
        "{actual} != {expected}"
    );
}
fn bubble(s: Atom) -> IntegralFamily {
    let gram = vec![vec![s]];
    IntegralFamily {
        name: "physical_invariant_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::new(), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1], Atom::new(), &gram).unwrap(),
        ],
        external_gram: gram,
        physical_propagators: 2,
        epsilon: symbol!("physical_derivative::eps"),
        dimension: 4,
    }
}

#[test]
fn nondiagonal_gram_vector_field_obeys_covariance_and_explicit_mass_chain_rule() {
    let s = Atom::var(symbol!("physical_derivative::s"));
    let t_symbol = symbol!("physical_derivative::t");
    let t = Atom::var(t_symbol);
    let u = Atom::var(symbol!("physical_derivative::u"));
    let gram = vec![vec![s.clone(), t.clone()], vec![t.clone(), u.clone()]];
    let family = IntegralFamily {
        name: "two_external_derivative".into(),
        loops: vec!["l".into()],
        external: vec!["p1".into(), "p2".into()],
        propagators: vec![
            Propagator::quadratic(&[1], &[0, 0], t.clone().pow(2), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1, 0], Atom::new(), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1, 2], t.clone().pow(3), &gram).unwrap(),
        ],
        external_gram: gram.clone(),
        physical_propagators: 3,
        epsilon: symbol!("physical_derivative::eps"),
        dimension: 4,
    };
    let derivative = KinematicDerivative::new(&family, t_symbol).unwrap();
    let determinant = &s * &u - t.clone().pow(2);
    let expected = [
        [
            -&t / (Atom::num(2) * &determinant),
            &s / (Atom::num(2) * &determinant),
        ],
        [
            &u / (Atom::num(2) * &determinant),
            -&t / (Atom::num(2) * &determinant),
        ],
    ];
    for (i, row) in derivative.external_vector_field().iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            equal(entry, &expected[i][j]);
            let covariance = (0..2).fold(Atom::new(), |sum, k| {
                sum + &row[k] * &gram[k][j]
                    + &gram[i][k] * &derivative.external_vector_field()[j][k]
            });
            equal(&covariance, &gram[i][j].derivative(t_symbol));
        }
    }
    let derivatives = derivative.denominator_derivatives();
    equal(&derivatives[0].constant, &(-Atom::num(2) * &t));
    equal(
        &derivatives[2].constant,
        &(Atom::num(4) - Atom::num(3) * t.clone().pow(2)),
    );
    equal(&derivatives[2].scalar_products[0], &Atom::new());
    equal(
        &derivatives[2].scalar_products[1],
        &((-&t + Atom::num(2) * &u) / &determinant),
    );
    equal(
        &derivatives[2].scalar_products[2],
        &((&s - Atom::num(2) * &t) / &determinant),
    );
    assert!(
        derivative
            .nonzero_conditions()
            .iter()
            .any(|c| (c - &determinant).together().cancel().is_zero())
    );
}

#[test]
fn explicit_routing_and_invariant_dependence_are_both_differentiated() {
    let x_symbol = symbol!("physical_derivative::routing_x");
    let x = Atom::var(x_symbol);
    let mut family = bubble(x.clone());
    family.propagators[1] = Propagator {
        constant: x.clone().pow(3) - x.clone().pow(2),
        scalar_products: vec![Atom::num(1), Atom::num(2) * &x],
    };
    let derivative = KinematicDerivative::new(&family, x_symbol).unwrap();
    let second = &derivative.denominator_derivatives()[1];
    // D2=(l+x*p)^2-x^2, p^2=x; d(l.p)/dx=(l.p)/(2*x).
    equal(
        &second.constant,
        &(Atom::num(3) * x.clone().pow(2) - Atom::num(2) * &x),
    );
    equal(&second.scalar_products[1], &Atom::num(3));
}

#[test]
fn negative_isp_powers_and_inherited_domains_survive_reduction_closure() {
    let s_symbol = symbol!("physical_derivative::isp_s");
    let s = Atom::var(s_symbol);
    let u_symbol = symbol!("physical_derivative::domain_u");
    let u = Atom::var(u_symbol);
    let mut family = bubble(s.clone());
    family.physical_propagators = 1;
    family.propagators[0] = Propagator {
        constant: Atom::num(-1),
        scalar_products: vec![Atom::num(1) / &u, Atom::new()],
    };
    family.propagators[1] = Propagator {
        constant: Atom::new(),
        scalar_products: vec![Atom::new(), Atom::num(1)],
    };
    let target = Integral(vec![2, -2]);
    let derivative = KinematicDerivative::new(&family, s_symbol).unwrap();
    let terms = derivative.integral(&target).unwrap();
    assert_eq!(terms.len(), 1);
    equal(&terms[&target], &(Atom::num(1) / &s));
    let backend = TableBackend {
        name: "physical_derivative_identity_table".into(),
        reduction: Reduction {
            residuals: vec![target.clone()],
            ..Default::default()
        },
    };
    let system = derivative
        .differential_system(&backend, &[target], 3, &RunContext::default())
        .unwrap();
    equal(&system.matrix[0][0], &(Atom::num(1) / s));
    for variable in [u_symbol, s_symbol] {
        let zero = BTreeMap::from([(Atom::var(variable), Atom::new())]);
        assert!(
            system.nonzero_conditions.iter().any(|condition| {
                !condition.derivative(variable).is_zero()
                    && family::substitute(condition, &zero)
                        .together()
                        .cancel()
                        .is_zero()
            }),
            "missing input/Gram guard for {variable:?}"
        );
    }
}

#[test]
fn native_massless_bubble_system_has_the_physical_scaling_derivative() {
    let s_symbol = symbol!("physical_derivative::bubble_s");
    let s = Atom::var(s_symbol);
    let family = bubble(s.clone());
    let target = Integral(vec![1, 1]);
    let derivative = KinematicDerivative::new(&family, s_symbol).unwrap();
    let terms = derivative.integral(&target).unwrap();
    let expected: LinearCombination = BTreeMap::from([
        (target.clone(), -Atom::num(1) / (Atom::num(2) * &s)),
        (Integral(vec![0, 2]), Atom::num(1) / (Atom::num(2) * &s)),
        (Integral(vec![1, 2]), Atom::num((-1, 2))),
    ]);
    assert_eq!(terms.len(), expected.len());
    for (integral, coefficient) in expected {
        equal(&terms[&integral], &coefficient);
    }
    let fixed = reduction::parameter_derivative(&family, &target, s_symbol).unwrap();
    assert_eq!(fixed, vec![(Integral(vec![1, 2]), Atom::num(-1))]);
    let system = derivative
        .differential_system(
            &RustRedBackend::default(),
            std::slice::from_ref(&target),
            8,
            &RunContext::default(),
        )
        .unwrap();
    assert_eq!(system.basis, vec![target]);
    equal(&system.matrix[0][0], &(-Atom::var(family.epsilon) / &s));
    assert!(system.nonzero_conditions.iter().any(|condition| {
        let zero = BTreeMap::from([(s.clone(), Atom::new())]);
        !condition.derivative(s_symbol).is_zero()
            && family::substitute(condition, &zero)
                .together()
                .cancel()
                .is_zero()
    }));
}

#[test]
fn constant_singular_gram_allows_mass_derivatives_but_a_varying_one_is_rejected() {
    let m_symbol = symbol!("physical_derivative::mass");
    let m = Atom::var(m_symbol);
    let mut family = bubble(Atom::new());
    family.physical_propagators = 1;
    family.propagators[0].constant = -m;
    family.propagators[1] = Propagator {
        constant: Atom::new(),
        scalar_products: vec![Atom::new(), Atom::num(1)],
    };
    let derivative = KinematicDerivative::new(&family, m_symbol).unwrap();
    assert!(derivative.external_vector_field()[0][0].is_zero());
    assert_eq!(
        derivative.integral(&Integral(vec![2, 0])).unwrap(),
        BTreeMap::from([(Integral(vec![3, 0]), Atom::num(2))])
    );

    let x = Atom::var(m_symbol);
    let gram = vec![vec![x.clone(), x.clone()], vec![x.clone(), x]];
    family.external = vec!["p1".into(), "p2".into()];
    family.external_gram = gram;
    family.propagators = (0..3)
        .map(|i| Propagator {
            constant: Atom::new(),
            scalar_products: (0..3).map(|j| Atom::num(i64::from(i == j))).collect(),
        })
        .collect();
    assert!(
        matches!(KinematicDerivative::new(&family, m_symbol), Err(Error::Unsupported(message)) if message.contains("singular external Gram"))
    );
    assert!(
        matches!(KinematicDerivative::new(&family, family.epsilon), Err(Error::Unsupported(message)) if message.contains("regulator"))
    );
}
