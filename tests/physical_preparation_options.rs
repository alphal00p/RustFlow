use std::collections::BTreeMap;
use std::sync::Mutex;
use symbolica::prelude::*;
use symbolica_amflow::{reduction::Reduction, transport_cache::*, *};

fn family(mass: Atom, scale: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "physical_preparation_options".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -&scale * mass,
            scalar_products: vec![scale],
        }],
        physical_propagators: 1,
        epsilon: symbol!("physical_options::eps"),
        dimension: 4,
    }
}
fn exact(a: &Atom, b: &Atom) {
    assert!((a - b).together().cancel().is_zero(), "{a} != {b}");
}
struct Recorded {
    table: TableBackend,
    requests: Mutex<Vec<Vec<Integral>>>,
    round_guards: Option<Symbol>,
}
impl ReductionBackend for Recorded {
    fn identity(&self) -> String {
        "physical-option-test-recorded-table".into()
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        let mut requests = self.requests.lock().unwrap();
        requests.push(targets.to_vec());
        let mut result = self.table.reduce(family, targets, context)?;
        if let Some(s) = self.round_guards {
            result.nonzero_conditions = vec![Atom::var(s) - Atom::num(requests.len() as i64 + 1)];
        }
        Ok(result)
    }
}
fn table(coefficient: Atom, third: Atom) -> Recorded {
    Recorded {
        table: TableBackend {
            name: "exact supplied coefficient relations".into(),
            reduction: Reduction {
                rules: BTreeMap::from([
                    (
                        Integral(vec![2]),
                        BTreeMap::from([(Integral(vec![1]), coefficient)]),
                    ),
                    (
                        Integral(vec![3]),
                        BTreeMap::from([(Integral(vec![1]), third)]),
                    ),
                ]),
                residuals: vec![Integral(vec![1])],
                ..Default::default()
            },
        },
        requests: Mutex::default(),
        round_guards: None,
    }
}

#[test]
fn physical_skip_reduces_derivatives_without_canceling_independent_coordinates() {
    let (s, t) = symbol!("physical_options::s", "physical_options::t");
    let mass = Atom::var(s) - Atom::var(t);
    let family = family(mass.clone(), Atom::one());
    let coefficient = (Atom::one() - Atom::var(family.epsilon)) / mass;
    let backend = table(coefficient.clone(), Atom::new());
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[s, t],
        &backend,
        &FlowOptions {
            skip_reduction: true,
            ..Default::default()
        },
        "regular mass sheet",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(
        *backend.requests.lock().unwrap(),
        vec![vec![Integral(vec![2])]]
    );
    exact(
        &prepared.flow().system().derivatives[&s][0][0],
        &coefficient,
    );
    exact(
        &prepared.flow().system().derivatives[&t][0][0],
        &(-coefficient),
    );
    assert_eq!(
        prepared.target_reductions()[0],
        BTreeMap::from([(Integral(vec![1]), Atom::one())])
    );
    assert!(prepared.basis_transformations().is_empty());
    assert!(prepared.basis_refinement().is_none());
}

#[test]
fn retained_constant_family_needs_no_reduction_call() {
    struct Forbidden;
    impl ReductionBackend for Forbidden {
        fn identity(&self) -> String {
            "unused".into()
        }
        fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
            panic!("a zero derivative needs no reduction")
        }
    }
    let s = symbol!("physical_options::s");
    let prepared = PreparedPhysicalFamily::new(
        &family(Atom::one(), Atom::one()),
        &[Integral(vec![2])],
        &[s],
        &Forbidden,
        &FlowOptions {
            skip_reduction: true,
            ..Default::default()
        },
        "constant sheet",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![2])]);
    assert!(prepared.flow().system().derivatives[&s][0][0].is_zero());
}

#[test]
fn automatic_seed_rejects_unproved_pole_bounds_before_reduction_or_cache_mutation() {
    struct Forbidden;
    impl ReductionBackend for Forbidden {
        fn identity(&self) -> String {
            "unused for pole admission".into()
        }
        fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
            panic!("unproved pole bounds must fail before seed reduction")
        }
    }
    let s = symbol!("physical_options::unused_s");
    let family = family(Atom::one(), Atom::var(symbol!("physical_options::eps")));
    let options = FlowOptions {
        skip_reduction: true,
        ..Default::default()
    };
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![4])],
        &[s],
        &Forbidden,
        &options,
        "epsilon-normalized tadpole",
        &RunContext::default(),
    )
    .unwrap();
    let mut cache = RustFlowCache::default();
    let result = prepared.seed_cache(
        &mut cache,
        &BTreeMap::from([(s, Atom::one())]),
        0,
        0,
        &options,
        &Forbidden,
        &RunContext::default(),
    );
    assert!(matches!(result, Err(Error::Unsupported(_))));
    assert!(cache.entries().is_empty());
}

#[test]
fn both_closure_modes_retain_conditions_from_every_round() {
    let s = symbol!("physical_options::s");
    let family = family(Atom::var(s), Atom::one());
    let eps = Atom::var(family.epsilon);
    for skip in [false, true] {
        let c = (Atom::one() - &eps) / Atom::var(s);
        let mut backend = table(c.clone(), -&eps * c / (Atom::num(2) * Atom::var(s)));
        backend.round_guards = Some(s);
        let prepared = PreparedPhysicalFamily::new(
            &family,
            &[Integral(vec![2])],
            &[s],
            &backend,
            &FlowOptions {
                skip_reduction: skip,
                ..Default::default()
            },
            "restricted supplied domain",
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(backend.requests.lock().unwrap().len(), 2);
        let p = Precision::decimal(60).unwrap();
        let entry = |point| CachedBoundary {
            identity: prepared.flow().identity().clone(),
            point: CachedPoint::Exact(BTreeMap::from([(s, Atom::num(point))])),
            kind: PointKind::Physical,
            range: EpsilonRange::new(0, 0).unwrap(),
            coefficients: vec![vec![p.zero(); prepared.basis().len()]],
            accuracy: BoundaryAccuracy::supplied(
                20,
                p.bits,
                vec![vec![p.real(0); prepared.basis().len()]],
                "exact zero fixture for cache domain admission",
            )
            .unwrap(),
        };
        for point in [2, 3] {
            assert!(
                RustFlowCache::default().insert(entry(point)).is_err(),
                "lost round guard at {point}, skip={skip}"
            );
        }
        RustFlowCache::default().insert(entry(4)).unwrap();
    }
}

#[test]
fn physical_refinement_uses_each_partial_derivative_and_exact_complex_coefficients() {
    let (s, t) = symbol!("physical_options::s", "physical_options::t");
    for complex in [false, true] {
        for dimension in [4, 6] {
            let phase = if complex {
                Atom::num(Complex::new(Rational::one(), Rational::one()))
            } else {
                Atom::one()
            };
            let eps = Atom::var(symbol!("physical_options::eps"));
            let mass = &phase * Atom::var(s);
            let scale = &phase * (&eps + Atom::var(t));
            let mut family = family(mass.clone(), scale.clone());
            family.dimension = dimension;
            let h = Atom::num(dimension / 2) - &eps;
            let c = (&h - Atom::one()) / (&mass * &scale);
            let backend = table(
                c.clone(),
                (&h - Atom::num(2)) * &c / (Atom::num(2) * &mass * &scale),
            );
            let prepared = PreparedPhysicalFamily::new(
                &family,
                &[Integral(vec![1])],
                &[s, t],
                &backend,
                &FlowOptions {
                    refine_basis: true,
                    dimension,
                    ..Default::default()
                },
                "fixed complex mass branch",
                &RunContext::default(),
            )
            .unwrap();
            assert_eq!(prepared.basis(), &[Integral(vec![2])]);
            let report = prepared.basis_refinement().unwrap();
            assert_eq!(report.basis_changes, 1);
            assert!(report.factorized);
            // I_a is proportional to c^-a (m²)^(D/2-a). Check its
            // independently derived logarithmic derivative in both variables.
            for variable in [s, t] {
                let expected = (&h - Atom::num(2)) * mass.derivative(variable) / &mass
                    - Atom::num(2) * scale.derivative(variable) / &scale;
                exact(
                    &prepared.flow().system().derivatives[&variable][0][0],
                    &expected,
                );
            }
            exact(
                &prepared.target_reductions()[0][&Integral(vec![2])],
                &(Atom::one() / &c),
            );
            exact(
                &prepared.basis_transformations()[0].matrix[0][0],
                &(Atom::one() / c),
            );
        }
    }
}

#[test]
fn redundant_retained_coordinates_preserve_the_physical_solution_after_refinement() {
    let (s, t) = symbol!("physical_options::s", "physical_options::t");
    let eps = Atom::var(symbol!("physical_options::eps"));
    let mass = Atom::var(s);
    let scale = &eps + Atom::var(t);
    let family = family(mass.clone(), scale.clone());
    let h = Atom::num(2) - &eps;
    let c = (&h - Atom::one()) / (&mass * &scale);
    let d = (&h - Atom::num(2)) * &c / (Atom::num(2) * &mass * &scale);
    let backend = table(c.clone(), d.clone());
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1]), Integral(vec![2])],
        &[s, t],
        &backend,
        &FlowOptions {
            skip_reduction: true,
            refine_basis: true,
            ..Default::default()
        },
        "dependent tadpole coordinates",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![3]), Integral(vec![2])]);
    // Factor out the common scalar tadpole value. These relative values and
    // its logarithmic derivative follow directly from its gamma formula.
    let relative = [d, c];
    for variable in [s, t] {
        let log_derivative = (&h - Atom::one()) * mass.derivative(variable) / &mass
            - scale.derivative(variable) / &scale;
        for (row, value) in prepared.flow().system().derivatives[&variable]
            .iter()
            .zip(&relative)
        {
            let action = row
                .iter()
                .zip(&relative)
                .fold(Atom::new(), |a, (m, y)| a + m * y);
            exact(
                &action,
                &(value.derivative(variable) + value * &log_derivative),
            );
        }
    }
}

#[test]
fn target_projection_respects_retained_basis_order_after_a_swap() {
    let s = symbol!("physical_options::s");
    let family = family(Atom::var(s), Atom::one());
    let eps = Atom::var(family.epsilon);
    // An explicit algebraic reduction fixture exercises the supplied-table
    // contract and projection. It is not a claim about this family's IBPs.
    let m = Atom::var(s) + &eps;
    let c = (Atom::one() - &eps) / &m;
    let backend = table(c.clone(), -&eps * &c / (Atom::num(2) * &m));
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1]), Integral(vec![2])],
        &[s],
        &backend,
        &FlowOptions {
            skip_reduction: true,
            refine_basis: true,
            ..Default::default()
        },
        "supplied algebraic fixture",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![3]), Integral(vec![2])]);
    let p = Precision::decimal(100).unwrap();
    let boundary = CachedBoundary {
        identity: prepared.flow().identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::num(3))])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(-2, 1).unwrap(),
        coefficients: vec![
            vec![p.zero(); 2],
            vec![p.zero(); 2],
            vec![p.i(2), p.i(7)],
            vec![p.zero(); 2],
        ],
        accuracy: BoundaryAccuracy::supplied(
            80,
            p.bits,
            vec![vec![p.real(0); 2]; 4],
            "exact polynomial values for a supplied algebraic projection fixture",
        )
        .unwrap(),
    };
    let projected = prepared
        .project_targets(&boundary, EpsilonRange::new(-1, 0).unwrap(), 20)
        .unwrap();
    assert!(p.close(&projected[0].coefficients[&-1], &p.i(-36), 50));
    assert!(p.close(&projected[0].coefficients[&0], &p.i(-60), 50));
    assert!(p.close(&projected[1].coefficients[&-1], &p.zero(), 50));
    assert!(p.close(&projected[1].coefficients[&0], &p.i(7), 50));
}
