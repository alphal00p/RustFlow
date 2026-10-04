use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn scalar_substitution_keeps_general_expression_rules() {
    use std::collections::BTreeMap;
    let expression = parse!("1+x^2+y");
    assert_eq!(
        family::substitute(&expression, &BTreeMap::new()),
        expression
    );
    let scalar = BTreeMap::from([(parse!("x"), parse!("z+1")), (parse!("y"), Atom::num(3))]);
    assert!(
        (family::substitute(&expression, &scalar) - parse!("4+(z+1)^2"))
            .together()
            .cancel()
            .is_zero()
    );
    let composite = BTreeMap::from([(parse!("x^2"), Atom::num(7)), (parse!("y"), Atom::num(3))]);
    assert_eq!(family::substitute(&expression, &composite), Atom::num(11));
}

#[test]
fn poles_of_large_integer_region_polynomial_stabilize_with_precision() {
    let denominator = Atom::parse(
        include_str!("../fixtures/regressions/region-pole-degree12.txt"),
        "pole_regression",
        Default::default(),
    )
    .unwrap();
    let system = DifferentialSystem {
        variable: symbol!("pole_regression::eta"),
        matrix: vec![vec![Atom::num(1) / denominator]],
    };
    let low = system
        .compile(Precision::decimal(60).unwrap(), &Default::default())
        .unwrap();
    let p = Precision::decimal(90).unwrap();
    let high = system.compile(p, &Default::default()).unwrap();
    assert_eq!(low.poles.len(), 12);
    assert_eq!(high.poles.len(), 12);
    for root in &high.poles {
        assert!(
            low.poles
                .iter()
                .any(|old| p.norm(&p.sub(root, old)) <= p.tolerance(40) * p.norm(root))
        );
    }
}

#[test]
fn pole_scaling_preserves_very_large_and_very_small_roots() {
    let p = Precision::decimal(60).unwrap();
    let x = symbol!("scaled_pole_x");
    for exponent in [-100, 100] {
        let radius = Atom::num(10).pow(exponent);
        let system = DifferentialSystem {
            variable: x,
            matrix: vec![vec![
                Atom::num(1) / (Atom::var(x).pow(2) - Atom::num(2) * radius.pow(2)),
            ]],
        };
        let roots = system.compile(p, &Default::default()).unwrap().poles;
        assert_eq!(roots.len(), 2);
        let expected = p.mul(
            &p.powi(&p.i(10), exponent),
            &ComplexFloat::new(p.real(2).sqrt(), p.real(0)),
        );
        for root in roots {
            let positive_error = p.norm(&p.sub(&root, &expected));
            let negative_error = p.norm(&p.add(&root, &expected));
            let error = if positive_error < negative_error {
                positive_error
            } else {
                negative_error
            };
            assert!(error <= p.tolerance(40) * p.norm(&expected));
        }
    }
}

fn tadpole() -> IntegralFamily {
    IntegralFamily {
        name: "tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: parse!("-m2"),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}

#[test]
fn exact_tadpole_reduction_and_derivative() {
    let family = tadpole();
    let result = RustRedBackend::default()
        .reduce(&family, &[Integral(vec![2])], &RunContext::default())
        .unwrap();
    let terms = result.expand(&Integral(vec![2])).unwrap();
    assert_eq!(terms.len(), 1);
    assert!(
        (&terms[&Integral(vec![1])] - parse!("(1-eps)/m2"))
            .together()
            .cancel()
            .is_zero()
    );
    let derivatives = family::eta_derivative(&Integral(vec![2, -1]), &[true, true]).unwrap();
    assert_eq!(
        derivatives,
        vec![
            (Integral(vec![3, -1]), Atom::num(2)),
            (Integral(vec![2, 0]), Atom::num(-1))
        ]
    );
}

#[test]
fn sampled_preparation_keeps_dimension_symbolic_for_basis_refinement() {
    let mut family = tadpole();
    let eps = Atom::var(family.epsilon);
    // This exact epsilon-dependent mass exposes a genuinely mixed denominator
    // in the raised tadpole's reduction, making an omitted swap observable.
    family.propagators[0].constant = -eps.clone();
    let target = Integral(vec![2]);
    let options = FlowOptions {
        refine_basis: true,
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let flow = PreparedFlow::new_at_epsilon(
        &family,
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
        &Rational::from((1, 100)),
    )
    .unwrap();
    assert_eq!(flow.basis_refinement.as_ref().unwrap().basis_changes, 1);
    assert!(flow.basis_refinement.as_ref().unwrap().factorized);
    assert_eq!(flow.reduced.basis, vec![target.clone()]);
    assert_eq!(flow.reduced.targets[0][&target], Atom::num(1));
    let mass = Atom::var(flow.system.variable) + &eps;
    assert!(
        (&flow.system.matrix[0][0] + &eps / &mass)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        (&flow.reduced.transformations[0].matrix[0][0] - &mass / (Atom::num(1) - &eps))
            .together()
            .cancel()
            .is_zero()
    );
    let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60).unwrap();
    for epsilon in [Rational::from((1, 100)), Rational::from((1, 200))] {
        let value = flow
            .evaluate(&epsilon, &options, &boundary, &context)
            .unwrap();
        let e = p.rational(&epsilon);
        let expected = p.mul(&p.gamma_real(&e.re).unwrap(), &p.pow(&e, &p.neg(&e)));
        assert!(p.close(&value[0], &expected, 20));
    }
}

#[test]
fn automatic_refinement_takes_precedence_over_sampled_reduction() {
    struct SymbolicOnly;
    impl ReductionBackend for SymbolicOnly {
        fn identity(&self) -> String {
            "refinement-symbolic-only".into()
        }
        fn reduce(
            &self,
            family: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<reduction::Reduction> {
            RustRedBackend::default().reduce(family, targets, context)
        }
        fn reduce_at_epsilon(
            &self,
            _: &IntegralFamily,
            _: &[Integral],
            _: &Rational,
            _: &RunContext,
        ) -> Result<reduction::Reduction> {
            Err(Error::Reduction(
                "dimension was specialized before refinement".into(),
            ))
        }
    }
    let family = IntegralFamily {
        name: "refined_product".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::num(2), &[]).unwrap(),
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::new(), Atom::num(1), Atom::new()],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("refined_product_eps"),
        dimension: 4,
    };
    let options = FlowOptions {
        refine_basis: true,
        digits: 10,
        series_order: 48,
        ..Default::default()
    };
    let values = solve_integrals(
        &family,
        &[Integral(vec![1, 1, 0])],
        &KinematicPoint::default(),
        0,
        &options,
        &SymbolicOnly,
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(values[0].verified_digits, Some(10));
    let p = Precision::decimal(50).unwrap();
    assert!(p.close(&values[0].coefficients[&-2], &p.i(2), 10));
}

#[test]
fn taylor_transport_exponential_and_branch() {
    let p = Precision::decimal(70).unwrap();
    let system = DifferentialSystem {
        variable: symbol!("x"),
        matrix: vec![vec![Atom::num(1)]],
    };
    let compiled = system.compile(p, &ahash::HashMap::default()).unwrap();
    let result = compiled
        .transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.i(1)],
            &FlowOptions::default(),
            &RunContext::default(),
        )
        .unwrap();
    assert!(p.close(&result.values[0], &p.exp(&p.i(1)), 28));
    let system = DifferentialSystem {
        variable: symbol!("x"),
        matrix: vec![vec![parse!("1/(2*x)")]],
    };
    let compiled = system.compile(p, &ahash::HashMap::default()).unwrap();
    let result = compiled
        .transport(
            &BoundaryData {
                point: p.i(1),
                values: vec![p.i(1)],
            },
            &[p.complex(0, -1), p.i(-1)],
            &FlowOptions::default(),
            &RunContext::default(),
        )
        .unwrap();
    assert!(p.close(&result.values[0], &p.complex(0, -1), 25));
}

#[test]
fn block_order_and_incomplete_table() {
    let sys = DifferentialSystem {
        variable: symbol!("x"),
        matrix: vec![
            vec![Atom::num(1), Atom::num(1)],
            vec![Atom::new(), Atom::num(2)],
        ],
    };
    assert_eq!(sys.blocks().unwrap(), vec![vec![1], vec![0]]);
    let table = TableBackend {
        name: "empty".into(),
        reduction: Default::default(),
    };
    assert!(matches!(
        table.reduce(&tadpole(), &[Integral(vec![1])], &RunContext::default()),
        Err(Error::IncompleteReduction(_))
    ));
}

#[test]
fn automatic_contour_tracks_the_selected_branch() {
    let x = symbol!("path_x");
    let z = Atom::var(x);
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::num(1) / (Atom::num(2) * (z - Atom::num(1)))]],
    };
    let p = Precision::decimal(60).unwrap();
    let compiled = system.compile(p, &ahash::HashMap::default()).unwrap();
    for side in [-1, 1] {
        let path = compiled.plan_path(&p.zero(), &p.i(2), side).unwrap();
        assert!(path.len() > 1);
        let result = compiled
            .transport(
                &BoundaryData {
                    point: p.zero(),
                    values: vec![p.i(1)],
                },
                &path,
                &FlowOptions::default(),
                &RunContext::default(),
            )
            .unwrap();
        assert!(p.close(&result.values[0], &p.complex(0, -side), 20));
    }
}

#[test]
fn unresolved_search_residuals_do_not_bypass_closure_and_cancellation() {
    struct Residuals;
    impl ReductionBackend for Residuals {
        fn identity(&self) -> String {
            "residual-only".into()
        }
        fn reduce(
            &self,
            _: &IntegralFamily,
            targets: &[Integral],
            _: &RunContext,
        ) -> Result<reduction::Reduction> {
            Ok(reduction::Reduction {
                residuals: targets.to_vec(),
                ..Default::default()
            })
        }
    }
    let family = tadpole();
    let (deformed, mask) = family.deform(symbol!("eta"), &MassMode::All).unwrap();
    assert!(matches!(
        reduction::differential_system(
            &Residuals,
            &deformed,
            &[Integral(vec![1])],
            &mask,
            3,
            &RunContext::default()
        ),
        Err(Error::IncompleteReduction(_))
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        PreparedFlow::new(
            &family,
            &[Integral(vec![1])],
            &KinematicPoint::default(),
            &RustRedBackend::default(),
            &FlowOptions::default(),
            &context
        ),
        Err(Error::Cancelled)
    ));
}

#[test]
fn scaleless_tadpole_is_removed_by_native_sector_analysis() {
    let mut family = tadpole();
    family.propagators[0].constant = Atom::new();
    let targets = [Integral(vec![1]), Integral(vec![2])];
    for factorized in [false, true] {
        let reduction = RustRedBackend {
            factorized,
            ..Default::default()
        }
        .reduce(&family, &targets, &RunContext::default())
        .unwrap();
        for target in &targets {
            assert!(reduction.expand(target).unwrap().is_empty());
        }
    }
}

#[test]
fn precision_refinement_also_tightens_taylor_truncation() {
    let p = Precision::decimal(80).unwrap();
    let system = DifferentialSystem {
        variable: symbol!("guard_x"),
        matrix: vec![vec![parse!("1/(2*guard_x)")]],
    }
    .compile(p, &ahash::HashMap::default())
    .unwrap();
    let expected = p.pow(&p.i(2), &p.rational(&Rational::from((1, 2))));
    let run = |guard_digits| {
        let options = FlowOptions {
            digits: 10,
            guard_digits,
            series_order: 40,
            ..Default::default()
        };
        let result = system
            .transport(
                &BoundaryData {
                    point: p.i(1),
                    values: vec![p.i(1)],
                },
                &[p.i(2)],
                &options,
                &RunContext::default(),
            )
            .unwrap();
        p.norm(&p.sub(&result.values[0], &expected))
    };
    let low = run(30);
    let high = run(60);
    assert!(high < p.tolerance(55));
    assert!(high < low * p.tolerance(15));
}

#[test]
fn sparse_ode_does_not_mistake_a_zero_taylor_tail_for_convergence() {
    let p = Precision::decimal(60).unwrap();
    let system = DifferentialSystem {
        variable: symbol!("sparse_x"),
        matrix: vec![vec![parse!("sparse_x^20")]],
    }
    .compile(p, &ahash::HashMap::default())
    .unwrap();
    let options = FlowOptions {
        digits: 20,
        guard_digits: 20,
        series_order: 16,
        max_steps: 10000,
        ..Default::default()
    };
    let result = system
        .transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.i(1)],
            &options,
            &RunContext::default(),
        )
        .unwrap();
    assert!(p.close(&result.values[0], &p.exp(&p.scale(&p.i(1), 1, 21)), 20));
    assert!(result.diagnostics.rejected_steps > 0);
    // Retain useful entire-system step proposals instead of retrying the
    // complete remaining interval at every center.
    assert!(result.diagnostics.steps + result.diagnostics.rejected_steps < 5000);
}
