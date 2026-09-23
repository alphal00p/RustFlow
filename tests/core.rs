use symbolica::prelude::*;
use symbolica_amflow::*;

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
