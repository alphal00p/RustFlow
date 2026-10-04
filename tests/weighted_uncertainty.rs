use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::diffexp::EpsilonSystem;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

#[test]
fn weighted_bound_accounts_for_physical_and_epsilon_indices() {
    let p = Precision::decimal(70).unwrap();
    let system = EpsilonSystem {
        variable: symbol!("weighted_bound::x"),
        matrices: vec![
            vec![
                vec![Atom::num(1), Atom::new()],
                vec![Atom::new(), Atom::num(2)],
            ],
            vec![
                vec![Atom::new(), Atom::new()],
                vec![Atom::num(3), Atom::new()],
            ],
        ],
    }
    .compile(p, &Default::default())
    .unwrap();
    let end = p.rational(&Rational::from((1, 100)));
    let weights = vec![p.real(100), p.real(1), p.real(1), p.real(2)];
    let weighted = system
        .error_amplification_weighted(&p.zero(), &end, &weights)
        .unwrap();
    // The largest augmented row is dY_(1,1)=2Y_(1,1)+3Y_(0,0).
    // Its weighted row sum is 2 + 3*100/2 = 152.
    assert!(p.close(
        &ComplexFloat::new(weighted, p.real(0)),
        &p.exp(&p.rational(&Rational::from((152, 100)))),
        60,
    ));
    let unweighted = system.error_amplification(&p.zero(), &end).unwrap();
    let unit_weights = system
        .error_amplification_weighted(&p.zero(), &end, &vec![p.real(1); 4])
        .unwrap();
    assert_eq!(unweighted, unit_weights);
    assert!(p.close(
        &ComplexFloat::new(unweighted, p.real(0)),
        &p.exp(&p.rational(&Rational::from((5, 100)))),
        60,
    ));
    assert!(matches!(
        system.error_amplification_weighted(&p.zero(), &end, &weights[..3]),
        Err(Error::InvalidInput(_))
    ));
    let mut invalid = weights;
    invalid[2] = p.real(0);
    assert!(matches!(
        system.error_amplification_weighted(&p.zero(), &end, &invalid),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn large_higher_coefficient_preserves_lower_order_accuracy_and_cache_reuse() {
    let epsilon = symbol!("weighted_cache::eps");
    let s = symbol!("weighted_cache::s");
    let engine = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon)]])]),
        },
        &[parse!("weighted_cache::I")],
        &Atom::num(1),
        Prescription::PlusI0,
        "entire scalar exponential",
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    let range = EpsilonRange::new(0, 4).unwrap();
    let mut cache = RustFlowCache::default();
    cache
        .insert(CachedBoundary {
            identity: engine.identity().clone(),
            point: CachedPoint::Exact(BTreeMap::from([(s, Atom::new())])),
            kind: PointKind::Physical,
            range,
            coefficients: vec![
                vec![p.zero()],
                vec![p.i(1)],
                vec![p.zero()],
                vec![p.zero()],
                vec![p.i(100_000_000)],
            ],
            accuracy: BoundaryAccuracy::supplied(
                24,
                p.bits,
                vec![vec![p.tolerance(30)]; 5],
                "24-digit boundary evidence, independent of its 80 working digits",
            )
            .unwrap(),
        })
        .unwrap();
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let options = FlowOptions::default();
    let first = engine
        .evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(1))]),
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert_eq!(first.boundary.accuracy.verified_digits(), 21);
    let second = engine
        .evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(2))]),
            range,
            &options,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    assert_eq!(second.boundary.accuracy.verified_digits(), 20);
    assert!(
        (&second.starting_point.restart_coordinates().unwrap()[&s] - Atom::num(1))
            .together()
            .cancel()
            .is_zero()
    );
    for (x, result) in [(1, first), (2, second)] {
        let p = Precision {
            bits: result.boundary.accuracy.working_bits(),
        };
        // Y=(epsilon+10^8 epsilon^4) exp(s epsilon), through epsilon^4.
        let expected = [
            p.zero(),
            p.i(1),
            p.i(x),
            p.scale(&p.i(x * x), 1, 2),
            p.add(&p.i(100_000_000), &p.scale(&p.i(x * x * x), 1, 6)),
        ];
        for (actual, expected) in result.boundary.coefficients.iter().flatten().zip(expected) {
            assert!(p.close(actual, &expected, 40));
        }
        let errors = result.boundary.accuracy.comparison_errors();
        assert!(errors[0][0] < p.tolerance(20));
        assert!(errors[4][0] > errors[0][0].clone() * p.real(10_000_000));
    }
}
