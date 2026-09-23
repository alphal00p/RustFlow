use symbolica::prelude::*;
use symbolica_amflow::*;
fn bubble() -> IntegralFamily {
    IntegralFamily {
        name: "ft_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-2)]],
        propagators: vec![
            Propagator {
                constant: Atom::num(-1),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(-5),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}
#[test]
fn gaussian_tensor_terminal_includes_shift_and_wick_contractions() {
    let mut family = bubble();
    family.propagators[1] = family.propagators[0].clone();
    family.propagators[0] = Propagator {
        constant: Atom::num(-3),
        scalar_products: vec![Atom::num(1), Atom::num(1)],
    };
    family.physical_propagators = 1;
    let eps = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let d = Rational::from((19, 5));
    let m = p.scale(&p.i(1), 5, 2);
    let tad = |n| vacuum::tadpole(n, &m, &d, p).unwrap();
    let first = gaussian::terminal(&family, &Integral(vec![3, -1]), &eps, p).unwrap();
    assert!(p.close(&first, &p.add(&tad(2), &tad(3)), 40));
    let second = gaussian::terminal(&family, &Integral(vec![3, -2]), &eps, p).unwrap();
    let expected = p.add(
        &p.add(&tad(1), &p.add(&p.scale(&tad(2), 2, 1), &tad(3))),
        &p.scale(&p.add(&tad(2), &p.mul(&m, &tad(3))), -10, 19),
    );
    assert!(p.close(&second, &expected, 40), "{second} != {expected}");
}
#[test]
fn recursive_ft_agrees_with_direct_bubble_parameter_system() {
    let family = bubble();
    let target = Integral(vec![1, 1]);
    let eps = Rational::from((1, 10));
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(60).unwrap();
    let direct = ft::evaluate_bubble(&family, &target, &eps, &options, &context).unwrap();
    for (skip_reduction, refine_basis) in
        [(false, false), (true, false), (false, true), (true, true)]
    {
        let options = FlowOptions {
            skip_reduction,
            refine_basis,
            ..options.clone()
        };
        let general = ft::FtEvaluator::new(&backend, &context)
            .evaluate(&family, &target, &eps, &options)
            .unwrap();
        assert!(p.close(&general, &direct, 20), "{general} != {direct}");
    }
}
#[test]
fn recursive_ft_two_loop_sunset_matches_schwinger_parameters() {
    let family = IntegralFamily {
        name: "ft_sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[1, 1], &[], Atom::new(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let eps = Rational::from((3, 4));
    let p = Precision::decimal(60).unwrap();
    let answer = ft::FtEvaluator::new(&backend, &context)
        .evaluate(&family, &Integral(vec![1, 1, 1]), &eps, &options)
        .unwrap();
    let gamma = |n, d| {
        p.gamma_real(&p.rational(&Rational::from((n, d))).re)
            .unwrap()
    };
    let expected = p.neg(&p.div(
        &p.mul(&gamma(1, 2), &p.powi(&gamma(3, 4), 2)),
        &p.scale(&gamma(3, 2), 1, 4),
    ));
    assert!(p.close(&answer, &expected, 20), "{answer} != {expected}");
}
