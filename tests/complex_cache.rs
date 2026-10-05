use symbolica::prelude::*;
use symbolica_amflow::{Error, Precision, Result, RunContext};
fn x() -> Symbol {
    symbol!("complex_path_probe::x")
}
fn parse(text: &str) -> Atom {
    Atom::parse(text, "complex_path_probe", Default::default()).unwrap()
}
use std::collections::BTreeMap;
use symbolica_amflow::algebraic::{
    AlgebraicKinematicSystem, AlgebraicSystem, RootSeed, SquareRoot,
};
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{DifferentialSystem, FlowOptions, Prescription, RustFlow};

fn symbols() -> (Symbol, Symbol, Symbol) {
    (
        x(),
        symbol!("complex_path_probe::eps"),
        symbol!("complex_path_probe::r"),
    )
}
fn root_germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(symbols().2, sheet)]),
    }
}
fn point(a: Atom, sheet: RootSheet) -> Result<CachedPoint> {
    CachedPoint::Exact(BTreeMap::from([(x(), a)])).with_root_germ(root_germ(sheet))
}
fn flow(radicand: Atom) -> Result<RustFlow<AlgebraicKinematicSystem>> {
    let (x, eps, r) = symbols();
    let system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([(
            x,
            vec![vec![
                Atom::var(eps) * radicand.derivative(x) / (2 * Atom::var(r)),
            ]],
        )]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand,
        }],
    };
    RustFlow::new_algebraic(
        system,
        &[parse("I")],
        &Atom::one(),
        Prescription::PlusI0,
        "explicit analytic exp(epsilon*(r-r_source)) branch on the registered root cover",
    )
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 35,
        series_order: 64,
        ..Default::default()
    }
}
fn boundary(
    flow: &RustFlow<AlgebraicKinematicSystem>,
    coordinate: Atom,
    sheet: RootSheet,
) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: point(coordinate, sheet)?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 3)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            45,
            p.bits,
            vec![vec![p.tolerance(60)]; 4],
            "exact unit analytic boundary with conservative numerical allowance",
        )?,
    })
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}
fn expected(
    p: Precision,
    delta: &symbolica_amflow::ComplexFloat,
) -> Vec<symbolica_amflow::ComplexFloat> {
    let mut out = vec![p.i(1)];
    for index in 1..4 {
        out.push(p.div(&p.mul(out.last().unwrap(), delta), &p.i(index)));
    }
    out
}
#[test]
fn complex_mass_cache_refines_and_restarts() -> Result<()> {
    let radicand = parse("x-1+𝑖/7");
    let flow = flow(radicand.clone())?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(&flow, Atom::new(), RootSheet::Principal)?)?;
    let p = Precision::decimal(80)?;
    let root0 = p.pow(
        &p.eval(&parse("-1+𝑖/7"), &Default::default())?,
        &p.rational(&Rational::from((1, 2))),
    );
    for end in [1, 2] {
        let result = flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(x(), Atom::num(end))]),
            &root_germ(RootSheet::Principal),
            EpsilonRange::new(0, 3)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        let root = p.pow(
            &p.eval(
                &symbolica_amflow::family::substitute(
                    &radicand,
                    &BTreeMap::from([(Atom::var(x()), Atom::num(end))]),
                ),
                &Default::default(),
            )?,
            &p.rational(&Rational::from((1, 2))),
        );
        for (actual, expected) in result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(expected(p, &p.sub(&root, &root0)))
        {
            assert!(p.close(actual, &expected, 28));
        }
        assert!(result.boundary.accuracy.verified_digits() >= 20);
        assert!(result.transport.is_some());
    }
    let directory =
        std::env::temp_dir().join(format!("rustflow-complex-cache-{}", std::process::id()));
    if directory.exists() {
        std::fs::remove_dir_all(&directory).unwrap();
    }
    cache.save(&directory)?;
    let mut cache = RustFlowCache::load(&directory)?;
    let hit = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(x(), Atom::num(2))]),
        &root_germ(RootSheet::Principal),
        EpsilonRange::new(0, 3)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(hit.transport.is_none());
    std::fs::remove_dir_all(directory).unwrap();
    Ok(())
}
#[test]
fn complex_straight_cut_crossing_changes_endpoint_germ() -> Result<()> {
    let flow = flow(Atom::var(x()))?;
    let mut cache = RustFlowCache::default();
    let start = parse("-1+𝑖");
    let end = parse("-1-𝑖");
    cache.insert(boundary(&flow, start.clone(), RootSheet::Principal)?)?;
    assert!(!flow.identity().conditions_admit_straight_path(
        &point(start.clone(), RootSheet::Principal)?,
        &point(end.clone(), RootSheet::Principal)?,
        Precision::decimal(60)?,
        20
    )?);
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(x(), end.clone())]),
        &root_germ(RootSheet::Opposite),
        EpsilonRange::new(0, 3)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let p = Precision::decimal(80)?;
    let half = p.rational(&Rational::from((1, 2)));
    let initial = p.pow(&p.eval(&start, &Default::default())?, &half);
    let final_root = p.neg(&p.pow(&p.eval(&end, &Default::default())?, &half));
    for (actual, expected) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(expected(p, &p.sub(&final_root, &initial)))
    {
        assert!(p.close(actual, &expected, 28));
    }
    assert_eq!(
        result.boundary.point.root_germ(),
        Some(&root_germ(RootSheet::Opposite))
    );
    Ok(())
}
#[test]
fn cache_multiple_complex_segments_wind_and_retain_sheets() -> Result<()> {
    let flow = flow(Atom::var(x()))?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(&flow, Atom::one(), RootSheet::Principal)?)?;
    let mut previous = point(Atom::one(), RootSheet::Principal)?;
    for (coordinate, sheet) in [
        (parse("𝑖"), RootSheet::Principal),
        (parse("-1"), RootSheet::Principal),
        (parse("-𝑖"), RootSheet::Opposite),
        (parse("1"), RootSheet::Opposite),
        (parse("𝑖"), RootSheet::Opposite),
        (parse("-1"), RootSheet::Opposite),
        (parse("-𝑖"), RootSheet::Principal),
        (parse("1"), RootSheet::Principal),
    ] {
        let selected = previous.clone();
        let strict_policy = ScaledDistance {
            scales: BTreeMap::new(),
            admissible: move |source: &CachedBoundary, _: &CachedPoint| {
                Ok(source.point.root_germ() == selected.root_germ()
                    && source.point.restart_coordinates()? == selected.restart_coordinates()?)
            },
        };
        let result = flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(x(), coordinate.clone())]),
            &root_germ(sheet),
            EpsilonRange::new(0, 3)?,
            &options(),
            &RunContext::default(),
            &strict_policy,
        )?;
        let p = Precision::decimal(80)?;
        let z = p.eval(&coordinate, &Default::default())?;
        let root = if z.im == p.real(0) && z.re < p.real(0) {
            p.complex(0, 1)
        } else {
            p.pow(&z, &p.rational(&Rational::from((1, 2))))
        };
        let root = if sheet == RootSheet::Opposite {
            p.neg(&root)
        } else {
            root
        };
        for (actual, expected) in result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(expected(p, &p.sub(&root, &p.i(1))))
        {
            assert!(p.close(actual, &expected, 27));
        }
        previous = result.boundary.point;
    }
    Ok(())
}
#[test]
fn exact_endpoint_principal_reference_ignores_rounding_cancellation() -> Result<()> {
    let (_, _, r) = symbols();
    let p = Precision::decimal(30)?;
    for (extra, sign) in [("", 1), ("+𝑖/10^100", 1), ("-𝑖/10^100", -1)] {
        let radicand = parse(&format!("-1+𝑖*(x^2/7+x/11-18/77){extra}"));
        let system = AlgebraicSystem::ordinary(
            DifferentialSystem {
                variable: x(),
                matrix: vec![vec![Atom::new()]],
            },
            vec![SquareRoot {
                symbol: r,
                radicand,
            }],
        );
        let compiled = system.compile(p)?;
        let state =
            compiled.branch_state_at(&p.i(1), &BTreeMap::from([(r, RootSeed::Principal)]))?;
        assert!(p.close(&state.roots[0].1, &p.complex(0, sign), 25));
    }
    Ok(())
}

#[test]
fn nearer_incompatible_sheet_is_skipped_before_integral_transport() -> Result<()> {
    let flow = flow(Atom::var(x()))?;
    let p = Precision::decimal(80)?;
    let half = p.rational(&Rational::from((1, 2)));
    let mut cache = RustFlowCache::default();
    for coordinate in [parse("-1+𝑖"), parse("2-𝑖")] {
        let mut data = boundary(&flow, coordinate.clone(), RootSheet::Principal)?;
        let root = p.pow(&p.eval(&coordinate, &Default::default())?, &half);
        data.coefficients = expected(p, &root).into_iter().map(|z| vec![z]).collect();
        cache.insert(data)?;
    }
    let target = parse("-1-𝑖");
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(x(), target.clone())]),
        &root_germ(RootSheet::Principal),
        EpsilonRange::new(0, 3)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert_eq!(
        result.starting_point.restart_coordinates()?[&x()],
        parse("2-𝑖")
    );
    let root = p.pow(&p.eval(&target, &Default::default())?, &half);
    for (actual, expected) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(expected(p, &root))
    {
        assert!(p.close(actual, &expected, 28));
    }
    Ok(())
}
#[test]
fn complex_admission_keeps_raw_guards_and_transactional_cancellation() -> Result<()> {
    let (x, eps, r) = symbols();
    let mut system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([(x, vec![vec![Atom::new()]])]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: parse("x+𝑖"),
        }],
    };
    let flow = RustFlow::with_algebraic_conditions(
        system.clone(),
        &[parse("I")],
        &Atom::one(),
        Prescription::PlusI0,
        "complex guarded zero DE",
        &[parse("x-1/2")],
    )?;
    assert!(!flow.identity().conditions_admit_straight_path(
        &point(Atom::new(), RootSheet::Principal)?,
        &point(Atom::one(), RootSheet::Principal)?,
        Precision::decimal(60)?,
        20
    )?);
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(&flow, Atom::new(), RootSheet::Principal)?)?;
    let initial = cache.len();
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(x, Atom::one())]),
            &root_germ(RootSheet::Principal),
            EpsilonRange::new(0, 3)?,
            &options(),
            &context,
            &policy()
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(cache.len(), initial);
    assert!(
        boundary(&flow, Atom::num((1, 2)), RootSheet::Principal)?
            .validate()
            .is_err()
    );
    system.roots[0].radicand = parse("x-1/2+𝑖*(x-1/2)");
    let flow = RustFlow::new_algebraic(
        system,
        &[parse("I")],
        &Atom::one(),
        Prescription::PlusI0,
        "complex branch zero",
    )?;
    assert!(!flow.identity().conditions_admit_straight_path(
        &point(Atom::new(), RootSheet::Principal)?,
        &point(Atom::one(), RootSheet::Principal)?,
        Precision::decimal(60)?,
        20
    )?);
    Ok(())
}

#[test]
fn canonical_complex_cut_cache_matches_letter_log_and_restarts() -> Result<()> {
    let (x, epsilon, r) = symbols();
    let system = symbolica_amflow::algebraic::CanonicalAlgebraicSystem::new(
        epsilon,
        &[x],
        &[Atom::one() + Atom::var(r)],
        &[vec![vec![Atom::one()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(x),
        }],
    )?;
    let flow = RustFlow::new_canonical(
        system,
        &[parse("I")],
        &Atom::one(),
        Prescription::PlusI0,
        "specified upper-root continuation; letter1+r does not cross its logarithm cut",
    )?;
    let p = Precision::decimal(80)?;
    let source = parse("-1+𝑖");
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: point(source.clone(), RootSheet::Principal)?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 3)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            45,
            p.bits,
            vec![vec![p.tolerance(60)]; 4],
            "unit analytic canonical boundary; independent conservative source allowance",
        )?,
    })?;
    let half = p.rational(&Rational::from((1, 2)));
    let root0 = p.pow(&p.eval(&source, &Default::default())?, &half);
    let mut last = None;
    for target in [parse("-1-𝑖"), parse("-2-𝑖")] {
        let result = flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(x, target.clone())]),
            Some(&root_germ(RootSheet::Opposite)),
            EpsilonRange::new(0, 3)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        let root = p.neg(&p.pow(&p.eval(&target, &Default::default())?, &half));
        let logarithm = p.log(&p.div(&p.add(&p.i(1), &root), &p.add(&p.i(1), &root0)));
        for (actual, expected) in result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(expected(p, &logarithm))
        {
            assert!(p.close(actual, &expected, 28));
        }
        assert!(result.transport.is_some());
        if let Some(previous) = last {
            assert_eq!(result.starting_point.restart_coordinates()?[&x], previous);
        }
        last = Some(target);
    }
    let directory = std::env::temp_dir().join(format!(
        "rustflow-canonical-complex-cache-{}",
        std::process::id()
    ));
    if directory.exists() {
        std::fs::remove_dir_all(&directory).unwrap();
    }
    cache.save(&directory)?;
    let mut cache = RustFlowCache::load(&directory)?;
    let repeated = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(x, parse("-2-𝑖"))]),
        Some(&root_germ(RootSheet::Opposite)),
        EpsilonRange::new(0, 3)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(repeated.transport.is_none());
    std::fs::remove_dir_all(directory).unwrap();
    Ok(())
}

#[test]
fn varying_complex_mass_coordinate_uses_both_partials_and_restarts() -> Result<()> {
    let (s, epsilon, r) = symbols();
    let mass = symbol!("complex_path_probe::m2");
    let derivative: Atom = Atom::var(epsilon) / (2 * Atom::var(r));
    let flow = RustFlow::new_algebraic(
        AlgebraicKinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![derivative.clone()]]),
                (mass, vec![vec![-derivative]]),
            ]),
            roots: vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(s) - Atom::var(mass),
            }],
        },
        &[parse("I")],
        &Atom::one(),
        Prescription::PlusI0,
        "upper-half radicand domain with independent complex mass and invariant coordinates",
    )?;
    let p = Precision::decimal(80)?;
    let source = BTreeMap::from([(s, Atom::new()), (mass, parse("1-𝑖/7"))]);
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(source.clone())
            .with_root_germ(root_germ(RootSheet::Principal))?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 3)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            45,
            p.bits,
            vec![vec![p.tolerance(60)]; 4],
            "exact unit boundary with explicit conservative allowance",
        )?,
    })?;
    let half = p.rational(&Rational::from((1, 2)));
    let root0 = p.pow(
        &p.eval(&(&source[&s] - &source[&mass]), &Default::default())?,
        &half,
    );
    let directory = std::env::temp_dir().join(format!(
        "rustflow-varying-complex-mass-{}",
        std::process::id()
    ));
    if directory.exists() {
        std::fs::remove_dir_all(&directory).unwrap();
    }
    for (index, destination) in [
        BTreeMap::from([(s, parse("1+𝑖/9")), (mass, parse("3/2-𝑖/5"))]),
        BTreeMap::from([(s, parse("2+𝑖/4")), (mass, parse("2-𝑖/3"))]),
    ]
    .into_iter()
    .enumerate()
    {
        let result = flow.evaluate_to(
            &mut cache,
            &destination,
            &root_germ(RootSheet::Principal),
            EpsilonRange::new(0, 3)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        let root = p.pow(
            &p.eval(
                &(&destination[&s] - &destination[&mass]),
                &Default::default(),
            )?,
            &half,
        );
        for (actual, expected) in result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(expected(p, &p.sub(&root, &root0)))
        {
            assert!(p.close(actual, &expected, 28));
        }
        assert_eq!(result.boundary.point.restart_coordinates()?, destination);
        if index == 0 {
            cache.save(&directory)?;
            cache = RustFlowCache::load(&directory)?;
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
    Ok(())
}
