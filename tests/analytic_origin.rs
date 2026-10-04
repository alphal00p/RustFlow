use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, AnalyticOriginOptions, RootSeed, SquareRoot};
use symbolica_amflow::diffexp;
use symbolica_amflow::{Error, FlowOptions, Precision, Result, RunContext};
fn system(matrix: Vec<Vec<Atom>>, count: usize) -> AlgebraicSystem {
    let n = matrix.len();
    let mut matrices = vec![vec![vec![Atom::new(); n]; n]; count];
    matrices[1] = matrix;
    AlgebraicSystem {
        system: diffexp::EpsilonSystem {
            variable: symbol!("api_test::x"),
            matrices,
        },
        roots: vec![],
        nonzero_conditions: vec![Atom::var(symbol!("api_test::x"))],
    }
}
fn options() -> AnalyticOriginOptions {
    AnalyticOriginOptions {
        order: 32,
        check_digits: 30,
        ..Default::default()
    }
}
#[test]
fn residue_coupled_analytic_solution_matches_exact_geometric_epsilon_series() -> Result<()> {
    let x = Atom::var(symbol!("api_test::x"));
    let sys = system(
        vec![
            vec![Atom::one() / x, Atom::one()],
            vec![Atom::new(), Atom::new()],
        ],
        5,
    );
    let mut boundary = vec![vec![Atom::new(); 2]; 5];
    boundary[0][1] = Atom::one();
    for denominator in [8, 16] {
        let p = Precision::decimal(60)?;
        let offset = p.rational(&Rational::from((1, denominator)));
        let seed = sys.analytic_origin(
            &boundary,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &offset,
            &options(),
        )?;
        for k in 0..5 {
            assert!(p.close(
                &seed.boundary.coefficients[k][0],
                &if k == 0 { p.zero() } else { offset.clone() },
                50
            ));
            assert_eq!(
                seed.boundary.coefficients[k][1],
                if k == 0 { p.i(1) } else { p.zero() }
            );
        }
        let answer = sys.compile(p)?.transport(
            &seed.boundary,
            &[p.i(1)],
            &seed.seeds,
            &FlowOptions {
                digits: 30,
                guard_digits: 30,
                series_order: 32,
                ..Default::default()
            },
            &RunContext::default(),
            false,
        )?;
        for k in 1..5 {
            assert!(p.close(&answer.solution.coefficients[k][0], &p.i(1), 45));
        }
    }
    Ok(())
}
#[test]
fn rejects_incompatible_residue_and_nonpure_connections() -> Result<()> {
    let p = Precision::decimal(50)?;
    let x = Atom::var(symbol!("api_test::x"));
    let mut s = system(vec![vec![Atom::one() / x]], 2);
    let b = vec![vec![Atom::one()], vec![Atom::new()]];
    assert!(matches!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &p.rational(&Rational::from((1, 8))),
            &options()
        ),
        Err(Error::InvalidInput(_))
    ));
    s.system.matrices[0][0][0] = Atom::one();
    assert!(matches!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &p.i(1),
            &options()
        ),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}
#[test]
fn rejects_higher_poles_and_origin_branch_points() -> Result<()> {
    let p = Precision::decimal(50)?;
    let x = Atom::var(symbol!("api_test::x"));
    let mut s = system(vec![vec![Atom::one() / (&x * &x)]], 2);
    let b = vec![vec![Atom::new()]; 2];
    assert!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &p.rational(&Rational::from((1, 8))),
            &options()
        )
        .is_err()
    );
    s.system.matrices[1][0][0] = Atom::one();
    s.roots = vec![SquareRoot {
        symbol: symbol!("api_test::r"),
        radicand: x,
    }];
    assert!(matches!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &p.i(1),
            &options()
        ),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}
#[test]
fn sheet_specific_residue_cancellation_is_conservatively_rejected() -> Result<()> {
    let p = Precision::decimal(50)?;
    let x = Atom::var(symbol!("api_test::x"));
    let r = symbol!("api_test::r");
    let mut s = system(vec![vec![(Atom::var(r) - Atom::one()) / &x]], 2);
    s.roots = vec![SquareRoot {
        symbol: r,
        radicand: Atom::one() + x,
    }];
    let b = vec![vec![Atom::one()], vec![Atom::new()]];
    assert!(matches!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::from([(r, RootSeed::Principal)]),
            p,
            &p.rational(&Rational::from((1, 8))),
            &options()
        ),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}
#[test]
fn validates_exact_constant_boundary_and_options() -> Result<()> {
    let p = Precision::decimal(50)?;
    let x = Atom::var(symbol!("api_test::x"));
    let s = system(vec![vec![Atom::new()]], 2);
    let offset = p.rational(&Rational::from((1, 8)));
    assert!(
        s.analytic_origin(
            &[vec![Atom::one()]],
            &Default::default(),
            &BTreeMap::new(),
            p,
            &offset,
            &options()
        )
        .is_err()
    );
    for a in [x, Atom::num(p.real(1))] {
        assert!(
            s.analytic_origin(
                &[vec![a], vec![Atom::new()]],
                &Default::default(),
                &BTreeMap::new(),
                p,
                &offset,
                &options()
            )
            .is_err()
        );
    }
    let b = vec![vec![Atom::one()], vec![Atom::new()]];
    let constants =
        ahash::HashMap::from_iter([(Atom::var(symbol!("api_test::unused")), p.parse("NaN", "0")?)]);
    assert!(
        s.analytic_origin(&b, &constants, &BTreeMap::new(), p, &offset, &options())
            .is_err()
    );
    assert!(
        s.analytic_origin(
            &b,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &offset,
            &AnalyticOriginOptions {
                order: 32,
                check_digits: 0,
                ..Default::default()
            }
        )
        .is_err()
    );
    Ok(())
}

#[test]
fn cancellation_is_typed_and_shared_with_the_callers_run_context() -> Result<()> {
    let context = RunContext::default();
    let options = AnalyticOriginOptions {
        cancellation: context.cancellation.clone(),
        ..Default::default()
    };
    let cloned = options.clone();
    let independent = AnalyticOriginOptions::default();
    context.cancellation.cancel();
    assert!(matches!(cloned.cancellation.check(), Err(Error::Cancelled)));
    assert!(independent.cancellation.check().is_ok());
    let p = Precision::decimal(50)?;
    let s = system(vec![vec![Atom::new()]], 2);
    let boundary = vec![vec![Atom::one()], vec![Atom::new()]];
    assert!(matches!(
        s.analytic_origin(
            &boundary,
            &Default::default(),
            &BTreeMap::new(),
            p,
            &p.rational(&Rational::from((1, 8))),
            &options
        ),
        Err(Error::Cancelled)
    ));
    Ok(())
}
