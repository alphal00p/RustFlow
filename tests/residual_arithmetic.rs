use symbolica::prelude::*;
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon};
use symbolica_amflow::*;

fn atom(a: &str) -> Atom {
    Atom::parse(a, "residual_arithmetic_public", Default::default()).unwrap()
}
fn symbol(a: &str) -> Symbol {
    match atom(a).as_view() {
        AtomView::Var(a) => a.get_symbol(),
        _ => unreachable!(),
    }
}

#[test]
fn default_is_explicit_ball_and_opt_in_retains_complex_transport() -> Result<()> {
    let p = Precision::decimal(70)?;
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let source = DifferentialSystem {
        variable: symbol("x"),
        matrix: vec![vec![i + atom("1/(2+x)")]],
    }
    .compile(p, &Default::default())?;
    let boundary = BoundaryData {
        point: p.zero(),
        values: vec![p.i(1)],
    };
    let end = p.parse("0.2", "0.1")?;
    let options = FlowOptions {
        digits: 25,
        series_order: 64,
        ..Default::default()
    };
    assert_eq!(options.residual_arithmetic, ResidualArithmetic::Ball);
    let default = source.transport(
        &boundary,
        std::slice::from_ref(&end),
        &options,
        &RunContext::default(),
    )?;
    let expected = p.mul(
        &p.exp(&p.mul(&p.complex(0, 1), &end)),
        &p.scale(&p.add(&p.i(2), &end), 1, 2),
    );
    for mode in [
        ResidualArithmetic::Ball,
        ResidualArithmetic::AdaptiveInteger,
    ] {
        assert_eq!(mode.as_str().parse::<ResidualArithmetic>()?, mode);
        let o = FlowOptions {
            residual_arithmetic: mode,
            ..options.clone()
        };
        for coordinate in [LocalCoordinate::Identity, LocalCoordinate::BalancedMobius] {
            let o = FlowOptions {
                local_coordinate: coordinate,
                ..o.clone()
            };
            let value = source.transport(
                &boundary,
                std::slice::from_ref(&end),
                &o,
                &RunContext::default(),
            )?;
            assert!(p.close(&value.values[0], &expected, 25));
            if mode == ResidualArithmetic::Ball && coordinate == LocalCoordinate::Identity {
                assert_eq!(value.values, default.values);
                assert_eq!(value.diagnostics.steps, default.diagnostics.steps);
            }
            let context = RunContext::default();
            context.cancellation.cancel();
            assert!(matches!(
                source.transport(&boundary, std::slice::from_ref(&end), &o, &context),
                Err(Error::Cancelled)
            ));
        }
    }
    assert!(matches!(
        "integer".parse::<ResidualArithmetic>(),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn opt_in_epsilon_refinement_keeps_independent_accuracy_and_saved_segments() -> Result<()> {
    let system = EpsilonSystem::from_differential_system(
        &DifferentialSystem {
            variable: symbol("x"),
            matrix: vec![vec![atom("eps/(1+x)")]],
        },
        symbol("eps"),
        4,
    )?;
    for mode in [
        ResidualArithmetic::Ball,
        ResidualArithmetic::AdaptiveInteger,
    ] {
        let result = transport_epsilon(
            &system,
            |p| {
                Ok(EpsilonBoundary {
                    point: p.zero(),
                    leading: 0,
                    coefficients: vec![
                        vec![p.i(1)],
                        vec![p.zero()],
                        vec![p.zero()],
                        vec![p.zero()],
                        vec![p.zero()],
                    ],
                })
            },
            &[Atom::one()],
            &FlowOptions {
                digits: 25,
                series_order: 40,
                residual_arithmetic: mode,
                ..Default::default()
            },
            &RunContext::default(),
            true,
        )?;
        assert_eq!(result.verified_digits, Some(25));
        assert!(!result.checkpoints.is_empty());
        assert!(!result.segments.is_empty());
        let p = Precision {
            bits: result.diagnostics.working_bits,
        };
        let mut expected = p.i(1);
        let logarithm = p.log(&p.i(2));
        for (k, row) in result.coefficients.iter().enumerate() {
            if k > 0 {
                expected = p.scale(&p.mul(&expected, &logarithm), 1, k as i64);
            }
            assert!(p.close(&row[0], &expected, 25));
        }
        let saved = result.segments.last().unwrap();
        let local = saved.coordinate.local_point(p, &saved.center, &saved.end)?;
        let replay = saved.evaluate_local(&local)?.0;
        for (value, row) in replay.iter().zip(&result.coefficients) {
            assert!(p.close(value, &row[0], 25));
        }
    }
    Ok(())
}
