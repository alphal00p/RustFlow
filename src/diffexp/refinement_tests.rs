use super::*;

// A real, exactly soluble transport keeps the supplied Precision on compiled
// arithmetic, diagnostics and retained segments, independently of FlowOptions.
fn constant_profile(
    p: Precision,
    options: &FlowOptions,
    value: C,
    leading: i32,
    save: bool,
) -> Result<EpsilonSolution> {
    EpsilonSystem {
        variable: symbol!("precision_recovery::x"),
        matrices: vec![vec![vec![Atom::zero()]]],
    }
    .compile(p, &Default::default())?
    .transport(
        &EpsilonBoundary {
            point: p.zero(),
            leading,
            coefficients: vec![vec![value]],
        },
        &[p.i(1)],
        options,
        &RunContext::default(),
        save,
    )
}

fn options() -> FlowOptions {
    FlowOptions {
        digits: 12,
        guard_digits: 48,
        series_order: 8,
        max_precision_attempts: 2,
        ..Default::default()
    }
}

#[test]
fn hinted_retry_keeps_local_budget_and_then_refines_after_success() -> Result<()> {
    let mut profiles = Vec::new();
    let minimum_bits = 451;
    let result = refine_epsilon_transport(
        &options(),
        &RunContext::default(),
        true,
        |p, local| {
            profiles.push((p.bits, local.digits, local.guard_digits, local.series_order));
            if p.bits < minimum_bits {
                return Err(Error::InsufficientPrecision {
                    minimum_bits,
                    context: "fixed local budget needs more arithmetic bits".into(),
                });
            }
            let delta = if profiles.len() == 2 { "1e-14" } else { "2e-14" };
            let value = p.add(&p.i(1), &p.parse(delta, "0")?);
            let result = constant_profile(p, local, value, 0, true)?;
            assert_eq!(result.diagnostics.working_bits, p.bits);
            assert!(result.segments.iter().all(|s| s.working_bits == p.bits));
            Ok(result)
        },
    )?;
    assert_eq!(profiles.len(), 3);
    assert_eq!(profiles.iter().map(|p| p.1).collect::<Vec<_>>(), [12, 12, 12]);
    assert_eq!(profiles.iter().map(|p| p.2).collect::<Vec<_>>(), [48, 48, 68]);
    assert_eq!(profiles.iter().map(|p| p.3).collect::<Vec<_>>(), [8, 40, 72]);
    assert_eq!(profiles[0].0, Precision::decimal(60)?.bits);
    assert!(profiles[1].0 >= minimum_bits && profiles[2].0 > profiles[1].0);
    // The retry adds arithmetic margin over the unchanged local profile.
    assert!(profiles[1].0 > Precision::decimal(profiles[1].1 + profiles[1].2)?.bits);
    assert!(profiles[2].0 > Precision::decimal(profiles[2].1 + profiles[2].2)?.bits);
    assert_eq!(result.verified_digits, Some(12));
    let p = Precision { bits: result.diagnostics.working_bits };
    assert!(result.comparison_errors[0][0] > p.real(0));
    assert!(result.comparison_errors[0][0] < p.tolerance(12));
    assert_eq!(result.checkpoints.len(), 1);
    assert_eq!(result.checkpoints[0].comparison_errors, result.comparison_errors);
    Ok(())
}

#[test]
fn failure_discards_previous_and_one_recovered_profile_is_not_verified() -> Result<()> {
    for max_precision_attempts in [2, 3] {
        let mut calls = 0;
        let mut guards = Vec::new();
        let opts = FlowOptions { max_precision_attempts, ..options() };
        let result = refine_epsilon_transport(&opts, &RunContext::default(), false, |p, local| {
            calls += 1;
            guards.push(local.guard_digits);
            if calls == 2 {
                return Err(Error::InsufficientPrecision {
                    minimum_bits: 451,
                    context: "intervening arithmetic failure".into(),
                });
            }
            // A stale first profile would also cause a spurious leading-power
            // mismatch. Only the two fresh leading=1 profiles may be compared.
            let leading = if calls == 1 { 0 } else { 1 };
            constant_profile(p, local, p.i(1), leading, false)
        });
        assert_eq!(calls, max_precision_attempts + 1);
        assert_eq!(guards, [48, 68, 68, 88][..calls]);
        if max_precision_attempts == 2 {
            assert!(matches!(result, Err(Error::Accuracy(_))));
        } else {
            let result = result?;
            assert_eq!(result.leading, 1);
            assert_eq!(result.verified_digits, Some(12));
            assert!(result.checkpoints.is_empty());
        }
    }
    Ok(())
}

#[test]
fn ordinary_success_still_requires_two_agreeing_independent_profiles() -> Result<()> {
    for max_precision_attempts in [1, 2] {
        let mut calls = 0;
        let opts = FlowOptions { max_precision_attempts, ..options() };
        let result = refine_epsilon_transport(&opts, &RunContext::default(), false, |p, local| {
            calls += 1;
            // With no arithmetic failures, retain the original +20/+32 profile.
            assert_eq!(p.bits, Precision::decimal(60 + 20 * (calls as u32 - 1))?.bits);
            assert_eq!(local.guard_digits, 48 + 20 * (calls as u32 - 1));
            assert_eq!(local.series_order, 8 + 32 * (calls - 1));
            constant_profile(p, local, p.i(if calls == 1 { 1 } else { 2 }), 0, false)
        });
        assert_eq!(calls, max_precision_attempts + 1);
        if max_precision_attempts == 1 {
            assert!(matches!(result, Err(Error::Accuracy(_))));
        } else {
            let result = result?;
            assert_eq!(result.verified_digits, Some(12));
            let p = Precision { bits: result.diagnostics.working_bits };
            assert_eq!(result.coefficients[0][0], p.i(2));
            assert_eq!(result.comparison_errors[0][0], p.real(0));
        }
    }
    Ok(())
}

#[test]
fn precision_exhaustion_retains_last_hint_and_never_verifies() {
    let mut calls = 0;
    let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |p, local| {
        calls += 1;
        assert_eq!(local.guard_digits, 48);
        Err(Error::InsufficientPrecision {
            minimum_bits: p.bits + 100,
            context: format!("unresolved profile {calls}"),
        })
    })
    .unwrap_err();
    assert_eq!(calls, 3);
    match error {
        Error::InsufficientPrecision { context, .. } => {
            assert!(context.contains("without two successful profiles"));
            assert!(context.ends_with("unresolved profile 3"));
        }
        other => panic!("unexpected error: {other}"),
    }
}

#[test]
fn ordinary_errors_and_changed_leading_power_are_still_rejected() {
    for error in [Error::Accuracy("local rejection".into()), Error::InvalidInput("bad source".into()), Error::Cancelled] {
        let expected = error.to_string();
        let mut failure = Some(error);
        let mut calls = 0;
        let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |_, _| {
            calls += 1;
            Err(failure.take().unwrap())
        })
        .unwrap_err();
        assert_eq!(error.to_string(), expected);
        assert_eq!(calls, 1);
    }
    let mut calls = 0;
    let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |p, local| {
        let leading = calls;
        calls += 1;
        constant_profile(p, local, p.i(1), leading, false)
    })
    .unwrap_err();
    assert_eq!(calls, 2);
    assert!(matches!(error, Error::InvalidInput(message) if message.contains("leading epsilon power")));
}
