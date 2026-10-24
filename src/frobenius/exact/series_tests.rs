use super::*;

#[test]
fn explicit_cell_caps_preserve_requested_admission_without_default_ceilings() {
    let defaults = ExactFrobeniusLimits::default();
    assert_eq!(
        (
            defaults.max_dimension,
            defaults.max_order,
            defaults.max_coefficient_bits,
            defaults.max_scalar_cells,
        ),
        (None, None, None, None)
    );
    for (cells, accepted) in [
        (0, false),
        (16_000_000, true),
        (16_000_001, true),
        (32_000_000, true),
        (32_000_001, true),
    ] {
        let limits = ExactFrobeniusLimits {
            max_dimension: Some(128),
            max_coefficient_bits: Some(1_048_576),
            max_scalar_cells: Some(cells),
            ..Default::default()
        };
        assert_eq!(limits.validate().is_ok(), accepted);
    }
    let context = RunContext::default();
    let old = ExactFrobeniusLimits {
        max_dimension: Some(128),
        max_coefficient_bits: Some(1_048_576),
        max_scalar_cells: Some(16_000_000),
        ..Default::default()
    };
    let larger = ExactFrobeniusLimits {
        max_scalar_cells: Some(32_000_000),
        ..old.clone()
    };
    // PR15: requested order7 plus the independently proved gauge span5.
    assert_eq!(127usize * 127 * 128 * (7 + 5 + 1), 26_838_656);
    let old_domain = ExactDomain {
        limits: &old,
        context: &context,
        elimination_mode: ExactEliminationMode::ObservedSparse,
    };
    let domain = ExactDomain {
        limits: &larger,
        context: &context,
        elimination_mode: ExactEliminationMode::ObservedSparse,
    };
    assert!(matches!(
        old_domain.preflight(127, 12),
        Err(Error::Limit(_))
    ));
    assert!(domain.preflight(127, 12).is_ok());
    assert!(matches!(domain.preflight(127, 15), Err(Error::Limit(_))));
    assert!(matches!(domain.preflight(129, 1), Err(Error::Limit(_))));
    let excess_bits = ExactFrobeniusLimits {
        max_coefficient_bits: Some(1_048_577),
        ..larger
    };
    assert!(excess_bits.validate().is_ok());
}

fn series(expression: &str) -> impl Iterator<Item = (Rational, Atom)> {
    let z = symbol!("exact_root_series_grid::z");
    Atom::parse(expression, "exact_root_series_grid", Default::default())
        .unwrap()
        .series(z, 0, 5)
        .unwrap()
        .terms()
        .map(|(power, coefficient)| (power, coefficient.clone()))
        .collect::<Vec<_>>()
        .into_iter()
}

#[test]
fn native_nonconstant_unit_root_grid_ignores_only_exact_zero_fractional_slots() {
    let terms = series("(1-z)^(-1/2)").collect::<Vec<_>>();
    assert!(
        terms
            .iter()
            .any(|(power, coefficient)| !power.is_integer() && coefficient.is_zero())
    );
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let coefficients = regular_root_coefficients(
        terms.iter().map(|(p, c)| (p.clone(), c)),
        4,
        &ExactDomain {
            limits: &limits,
            context: &context,
            elimination_mode: ExactEliminationMode::ConservativeDense,
        },
    )
    .unwrap();
    for (value, expected) in coefficients.iter().zip([
        Rational::one(),
        Rational::from((1, 2)),
        Rational::from((3, 8)),
        Rational::from((5, 16)),
        Rational::from((35, 128)),
    ]) {
        assert_eq!(*value, Gaussian::new(expected, Rational::zero()));
    }
}

#[test]
fn nonzero_fractional_power_is_not_discarded_as_a_ramification_slot() {
    let terms = series("z^(1/2)+(1-z)^(-1/2)").collect::<Vec<_>>();
    assert!(
        terms
            .iter()
            .any(|(power, coefficient)| !power.is_integer() && !coefficient.is_zero())
    );
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    assert!(matches!(
        regular_root_coefficients(
            terms.iter().map(|(p, c)| (p.clone(), c)),
            4,
            &ExactDomain {
                limits: &limits,
                context: &context,
                elimination_mode: ExactEliminationMode::ConservativeDense
            }
        ),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn disabled_height_estimates_preserve_compact_native_integer_power_series() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let domain = ExactDomain::dense(&limits, &context);
    let x = symbol!("exact_compact_height::x");
    let expression = (Atom::one() + Atom::var(x)).pow(i64::MAX);
    let coefficients = domain.series(&expression, x, 1).unwrap();
    assert_eq!(
        coefficients,
        vec![domain.integer(1), domain.integer(i64::MAX)]
    );
    let diagnostic = ExactFrobeniusLimits {
        max_coefficient_bits: Some(u64::MAX),
        ..limits.clone()
    };
    assert!(matches!(
        ExactDomain::dense(&diagnostic, &context).series(&expression, x, 1),
        Err(Error::Limit(_))
    ));
    let beyond_native = (Atom::one() + Atom::var(x)).pow(Atom::num(Integer::from(i64::MAX) + 1));
    assert!(
        matches!(domain.preflight_series(&beyond_native, x, 1), Err(Error::Limit(message)) if message.contains("i64 representation"))
    );
    let budget = crate::ResourceBudget::with_memory_monitor(
        crate::ResourceLimits {
            max_process_rss_bytes: Some(1),
            ..Default::default()
        },
        || Ok(0),
    );
    assert!(matches!(
        budget.run(|| domain.series(&expression, x, 1)),
        Err(Error::Resource(crate::ResourceLimitError::Memory { .. }))
    ));
}
