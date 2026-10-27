use super::*;
use crate::frobenius::exact::{ExactFrobeniusLimits, multiply};
use crate::{Error, RunContext};
use symbolica::prelude::*;

fn g(value: i64) -> Gaussian {
    Gaussian::new(Rational::from(value), Rational::from(0))
}

#[test]
fn shared_large_denominators_preserve_the_original_coefficient_ceiling() {
    let context = RunContext::default();
    let limits = ExactFrobeniusLimits {
        max_dimension: 127,
        max_coefficient_bits: 1_048_576,
        ..Default::default()
    };
    let q = Integer::from(2).pow(8192);
    let x = Rational::from((Integer::from(1), q));
    let left = vec![vec![Gaussian::new(x.clone(), x.clone()); 127]];
    let right = vec![vec![Gaussian::new(x.clone(), -x)]; 127];
    let result = multiply(left, right, &limits, &context).unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].len(), 1);
    assert_eq!(
        result[0][0],
        Gaussian::new(
            Rational::from((Integer::from(254), Integer::from(2).pow(16384))),
            Rational::from(0),
        )
    );
    assert!(height(&result[0][0]) <= limits.max_coefficient_bits);
}

#[test]
fn cheap_integer_product_keeps_native_exact_results() {
    let result = multiply(
        vec![vec![g(1), g(2)], vec![g(3), g(4)]],
        vec![vec![g(5), g(6)], vec![g(7), g(8)]],
        &ExactFrobeniusLimits::default(),
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(result, vec![vec![g(19), g(22)], vec![g(43), g(50)]]);
}

#[test]
fn denominator_fallback_handles_complex_signs_and_exact_zeros() {
    let value = |re, im| Gaussian::new(Rational::from((re, 257i64)), Rational::from((im, 257i64)));
    let result = multiply(
        vec![vec![value(1, 1), value(2, 0)], vec![g(0), value(1, -1)]],
        vec![vec![value(1, -1), g(0)], vec![value(3, 0), value(1, 1)]],
        &ExactFrobeniusLimits {
            max_coefficient_bits: 64,
            ..Default::default()
        },
        &RunContext::default(),
    )
    .unwrap();
    let expected = |re, im| {
        Gaussian::new(
            Rational::from((re, 66049i64)),
            Rational::from((im, 66049i64)),
        )
    };
    assert_eq!(
        result,
        vec![
            vec![expected(8, 0), expected(2, 2)],
            vec![expected(3, -3), expected(2, 0)]
        ]
    );
}

#[test]
fn genuinely_oversized_product_remains_refused() {
    let value = Gaussian::new(
        Rational::from(Integer::from(2).pow(524_288)),
        Rational::from(0),
    );
    assert!(matches!(
        multiply(
            vec![vec![value.clone()]],
            vec![vec![value]],
            &ExactFrobeniusLimits {
                max_coefficient_bits: 1_048_576,
                ..Default::default()
            },
            &RunContext::default(),
        ),
        Err(Error::Limit(_))
    ));
}

#[test]
fn lcm_growth_is_refused_before_crossing_the_same_cap() {
    let common = Integer::from(257);
    let denominator = Integer::from(263);
    assert!(common.lcm(&denominator).significant_bits() > 16);
    assert!(
        bounded_lcm(&common, &denominator, 16, &mut || Ok(()))
            .unwrap()
            .is_none()
    );
    let value = |denominator| Gaussian::new(Rational::from((1, denominator)), Rational::from(0));
    assert!(matches!(
        multiply(
            vec![vec![value(257), value(263)]],
            vec![vec![g(1)], vec![g(1)]],
            &ExactFrobeniusLimits {
                max_coefficient_bits: 16,
                ..Default::default()
            },
            &RunContext::default(),
        ),
        Err(Error::Limit(_))
    ));
}

#[test]
fn dimensions_work_and_cells_are_still_checked_before_height_fallback() {
    let context = RunContext::default();
    assert!(matches!(
        multiply(
            vec![vec![g(1)]],
            vec![],
            &ExactFrobeniusLimits::default(),
            &context
        ),
        Err(Error::InvalidInput(_))
    ));
    let small = ExactFrobeniusLimits {
        max_scalar_cells: 1,
        ..Default::default()
    };
    assert!(matches!(
        multiply(vec![vec![g(1); 17]], vec![vec![g(1)]; 17], &small, &context),
        Err(Error::Limit(_))
    ));
    let small = ExactFrobeniusLimits {
        max_scalar_cells: 3,
        ..Default::default()
    };
    assert!(matches!(
        multiply(
            vec![vec![g(1)], vec![g(1)]],
            vec![vec![g(1), g(1)]],
            &small,
            &context
        ),
        Err(Error::Limit(_))
    ));
}

#[test]
fn cancellation_propagates_at_entry_and_during_new_denominator_work() {
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        multiply(
            vec![vec![g(1)]],
            vec![vec![g(1)]],
            &ExactFrobeniusLimits::default(),
            &context
        ),
        Err(Error::Cancelled)
    ));
    let context = RunContext::default();
    let mut checkpoints = 0;
    let mut check = || {
        checkpoints += 1;
        if checkpoints > 10 {
            context.cancellation.cancel();
        }
        context.cancellation.check()
    };
    assert!(matches!(
        bound(
            &[vec![g(1), g(2)]],
            &[vec![g(3)], vec![g(4)]],
            64,
            &mut check
        ),
        Err(Error::Cancelled)
    ));
    assert!(checkpoints > 10);
}

#[test]
fn an_exact_zero_product_is_not_a_small_residual_guess() {
    let result = multiply(
        vec![vec![g(0); 127]],
        vec![vec![g(17)]; 127],
        &ExactFrobeniusLimits {
            max_coefficient_bits: 64,
            ..Default::default()
        },
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(result, vec![vec![g(0)]]);
}
