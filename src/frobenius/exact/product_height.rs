//! A denominator-aware admission bound; all arithmetic and products stay native.
use super::{Gaussian, height};
use crate::Result;
use symbolica::prelude::Integer;

fn product_bits(left: &Integer, right: &Integer) -> u64 {
    let a = left.significant_bits();
    let b = right.significant_bits();
    if a == 0 || b == 0 {
        0
    } else if a == 1 {
        b // A signed unit.
    } else if b == 1 {
        a
    } else {
        a.saturating_add(b)
    }
}

fn bounded_lcm(
    common: &Integer,
    denominator: &Integer,
    cap: u64,
    check: &mut impl FnMut() -> Result<()>,
) -> Result<Option<Integer>> {
    check()?;
    let denominator = denominator.abs();
    if denominator.significant_bits() == 0 || denominator.significant_bits() > cap {
        return Ok(None);
    }
    let gcd = common.gcd(&denominator);
    check()?;
    let cofactor = common / &gcd;
    if product_bits(&cofactor, &denominator) > cap {
        return Ok(None);
    }
    check()?;
    let next = common.lcm(&denominator);
    Ok((next.significant_bits() <= cap).then_some(next))
}

struct Summary {
    denominator_bits: u64,
    numerator_bits: u64,
    zero: bool,
}

fn summary<'a>(
    values: impl Iterator<Item = &'a Gaussian> + Clone,
    cap: u64,
    check: &mut impl FnMut() -> Result<()>,
) -> Result<Option<Summary>> {
    let mut common = Integer::from(1);
    for value in values.clone() {
        for part in [&value.re, &value.im] {
            let Some(next) = bounded_lcm(&common, part.denominator_ref(), cap, check)? else {
                return Ok(None);
            };
            common = next;
        }
    }
    let mut numerator_bits = 0;
    let mut zero = true;
    for value in values {
        for part in [&value.re, &value.im] {
            check()?;
            if part.numerator_ref().significant_bits() == 0 {
                continue;
            }
            zero = false;
            let denominator = part.denominator_ref().abs();
            let scale = &common / &denominator;
            // Bound the scaled numerator without constructing its product.
            numerator_bits = numerator_bits.max(product_bits(part.numerator_ref(), &scale));
        }
    }
    Ok(Some(Summary {
        denominator_bits: common.significant_bits(),
        numerator_bits,
        zero,
    }))
}

/// Called only after the original dimension/work/cell checks and a failed
/// conservative estimate. At most one row/column common denominator is held.
///
/// If D,E clear a row and column, their Gaussian dot numerator contains at most
/// 2m signed integer products. Scaled component bit bounds a,b therefore give
/// numerator bits <= a+b+ceil(log2(2m)); the reduced denominator divides D*E.
/// Global row/column maxima remain conservative. The normal native Matrix
/// multiplication and its observed result-height check are never replaced.
pub(super) fn bound(
    left: &[Vec<Gaussian>],
    right: &[Vec<Gaussian>],
    cap: u64,
    check: &mut impl FnMut() -> Result<()>,
) -> Result<Option<u64>> {
    check()?;
    for value in left.iter().chain(right).flatten() {
        check()?;
        if height(value) > cap {
            return Ok(None);
        }
    }
    let mut left_numerator_bits = 0;
    let mut left_denominator_bits = 0;
    let mut left_zero = true;
    for row in left {
        let Some(value) = summary(row.iter(), cap, check)? else {
            return Ok(None);
        };
        left_numerator_bits = left_numerator_bits.max(value.numerator_bits);
        left_denominator_bits = left_denominator_bits.max(value.denominator_bits);
        left_zero &= value.zero;
    }
    let mut right_numerator_bits = 0;
    let mut right_denominator_bits = 0;
    let mut right_zero = true;
    for column in 0..right[0].len() {
        let Some(value) = summary(right.iter().map(|row| &row[column]), cap, check)? else {
            return Ok(None);
        };
        right_numerator_bits = right_numerator_bits.max(value.numerator_bits);
        right_denominator_bits = right_denominator_bits.max(value.denominator_bits);
        right_zero &= value.zero;
    }
    if left_zero || right_zero {
        return Ok(Some(1));
    }
    let sum_bits = u64::from(usize::BITS - (right.len() - 1).leading_zeros()) + 1;
    let bits = left_numerator_bits
        .saturating_add(right_numerator_bits)
        .saturating_add(sum_bits)
        .max(left_denominator_bits.saturating_add(right_denominator_bits));
    Ok((bits <= cap).then_some(bits))
}

#[cfg(test)]
#[path = "product_height_tests.rs"]
mod tests;
