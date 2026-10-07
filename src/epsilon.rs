//! Numerical Laurent reconstruction from exact, nonzero regulator samples.
use crate::numeric::solve;
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct LaurentExpansion {
    pub coefficients: BTreeMap<i32, C>,
    /// Set only after independent precision/sample refinement. Relative digits
    /// for resolved nonzero coefficients; absolute digits near numerical zero.
    pub verified_digits: Option<u32>,
    pub working_bits: u32,
    pub samples: usize,
    pub validation_samples: usize,
    pub refinements: usize,
    /// Absolute coefficient changes between the last independent fits.
    pub comparison_errors: BTreeMap<i32, Float>,
}

/// Fit using all samples, returning powers `leading..=last`.
/// This operation alone does not certify accuracy; use independent refinements.
pub fn fit_epsilon(
    samples: &[Rational],
    values: &[C],
    leading: i32,
    last: i32,
    p: Precision,
) -> Result<LaurentExpansion> {
    if last < leading
        || samples.len() != values.len()
        || samples.len() < (last as i64 - leading as i64 + 1) as usize
    {
        return Err(Error::InvalidInput(
            "epsilon fit dimensions or power range".into(),
        ));
    }
    if samples.iter().any(|v| v.is_zero()) {
        return Err(Error::InvalidInput(
            "epsilon samples must be nonzero".into(),
        ));
    }
    let mut sorted = samples.to_vec();
    sorted.sort();
    sorted.dedup();
    if sorted.len() != samples.len() {
        return Err(Error::InvalidInput(
            "epsilon samples must be distinct".into(),
        ));
    }
    if values.iter().any(|v| !p.finite(v)) {
        return Err(Error::InvalidInput("nonfinite epsilon sample".into()));
    }
    let radius = samples
        .iter()
        .map(|r| p.norm(&p.rational(r)))
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let radius = C::new(radius, p.real(0));
    let n = samples.len();
    let mut matrix = vec![vec![p.zero(); n]; n];
    let mut rhs = Vec::new();
    for (i, (eps, value)) in samples.iter().zip(values).enumerate() {
        let eps = p.rational(eps);
        let t = p.div(&eps, &radius);
        let mut power = p.i(1);
        for entry in &mut matrix[i] {
            *entry = power.clone();
            power = p.mul(&power, &t);
        }
        rhs.push(p.mul(value, &p.powi(&eps, -(leading as i64))));
    }
    let fitted = solve(p, matrix, rhs)?;
    let mut coefficients = BTreeMap::new();
    for k in leading..=last {
        coefficients.insert(
            k,
            p.div(
                &fitted[(k - leading) as usize],
                &p.powi(&radius, (k - leading) as i64),
            ),
        );
    }
    Ok(LaurentExpansion {
        coefficients,
        verified_digits: None,
        working_bits: p.bits,
        samples: n,
        validation_samples: 0,
        refinements: 0,
        comparison_errors: BTreeMap::new(),
    })
}

pub fn epsilon_samples(count: usize, denominator: i64) -> Result<Vec<Rational>> {
    if count == 0 || count > 10000 || denominator <= 0 {
        return Err(Error::InvalidInput("epsilon sampling parameters".into()));
    }
    let total = denominator
        .checked_mul(count as i64)
        .ok_or_else(|| Error::InvalidInput("epsilon denominator overflow".into()))?;
    Ok((1..=count)
        .map(|i| Rational::from((i as i64, total)))
        .collect())
}

/// Balanced real nodes for automatic Laurent reconstruction. Both signs are
/// analytic regulator samples, not different physical i0 prescriptions. Avoid
/// zero and retain exact rational coordinates for reduction and independent fits.
#[cfg(feature = "automatic")]
pub(crate) fn symmetric_epsilon_samples(count: usize, denominator: i64) -> Result<Vec<Rational>> {
    if count == 0 || count > 10000 || denominator <= 0 {
        return Err(Error::InvalidInput("epsilon sampling parameters".into()));
    }
    let half = count.div_ceil(2) as i64;
    let total = denominator
        .checked_mul(half)
        .ok_or_else(|| Error::InvalidInput("epsilon denominator overflow".into()))?;
    Ok((0..count)
        .map(|index| {
            let magnitude = (index / 2 + 1) as i64;
            Rational::from((
                if index % 2 == 0 {
                    magnitude
                } else {
                    -magnitude
                },
                total,
            ))
        })
        .collect())
}
