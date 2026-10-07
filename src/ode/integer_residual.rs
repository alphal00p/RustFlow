//! Bounded native integer products for exact stored-dyadic source defects.
//!
//! The old ball constructor remains the general fallback. Integer admission
//! uses resource estimates, not hard allocator quotas. A failed error-budget
//! check consults that same chart's lazily retained ball calculation so tighter
//! exact arithmetic cannot erase its useful working-precision retry signal.
use super::source::{
    self, BallPolynomial, DenominatorEnclosure, ExactPolynomialRow, ExactSourceResidual, Gaussian,
    exact_ball,
};
use crate::{CancellationToken, ComplexFloat as C, Error, Precision, Result};
use std::sync::{Arc, Mutex};
use symbolica::domains::float::{FloatField, RoundingDirection};
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

type ExactPolynomial = UnivariatePolynomial<FloatField<Gaussian>>;
type IntegerPolynomial = UnivariatePolynomial<Z>;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Limits {
    input_bits: u64,
    integer_bits: u64,
    work: u64,
    bytes: u64,
    shifted_degree: usize,
    degree: usize,
    components: usize,
    order: usize,
}

#[cfg(test)]
#[path = "integer_residual_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "integer_residual_bench.rs"]
mod bench;
impl Default for Limits {
    fn default() -> Self {
        Self {
            input_bits: 4096,
            integer_bits: 131072,
            work: 100_000_000,
            bytes: 256 * 1024 * 1024,
            shifted_degree: 24,
            degree: 1024,
            components: 4096,
            order: 1024,
        }
    }
}
fn limit(message: &str) -> Error {
    Error::Limit(format!("integer residual trial: {message}"))
}
fn height(a: &Gaussian) -> u64 {
    [&a.re, &a.im]
        .into_iter()
        .map(|a| {
            a.numerator_ref()
                .significant_bits()
                .max(a.denominator_ref().significant_bits())
        })
        .max()
        .unwrap_or(0)
}
fn polynomial_height(a: &ExactPolynomial) -> u64 {
    a.coefficients().iter().map(height).max().unwrap_or(0)
}
fn float_height(a: &Float) -> u64 {
    u64::from(a.prec()).saturating_add(u64::from(a.as_raw().get_exp().unwrap_or(0).unsigned_abs()))
}
fn gaussian(a: &C) -> Gaussian {
    Gaussian::new(a.re.to_rational(), a.im.to_rational())
}

struct Budget<'a> {
    limits: Limits,
    work: u64,
    bytes: u64,
    cancellation: &'a CancellationToken,
}
impl Budget<'_> {
    fn charge(&mut self, terms: usize, bits: u64, operations: u64) -> Result<()> {
        self.cancellation.check()?;
        if bits > self.limits.integer_bits {
            return Err(limit("intermediate coefficient height"));
        }
        let words = bits.div_ceil(64).max(1);
        self.work = self.work.saturating_add(operations.saturating_mul(words));
        // Include coefficient headers and transient copies, as well as limbs.
        self.bytes = self
            .bytes
            .saturating_add((terms as u64).saturating_mul(bits.div_ceil(8).saturating_add(64)));
        if self.work > self.limits.work {
            return Err(limit("cumulative native work estimate"));
        }
        if self.bytes > self.limits.bytes {
            return Err(limit("cumulative allocation estimate"));
        }
        Ok(())
    }
    fn check(&self, bits: u64) -> Result<()> {
        self.cancellation.check()?;
        if bits > self.limits.integer_bits {
            Err(limit("observed coefficient height"))
        } else {
            Ok(())
        }
    }
    fn denominator<'a>(&mut self, values: impl Iterator<Item = &'a Gaussian>) -> Result<Integer> {
        let mut denominator = Integer::one();
        for a in values {
            for q in [&a.re, &a.im] {
                let next = q.denominator_ref();
                if *next <= Integer::zero() {
                    return Err(Error::InvalidInput(
                        "nonpositive exact coefficient denominator".into(),
                    ));
                }
                let bits = denominator
                    .significant_bits()
                    .saturating_add(next.significant_bits());
                self.charge(4, bits, bits.div_ceil(64).max(1))?;
                denominator = denominator.lcm(next);
                self.check(denominator.significant_bits())?;
            }
        }
        Ok(denominator)
    }
    fn shift(&mut self, a: &[Gaussian], center: &Gaussian) -> Result<ExactPolynomial> {
        let ring = FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero()));
        let degree = a.len().saturating_sub(1);
        let bits = a.iter().map(height).max().unwrap_or(0);
        if center.is_zero() {
            self.charge(a.len().saturating_mul(4), bits, a.len() as u64)?;
            return Ok(UnivariatePolynomial::from_coefficients(
                &ring,
                a.to_vec(),
                Arc::new(PolyVariable::Temporary(0)),
            ));
        }
        let predicted = bits
            .saturating_add(height(center))
            .saturating_mul(a.len() as u64)
            .saturating_mul(4)
            .saturating_add(16);
        self.charge(
            a.len().saturating_mul(8),
            predicted,
            (degree as u64)
                .saturating_mul(degree as u64)
                .saturating_mul(4),
        )?;
        let polynomial = UnivariatePolynomial::from_coefficients(
            &ring,
            a.to_vec(),
            Arc::new(PolyVariable::Temporary(0)),
        )
        .shift_var(center);
        self.check(polynomial_height(&polynomial))?;
        Ok(polynomial)
    }
}

struct IntegerPair {
    re: IntegerPolynomial,
    im: IntegerPolynomial,
}
impl IntegerPair {
    fn height(&self) -> u64 {
        self.re
            .coefficients()
            .iter()
            .chain(self.im.coefficients())
            .map(Integer::significant_bits)
            .max()
            .unwrap_or(0)
    }
    fn len(&self) -> usize {
        self.re
            .coefficients()
            .len()
            .max(self.im.coefficients().len())
    }
    fn clear(a: &ExactPolynomial, denominator: &Integer, budget: &mut Budget<'_>) -> Result<Self> {
        let bits = polynomial_height(a).saturating_add(denominator.significant_bits());
        budget.charge(
            a.coefficients().len().saturating_mul(4),
            bits,
            (a.coefficients().len() as u64)
                .saturating_mul(bits.div_ceil(64))
                .saturating_mul(4),
        )?;
        let component = |imaginary: bool| -> Result<IntegerPolynomial> {
            let coefficients = a
                .coefficients()
                .iter()
                .map(|a| {
                    budget.cancellation.check()?;
                    let a = if imaginary { &a.im } else { &a.re };
                    let factor = Z.try_div(denominator, a.denominator_ref()).ok_or_else(|| {
                        Error::Numerical(
                            "exact source scale is not divisible by its coefficient denominator"
                                .into(),
                        )
                    })?;
                    Ok(a.numerator_ref() * factor)
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(UnivariatePolynomial::from_coefficients(
                &Z,
                coefficients,
                Arc::new(PolyVariable::Temporary(0)),
            ))
        };
        let result = Self {
            re: component(false)?,
            im: component(true)?,
        };
        budget.check(result.height())?;
        Ok(result)
    }
    fn derivative(&self, budget: &mut Budget<'_>) -> Result<Self> {
        let bits = self
            .height()
            .saturating_add(u64::from(usize::BITS - self.len().max(1).leading_zeros()));
        budget.charge(self.len().saturating_mul(2), bits, self.len() as u64)?;
        let result = Self {
            re: self.re.derivative(),
            im: self.im.derivative(),
        };
        budget.check(result.height())?;
        Ok(result)
    }
    fn multiply(&self, b: &Self, budget: &mut Budget<'_>) -> Result<Self> {
        let terms = self.len().saturating_add(b.len()).max(1);
        let bits = self
            .height()
            .saturating_add(b.height())
            .saturating_add(u64::from(usize::BITS - terms.leading_zeros()))
            .saturating_add(2);
        budget.charge(
            terms.saturating_mul(8),
            bits,
            (self.len() as u64)
                .saturating_mul(b.len() as u64)
                .saturating_mul(4),
        )?;
        let result = Self {
            re: &(&self.re * &b.re) - &(&self.im * &b.im),
            im: &(&self.re * &b.im) + &(&self.im * &b.re),
        };
        budget.check(result.height())?;
        Ok(result)
    }
    fn subtract(&self, b: &Self, budget: &mut Budget<'_>) -> Result<Self> {
        let terms = self.len().max(b.len());
        let bits = self.height().max(b.height()).saturating_add(1);
        budget.charge(
            terms.saturating_mul(2),
            bits,
            (terms as u64).saturating_mul(2),
        )?;
        let result = Self {
            re: &self.re - &b.re,
            im: &self.im - &b.im,
        };
        budget.check(result.height())?;
        Ok(result)
    }
    fn rational(&self, denominator: &Integer, budget: &mut Budget<'_>) -> Result<ExactPolynomial> {
        let bits = self.height().max(denominator.significant_bits());
        budget.charge(
            self.len().saturating_mul(8),
            bits,
            (self.len() as u64)
                .saturating_mul(bits.div_ceil(64))
                .saturating_mul(2),
        )?;
        let zero = Integer::zero();
        let coefficients = (0..self.len())
            .map(|i| {
                budget.cancellation.check()?;
                Ok(Gaussian::new(
                    Rational::from((
                        self.re.coefficients().get(i).unwrap_or(&zero).clone(),
                        denominator.clone(),
                    )),
                    Rational::from((
                        self.im.coefficients().get(i).unwrap_or(&zero).clone(),
                        denominator.clone(),
                    )),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(UnivariatePolynomial::from_coefficients(
            &FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero())),
            coefficients,
            Arc::new(PolyVariable::Temporary(0)),
        ))
    }
}

fn selected(
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    size: usize,
    limits: Limits,
    cancellation: &CancellationToken,
) -> Result<()> {
    cancellation.check()?;
    if coefficients.len() > limits.order || size > limits.components {
        return Err(limit("solution dimensions"));
    }
    for value in std::iter::once(center).chain(coefficients.iter().flatten()) {
        cancellation.check()?;
        if [&value.re, &value.im]
            .into_iter()
            .any(|a| float_height(a) > limits.input_bits)
        {
            return Err(limit("stored dyadic representation"));
        }
    }
    let shifted = !(center.re.is_zero() && center.im.is_zero());
    for row in rows {
        for a in std::iter::once(&row.denominator)
            .chain(row.entries.iter().map(|(_, a)| a))
            .chain(row.denominator_factors.iter().map(|(a, _)| a))
        {
            cancellation.check()?;
            if a.len().saturating_sub(1)
                > if shifted {
                    limits.shifted_degree
                } else {
                    limits.degree
                }
            {
                return Err(limit("source shift degree"));
            }
            if a.iter().any(|a| height(a) > limits.input_bits) {
                return Err(limit("exact source coefficient height"));
            }
        }
    }
    Ok(())
}

fn construct(
    p: Precision,
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
    budget: &mut Budget<'_>,
) -> Result<ExactSourceResidual> {
    construct_observed(p, rows, center, coefficients, channels, budget, |_| {})
}

fn construct_observed(
    p: Precision,
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
    budget: &mut Budget<'_>,
    mut observe: impl FnMut(&ExactPolynomial),
) -> Result<ExactSourceResidual> {
    let n = rows.len();
    let size = n * channels;
    // Preflight both conversion and retained fallback inputs. These estimates
    // intentionally count copies that may have shorter lifetimes in practice.
    for a in std::iter::once(center).chain(coefficients.iter().flatten()) {
        let bits = float_height(&a.re).max(float_height(&a.im));
        budget.charge(8, bits, 4)?;
    }
    for row in rows {
        for a in std::iter::once(&row.denominator)
            .chain(row.entries.iter().map(|(_, a)| a))
            .chain(row.denominator_factors.iter().map(|(a, _)| a))
        {
            budget.charge(
                a.len().saturating_mul(8),
                a.iter().map(height).max().unwrap_or(0),
                a.len() as u64,
            )?;
        }
    }
    let ring = FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero()));
    let variable = Arc::new(PolyVariable::Temporary(0));
    let polynomial = |a| ExactPolynomial::from_coefficients(&ring, a, variable.clone());
    let solutions = (0..size)
        .map(|i| {
            budget.cancellation.check()?;
            Ok(polynomial(
                coefficients
                    .iter()
                    .map(|r| {
                        budget.cancellation.check()?;
                        Ok(gaussian(&r[i]))
                    })
                    .collect::<Result<Vec<_>>>()?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let solution_scale = budget.denominator(solutions.iter().flat_map(|a| a.coefficients()))?;
    let solutions = solutions
        .iter()
        .map(|a| IntegerPair::clear(a, &solution_scale, budget))
        .collect::<Result<Vec<_>>>()?;
    let derivatives = solutions
        .iter()
        .map(|a| a.derivative(budget))
        .collect::<Result<Vec<_>>>()?;
    let ball_ring = FloatField::from_rep(exact_ball(&ring.zero(), p));
    let enclose = |a: &ExactPolynomial, budget: &mut Budget<'_>| -> Result<BallPolynomial> {
        budget.charge(
            a.coefficients().len().saturating_mul(8),
            u64::from(p.bits),
            a.coefficients().len() as u64,
        )?;
        Ok(BallPolynomial::from_coefficients(
            &ball_ring,
            a.coefficients()
                .iter()
                .map(|a| {
                    budget.cancellation.check()?;
                    Ok(exact_ball(a, p))
                })
                .collect::<Result<Vec<_>>>()?,
            variable.clone(),
        ))
    };
    let mut residuals = vec![enclose(&polynomial(vec![ring.zero()]), budget)?; size];
    let mut denominators = Vec::with_capacity(n);
    let center = gaussian(center);
    for (i, row) in rows.iter().enumerate() {
        let denominator = budget.shift(&row.denominator, &center)?;
        let entries = row
            .entries
            .iter()
            .map(|(j, a)| Ok((*j, budget.shift(a, &center)?)))
            .collect::<Result<Vec<_>>>()?;
        let source_scale = budget.denominator(
            denominator
                .coefficients()
                .iter()
                .chain(entries.iter().flat_map(|(_, a)| a.coefficients())),
        )?;
        let source_denominator = IntegerPair::clear(&denominator, &source_scale, budget)?;
        let entries = entries
            .iter()
            .map(|(j, a)| Ok((*j, IntegerPair::clear(a, &source_scale, budget)?)))
            .collect::<Result<Vec<_>>>()?;
        let bits = source_scale
            .significant_bits()
            .saturating_add(solution_scale.significant_bits());
        budget.charge(2, bits, bits.div_ceil(64))?;
        let common_scale = &source_scale * &solution_scale;
        budget.check(common_scale.significant_bits())?;
        for channel in 0..channels {
            let index = channel * n + i;
            let mut numerator = source_denominator.multiply(&derivatives[index], budget)?;
            for (column, a) in &entries {
                let shift = column / n;
                if shift <= channel {
                    let source = (channel - shift) * n + column % n;
                    let term = a.multiply(&solutions[source], budget)?;
                    numerator = numerator.subtract(&term, budget)?;
                }
            }
            let numerator = numerator.rational(&common_scale, budget)?;
            observe(&numerator);
            residuals[index] = enclose(&numerator, budget)?;
        }
        observe(&denominator);
        let expanded = enclose(&denominator, budget)?;
        let factors = row
            .denominator_factors
            .iter()
            .map(|(a, k)| {
                let factor = budget.shift(a, &center)?;
                observe(&factor);
                Ok((enclose(&factor, budget)?, *k))
            })
            .collect::<Result<Vec<_>>>()?;
        denominators.push(DenominatorEnclosure { expanded, factors });
    }
    Ok(ExactSourceResidual::from_polynomials(
        p,
        residuals,
        denominators,
    ))
}

#[derive(Debug)]
pub(crate) struct LazyBall {
    p: Precision,
    rows: Vec<ExactPolynomialRow>,
    center: C,
    coefficients: Vec<Vec<C>>,
    channels: usize,
    cancellation: CancellationToken,
    residual: Mutex<Option<ExactSourceResidual>>,
}
impl LazyBall {
    fn bounds(&self, step: &C, budget: Option<(&[C], &Float)>) -> Result<Vec<Float>> {
        self.cancellation.check()?;
        let mut retained = self
            .residual
            .lock()
            .map_err(|_| Error::Numerical("poisoned lazy source enclosure".into()))?;
        self.cancellation.check()?;
        if retained.is_none() {
            let source = ExactSourceResidual::new(
                self.p,
                &self.rows,
                &self.center,
                &self.coefficients,
                self.channels,
            )?;
            self.cancellation.check()?;
            *retained = Some(source);
        }
        let result = retained
            .as_ref()
            .unwrap()
            .defect_bounds(self.p, step, budget);
        self.cancellation.check()?;
        result
    }
}

#[derive(Clone, Debug)]
pub(crate) enum AdaptiveResidual {
    Ball {
        source: ExactSourceResidual,
        cancellation: CancellationToken,
    },
    Integer {
        source: ExactSourceResidual,
        fallback: Arc<LazyBall>,
    },
}
impl AdaptiveResidual {
    pub(crate) fn ball(
        p: Precision,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        cancellation: &CancellationToken,
    ) -> Result<Self> {
        cancellation.check()?;
        let source = ExactSourceResidual::new(p, rows, center, coefficients, channels)?;
        cancellation.check()?;
        Ok(Self::Ball {
            source,
            cancellation: cancellation.clone(),
        })
    }
    pub(crate) fn new(
        p: Precision,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        cancellation: &CancellationToken,
    ) -> Result<Self> {
        Self::with_limits(
            p,
            rows,
            center,
            coefficients,
            channels,
            cancellation,
            Limits::default(),
        )
    }
    fn with_limits(
        p: Precision,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        cancellation: &CancellationToken,
        limits: Limits,
    ) -> Result<Self> {
        cancellation.check()?;
        let (_, size) = source::validate_inputs(p, rows, center, coefficients, channels)?;
        let integer =
            selected(rows, center, coefficients, size, limits, cancellation).and_then(|()| {
                construct(
                    p,
                    rows,
                    center,
                    coefficients,
                    channels,
                    &mut Budget {
                        limits,
                        work: 0,
                        bytes: 0,
                        cancellation,
                    },
                )
            });
        match integer {
            Ok(source) => {
                cancellation.check()?;
                Ok(Self::Integer {
                    source,
                    fallback: Arc::new(LazyBall {
                        p,
                        rows: rows.to_vec(),
                        center: center.clone(),
                        coefficients: coefficients.to_vec(),
                        channels,
                        cancellation: cancellation.clone(),
                        residual: Mutex::new(None),
                    }),
                })
            }
            Err(Error::Limit(_)) => {
                Self::ball(p, rows, center, coefficients, channels, cancellation)
            }
            Err(error) => Err(error),
        }
    }
    pub(crate) fn defect_bounds(
        &self,
        p: Precision,
        step: &C,
        budget: Option<(&[C], &Float)>,
    ) -> Result<Vec<Float>> {
        match self {
            Self::Ball {
                source,
                cancellation,
            } => {
                cancellation.check()?;
                let result = source.defect_bounds(p, step, budget);
                cancellation.check()?;
                result
            }
            Self::Integer { source, fallback } => {
                fallback.cancellation.check()?;
                let result = source.defect_bounds(p, step, budget);
                fallback.cancellation.check()?;
                let Some((values, tolerance)) = budget else {
                    return result;
                };
                match result {
                    Ok(bounds)
                        if bounds.iter().zip(values).all(|(bound, value)| {
                            let re = value.re.norm();
                            let im = value.im.norm();
                            let magnitude = if re > im { re } else { im };
                            let scale = if magnitude > p.real(1) {
                                magnitude
                            } else {
                                p.real(1)
                            };
                            *bound <= tolerance.mul_round(&scale, p.bits, RoundingDirection::Down)
                        }) =>
                    {
                        Ok(bounds)
                    }
                    Err(error @ (Error::InvalidInput(_) | Error::Cancelled)) => Err(error),
                    _ => fallback.bounds(step, budget),
                }
            }
        }
    }
}
