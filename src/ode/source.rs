//! Exact-source coefficient enclosures for ordinary rational Taylor charts.
//!
//! The stored Taylor coefficients define an exact dyadic polynomial. Native
//! ball polynomial arithmetic encloses its differential defect against the
//! original rational source; it does not duplicate the Taylor recurrence.
//! The resulting local defect bound does not include global amplification,
//! boundary uncertainty, coordinate evaluation or singular matching errors.
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::sync::Arc;
use symbolica::domains::float::{ComplexBall, FloatField, RealBall, RoundingDirection};
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

pub(crate) type Gaussian = Complex<Rational>;
type BallPolynomial = UnivariatePolynomial<FloatField<ComplexBall>>;

/// Exact coefficients after specializing numerical parameters as their stored
/// dyadics. Any uncertainty in those parameters remains a separate contract.
#[derive(Clone, Debug)]
pub(crate) struct ExactPolynomialRow {
    pub(crate) denominator: Vec<Gaussian>,
    pub(crate) entries: Vec<(usize, Vec<Gaussian>)>,
}

pub(crate) struct ExactSpecialization {
    values: ahash::HashMap<Atom, Gaussian>,
    ring: FloatField<Gaussian>,
}

impl ExactSpecialization {
    pub(crate) fn new(p: Precision, values: &ahash::HashMap<Atom, C>) -> Result<Self> {
        let mut exact = ahash::HashMap::default();
        for (key, value) in values {
            if !p.finite(value) {
                return Err(Error::InvalidInput(
                    "nonfinite source specialization".into(),
                ));
            }
            exact.insert(
                key.clone(),
                Gaussian::new(value.re.to_rational(), value.im.to_rational()),
            );
        }
        exact.insert(
            Atom::var(crate::family::imaginary_parameter()),
            Gaussian::new(Rational::zero(), Rational::one()),
        );
        Ok(Self {
            values: exact,
            ring: FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero())),
        })
    }

    pub(crate) fn polynomial(&self, a: &Atom, variable: Symbol) -> Result<Vec<Gaussian>> {
        let mut result = vec![self.ring.zero()];
        for (degree, coefficient) in super::polynomial_coefficient_terms(a, variable)? {
            let value = coefficient
                .evaluate_in_ring(&self.values, &self.ring)
                .map_err(|e| {
                    Error::Unsupported(format!("exact rational source specialization: {e}"))
                })?;
            result.resize(result.len().max(degree + 1), self.ring.zero());
            result[degree] += value;
        }
        while result.len() > 1 && result.last().is_some_and(SingleFloat::is_zero) {
            result.pop();
        }
        Ok(result)
    }
}

fn dyadic_ball(p: Precision, value: &C) -> ComplexBall {
    let mut ball = ComplexBall::new(
        RealBall::exact(value.re.clone()),
        RealBall::exact(value.im.clone()),
    );
    // These are exact stored dyadics, not estimates of source accuracy. Native
    // outward rescaling preserves their value/enclosure at the working bits.
    ball.set_precision(p.bits);
    ball
}

fn exact_ball(value: &Gaussian, p: Precision) -> ComplexBall {
    ComplexBall::from_rational_ball(value, &Rational::zero(), p.bits)
}

#[derive(Clone, Debug)]
pub(crate) struct ExactSourceResidual {
    residuals: Vec<(BallPolynomial, BallPolynomial)>,
}

impl ExactSourceResidual {
    pub(crate) fn new(
        p: Precision,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
    ) -> Result<Self> {
        let n = rows.len();
        let size = n
            .checked_mul(channels)
            .ok_or_else(|| Error::Limit("exact-source residual dimensions overflow".into()))?;
        if n == 0
            || channels == 0
            || coefficients.is_empty()
            || !p.finite(center)
            || coefficients
                .iter()
                .any(|r| r.len() != size || r.iter().any(|a| !p.finite(a)))
            || rows.iter().any(|r| {
                r.denominator.is_empty()
                    || r.entries.iter().any(|(j, a)| *j >= size || a.is_empty())
            })
        {
            return Err(Error::InvalidInput(
                "exact-source residual dimensions or coefficients".into(),
            ));
        }
        let ring = FloatField::from_rep(dyadic_ball(p, &p.zero()));
        let variable = Arc::new(PolyVariable::Temporary(0));
        let polynomial = |a| UnivariatePolynomial::from_coefficients(&ring, a, variable.clone());
        let exact_polynomial =
            |a: &[Gaussian]| polynomial(a.iter().map(|a| exact_ball(a, p)).collect());
        let solutions = (0..size)
            .map(|i| polynomial(coefficients.iter().map(|r| dyadic_ball(p, &r[i])).collect()))
            .collect::<Vec<_>>();
        let mut residuals =
            vec![(polynomial(vec![ring.zero()]), polynomial(vec![ring.one()])); size];
        let center = dyadic_ball(p, center);
        for (i, row) in rows.iter().enumerate() {
            let denominator = exact_polynomial(&row.denominator).shift_var(&center);
            let entries = row
                .entries
                .iter()
                .map(|(j, a)| (*j, exact_polynomial(a).shift_var(&center)))
                .collect::<Vec<_>>();
            for channel in 0..channels {
                let index = channel * n + i;
                let mut numerator = &denominator * &solutions[index].derivative();
                for (column, a) in &entries {
                    let shift = column / n;
                    if shift <= channel {
                        let source = (channel - shift) * n + column % n;
                        numerator = &numerator - &(a * &solutions[source]);
                    }
                }
                residuals[index] = (numerator, denominator.clone());
            }
        }
        Ok(Self { residuals })
    }

    /// Bound |h| sup |P'-AP| on a disk containing the entire straight step.
    /// Directed arithmetic certifies this local defect enclosure only.
    pub(crate) fn defect_bounds(
        &self,
        p: Precision,
        step: &C,
        budget: Option<(&[C], &Float)>,
    ) -> Result<Vec<Float>> {
        if !p.finite(step)
            || budget.is_some_and(|(values, tolerance)| {
                values.len() != self.residuals.len()
                    || values.iter().any(|a| !p.finite(a))
                    || !tolerance.is_finite()
                    || *tolerance <= p.real(0)
            })
        {
            return Err(Error::InvalidInput(
                "exact-source residual step or budget".into(),
            ));
        }
        // A directed component one-norm is a certified upper complex radius;
        // no uncertified square root or transcendental ball operation is used.
        let radius = step
            .re
            .norm()
            .add_round(&step.im.norm(), p.bits, RoundingDirection::Up);
        let mut bounds = Vec::with_capacity(self.residuals.len());
        for (i, (numerator, denominator)) in self.residuals.iter().enumerate() {
            let mut lower = complex_lower(&denominator.get_constant(), p);
            let mut power = radius.clone();
            for a in denominator.coefficients().iter().skip(1) {
                lower = lower.sub_round(
                    &complex_upper(a, p).mul_round(&power, p.bits, RoundingDirection::Up),
                    p.bits,
                    RoundingDirection::Down,
                );
                power = power.mul_round(&radius, p.bits, RoundingDirection::Up);
            }
            if !lower.is_finite() || lower <= p.real(0) {
                return Err(Error::Accuracy(
                    "exact-source denominator enclosure is inconclusive on the trial disk".into(),
                ));
            }
            let mut upper = p.real(0);
            let mut width = p.real(0);
            let mut central = p.real(0);
            let mut power = p.real(1);
            for a in numerator.coefficients() {
                upper = upper.add_round(
                    &complex_upper(a, p).mul_round(&power, p.bits, RoundingDirection::Up),
                    p.bits,
                    RoundingDirection::Up,
                );
                let central_upper = a.re.center.norm().add_round(
                    &a.im.center.norm(),
                    p.bits,
                    RoundingDirection::Up,
                );
                central = central.add_round(
                    &central_upper.mul_round(&power, p.bits, RoundingDirection::Up),
                    p.bits,
                    RoundingDirection::Up,
                );
                let coefficient_width =
                    a.re.radius
                        .add_round(&a.im.radius, p.bits, RoundingDirection::Up);
                width = width.add_round(
                    &coefficient_width.mul_round(&power, p.bits, RoundingDirection::Up),
                    p.bits,
                    RoundingDirection::Up,
                );
                power = power.mul_round(&radius, p.bits, RoundingDirection::Up);
            }
            let bound = upper
                .div_round(&lower, p.bits, RoundingDirection::Up)
                .mul_round(&radius, p.bits, RoundingDirection::Up);
            let width = width
                .div_round(&lower, p.bits, RoundingDirection::Up)
                .mul_round(&radius, p.bits, RoundingDirection::Up);
            let central = central
                .div_round(&lower, p.bits, RoundingDirection::Up)
                .mul_round(&radius, p.bits, RoundingDirection::Up);
            if !bound.is_finite() || !width.is_finite() {
                return Err(Error::InsufficientPrecision {
                    minimum_bits: p.bits.saturating_mul(2),
                    context: "nonfinite exact-source residual enclosure".into(),
                });
            }
            if let Some((values, tolerance)) = budget {
                // Use a lower component modulus for the mixed-scale budget so
                // rounding cannot enlarge the allowed local error.
                let re = values[i].re.norm();
                let im = values[i].im.norm();
                let magnitude = if re > im { re } else { im };
                let scale = if magnitude > p.real(1) {
                    magnitude
                } else {
                    p.real(1)
                };
                let threshold = tolerance.mul_round(&scale, p.bits, RoundingDirection::Down);
                if threshold.is_zero() {
                    return Err(Error::InsufficientPrecision {
                        minimum_bits: p.bits.saturating_mul(2),
                        context: "exact-source local error budget underflows".into(),
                    });
                }
                if width > threshold && width >= central {
                    let ratio = width.to_rational() / threshold.to_rational();
                    let extra = ratio
                        .numerator_ref()
                        .significant_bits()
                        .saturating_sub(ratio.denominator_ref().significant_bits())
                        .saturating_add(1);
                    let minimum_bits = u32::try_from(extra)
                        .ok()
                        .and_then(|n| p.bits.checked_add(n))
                        .and_then(|n| n.checked_add(32))
                        .ok_or_else(|| {
                            Error::Limit("exact-source precision hint overflow".into())
                        })?;
                    return Err(Error::InsufficientPrecision {
                        minimum_bits,
                        context: format!(
                            "exact-source residual enclosure in component {i} exceeds the local error budget"
                        ),
                    });
                }
            }
            bounds.push(bound);
        }
        Ok(bounds)
    }
}

fn real_upper(value: &RealBall) -> Float {
    let lower = value.lower_bound().norm();
    let upper = value.upper_bound().norm();
    if lower > upper { lower } else { upper }
}

fn complex_upper(value: &ComplexBall, p: Precision) -> Float {
    real_upper(&value.re).add_round(&real_upper(&value.im), p.bits, RoundingDirection::Up)
}

fn real_lower(value: &RealBall, p: Precision) -> Float {
    let lower = value.lower_bound();
    let upper = value.upper_bound();
    if lower > p.real(0) {
        lower
    } else if upper < p.real(0) {
        -upper
    } else {
        p.real(0)
    }
}

fn complex_lower(value: &ComplexBall, p: Precision) -> Float {
    let re = real_lower(&value.re, p);
    let im = real_lower(&value.im, p);
    if re > im { re } else { im }
}
