//! Arithmetic-conditioning diagnostic, not a full Taylor-coefficient error certificate.
//! Native directed balls enclose the polynomial on explicitly stated one-unit
//! working-precision perturbations. Coefficient-generation and inherited source
//! errors remain separate from this endpoint cancellation diagnostic.
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::sync::Arc;
use symbolica::domains::float::{ComplexBall, FloatField, RealBall};
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

type BallPolynomial = UnivariatePolynomial<FloatField<ComplexBall>>;
type RealPolynomial = UnivariatePolynomial<FloatField<RealBall>>;
#[derive(Clone, Debug)]
enum ConditioningPolynomials {
    Real(Vec<RealPolynomial>),
    Complex(Vec<BallPolynomial>),
}
#[derive(Clone, Debug)]
pub(crate) struct ConditioningChart {
    polynomials: ConditioningPolynomials,
}

fn perturbation_real(p: Precision, value: &Float) -> RealBall {
    let bits = p.bits.min(value.prec());
    let epsilon = Rational::from((Integer::one(), Integer::from(2).pow(u64::from(bits))));
    let center = value.to_rational();
    // Relative perturbations of the stored value measure conditioning
    // independently of coordinate scale. An absolute floor would invent
    // high-degree uncertainty even for exact zero Taylor coefficients.
    let radius = epsilon * center.clone().abs();
    RealBall::from_rational_ball(&center, &radius, p.bits)
}

pub(crate) fn perturbation_ball(p: Precision, value: &C) -> ComplexBall {
    ComplexBall::new(
        perturbation_real(p, &value.re),
        perturbation_real(p, &value.im),
    )
}

impl ConditioningChart {
    pub(crate) fn new(p: Precision, coefficients: &[Vec<C>]) -> Result<Self> {
        let n = coefficients.first().map_or(0, Vec::len);
        if n == 0
            || coefficients
                .iter()
                .any(|r| r.len() != n || r.iter().any(|x| !p.finite(x)))
        {
            return Err(Error::InvalidInput(
                "conditioning polynomial dimensions or coefficients".into(),
            ));
        }
        let variable = Arc::new(PolyVariable::Temporary(0));
        if coefficients
            .iter()
            .flatten()
            .all(|value| value.im.is_zero())
        {
            let ring = FloatField::from_rep(RealBall::exact(p.real(0)));
            let polynomials = (0..n)
                .map(|i| {
                    RealPolynomial::from_coefficients(
                        &ring,
                        coefficients
                            .iter()
                            .map(|r| perturbation_real(p, &r[i].re))
                            .collect(),
                        variable.clone(),
                    )
                })
                .collect();
            return Ok(Self {
                polynomials: ConditioningPolynomials::Real(polynomials),
            });
        }
        let ring = FloatField::from_rep(ComplexBall::new(
            RealBall::exact(p.real(0)),
            RealBall::exact(p.real(0)),
        ));
        let polynomials = (0..n)
            .map(|i| {
                UnivariatePolynomial::from_coefficients(
                    &ring,
                    coefficients
                        .iter()
                        .map(|r| perturbation_ball(p, &r[i]))
                        .collect(),
                    variable.clone(),
                )
            })
            .collect();
        Ok(Self {
            polynomials: ConditioningPolynomials::Complex(polynomials),
        })
    }
    /// Return the strongest decimal mixed-scale tolerance checked for every
    /// component, or request re-evaluation from fresh sources at higher bits.
    pub(crate) fn check(&self, p: Precision, step: &C, values: &[C], digits: u32) -> Result<u32> {
        let n = match &self.polynomials {
            ConditioningPolynomials::Real(a) => a.len(),
            ConditioningPolynomials::Complex(a) => a.len(),
        };
        if values.len() != n || digits == 0 || !p.finite(step) {
            return Err(Error::InvalidInput(
                "conditioning evaluation dimensions".into(),
            ));
        }
        match &self.polynomials {
            ConditioningPolynomials::Real(polynomials) if step.im.is_zero() => {
                let argument = perturbation_real(p, &step.re);
                check_balls(
                    p,
                    polynomials.iter().map(|a| {
                        ComplexBall::new(a.evaluate(&argument), RealBall::exact(p.real(0)))
                    }),
                    values,
                    digits,
                )
            }
            ConditioningPolynomials::Real(polynomials) => {
                // A real polynomial may still be evaluated on a complex path.
                // Promote its complete balls without discarding their widths.
                let ring = FloatField::from_rep(perturbation_ball(p, &p.zero()));
                let variable = Arc::new(PolyVariable::Temporary(0));
                let argument = perturbation_ball(p, step);
                check_balls(
                    p,
                    polynomials.iter().map(|a| {
                        BallPolynomial::from_coefficients(
                            &ring,
                            a.coefficients()
                                .iter()
                                .map(|r| ComplexBall::new(r.clone(), RealBall::exact(p.real(0))))
                                .collect(),
                            variable.clone(),
                        )
                        .evaluate(&argument)
                    }),
                    values,
                    digits,
                )
            }
            ConditioningPolynomials::Complex(polynomials) => {
                let argument = perturbation_ball(p, step);
                check_balls(
                    p,
                    polynomials.iter().map(|a| a.evaluate(&argument)),
                    values,
                    digits,
                )
            }
        }
    }
}

/// Shared diagnostic for explicitly enclosed candidate evaluations.
pub(crate) fn check_balls(
    p: Precision,
    balls: impl Iterator<Item = ComplexBall>,
    values: &[C],
    digits: u32,
) -> Result<u32> {
    let ceiling = p.bits / 4; // conservative decimal ceiling; no f64 conversion
    let mut checked = ceiling.max(digits);
    let mut minimum_bits = p.bits;
    let mut witness = 0;
    for (index, (ball, value)) in balls.zip(values).enumerate() {
        let center = C::new(ball.re.center.clone(), ball.im.center.clone());
        let delta = p.sub(&center, value);
        let radius = ball.re.radius + ball.im.radius + delta.re.norm() + delta.im.norm();
        if !radius.is_finite() {
            return Err(Error::InsufficientPrecision {
                minimum_bits: p.bits.saturating_mul(2),
                context: "nonfinite endpoint conditioning estimate".into(),
            });
        }
        let magnitude = p.norm(value);
        let scale = if magnitude > p.real(1) {
            magnitude
        } else {
            p.real(1)
        };
        let ratio = radius.clone() / (p.tolerance(digits) * &scale);
        if ratio > p.real(1) {
            let q = ratio.to_rational();
            let extra = q
                .numerator_ref()
                .significant_bits()
                .saturating_sub(q.denominator_ref().significant_bits())
                .saturating_add(1);
            let extra = u32::try_from(extra)
                .map_err(|_| Error::Limit("conditioning precision hint overflow".into()))?;
            let need = p
                .bits
                .checked_add(extra)
                .and_then(|n| n.checked_add(32))
                .ok_or_else(|| Error::Limit("conditioning precision hint overflow".into()))?;
            if need > minimum_bits {
                minimum_bits = need;
                witness = index;
            }
        } else {
            checked = checked_decimal_tolerance(p, &radius, &scale, digits, checked);
        }
    }
    if minimum_bits > p.bits {
        Err(Error::InsufficientPrecision {
            minimum_bits,
            context: format!(
                "endpoint arithmetic cancellation in component {witness}; {digits} requested decimal digits are unresolved"
            ),
        })
    } else {
        Ok(checked)
    }
}

/// The decimal thresholds decrease monotonically. Usually the working
/// precision already supports the ceiling; otherwise find the last passing
/// threshold without rebuilding every intervening arbitrary-precision power.
fn checked_decimal_tolerance(
    p: Precision,
    radius: &Float,
    scale: &Float,
    minimum: u32,
    ceiling: u32,
) -> u32 {
    if *radius <= p.tolerance(ceiling) * scale {
        return ceiling;
    }
    let mut low = minimum;
    let mut high = ceiling;
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if *radius <= p.tolerance(middle) * scale {
            low = middle;
        } else {
            high = middle;
        }
    }
    low
}

#[cfg(test)]
mod decimal_search_tests {
    use super::*;

    #[test]
    fn binary_threshold_search_preserves_the_linear_scan() {
        for bits in [53, 106, 201, 400] {
            let p = Precision { bits };
            for minimum in [1, 8, 20] {
                let ceiling = (bits / 4).max(minimum);
                for scale in [p.real(1), p.real(17), p.real(1024)] {
                    for digits in 0..=ceiling + 2 {
                        for factor in [p.real(0), p.real(1), p.real(2), p.real(3) / 2] {
                            let radius = p.tolerance(digits) * &scale * factor;
                            let mut expected = minimum;
                            for candidate in minimum + 1..=ceiling {
                                if radius > p.tolerance(candidate) * &scale {
                                    break;
                                }
                                expected = candidate;
                            }
                            assert_eq!(
                                checked_decimal_tolerance(p, &radius, &scale, minimum, ceiling),
                                expected
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Exact waypoints still identify requested displacements even if conversion
/// rounds two of them together. The first boundary point is only a Float:
/// identity is proven against that exact dyadic, not an invented exact input.
/// The numeric low-level transport deliberately has no such preflight.
pub(crate) fn exact_waypoint_conditioning(
    p: Precision,
    boundary: &C,
    exact_waypoints: &[Atom],
    points: &[C],
    digits: u32,
) -> Result<u32> {
    if exact_waypoints.len() != points.len() {
        return Err(Error::InvalidInput(
            "coordinate conditioning dimensions".into(),
        ));
    }
    let mut exact_previous = Atom::num(Complex::new(
        boundary.re.to_rational(),
        boundary.im.to_rational(),
    ));
    let mut previous = boundary.clone();
    let mut checked = (p.bits / 4).max(digits);
    for (exact, point) in exact_waypoints.iter().zip(points) {
        if !(exact - &exact_previous).together().cancel().is_zero() {
            let chart = ConditioningChart::new(p, &[vec![p.neg(&previous)], vec![point.clone()]])?;
            checked = checked.min(
                chart
                    .check(p, &p.i(1), &[p.sub(point, &previous)], digits)
                    .map_err(|error| match error {
                        Error::InsufficientPrecision {
                            minimum_bits,
                            context,
                        } => Error::InsufficientPrecision {
                            minimum_bits,
                            context: format!("exact waypoint separation: {context}"),
                        },
                        other => other,
                    })?,
            );
        }
        exact_previous = exact.clone();
        previous = point.clone();
    }
    Ok(checked)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_perturbations_preserve_coordinate_scaling_and_exact_zero_tails() -> Result<()> {
        let p = Precision::decimal(60)?;
        let radius = p.powi(&p.i(2), 100);
        for scale in [p.i(1), radius.clone()] {
            let mut coefficients = vec![vec![p.zero()]; 17];
            coefficients[0][0] = p.i(1);
            coefficients[16][0] = p.powi(&scale, -16);
            let chart = ConditioningChart::new(p, &coefficients)?;
            assert!(chart.check(p, &scale, &[p.i(2)], 20)? >= 20);
        }
        let mut padded = vec![vec![p.zero()]; 65];
        padded[0][0] = p.i(3);
        let chart = ConditioningChart::new(p, &padded)?;
        assert!(chart.check(p, &radius, &[p.i(3)], 20)? >= 20);
        assert!(chart.check(p, &p.div(&p.i(1), &radius), &[p.i(3)], 20)? >= 20);
        Ok(())
    }

    #[test]
    fn exact_zero_component_does_not_limit_other_component_precision() -> Result<()> {
        let p = Precision::decimal(80)?;
        let value = C::new(p.real(3), Float::with_val(16, 0));
        let chart = ConditioningChart::new(p, &[vec![value.clone()]])?;
        assert!(chart.check(p, &p.i(1), &[value], 40)? >= 40);
        Ok(())
    }
}
