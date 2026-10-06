//! Opt-in native Padé candidates and exact-source local defect certificates.
//!
//! Candidates represent the exact stored Taylor dyadics, not exact solutions.
//! The disk defect and arithmetic conditioning are local checks; inherited
//! input uncertainty and global amplification remain separate contracts.
use super::{
    conditioning,
    source::{
        BallPolynomial, DenominatorEnclosure, ExactPolynomialRow, ExactSourceResidual, Gaussian,
        exact_ball,
    },
};
use crate::{ComplexFloat as C, Error, PadeOptions, Precision, Result, RunContext};
use std::sync::Arc;
use symbolica::domains::float::{ComplexBall, FloatField};
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

type Polynomial = UnivariatePolynomial<FloatField<Gaussian>>;

/// Saved rational dense output. Construction is private to the checked solver;
/// arbitrary external fractions cannot be mistaken for accepted trajectories.
#[derive(Clone, Debug)]
pub struct RationalCandidate {
    numerators: Vec<Polynomial>,
    denominator: Polynomial,
    numerator_derivatives: Vec<Polynomial>,
    denominator_derivative: Polynomial,
    residual: ExactSourceResidual,
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
fn polynomial_height(a: &Polynomial) -> u64 {
    a.coefficients().iter().map(height).max().unwrap_or(0)
}
fn limit(message: &str) -> Error {
    Error::Limit(format!("Padé construction: {message}"))
}

/// Bound the exact representation before materializing a dyadic numerator or
/// denominator. Stored precision and binary exponent are both relevant.
fn check_float(a: &Float, options: &PadeOptions) -> Result<()> {
    if !a.is_finite() {
        return Err(Error::InvalidInput("nonfinite Padé coefficient".into()));
    }
    if a.is_zero() {
        return Ok(());
    }
    let exponent = a
        .as_raw()
        .get_exp()
        .ok_or_else(|| limit("missing binary exponent"))?;
    let bits = u64::from(exponent.unsigned_abs()).saturating_add(u64::from(a.prec()));
    if bits > options.max_input_bits {
        return Err(limit("stored dyadic exceeds input bit limit"));
    }
    Ok(())
}
fn gaussian(a: &C) -> Gaussian {
    Gaussian::new(a.re.to_rational(), a.im.to_rational())
}

/// Finite estimates limit each requested native operation before its work begins.
/// This is intentionally conservative: rejection selects the existing Taylor path.
struct Budget<'a> {
    options: &'a PadeOptions,
    work: u64,
    context: &'a RunContext,
}
impl Budget<'_> {
    fn charge(&mut self, degree: usize, bits: u64) -> Result<()> {
        self.context.cancellation.check()?;
        if bits > self.options.max_coefficient_bits {
            return Err(limit("predicted coefficient growth exceeds bit limit"));
        }
        self.work = self
            .work
            .saturating_add((degree as u64 + 1).saturating_mul(bits.max(1)));
        if self.work > self.options.max_work {
            return Err(limit("native polynomial work estimate exhausted"));
        }
        Ok(())
    }
    fn check(&self, a: &Polynomial) -> Result<()> {
        if polynomial_height(a) > self.options.max_coefficient_bits {
            return Err(limit("coefficient bit limit"));
        }
        Ok(())
    }
    fn multiply(&mut self, a: &Polynomial, b: &Polynomial) -> Result<Polynomial> {
        let degree = a.degree().saturating_add(b.degree());
        let bits = polynomial_height(a)
            .saturating_mul(a.coefficients().len() as u64)
            .saturating_add(polynomial_height(b).saturating_mul(b.coefficients().len() as u64))
            .saturating_mul(4)
            .saturating_add(16);
        self.charge(degree, bits)?;
        let result = a * b;
        self.check(&result)?;
        Ok(result)
    }
    fn subtract(&mut self, a: &Polynomial, b: &Polynomial) -> Result<Polynomial> {
        self.charge(
            a.degree().max(b.degree()),
            polynomial_height(a)
                .saturating_add(polynomial_height(b))
                .saturating_mul(4)
                .saturating_add(16),
        )?;
        let result = a - b;
        self.check(&result)?;
        Ok(result)
    }
    fn gcd(&mut self, a: &Polynomial, b: &Polynomial) -> Result<Polynomial> {
        let degree = a.degree().max(b.degree());
        self.charge(
            degree,
            polynomial_height(a)
                .max(polynomial_height(b))
                .saturating_mul((degree as u64 + 1).saturating_pow(2))
                .saturating_mul(8),
        )?;
        let result = a.gcd_euclidean(b);
        self.check(&result)?;
        Ok(result)
    }
    fn quotient(&mut self, a: &Polynomial, b: &Polynomial) -> Result<Polynomial> {
        let degree = a.degree().max(b.degree());
        self.charge(
            degree,
            polynomial_height(a)
                .max(polynomial_height(b))
                .saturating_mul((degree as u64 + 1).saturating_pow(2))
                .saturating_mul(8),
        )?;
        let result = a
            .try_div(b)
            .ok_or_else(|| Error::Accuracy("native Padé exact quotient failed".into()))?;
        self.check(&result)?;
        Ok(result)
    }
    fn normalize(&mut self, a: Polynomial, constant: &Gaussian) -> Result<Polynomial> {
        self.charge(
            a.degree(),
            polynomial_height(&a)
                .saturating_add(height(constant))
                .saturating_mul(8)
                .saturating_add(16),
        )?;
        if constant.is_zero() {
            return Err(Error::Accuracy("zero Padé normalization".into()));
        }
        let inverse = a.coefficient_ring().inv(constant);
        let result = a.mul_coeff(&inverse);
        self.check(&result)?;
        Ok(result)
    }
    fn derivative(&mut self, a: &Polynomial) -> Result<Polynomial> {
        self.charge(
            a.degree(),
            polynomial_height(a).saturating_add(usize::BITS as u64),
        )?;
        let result = a.derivative();
        self.check(&result)?;
        Ok(result)
    }
    fn shift(&mut self, a: &Polynomial, center: &Gaussian) -> Result<Polynomial> {
        let bits = polynomial_height(a)
            .saturating_mul(a.coefficients().len() as u64)
            .saturating_add(height(center).saturating_mul(a.degree() as u64))
            .saturating_mul(4)
            .saturating_add(16);
        self.charge(a.degree(), bits)?;
        let result = a.shift_var(center);
        self.check(&result)?;
        Ok(result)
    }
}

impl RationalCandidate {
    pub(crate) fn build(
        p: Precision,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        options: &PadeOptions,
        context: &RunContext,
    ) -> Result<Self> {
        context.cancellation.check()?;
        let n = rows.len();
        let size = n
            .checked_mul(channels)
            .ok_or_else(|| limit("dimensions overflow"))?;
        if size == 0 || coefficients.len() < 3 || coefficients.iter().any(|r| r.len() != size) {
            return Err(Error::InvalidInput(
                "Padé source or series dimensions".into(),
            ));
        }
        let degree = options.degree.min((coefficients.len() - 1) / 2);
        let count = 2 * degree + 1;
        for a in std::iter::once(center).chain(coefficients.iter().take(count).flatten()) {
            check_float(&a.re, options)?;
            check_float(&a.im, options)?;
        }
        // Scan source representation before shifts and products. The exact owner
        // retains every sparse high-degree term and original denominator factor.
        for row in rows {
            for a in std::iter::once(&row.denominator)
                .chain(row.denominator_factors.iter().map(|(a, _)| a))
                .chain(row.entries.iter().map(|(_, a)| a))
            {
                if a.is_empty()
                    || a.len() > options.max_work as usize
                    || a.iter().any(|a| height(a) > options.max_input_bits)
                {
                    return Err(limit("exact source representation exceeds input limits"));
                }
            }
            if row.entries.iter().any(|(j, _)| *j >= size) {
                return Err(Error::InvalidInput("Padé source column".into()));
            }
        }
        let ring = FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero()));
        let variable = Arc::new(PolyVariable::Temporary(0));
        let poly = |a| Polynomial::from_coefficients(&ring, a, variable.clone());
        let one = poly(vec![ring.one()]);
        let mut budget = Budget {
            options,
            work: 0,
            context,
        };
        budget.charge(
            size.checked_mul(count)
                .ok_or_else(|| limit("conversion dimensions overflow"))?,
            options.max_input_bits,
        )?;
        let mut fractions = Vec::with_capacity(size);
        for i in 0..size {
            let series = poly(
                coefficients
                    .iter()
                    .take(count)
                    .map(|r| gaussian(&r[i]))
                    .collect(),
            );
            let (mut numerator, mut denominator) = if series.degree() == 0 {
                (series.clone(), one.clone())
            } else {
                budget.charge(
                    degree,
                    polynomial_height(&series)
                        .saturating_mul((degree as u64 + 1).saturating_pow(2))
                        .saturating_mul(8),
                )?;
                let (a, b) = series
                    .clone()
                    .to_multivariate::<u16>()
                    .rational_approximant_univariate(degree as u32, degree as u32)
                    .ok_or_else(|| {
                        Error::Accuracy("native Padé approximant is defective".into())
                    })?;
                (
                    a.to_univariate_from_univariate(0),
                    b.to_univariate_from_univariate(0),
                )
            };
            if denominator.is_zero() {
                return Err(Error::Accuracy("zero Padé denominator".into()));
            }
            let gcd = budget.gcd(&numerator, &denominator)?;
            numerator = budget.quotient(&numerator, &gcd)?;
            denominator = budget.quotient(&denominator, &gcd)?;
            if denominator.get_constant().is_zero()
                || numerator.degree() > degree
                || denominator.degree() > degree
            {
                return Err(Error::Accuracy(
                    "native Padé degree or origin contract failed".into(),
                ));
            }
            let constant = denominator.get_constant();
            numerator = budget.normalize(numerator, &constant)?;
            denominator = budget.normalize(denominator, &constant)?;
            let matched = budget.multiply(&denominator, &series)?;
            let remainder = budget.subtract(&matched, &numerator)?;
            if remainder
                .coefficients()
                .iter()
                .take(count)
                .any(|a| !a.is_zero())
            {
                return Err(Error::Accuracy(
                    "native Padé series matching contract failed".into(),
                ));
            }
            fractions.push((numerator, denominator));
        }
        let mut denominator = one;
        let mut denominator_factors = Vec::new();
        for (_, q) in &fractions {
            let gcd = budget.gcd(&denominator, q)?;
            let factor = budget.quotient(q, &gcd)?;
            let constant = factor.get_constant();
            let factor = budget.normalize(factor, &constant)?;
            if denominator.degree() + factor.degree() > options.max_common_degree {
                return Err(limit("common denominator degree limit"));
            }
            denominator = budget.multiply(&denominator, &factor)?;
            if factor.degree() > 0 {
                denominator_factors.push(factor);
            }
        }
        let numerators = fractions
            .into_iter()
            .map(|(a, b)| {
                let factor = budget.quotient(&denominator, &b)?;
                budget.multiply(&a, &factor)
            })
            .collect::<Result<Vec<_>>>()?;
        let ball_ring = FloatField::from_rep(exact_ball(&ring.zero(), p));
        let enclose = |a: &Polynomial| {
            BallPolynomial::from_coefficients(
                &ball_ring,
                a.coefficients().iter().map(|a| exact_ball(a, p)).collect(),
                variable.clone(),
            )
        };
        let center = gaussian(center);
        let dq = budget.derivative(&denominator)?;
        let numerator_derivatives = numerators
            .iter()
            .map(|a| budget.derivative(a))
            .collect::<Result<Vec<_>>>()?;
        let q2 = budget.multiply(&denominator, &denominator)?;
        let mut residuals = vec![enclose(&poly(vec![ring.zero()])); size];
        let mut denominators = Vec::with_capacity(n);
        for (i, row) in rows.iter().enumerate() {
            let d = budget.shift(&poly(row.denominator.clone()), &center)?;
            let entries = row
                .entries
                .iter()
                .map(|(j, a)| Ok((*j, budget.shift(&poly(a.clone()), &center)?)))
                .collect::<Result<Vec<_>>>()?;
            for channel in 0..channels {
                let index = channel * n + i;
                let first = budget.multiply(&numerator_derivatives[index], &denominator)?;
                let second = budget.multiply(&numerators[index], &dq)?;
                let derivative = budget.subtract(&first, &second)?;
                let mut residual = budget.multiply(&d, &derivative)?;
                for (column, a) in &entries {
                    let shift = column / n;
                    if shift <= channel {
                        let source = (channel - shift) * n + column % n;
                        let term = budget.multiply(a, &numerators[source])?;
                        let term = budget.multiply(&denominator, &term)?;
                        residual = budget.subtract(&residual, &term)?;
                    }
                }
                residuals[index] = enclose(&residual);
            }
            let expanded = enclose(&budget.multiply(&d, &q2)?);
            let mut factors = row
                .denominator_factors
                .iter()
                .map(|(a, k)| Ok((enclose(&budget.shift(&poly(a.clone()), &center)?), *k)))
                .collect::<Result<Vec<_>>>()?;
            // Include D itself when the owner supplied no nonconstant factors.
            if factors.is_empty() {
                factors.push((enclose(&d), 1));
            }
            factors.extend(denominator_factors.iter().map(|a| (enclose(a), 2)));
            denominators.push(DenominatorEnclosure { expanded, factors });
        }
        Ok(Self {
            numerators,
            denominator,
            numerator_derivatives,
            denominator_derivative: dq,
            residual: ExactSourceResidual {
                residuals,
                denominators,
            },
        })
    }

    /// Evaluate the same rational function admitted by the transport controller.
    pub fn evaluate(&self, p: Precision, point: &C) -> Result<(Vec<C>, Vec<C>)> {
        let ring = crate::fixed_series::FixedComplexRing::new(p);
        let numerical = |a: &Polynomial| {
            UnivariatePolynomial::from_coefficients(
                &ring,
                a.coefficients()
                    .iter()
                    .map(|a| C::new(p.rational(&a.re).re, p.rational(&a.im).re))
                    .collect(),
                Arc::new(PolyVariable::Temporary(0)),
            )
        };
        let q = numerical(&self.denominator).evaluate(point);
        if !p.finite(&q) || q == p.zero() {
            return Err(Error::Accuracy(
                "Padé evaluation denominator vanishes".into(),
            ));
        }
        let dq = numerical(&self.denominator_derivative).evaluate(point);
        let mut values = Vec::with_capacity(self.numerators.len());
        let mut derivatives = Vec::with_capacity(self.numerators.len());
        for (a, derivative) in self.numerators.iter().zip(&self.numerator_derivatives) {
            let value = p.div(&numerical(a).evaluate(point), &q);
            derivatives.push(p.div(
                &p.sub(&numerical(derivative).evaluate(point), &p.mul(&value, &dq)),
                &q,
            ));
            values.push(value);
        }
        if values.iter().chain(&derivatives).any(|a| !p.finite(a)) {
            return Err(Error::Accuracy("nonfinite Padé evaluation".into()));
        }
        Ok((values, derivatives))
    }
    pub(crate) fn defect_bounds(
        &self,
        p: Precision,
        step: &C,
        values: &[C],
        tolerance: &Float,
    ) -> Result<Vec<Float>> {
        self.residual
            .defect_bounds(p, step, Some((values, tolerance)))
    }
    pub(crate) fn check_conditioning(
        &self,
        p: Precision,
        point: &C,
        values: &[C],
        digits: u32,
    ) -> Result<u32> {
        let ring = FloatField::from_rep(exact_ball(
            &Gaussian::new(Rational::zero(), Rational::zero()),
            p,
        ));
        let poly = |a: &Polynomial| {
            UnivariatePolynomial::from_coefficients(
                &ring,
                a.coefficients()
                    .iter()
                    .map(|a| {
                        // Outward exact conversion plus one working unit perturbation.
                        let radius = Rational::from((
                            Integer::one(),
                            Integer::from(2).pow(u64::from(p.bits)),
                        ));
                        let re = radius.clone() * a.re.clone().abs();
                        let im = radius * a.im.clone().abs();
                        ComplexBall::new(
                            symbolica::domains::float::RealBall::from_rational_ball(
                                &a.re, &re, p.bits,
                            ),
                            symbolica::domains::float::RealBall::from_rational_ball(
                                &a.im, &im, p.bits,
                            ),
                        )
                    })
                    .collect(),
                Arc::new(PolyVariable::Temporary(0)),
            )
        };
        let argument = conditioning::perturbation_ball(p, point);
        let denominator = poly(&self.denominator).evaluate(&argument);
        conditioning::check_balls(
            p,
            self.numerators
                .iter()
                .map(|a| poly(a).evaluate(&argument) / &denominator),
            values,
            digits,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DifferentialSystem;
    fn row(entry: Atom, p: Precision) -> Result<super::super::CompiledSystem> {
        DifferentialSystem {
            variable: symbol!("pade_unit::x"),
            matrix: vec![vec![entry]],
        }
        .compile(p, &Default::default())
    }
    #[test]
    fn native_exact_rational_zero_and_constant_contracts() -> Result<()> {
        let p = Precision::decimal(60)?;
        let x = Atom::var(symbol!("pade_unit::x"));
        for (entry, value) in [
            (-Atom::one() / (Atom::one() + x), p.i(1)),
            (Atom::zero(), p.zero()),
            (Atom::zero(), p.i(3)),
        ] {
            let system = row(entry, p)?;
            let coefficients = system.taylor(&p.zero(), std::slice::from_ref(&value), 16)?;
            let candidate = RationalCandidate::build(
                p,
                &system.exact_source_rows,
                &p.zero(),
                &coefficients,
                1,
                &PadeOptions {
                    degree: 8,
                    ..Default::default()
                },
                &RunContext::default(),
            )?;
            let step = p.rational(&Rational::from((1, 4)));
            let expected = if system.matrix[0][0].numerator.iter().all(|a| *a == p.zero()) {
                value
            } else {
                p.rational(&Rational::from((4, 5)))
            };
            let values = candidate.evaluate(p, &step)?.0;
            assert!(p.close(&values[0], &expected, 50));
            assert!(
                candidate.defect_bounds(p, &step, &values, &p.tolerance(40))?[0] <= p.tolerance(50)
            );
            assert!(candidate.check_conditioning(p, &step, &values, 30)? >= 30);
        }
        Ok(())
    }
    #[test]
    fn interior_pole_with_tiny_residue_and_sparse_source_are_rejected() -> Result<()> {
        let p = Precision::decimal(180)?;
        let system = row(Atom::zero(), p)?;
        // R=1+delta*z/(z-1/3), whose sampled defects are tiny but which
        // has an interior pole. Exact powers of two preserve stored dyadics.
        let delta = p.powi(&p.i(2), -400);
        let mut coefficients = vec![vec![p.i(1)]];
        for n in 1..=16 {
            coefficients.push(vec![p.neg(&p.mul(&delta, &p.powi(&p.i(3), n)))]);
        }
        let candidate = RationalCandidate::build(
            p,
            &system.exact_source_rows,
            &p.zero(),
            &coefficients,
            1,
            &PadeOptions {
                degree: 1,
                ..Default::default()
            },
            &RunContext::default(),
        )?;
        let values = candidate.evaluate(p, &p.i(1))?.0;
        assert!(matches!(
            candidate.defect_bounds(p, &p.i(1), &values, &p.tolerance(20)),
            Err(Error::Accuracy(_))
        ));
        let x = Atom::var(symbol!("pade_unit::x"));
        let sparse = row(
            x.clone().pow(100) * (x.clone() - Atom::num((1, 2))) * (x - Atom::one()),
            p,
        )?;
        let coefficients = sparse.taylor(&p.zero(), &[p.i(1)], 16)?;
        let candidate = RationalCandidate::build(
            p,
            &sparse.exact_source_rows,
            &p.zero(),
            &coefficients,
            1,
            &PadeOptions {
                degree: 1,
                ..Default::default()
            },
            &RunContext::default(),
        )?;
        assert!(candidate.defect_bounds(p, &p.i(1), &[p.i(1)], &p.tolerance(20))?[0] > p.real(1));
        Ok(())
    }
    #[test]
    fn defective_series_and_preconversion_limits_are_explicit() -> Result<()> {
        let p = Precision::decimal(60)?;
        let system = row(Atom::zero(), p)?;
        // No [1/1] approximant with Q(0)!=0 matches z^2.
        let bad = vec![vec![p.zero()], vec![p.zero()], vec![p.i(1)]];
        assert!(matches!(
            RationalCandidate::build(
                p,
                &system.exact_source_rows,
                &p.zero(),
                &bad,
                1,
                &PadeOptions {
                    degree: 1,
                    ..Default::default()
                },
                &RunContext::default()
            ),
            Err(Error::Accuracy(_))
        ));
        let huge = p.powi(&p.i(2), 100_000);
        let coefficients = vec![vec![huge], vec![p.zero()], vec![p.zero()]];
        assert!(matches!(
            RationalCandidate::build(
                p,
                &system.exact_source_rows,
                &p.zero(),
                &coefficients,
                1,
                &PadeOptions::default(),
                &RunContext::default()
            ),
            Err(Error::Limit(_))
        ));
        let tiny_budget = PadeOptions {
            max_work: 1,
            ..Default::default()
        };
        let coefficients = vec![vec![p.i(1)], vec![p.i(-1)], vec![p.i(1)]];
        assert!(matches!(
            RationalCandidate::build(
                p,
                &system.exact_source_rows,
                &p.zero(),
                &coefficients,
                1,
                &tiny_budget,
                &RunContext::default()
            ),
            Err(Error::Limit(_))
        ));
        let context = RunContext::default();
        context.cancellation.cancel();
        assert!(matches!(
            RationalCandidate::build(
                p,
                &system.exact_source_rows,
                &p.zero(),
                &coefficients,
                1,
                &PadeOptions::default(),
                &context
            ),
            Err(Error::Cancelled)
        ));
        Ok(())
    }
    #[test]
    fn independent_profiles_change_actual_rational_approximants_and_cap_uses_taylor() -> Result<()>
    {
        use crate::diffexp::{EpsilonBoundary, EpsilonSystem, refine_epsilon_transport};
        let x = symbol!("pade_unit::x");
        let system = EpsilonSystem {
            variable: x,
            matrices: vec![vec![vec![Atom::one()]]],
        };
        let options = crate::FlowOptions {
            digits: 2,
            guard_digits: 14,
            series_order: 8,
            pade: Some(PadeOptions {
                degree: 1,
                ..Default::default()
            }),
            ..Default::default()
        };
        let mut observed = Vec::new();
        let result =
            refine_epsilon_transport(&options, &RunContext::default(), true, |p, opts| {
                let ordinary = row(Atom::one(), p)?;
                let coefficients = ordinary.taylor(&p.zero(), &[p.i(1)], opts.series_order)?;
                let candidate = RationalCandidate::build(
                    p,
                    &ordinary.exact_source_rows,
                    &p.zero(),
                    &coefficients,
                    1,
                    opts.pade.as_ref().unwrap(),
                    &RunContext::default(),
                )?;
                observed.push((
                    candidate.denominator.degree(),
                    candidate
                        .evaluate(p, &p.rational(&Rational::from((1, 4))))?
                        .0[0]
                        .clone(),
                ));
                system.compile(p, &Default::default())?.transport(
                    &EpsilonBoundary {
                        point: p.zero(),
                        leading: 0,
                        coefficients: vec![vec![p.i(1)]],
                    },
                    &[p.rational(&Rational::from((1, 100_000_000)))],
                    opts,
                    &RunContext::default(),
                    true,
                )
            })?;
        assert_eq!(observed.iter().map(|a| a.0).collect::<Vec<_>>(), vec![1, 2]);
        let p = Precision {
            bits: result.diagnostics.working_bits,
        };
        assert!(!p.close(&observed[0].1, &observed[1].1, 6));
        assert!(result.diagnostics.pade_steps > 0);
        assert_eq!(result.verified_digits, Some(2));
        let options = crate::FlowOptions {
            digits: 2,
            guard_digits: 14,
            series_order: 80,
            pade: Some(PadeOptions {
                degree: 32,
                ..Default::default()
            }),
            ..Default::default()
        };
        let constant = EpsilonSystem {
            variable: x,
            matrices: vec![vec![vec![Atom::zero()]]],
        };
        let mut profiles = Vec::new();
        let result =
            refine_epsilon_transport(&options, &RunContext::default(), true, |p, opts| {
                profiles.push(opts.pade.as_ref().map(|p| p.degree));
                constant.compile(p, &Default::default())?.transport(
                    &EpsilonBoundary {
                        point: p.zero(),
                        leading: 0,
                        coefficients: vec![vec![p.i(1)]],
                    },
                    &[p.i(1)],
                    opts,
                    &RunContext::default(),
                    true,
                )
            })?;
        assert_eq!(profiles, vec![Some(32), None]);
        assert!(result.segments.iter().all(|s| s.pade.is_none()));
        assert_eq!(result.verified_digits, Some(2));
        Ok(())
    }
}
