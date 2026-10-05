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
    /// Certified native factorization, including constant content.
    pub(crate) denominator_factors: Vec<(Vec<Gaussian>, usize)>,
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
    residuals: Vec<BallPolynomial>,
    /// One denominator per physical row, shared by every epsilon channel.
    denominators: Vec<DenominatorEnclosure>,
}

#[derive(Clone, Debug)]
struct DenominatorEnclosure {
    expanded: BallPolynomial,
    factors: Vec<(BallPolynomial, usize)>,
}

/// Triangle inequality on a disk; a nonpositive result is inconclusive.
fn polynomial_lower(a: &BallPolynomial, radius: &Float, p: Precision) -> Float {
    let mut lower = complex_lower(&a.get_constant(), p);
    let mut power = radius.clone();
    for coefficient in a.coefficients().iter().skip(1) {
        lower = lower.sub_round(
            &complex_upper(coefficient, p).mul_round(&power, p.bits, RoundingDirection::Up),
            p.bits,
            RoundingDirection::Down,
        );
        power = power.mul_round(radius, p.bits, RoundingDirection::Up);
    }
    if lower.is_finite() && lower > p.real(0) {
        lower
    } else {
        p.real(0)
    }
}

impl DenominatorEnclosure {
    fn lower(&self, radius: &Float, p: Precision) -> Float {
        let expanded = polynomial_lower(&self.expanded, radius, p);
        let mut factored = p.real(1);
        for (factor, multiplicity) in &self.factors {
            let lower = polynomial_lower(factor, radius, p);
            // Inconclusive factors cannot certify a product. In particular,
            // two negative triangle bounds must never become positive.
            if lower.is_zero() {
                return expanded;
            }
            for _ in 0..*multiplicity {
                factored = factored.mul_round(&lower, p.bits, RoundingDirection::Down);
            }
            if !factored.is_finite() || factored <= p.real(0) {
                return expanded;
            }
        }
        if factored > expanded {
            factored
        } else {
            expanded
        }
    }
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
                    || (r.denominator_factors.is_empty()
                        && !(r.denominator.len() == 1 && r.denominator[0].is_one()))
                    || r.denominator_factors
                        .iter()
                        .any(|(a, multiplicity)| a.is_empty() || *multiplicity == 0)
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
        let mut residuals = vec![polynomial(vec![ring.zero()]); size];
        let mut denominators = Vec::with_capacity(n);
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
                residuals[index] = numerator;
            }
            denominators.push(DenominatorEnclosure {
                expanded: denominator,
                factors: row
                    .denominator_factors
                    .iter()
                    .map(|(a, multiplicity)| {
                        (exact_polynomial(a).shift_var(&center), *multiplicity)
                    })
                    .collect(),
            });
        }
        Ok(Self {
            residuals,
            denominators,
        })
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
        let denominator_lowers = self
            .denominators
            .iter()
            .map(|denominator| denominator.lower(&radius, p))
            .collect::<Vec<_>>();
        let mut bounds = Vec::with_capacity(self.residuals.len());
        for (i, numerator) in self.residuals.iter().enumerate() {
            let lower = &denominator_lowers[i % self.denominators.len()];
            if !lower.is_finite() || *lower <= p.real(0) {
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
                .div_round(lower, p.bits, RoundingDirection::Up)
                .mul_round(&radius, p.bits, RoundingDirection::Up);
            let width = width
                .div_round(lower, p.bits, RoundingDirection::Up)
                .mul_round(&radius, p.bits, RoundingDirection::Up);
            let central = central
                .div_round(lower, p.bits, RoundingDirection::Up)
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

#[cfg(test)]
mod factor_tests {
    use super::*;
    use crate::ode::compile_rows;

    fn atom(input: &str) -> Atom {
        Atom::parse(input, "factored_source_disk", Default::default()).unwrap()
    }

    fn variable() -> Symbol {
        match atom("x").as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        }
    }

    fn chart(
        p: Precision,
        denominator: &str,
        center: &C,
        parameters: &ahash::HashMap<Atom, C>,
    ) -> ExactSourceResidual {
        let rows = compile_rows(
            variable(),
            &[vec![Atom::one() / atom(denominator)]],
            p,
            parameters,
        )
        .unwrap();
        ExactSourceResidual::new(p, &rows.exact_source_rows, center, &[vec![p.i(1)]], 1).unwrap()
    }

    #[test]
    fn repeated_pole_disk_gets_a_strictly_positive_directed_lower_bound() {
        let p = Precision::decimal(70).unwrap();
        let chart = chart(p, "(1-x)^10", &p.zero(), &Default::default());
        let radius = p.rational(&Rational::from((1, 3))).re;
        let denominator = &chart.denominators[0];
        assert_eq!(
            polynomial_lower(&denominator.expanded, &radius, p),
            p.real(0)
        );
        let lower = denominator.lower(&radius, p);
        let exact = (Rational::one() - radius.to_rational()).pow(10);
        assert!(lower > p.real(0));
        assert!(lower.to_rational() <= exact);
        assert!(p.close(&C::new(lower, p.real(0)), &p.rational(&exact), 60));
        assert!(
            chart
                .defect_bounds(p, &C::new(radius, p.real(0)), None)
                .is_ok()
        );
    }

    #[test]
    fn expanded_bound_is_retained_when_factorization_is_weaker() {
        let p = Precision::decimal(60).unwrap();
        let chart = chart(p, "(1-x)*(1+x)", &p.zero(), &Default::default());
        let radius = p.rational(&Rational::from((1, 2))).re;
        let denominator = &chart.denominators[0];
        assert_eq!(
            denominator.lower(&radius, p),
            p.rational(&Rational::from((3, 4))).re
        );
    }

    #[test]
    fn two_inconclusive_factors_cannot_certify_a_disk_containing_poles() {
        let p = Precision::decimal(60).unwrap();
        for denominator in ["(1-x)*(1-2*x)", "(1-x)^2"] {
            let chart = chart(p, denominator, &p.zero(), &Default::default());
            assert_eq!(chart.denominators[0].lower(&p.real(2), p), p.real(0));
            assert!(matches!(
                chart.defect_bounds(p, &p.i(2), None),
                Err(Error::Accuracy(_))
            ));
        }
    }

    #[test]
    fn exact_content_multiplicity_and_complex_specialization_survive_shifting() {
        let p = Precision::decimal(70).unwrap();
        let parameters = ahash::HashMap::from_iter([
            (atom("a"), p.complex(3, 1)),
            (atom("m"), p.complex(2, -1)),
        ]);
        let chart = chart(p, "6*a*(m-x)^5", &p.complex(0, -1), &parameters);
        let denominator = &chart.denominators[0];
        let radius = p.rational(&Rational::from((1, 2))).re;
        let lower = denominator.lower(&radius, p);
        // The constant content gives lower |6(3+i)| >= 18; each shifted
        // linear factor gives 2 - 1/2. All inputs here are exact dyadics.
        let expected = p
            .rational(&(Rational::from(18) * Rational::from((3, 2)).pow(5)))
            .re;
        assert_eq!(lower, expected);
        assert!(
            denominator
                .factors
                .iter()
                .any(|(_, multiplicity)| *multiplicity == 5)
        );
        assert!(
            chart
                .defect_bounds(p, &C::new(radius, p.real(0)), None)
                .is_ok()
        );
    }

    #[test]
    fn factors_becoming_constants_preserve_the_source_domain() {
        let p = Precision::decimal(60).unwrap();
        let parameters = ahash::HashMap::from_iter([(atom("a"), p.zero())]);
        let chart = chart(p, "3*(1+a*x)^4", &p.i(2), &parameters);
        assert_eq!(chart.denominators[0].lower(&p.real(100), p), p.real(3));
        assert!(matches!(
            compile_rows(variable(), &[vec![atom("1/(a*x)")]], p, &parameters),
            Err(Error::Numerical(_))
        ));
    }

    #[test]
    fn polynomial_ode_accepts_the_certified_empty_product_for_one() {
        let p = Precision::decimal(60).unwrap();
        let mut rows =
            compile_rows(variable(), &[vec![atom("2*x")]], p, &Default::default()).unwrap();
        let initial =
            ExactSourceResidual::new(p, &rows.exact_source_rows, &p.zero(), &[vec![p.i(1)]], 1)
                .unwrap();
        let expected = initial.defect_bounds(p, &p.i(1), None).unwrap();
        assert_eq!(expected, vec![p.real(2)]);
        // Native factorization may represent 1 by the empty product.
        rows.exact_source_rows[0].denominator_factors.clear();
        let empty =
            ExactSourceResidual::new(p, &rows.exact_source_rows, &p.zero(), &[vec![p.i(1)]], 1)
                .unwrap();
        assert_eq!(empty.defect_bounds(p, &p.i(1), None).unwrap(), expected);
        rows.exact_source_rows[0].denominator[0] =
            Gaussian::new(Rational::from(2), Rational::zero());
        assert!(matches!(
            ExactSourceResidual::new(p, &rows.exact_source_rows, &p.zero(), &[vec![p.i(1)]], 1,),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn coupled_epsilon_channels_share_denominator_enclosures() {
        let p = Precision::decimal(60).unwrap();
        let rows = compile_rows(
            variable(),
            &[vec![atom("0"), atom("1/(1-x)^10")]],
            p,
            &Default::default(),
        )
        .unwrap();
        let chart = ExactSourceResidual::new(
            p,
            &rows.exact_source_rows,
            &p.zero(),
            &[vec![p.i(2), p.i(3), p.i(5)]],
            3,
        )
        .unwrap();
        assert_eq!(chart.denominators.len(), 1);
        let bounds = chart
            .defect_bounds(p, &p.rational(&Rational::from((1, 4))), None)
            .unwrap();
        assert_eq!(bounds[0], p.real(0));
        assert!(bounds[1] > p.real(0));
        assert!(p.close(
            &C::new(bounds[2].clone(), p.real(0)),
            &p.scale(&C::new(bounds[1].clone(), p.real(0)), 3, 2),
            50
        ));
    }

    #[test]
    fn distinct_rows_and_cross_channel_coupling_match_an_ordinary_system() {
        let p = Precision::decimal(60).unwrap();
        let a = atom("1/(1-x)^10");
        let b = atom("1/(1+x)^4");
        let z = Atom::zero();
        let twice_a = &a * Atom::num(2);
        let thrice_a = &a * Atom::num(3);
        let four_b = &b * Atom::num(4);
        let epsilon_rows = compile_rows(
            variable(),
            &[
                vec![a.clone(), twice_a.clone(), z.clone(), thrice_a.clone()],
                vec![z.clone(), b.clone(), four_b.clone(), z.clone()],
            ],
            p,
            &Default::default(),
        )
        .unwrap();
        // Write the independent ordinary six-component system explicitly,
        // including the lower-channel source terms and distinct row poles.
        let ordinary_rows = compile_rows(
            variable(),
            &[
                vec![
                    a.clone(),
                    twice_a.clone(),
                    z.clone(),
                    z.clone(),
                    z.clone(),
                    z.clone(),
                ],
                vec![
                    z.clone(),
                    b.clone(),
                    z.clone(),
                    z.clone(),
                    z.clone(),
                    z.clone(),
                ],
                vec![
                    z.clone(),
                    thrice_a.clone(),
                    a.clone(),
                    twice_a.clone(),
                    z.clone(),
                    z.clone(),
                ],
                vec![
                    four_b.clone(),
                    z.clone(),
                    z.clone(),
                    b.clone(),
                    z.clone(),
                    z.clone(),
                ],
                vec![z.clone(), z.clone(), z.clone(), thrice_a, a, twice_a],
                vec![z.clone(), z.clone(), four_b, z.clone(), z, b],
            ],
            p,
            &Default::default(),
        )
        .unwrap();
        let values = vec![vec![p.i(2), p.i(3), p.i(5), p.i(7), p.i(11), p.i(13)]];
        let coupled =
            ExactSourceResidual::new(p, &epsilon_rows.exact_source_rows, &p.zero(), &values, 3)
                .unwrap();
        let ordinary =
            ExactSourceResidual::new(p, &ordinary_rows.exact_source_rows, &p.zero(), &values, 1)
                .unwrap();
        assert_eq!(coupled.denominators.len(), 2);
        assert_eq!(ordinary.denominators.len(), 6);
        for step in [
            p.rational(&Rational::from((1, 4))),
            p.parse("0", "0.25").unwrap(),
        ] {
            assert_eq!(
                coupled.defect_bounds(p, &step, None).unwrap(),
                ordinary.defect_bounds(p, &step, None).unwrap()
            );
        }
    }
}
