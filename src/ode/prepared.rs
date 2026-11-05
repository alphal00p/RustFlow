//! Immutable native rational-row algebra, before any numerical specialization.
use super::*;
use symbolica::domains::float::{
    ComplexBall,
    RoundingDirection::{Down, Up},
};

type ExactPolePolynomial = UnivariatePolynomial<FloatField<source::Gaussian>>;

/// Only the pole scan discards multiplicities. Row denominators and their
/// enclosure certificates continue to use the original factors and content.
fn square_free_pole_polynomial(
    coefficients: Vec<source::Gaussian>,
    variable: Symbol,
) -> Result<ExactPolePolynomial> {
    let ring = FloatField::from_rep(source::Gaussian::new(Rational::zero(), Rational::zero()));
    let polynomial = ExactPolePolynomial::from_coefficients(
        &ring,
        coefficients,
        Arc::new(PolyVariable::Symbol(variable)),
    );
    // A vanished factor is a domain failure, not a constant with no poles.
    if polynomial.is_zero() {
        return Err(Error::Numerical(
            "identically zero specialized pole factor".into(),
        ));
    }
    let content = polynomial.lcoeff();
    let monic = polynomial.clone().div_coeff(&content);
    let derivative = monic.derivative();
    // FloatField's blanket PolynomialGCD implementation returns one, even for
    // exact Gaussian rationals. Explicit Euclidean GCD stays in this field and
    // avoids the u16 exponent conversion performed by ordinary .gcd().
    let repeated = monic.gcd_euclidean(&derivative);
    let square_free = monic
        .try_div(&repeated)
        .ok_or_else(|| Error::Numerical("exact square-free pole division failed".into()))?;
    if monic.clone().mul_coeff(&content) != polynomial
        || &repeated * &square_free != monic
        || derivative.try_div(&repeated).is_none()
        || repeated.degree().checked_add(square_free.degree()) != Some(monic.degree())
        || !square_free
            .gcd_euclidean(&square_free.derivative())
            .is_one()
    {
        return Err(Error::Numerical(
            "exact square-free pole verification failed".into(),
        ));
    }
    Ok(square_free)
}

/// Remove exactly shared roots before numerical enumeration. Previously solved
/// pieces are monic, square-free and pairwise coprime; no numerical GCD or
/// assumption about irreducibility is involved.
fn unscanned_pole_polynomial(
    original: &ExactPolePolynomial,
    previous: &[(ExactPolePolynomial, Vec<C>)],
) -> Result<ExactPolePolynomial> {
    let mut remaining = original.clone();
    let mut common_product = original.one();
    for (scanned, _) in previous {
        let common = remaining.gcd_euclidean(scanned);
        if common.is_zero() {
            return Err(Error::Numerical("zero common pole factor".into()));
        }
        // Ordinary quot_rem handles constants and degree collapse; never use
        // quot_rem_field, whose degree assumptions do not cover these cases.
        let (next, remainder) = remaining.quot_rem(&common);
        let scanned_quotient = scanned
            .try_div(&common)
            .ok_or_else(|| Error::Numerical("exact common pole factor division failed".into()))?;
        if !remainder.is_zero()
            || &common * &next != remaining
            || &common * &scanned_quotient != *scanned
            || common.degree().checked_add(next.degree()) != Some(remaining.degree())
            || common.degree().checked_add(scanned_quotient.degree()) != Some(scanned.degree())
        {
            return Err(Error::Numerical(
                "exact common pole factor verification failed".into(),
            ));
        }
        common_product = &common_product * &common;
        remaining = next;
    }
    if remaining.is_zero()
        || !remaining.lcoeff().is_one()
        || &common_product * &remaining != *original
        || common_product.degree().checked_add(remaining.degree()) != Some(original.degree())
        || !remaining.gcd_euclidean(&remaining.derivative()).is_one()
        || previous
            .iter()
            .any(|(scanned, _)| !remaining.gcd_euclidean(scanned).is_one())
    {
        return Err(Error::Numerical(
            "exact unscanned pole verification failed".into(),
        ));
    }
    Ok(remaining)
}

fn rounded_pole_coefficients(p: Precision, polynomial: &ExactPolePolynomial) -> Vec<C> {
    polynomial
        .coefficients()
        .iter()
        .map(|coefficient| {
            let mut rounded = p.rational(&coefficient.re);
            rounded.im = p.rational(&coefficient.im).re;
            rounded
        })
        .collect()
}

/// Lower modulus using directed component intervals, without a square root.
fn pole_modulus_lower(value: &ComplexBall, p: Precision) -> Float {
    let lower = |component: &symbolica::domains::float::RealBall| {
        let lo = component.lower_bound();
        let hi = component.upper_bound();
        if lo > p.real(0) {
            lo
        } else if hi < p.real(0) {
            -hi
        } else {
            p.real(0)
        }
    };
    let re = lower(&value.re);
    let im = lower(&value.im);
    if re > im { re } else { im }
}

/// P(c+w)=sum b_k w^k, with every b_k enclosed from the exact Gaussian
/// polynomial. Rouché's strict inequality |b_1|r > |b_0|+sum_{k>=2}|b_k|r^k
/// proves exactly one root in this disk. Neither a Newton correction alone nor
/// the precision attached to the center is used as an error certificate.
fn isolated_pole_radius(
    p: Precision,
    polynomial: &source::BallPolynomial,
    root: &C,
    zero_root_exponent: i32,
) -> Option<Float> {
    if !p.finite(root) {
        return None;
    }
    let shifted = polynomial.shift_var(&source::dyadic_ball(p, root));
    let coefficients = shifted.coefficients();
    let derivative = pole_modulus_lower(&coefficients[1], p);
    if !derivative.is_finite() || derivative <= p.real(0) {
        return None;
    }
    let residual = source::complex_upper(&coefficients[0], p);
    // A small positive disk even for exactly representable roots. The radius
    // is a proposal only; the outward-rounded inequality below is mandatory.
    let exponent = [root.re.as_raw().get_exp(), root.im.as_raw().get_exp()]
        .into_iter()
        .flatten()
        .max()
        .unwrap_or(zero_root_exponent);
    let ulp = p.powi(&p.i(2), i64::from(exponent) - i64::from(p.bits)).re;
    let correction = residual.div_round(&derivative, p.bits, Up);
    let radius = if correction > ulp { correction } else { ulp }.mul_round(&p.real(4), p.bits, Up);
    let mut tail = p.real(0);
    for coefficient in coefficients[2..].iter().rev() {
        tail = tail.mul_round(&radius, p.bits, Up).add_round(
            &source::complex_upper(coefficient, p),
            p.bits,
            Up,
        );
    }
    let remainder = tail
        .mul_round(&radius, p.bits, Up)
        .mul_round(&radius, p.bits, Up)
        .add_round(&residual, p.bits, Up);
    let linear = derivative.mul_round(&radius, p.bits, Down);
    (radius.is_finite()
        && radius > p.real(0)
        && remainder.is_finite()
        && linear.is_finite()
        && linear > remainder)
        .then_some(radius)
}

/// Exact common roots have already been removed. Certify the remaining union
/// with disjoint local isolation disks, refining only pieces whose disks are
/// inconclusive or overlap. All retry coefficients come from the retained exact
/// specialization, never a lift of the p-bit row arrays. Root-solver cluster
/// checks remain mandatory, including on each higher-precision solve.
fn resolve_coprime_pole_geometry(
    p: Precision,
    pieces: &mut [(ExactPolePolynomial, Vec<C>)],
    variable: Symbol,
    context: &RunContext,
) -> Result<()> {
    let mut bits = vec![p.bits; pieces.len()];
    let mut retry = vec![false; pieces.len()];
    let mut reason = "inconclusive exact-coefficient isolation disk";
    // At most six additional solves per affected piece. polynomial_roots also
    // has its own bounded scaling guard; none of these bits enter ODE arithmetic.
    for extra in [0_u32, 128, 256, 512, 1024, 2048, 4096] {
        context.cancellation.check()?;
        let geometry = Precision {
            bits: p
                .bits
                .checked_add(extra)
                .ok_or_else(|| Error::Limit("pole geometry precision overflow".into()))?,
        };
        for (i, (polynomial, roots)) in pieces.iter_mut().enumerate() {
            if retry[i] {
                context.cancellation.check()?;
                let coefficients = rounded_pole_coefficients(geometry, polynomial);
                *roots = polynomial_roots(geometry, &coefficients, variable).map_err(|error| {
                    // Only static numerical reasons, degree and precision; no
                    // source expressions, parameters or master values.
                    Error::Numerical(format!(
                        "degree-{} pole geometry at {} bits: {error}",
                        polynomial.degree(),
                        geometry.bits
                    ))
                })?;
                bits[i] = geometry.bits;
                context.cancellation.check()?;
            }
        }
        retry.fill(false);
        // An exact zero has no exponent of its own. Use the smallest
        // nonzero root scale for its radius proposal, so a tiny neighboring
        // pole does not introduce an artificial absolute precision demand.
        // The directed Rouché inequality still verifies every proposed disk.
        let zero_root_exponent = pieces
            .iter()
            .flat_map(|(_, roots)| roots)
            .filter_map(|root| {
                [root.re.as_raw().get_exp(), root.im.as_raw().get_exp()]
                    .into_iter()
                    .flatten()
                    .max()
            })
            .min()
            .unwrap_or(0);
        let mut disks = Vec::new();
        for (i, (polynomial, roots)) in pieces.iter().enumerate() {
            context.cancellation.check()?;
            let q = Precision { bits: bits[i] };
            let ball_polynomial = source::BallPolynomial::from_coefficients(
                &FloatField::from_rep(source::dyadic_ball(q, &q.zero())),
                polynomial
                    .coefficients()
                    .iter()
                    .map(|c| source::exact_ball(c, q))
                    .collect(),
                Arc::new(PolyVariable::Symbol(variable)),
            );
            for root in roots {
                context.cancellation.check()?;
                if let Some(radius) =
                    isolated_pole_radius(q, &ball_polynomial, root, zero_root_exponent)
                {
                    disks.push((i, root, radius));
                } else {
                    retry[i] = true;
                    reason = "inconclusive exact-coefficient isolation disk";
                }
            }
        }
        for (j, (i, root, radius)) in disks.iter().enumerate() {
            context.cancellation.check()?;
            for (other_i, other, other_radius) in &disks[..j] {
                let distance = pole_modulus_lower(
                    &(source::dyadic_ball(geometry, root) - source::dyadic_ball(geometry, other)),
                    geometry,
                );
                let required = radius.add_round(other_radius, geometry.bits, Up);
                if !distance.is_finite() || distance <= required {
                    retry[*i] = true;
                    retry[*other_i] = true;
                    reason = "overlapping exact-coefficient isolation disks";
                }
            }
        }
        if !retry.iter().any(|r| *r) {
            return Ok(());
        }
    }
    let degree = pieces
        .iter()
        .zip(&retry)
        .filter(|(_, retry)| **retry)
        .map(|((polynomial, _), _)| polynomial.degree())
        .max()
        .unwrap_or(0);
    Err(Error::Numerical(format!(
        "ambiguous nearby poles: degree-{degree} geometry at {} bits: {reason}; bounded extra precision exhausted",
        p.bits.saturating_add(4096)
    )))
}

#[derive(Clone, Debug)]
struct PreparedRational {
    numerator: Atom,
    denominator: Atom,
    /// First occurrence of each exact factor, in the original scan order.
    new_pole_factors: Vec<Atom>,
}

#[derive(Clone, Debug)]
struct PreparedRow {
    matrix: Vec<PreparedRational>,
    denominator: Atom,
    denominator_factors: Vec<(Atom, usize)>,
    cleared: Vec<(usize, Atom)>,
}

/// Exact source snapshot and denominator clearing shared by precision profiles.
/// Numeric parameter values are supplied to each compilation, never retained here.
#[derive(Clone, Debug)]
pub(crate) struct PreparedRows {
    variable: Symbol,
    source_rows: Vec<Vec<Atom>>,
    rows: Vec<PreparedRow>,
}

impl PreparedRows {
    /// Number of dense coefficient slots allocated by rows-only specialization.
    /// Count exact sparse polynomial degrees before allocating numerical arrays.
    pub(crate) fn polynomial_storage_values(&self, context: &RunContext) -> Result<u128> {
        let mut count = 0_u128;
        for row in &self.rows {
            for expression in std::iter::once(&row.denominator)
                .chain(row.cleared.iter().map(|(_, expression)| expression))
                .chain(
                    row.denominator_factors
                        .iter()
                        .map(|(expression, _)| expression),
                )
            {
                context.cancellation.check()?;
                let length = polynomial_coefficient_terms(expression, self.variable)?
                    .iter()
                    .map(|(degree, _)| *degree as u128 + 1)
                    .max()
                    .unwrap_or(1);
                count = count.saturating_add(length);
            }
        }
        Ok(count)
    }

    pub(crate) fn new(variable: Symbol, rows: &[Vec<Atom>], context: &RunContext) -> Result<Self> {
        context.cancellation.check()?;
        validate_ode_variable(variable)?;
        if rows.is_empty() || rows.iter().any(Vec::is_empty) {
            return Err(Error::InvalidInput(
                "nonempty rational compilation rows required".into(),
            ));
        }
        let mut prepared = Vec::with_capacity(rows.len());
        let mut seen_factors = ahash::HashSet::default();
        let mut factorizations = ahash::HashMap::default();
        for row in rows {
            context.cancellation.check()?;
            let mut matrix = Vec::with_capacity(row.len());
            let mut exact = Vec::with_capacity(row.len());
            let mut common_denominator: Option<ExactPolynomial> = None;
            for value in row {
                context.cancellation.check()?;
                let rational: RationalPolynomial<IntegerRing, u16> =
                    crate::family::encode_complex(value)
                        .try_to_rational_polynomial(&Q, &Z, None)
                        .map_err(|e| Error::InvalidInput(e.to_string()))?;
                context.cancellation.check()?;
                let factors = denominator_factors(&rational.denominator, &mut factorizations)?;
                context.cancellation.check()?;
                let new_pole_factors = factors
                    .iter()
                    .filter_map(|(factor, _)| {
                        let expression = factor.to_expression();
                        seen_factors
                            .insert(expression.clone())
                            .then_some(expression)
                    })
                    .collect();
                matrix.push(PreparedRational {
                    numerator: rational.numerator.to_expression(),
                    denominator: rational.denominator.to_expression(),
                    new_pole_factors,
                });
                common_denominator = Some(if let Some(previous) = common_denominator {
                    let quotient = previous
                        .try_div(&previous.gcd(&rational.denominator))
                        .ok_or_else(|| {
                            Error::Numerical("exact denominator LCM division failed".into())
                        })?;
                    multiply_polynomials(&quotient, &rational.denominator)?
                } else {
                    rational.denominator.clone()
                });
                context.cancellation.check()?;
                exact.push(rational);
            }
            let common_denominator = common_denominator.unwrap();
            let denominator_factors =
                denominator_factors(&common_denominator, &mut factorizations)?
                    .iter()
                    .map(|(factor, multiplicity)| (factor.to_expression(), *multiplicity))
                    .collect();
            context.cancellation.check()?;
            let mut cleared = Vec::new();
            for (column, rational) in exact.iter().enumerate() {
                context.cancellation.check()?;
                if rational.numerator.is_zero() {
                    continue;
                }
                let multiplier = common_denominator
                    .try_div(&rational.denominator)
                    .ok_or_else(|| Error::Numerical("exact denominator clearing failed".into()))?;
                let polynomial = multiply_polynomials(&rational.numerator, &multiplier)?;
                context.cancellation.check()?;
                cleared.push((column, polynomial.to_expression()));
            }
            prepared.push(PreparedRow {
                matrix,
                denominator: common_denominator.to_expression(),
                denominator_factors,
                cleared,
            });
        }
        context.cancellation.check()?;
        Ok(Self {
            variable,
            source_rows: rows.to_vec(),
            rows: prepared,
        })
    }

    /// Preserve row arithmetic and fresh stored-dyadic specialization. Only the
    /// pole scan uses an exactly specialized square-free polynomial.
    pub(crate) fn compile(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        context: &RunContext,
    ) -> Result<CompiledSystem> {
        context.cancellation.check()?;
        let variable = self.variable;
        let exact_values = source::ExactSpecialization::new(p, values)?;
        let mut matrix = Vec::with_capacity(self.rows.len());
        let mut scanned_poles = Vec::new();
        let mut pole_polynomials = Vec::new();
        let mut polynomial_rows = Vec::with_capacity(self.rows.len());
        let mut exact_source_rows = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            context.cancellation.check()?;
            let mut out = Vec::with_capacity(row.matrix.len());
            for rational in &row.matrix {
                context.cancellation.check()?;
                let numerator = polynomial_coefficients(&rational.numerator, variable, p, values)?;
                let denominator =
                    polynomial_coefficients(&rational.denominator, variable, p, values)?;
                if denominator.iter().all(|v| *v == p.zero()) {
                    return Err(Error::Numerical(
                        "identically zero specialized denominator".into(),
                    ));
                }
                for factor in &rational.new_pole_factors {
                    context.cancellation.check()?;
                    // Specialize J=i and stored parameter dyadics before any
                    // rounding: specialization can create exact multiplicities
                    // invisible to the formal multivariate factorization.
                    let square_free = square_free_pole_polynomial(
                        exact_values.polynomial(factor, variable)?,
                        variable,
                    )?;
                    context.cancellation.check()?;
                    if square_free.is_constant() {
                        continue;
                    }
                    // Retain every original nonconstant factor for source-domain
                    // recompilation, even when its poles were already scanned.
                    pole_polynomials.push(factor.clone());
                    let remaining = unscanned_pole_polynomial(&square_free, &scanned_poles)?;
                    context.cancellation.check()?;
                    if remaining.is_constant() {
                        continue;
                    }
                    let coefficients = rounded_pole_coefficients(p, &remaining);
                    if coefficients.iter().any(|c| !p.finite(c)) {
                        return Err(Error::Numerical(format!(
                            "degree-{} pole geometry at {} bits: nonfinite pole polynomial",
                            remaining.degree(),
                            p.bits
                        )));
                    }
                    let roots = polynomial_roots(p, &coefficients, variable).map_err(|error| {
                        Error::Numerical(format!(
                            "degree-{} pole geometry at {} bits: {error}",
                            remaining.degree(),
                            p.bits
                        ))
                    })?;
                    context.cancellation.check()?;
                    scanned_poles.push((remaining, roots));
                }
                out.push(NumericRational {
                    numerator,
                    denominator,
                });
            }
            let (polynomial, exact_source) =
                self.compile_polynomial_row(row, p, values, &exact_values, context)?;
            polynomial_rows.push(polynomial);
            exact_source_rows.push(exact_source);
            matrix.push(out);
        }
        context.cancellation.check()?;
        resolve_coprime_pole_geometry(p, &mut scanned_poles, variable, context)?;
        // Keep geometrically refined centers: rounding them back to p can
        // collapse distinct exact poles. SeriesSystem still returns p for all
        // row/chart arithmetic and for source/boundary accuracy accounting.
        let poles = scanned_poles
            .into_iter()
            .flat_map(|(_, roots)| roots)
            .collect();
        Ok(CompiledSystem {
            source: Arc::new(ExactRows {
                variable,
                rows: self.source_rows.clone(),
                values: values.clone(),
            }),
            p,
            matrix,
            poles,
            pole_polynomials,
            polynomial_rows,
            exact_source_rows,
        })
    }

    /// Specialize the shared denominator-cleared rows without finding poles.
    /// Frobenius recurrences need these finite polynomials, but do not need a
    /// numerical root solve for each off-center denominator factor. In
    /// particular, specialization may create repeated off-center roots without
    /// obstructing a local series at the chosen regular-singular point.
    pub(crate) fn compile_polynomial_rows(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        context: &RunContext,
    ) -> Result<Vec<PolynomialRow>> {
        context.cancellation.check()?;
        let exact_values = source::ExactSpecialization::new(p, values)?;
        self.rows
            .iter()
            .map(|row| {
                context.cancellation.check()?;
                // Preserve the entry-wise numerical domain checks made by the
                // full compiler, in addition to the shared exact factor checks.
                for rational in &row.matrix {
                    context.cancellation.check()?;
                    let denominator =
                        polynomial_coefficients(&rational.denominator, self.variable, p, values)?;
                    if denominator.iter().all(|v| *v == p.zero()) {
                        return Err(Error::Numerical(
                            "identically zero specialized denominator".into(),
                        ));
                    }
                }
                self.compile_polynomial_row(row, p, values, &exact_values, context)
                    .map(|(polynomial, _)| polynomial)
            })
            .collect()
    }

    /// One coefficient-extraction owner for full transport and rows-only
    /// consumers. Exact stored-dyadic source values and sparse-zero decisions
    /// remain identical to the full compiler's original operation order.
    fn compile_polynomial_row(
        &self,
        row: &PreparedRow,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        exact_values: &source::ExactSpecialization,
        context: &RunContext,
    ) -> Result<(PolynomialRow, source::ExactPolynomialRow)> {
        let variable = self.variable;
        let denominator = polynomial_coefficients(&row.denominator, variable, p, values)?;
        let exact_denominator = exact_values.polynomial(&row.denominator, variable)?;
        let exact_denominator_factors = row
            .denominator_factors
            .iter()
            .map(|(factor, multiplicity)| {
                context.cancellation.check()?;
                Ok((exact_values.polynomial(factor, variable)?, *multiplicity))
            })
            .collect::<Result<Vec<_>>>()?;
        if exact_denominator_factors
            .iter()
            .any(|(factor, _)| factor.iter().all(SingleFloat::is_zero))
        {
            return Err(Error::Numerical(
                "identically zero specialized denominator factor".into(),
            ));
        }
        let mut exact_entries = Vec::new();
        let mut entries = Vec::new();
        for (column, cleared) in &row.cleared {
            context.cancellation.check()?;
            let coefficients = polynomial_coefficients(cleared, variable, p, values)?;
            let exact_coefficients = exact_values.polynomial(cleared, variable)?;
            if exact_coefficients.iter().any(|c| !c.is_zero()) {
                exact_entries.push((*column, exact_coefficients));
            }
            if coefficients.iter().any(|c| *c != p.zero()) {
                entries.push((*column, coefficients));
            }
        }
        Ok((
            PolynomialRow {
                denominator,
                entries,
            },
            source::ExactPolynomialRow {
                denominator: exact_denominator,
                denominator_factors: exact_denominator_factors,
                entries: exact_entries,
            },
        ))
    }
}

#[cfg(test)]
mod polynomial_rows_tests {
    use super::*;

    #[test]
    fn rows_only_specialization_matches_full_compiler() -> Result<()> {
        let context = RunContext::default();
        let p = Precision::decimal(80)?;
        let rows = PreparedRows::new(
            symbol!("x"),
            &[vec![parse!("1/(1-x)"), parse!("x/(2-x)")]],
            &context,
        )?;
        let ordinary = rows.compile(p, &Default::default(), &context)?;
        let finite = rows.compile_polynomial_rows(p, &Default::default(), &context)?;
        for (actual, expected) in finite.iter().zip(&ordinary.polynomial_rows) {
            assert_eq!(actual.denominator, expected.denominator);
            assert_eq!(actual.entries, expected.entries);
        }
        Ok(())
    }

    #[test]
    fn rows_only_accepts_specialized_repeated_roots_and_keeps_domain_checks() -> Result<()> {
        let context = RunContext::default();
        let p = Precision::decimal(80)?;
        let values = ahash::HashMap::from_iter([(parse!("a"), p.zero())]);
        let rows = PreparedRows::new(symbol!("x"), &[vec![parse!("1/((1-x)^2+a)")]], &context)?;
        let finite = rows.compile_polynomial_rows(p, &values, &context)?;
        assert_eq!(finite[0].denominator, vec![p.i(1), p.i(-2), p.i(1)]);
        let undefined = PreparedRows::new(symbol!("x"), &[vec![parse!("1/(a*(1-x))")]], &context)?;
        assert!(matches!(
            undefined.compile_polynomial_rows(p, &values, &context),
            Err(Error::Numerical(_))
        ));
        context.cancellation.cancel();
        assert!(rows.compile_polynomial_rows(p, &values, &context).is_err());
        Ok(())
    }
}

#[cfg(test)]
#[path = "repeated_poles_tests.rs"]
mod repeated_poles_tests;
