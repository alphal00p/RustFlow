//! Immutable native rational-row algebra, before any numerical specialization.
use super::*;

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

/// Exact common factors have already been removed. Numerical overlap between
/// these coprime pieces is unresolved, so reject it instead of merging poles.
fn check_coprime_pole_separation(
    p: Precision,
    roots: &[C],
    previous: &[(ExactPolePolynomial, Vec<C>)],
) -> Result<()> {
    for (_, other_roots) in previous {
        for root in roots {
            for other in other_roots {
                let a = p.norm(root);
                let b = p.norm(other);
                let scale = if a < b { a } else { b };
                if p.norm(&p.sub(root, other)) <= p.tolerance(p.bits / 5) * scale {
                    return Err(Error::Numerical(
                        "ambiguous nearby poles from coprime specialized factors".into(),
                    ));
                }
            }
        }
    }
    Ok(())
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
        let mut poles = Vec::new();
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
                    let coefficients = remaining
                        .coefficients()
                        .iter()
                        .map(|coefficient| {
                            let mut rounded = p.rational(&coefficient.re);
                            rounded.im = p.rational(&coefficient.im).re;
                            rounded
                        })
                        .collect::<Vec<_>>();
                    if coefficients.iter().any(|c| !p.finite(c)) {
                        return Err(Error::Numerical("nonfinite pole polynomial".into()));
                    }
                    let roots = polynomial_roots(p, &coefficients, variable)?;
                    context.cancellation.check()?;
                    check_coprime_pole_separation(p, &roots, &scanned_poles)?;
                    // Exact coprimality and numerical separation were checked;
                    // no tolerance-based collector is needed for these pieces.
                    poles.extend(roots.iter().cloned());
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
