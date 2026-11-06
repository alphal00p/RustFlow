//! Immutable native rational-row algebra, before any numerical specialization.
use super::*;
use crate::physical_conditions::PhysicalPhaseTiming;
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
    // Helper indices carry degrees/counts; the enclosing compile scope carries the
    // pole index. Timings never include expressions or coefficients.
    let repeated = {
        let _timing = PhysicalPhaseTiming::start(
            "ode.prepared.pole.square_free.gcd",
            monic.degree(),
            derivative.coefficients().len(),
            None,
        );
        monic.gcd_euclidean(&derivative)
    };
    let square_free = {
        let _timing = PhysicalPhaseTiming::start(
            "ode.prepared.pole.square_free.division",
            monic.degree(),
            repeated.degree(),
            None,
        );
        monic
            .try_div(&repeated)
            .ok_or_else(|| Error::Numerical("exact square-free pole division failed".into()))?
    };
    if {
        let _timing = PhysicalPhaseTiming::start(
            "ode.prepared.pole.square_free.verification",
            monic.degree(),
            square_free.degree(),
            None,
        );
        monic.clone().mul_coeff(&content) != polynomial
            || &repeated * &square_free != monic
            || derivative.try_div(&repeated).is_none()
            || repeated.degree().checked_add(square_free.degree()) != Some(monic.degree())
            || !square_free
                .gcd_euclidean(&square_free.derivative())
                .is_one()
    } {
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
    for (index, (scanned, _)) in previous.iter().enumerate() {
        let common = {
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.pole.shared.gcd",
                index,
                remaining.degree(),
                None,
            );
            remaining.gcd_euclidean(scanned)
        };
        if common.is_zero() {
            return Err(Error::Numerical("zero common pole factor".into()));
        }
        // Ordinary quot_rem handles constants and degree collapse; never use
        // quot_rem_field, whose degree assumptions do not cover these cases.
        let (next, remainder, scanned_quotient) = {
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.pole.shared.removal",
                index,
                common.degree(),
                None,
            );
            let (next, remainder) = remaining.quot_rem(&common);
            let scanned_quotient = scanned.try_div(&common).ok_or_else(|| {
                Error::Numerical("exact common pole factor division failed".into())
            })?;
            (next, remainder, scanned_quotient)
        };
        if {
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.pole.shared.verification",
                index,
                common.degree(),
                None,
            );
            !remainder.is_zero()
                || &common * &next != remaining
                || &common * &scanned_quotient != *scanned
                || common.degree().checked_add(next.degree()) != Some(remaining.degree())
                || common.degree().checked_add(scanned_quotient.degree()) != Some(scanned.degree())
        } {
            return Err(Error::Numerical(
                "exact common pole factor verification failed".into(),
            ));
        }
        common_product = &common_product * &common;
        remaining = next;
    }
    if {
        let _timing = PhysicalPhaseTiming::start(
            "ode.prepared.pole.shared.final_verification",
            original.degree(),
            previous.len(),
            None,
        );
        remaining.is_zero()
            || !remaining.lcoeff().is_one()
            || &common_product * &remaining != *original
            || common_product.degree().checked_add(remaining.degree()) != Some(original.degree())
            || !remaining.gcd_euclidean(&remaining.derivative()).is_one()
            || previous.iter().enumerate().any(|(index, (scanned, _))| {
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.pole.shared.verify_coprime",
                    index,
                    scanned.degree(),
                    None,
                );
                !remaining.gcd_euclidean(scanned).is_one()
            })
    } {
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
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.pole.root_retry",
                    i,
                    polynomial.degree(),
                    Some(geometry.bits),
                );
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

/// Certify the actual LCM using already checked entry factors. Exact division
/// consumes associates (including factors with different variable maps) only
/// once. Integer content is the final residual, never a maximum of entry content.
fn reuse_row_factor_certificate(
    denominator: &ExactPolynomial,
    entries: &[DenominatorFactors],
    context: &RunContext,
) -> Result<Option<DenominatorFactors>> {
    if denominator.is_zero() {
        return Err(Error::Numerical("zero row denominator".into()));
    }
    let mut seen = ahash::HashSet::default();
    let mut candidates = Vec::new();
    for entry in entries {
        context.cancellation.check()?;
        for (factor, _) in entry.iter() {
            if !factor.is_constant() && seen.insert(factor.to_expression()) {
                candidates.push(factor);
            }
        }
    }
    let mut remaining = denominator.clone();
    let mut factors = Vec::new();
    for factor in candidates {
        context.cancellation.check()?;
        // A positive-degree witness proves progress on every successful divide,
        // independently of associates and variable ordering. Constants must not
        // enter this loop: units would never terminate, and 2/4 cannot max-merge.
        let (variable, degree) = factor
            .variables()
            .iter()
            .enumerate()
            .find_map(|(i, variable)| {
                let degree = factor.degree(i);
                (degree > 0).then_some((variable, degree))
            })
            .ok_or_else(|| Error::Numerical("nonconstant row factor has no degree".into()))?;
        let degree_in = |polynomial: &ExactPolynomial| {
            polynomial
                .variables()
                .iter()
                .position(|v| v == variable)
                .map_or(0, |i| polynomial.degree(i))
        };
        let mut multiplicity = 0_usize;
        loop {
            context.cancellation.check()?;
            let Some(quotient) = remaining.try_div(factor) else {
                break;
            };
            if degree_in(&remaining).checked_sub(degree) != Some(degree_in(&quotient)) {
                return Err(Error::Numerical(
                    "exact row factor division degree failed".into(),
                ));
            }
            multiplicity = multiplicity
                .checked_add(1)
                .ok_or_else(|| Error::Limit("row factor multiplicity overflow".into()))?;
            remaining = quotient;
        }
        if multiplicity != 0 {
            factors.push((factor.clone(), multiplicity));
        }
    }
    if !remaining.is_constant() {
        // Incomplete candidates are not a pole/domain certificate.
        return Ok(None);
    }
    if !remaining.is_one() {
        factors.push((remaining, 1));
    }
    let mut product = denominator.one();
    for (factor, multiplicity) in &factors {
        context.cancellation.check()?;
        if *multiplicity == 0
            || factor.variables().iter().enumerate().any(|(i, _)| {
                usize::from(factor.degree(i))
                    .checked_mul(*multiplicity)
                    .is_none_or(|d| d > usize::from(u16::MAX))
            })
        {
            return Err(Error::Limit("denominator factor degree overflow".into()));
        }
        product = multiply_polynomials(&product, &factor.pow(*multiplicity))?;
    }
    // Polynomial equality is structural in Symbolica: unify maps for the exact
    // comparison without changing the LCM used by denominator clearing.
    let mut expected = denominator.clone();
    expected.unify_variables(&mut product);
    if product != expected {
        return Err(Error::Numerical(
            "exact row denominator factor reconstruction failed".into(),
        ));
    }
    context.cancellation.check()?;
    Ok(Some(Arc::new(factors)))
}

fn row_denominator_factors(
    denominator: &ExactPolynomial,
    entries: &[DenominatorFactors],
    cache: &mut ahash::HashMap<ExactPolynomial, DenominatorFactors>,
    context: &RunContext,
    row_index: usize,
) -> Result<DenominatorFactors> {
    context.cancellation.check()?;
    if let Some(factors) = cache.get(denominator) {
        return Ok(factors.clone());
    }
    let reused = {
        let _timing = PhysicalPhaseTiming::start(
            "ode.prepared.new.row_factor_reuse_certification",
            row_index,
            entries.len(),
            None,
        );
        reuse_row_factor_certificate(denominator, entries, context)?
    };
    if let Some(factors) = reused {
        // Only a complete, exactly reconstructed certificate enters the cache.
        cache.insert(denominator.clone(), factors.clone());
        return Ok(factors);
    }
    context.cancellation.check()?;
    let _timing = PhysicalPhaseTiming::start(
        "ode.prepared.new.row_factor_reuse_fallback",
        row_index,
        entries.len(),
        None,
    );
    denominator_factors(denominator, cache)
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
        for (row_index, row) in rows.iter().enumerate() {
            let _row_timing =
                PhysicalPhaseTiming::start("ode.prepared.new.row", row_index, row.len(), None);
            context.cancellation.check()?;
            let mut matrix = Vec::with_capacity(row.len());
            let mut exact = Vec::with_capacity(row.len());
            let mut entry_factors = Vec::with_capacity(row.len());
            let mut common_denominator: Option<ExactPolynomial> = None;
            for (column, value) in row.iter().enumerate() {
                context.cancellation.check()?;
                let rational: RationalPolynomial<IntegerRing, u16> = {
                    let _timing = PhysicalPhaseTiming::start(
                        "ode.prepared.new.rational_conversion",
                        row_index,
                        column,
                        None,
                    );
                    crate::family::encode_complex(value)
                        .try_to_rational_polynomial(&Q, &Z, None)
                        .map_err(|e| Error::InvalidInput(e.to_string()))?
                };
                context.cancellation.check()?;
                let factors = {
                    let _timing = PhysicalPhaseTiming::start(
                        "ode.prepared.new.denominator_factorization",
                        row_index,
                        column,
                        None,
                    );
                    denominator_factors(&rational.denominator, &mut factorizations)?
                };
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
                entry_factors.push(factors);
                common_denominator = Some({
                    let _timing = PhysicalPhaseTiming::start(
                        "ode.prepared.new.row_lcm",
                        row_index,
                        column,
                        None,
                    );
                    if let Some(previous) = common_denominator {
                        let quotient = previous
                            .try_div(&previous.gcd(&rational.denominator))
                            .ok_or_else(|| {
                                Error::Numerical("exact denominator LCM division failed".into())
                            })?;
                        multiply_polynomials(&quotient, &rational.denominator)?
                    } else {
                        rational.denominator.clone()
                    }
                });
                context.cancellation.check()?;
                exact.push(rational);
            }
            let common_denominator = common_denominator.unwrap();
            let denominator_factors = {
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.new.row_denominator_factorization",
                    row_index,
                    row.len(),
                    None,
                );
                row_denominator_factors(
                    &common_denominator,
                    &entry_factors,
                    &mut factorizations,
                    context,
                    row_index,
                )?
                .iter()
                .map(|(factor, multiplicity)| (factor.to_expression(), *multiplicity))
                .collect()
            };
            context.cancellation.check()?;
            let mut cleared = Vec::new();
            for (column, rational) in exact.iter().enumerate() {
                context.cancellation.check()?;
                if rational.numerator.is_zero() {
                    continue;
                }
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.new.denominator_clearing",
                    row_index,
                    column,
                    None,
                );
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
        for (row_index, row) in self.rows.iter().enumerate() {
            let _row_timing = PhysicalPhaseTiming::start(
                "ode.prepared.compile.row",
                row_index,
                row.matrix.len(),
                Some(p.bits),
            );
            context.cancellation.check()?;
            let mut out = Vec::with_capacity(row.matrix.len());
            for (column, rational) in row.matrix.iter().enumerate() {
                let _entry_timing = PhysicalPhaseTiming::start(
                    "ode.prepared.compile.entry",
                    row_index,
                    column,
                    Some(p.bits),
                );
                context.cancellation.check()?;
                let (numerator, denominator) = {
                    let _timing = PhysicalPhaseTiming::start(
                        "ode.prepared.compile.rational_coefficients",
                        row_index,
                        column,
                        Some(p.bits),
                    );
                    let numerator =
                        polynomial_coefficients(&rational.numerator, variable, p, values)?;
                    let denominator =
                        polynomial_coefficients(&rational.denominator, variable, p, values)?;
                    (numerator, denominator)
                };
                if denominator.iter().all(|v| *v == p.zero()) {
                    return Err(Error::Numerical(
                        "identically zero specialized denominator".into(),
                    ));
                }
                for (pole_index, factor) in rational.new_pole_factors.iter().enumerate() {
                    // Pole indices are local to the enclosing row/column entry.
                    let _pole_timing = PhysicalPhaseTiming::start(
                        "ode.prepared.compile.pole",
                        pole_index,
                        rational.new_pole_factors.len(),
                        Some(p.bits),
                    );
                    context.cancellation.check()?;
                    // Specialize J=i and stored parameter dyadics before any
                    // rounding: specialization can create exact multiplicities
                    // invisible to the formal multivariate factorization.
                    let coefficients = {
                        let _timing = PhysicalPhaseTiming::start(
                            "ode.prepared.compile.pole_specialization",
                            pole_index,
                            0,
                            Some(p.bits),
                        );
                        exact_values.polynomial(factor, variable)?
                    };
                    let square_free = {
                        let _timing = PhysicalPhaseTiming::start(
                            "ode.prepared.compile.pole_square_free",
                            pole_index,
                            coefficients.len(),
                            None,
                        );
                        square_free_pole_polynomial(coefficients, variable)?
                    };
                    context.cancellation.check()?;
                    if square_free.is_constant() {
                        continue;
                    }
                    // Retain every original nonconstant factor for source-domain
                    // recompilation, even when its poles were already scanned.
                    pole_polynomials.push(factor.clone());
                    let remaining = {
                        let _timing = PhysicalPhaseTiming::start(
                            "ode.prepared.compile.pole_shared_roots",
                            pole_index,
                            scanned_poles.len(),
                            None,
                        );
                        unscanned_pole_polynomial(&square_free, &scanned_poles)?
                    };
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
                    let roots = {
                        let _timing = PhysicalPhaseTiming::start(
                            "ode.prepared.compile.pole_roots",
                            pole_index,
                            remaining.degree(),
                            Some(p.bits),
                        );
                        polynomial_roots(p, &coefficients, variable).map_err(|error| {
                            Error::Numerical(format!(
                                "degree-{} pole geometry at {} bits: {error}",
                                remaining.degree(),
                                p.bits
                            ))
                        })?
                    };
                    context.cancellation.check()?;
                    scanned_poles.push((remaining, roots));
                }
                out.push(NumericRational {
                    numerator,
                    denominator,
                });
            }
            let (polynomial, exact_source) = {
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.compile.polynomial_row",
                    row_index,
                    row.cleared.len(),
                    Some(p.bits),
                );
                self.compile_polynomial_row(row, p, values, &exact_values, context)?
            };
            polynomial_rows.push(polynomial);
            exact_source_rows.push(exact_source);
            matrix.push(out);
        }
        context.cancellation.check()?;
        {
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.compile.pole_geometry",
                scanned_poles.len(),
                0,
                Some(p.bits),
            );
            resolve_coprime_pole_geometry(p, &mut scanned_poles, variable, context)?;
        }
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
            .enumerate()
            .map(|(row_index, row)| {
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.compile.polynomial_only_row",
                    row_index,
                    row.matrix.len(),
                    Some(p.bits),
                );
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
                {
                    let _timing = PhysicalPhaseTiming::start(
                        "ode.prepared.compile.polynomial_row",
                        row_index,
                        row.cleared.len(),
                        Some(p.bits),
                    );
                    self.compile_polynomial_row(row, p, values, &exact_values, context)
                        .map(|(polynomial, _)| polynomial)
                }
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
        let (denominator, exact_denominator) = {
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.compile.row_denominator",
                0,
                0,
                Some(p.bits),
            );
            let denominator = polynomial_coefficients(&row.denominator, variable, p, values)?;
            let exact_denominator = exact_values.polynomial(&row.denominator, variable)?;
            (denominator, exact_denominator)
        };
        let exact_denominator_factors = row
            .denominator_factors
            .iter()
            .enumerate()
            .map(|(index, (factor, multiplicity))| {
                let _timing = PhysicalPhaseTiming::start(
                    "ode.prepared.compile.row_denominator_factor",
                    index,
                    row.denominator_factors.len(),
                    Some(p.bits),
                );
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
            let _timing = PhysicalPhaseTiming::start(
                "ode.prepared.compile.cleared_entry",
                *column,
                row.cleared.len(),
                Some(p.bits),
            );
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
mod row_factor_reuse_tests {
    use super::*;

    fn polynomial(expression: &Atom) -> ExactPolynomial {
        let rational: RationalPolynomial<IntegerRing, u16> =
            crate::family::encode_complex(expression)
                .try_to_rational_polynomial(&Q, &Z, None)
                .unwrap();
        assert!(rational.denominator.is_one());
        rational.numerator
    }

    fn assert_same_polynomial(mut actual: ExactPolynomial, mut expected: ExactPolynomial) {
        actual.unify_variables(&mut expected);
        assert_eq!(actual, expected);
    }

    fn assert_factors(factors: &[(ExactPolynomial, usize)], expected: &Atom) {
        let expected = polynomial(expected);
        let mut product = expected.one();
        for (factor, multiplicity) in factors {
            assert!(*multiplicity > 0);
            product = multiply_polynomials(&product, &factor.pow(*multiplicity)).unwrap();
        }
        assert_same_polynomial(product, expected);
    }

    fn assert_row(
        entries: Vec<Atom>,
        denominator: Atom,
        cleared: Vec<(usize, Atom)>,
    ) -> Result<PreparedRows> {
        let prepared = PreparedRows::new(
            symbol!("row_factor_reuse_x"),
            &[entries.clone()],
            &RunContext::default(),
        )?;
        assert_eq!(prepared.source_rows, vec![entries]);
        let row = &prepared.rows[0];
        assert_same_polynomial(polynomial(&row.denominator), polynomial(&denominator));
        assert_eq!(row.cleared.len(), cleared.len());
        for ((column, actual), (expected_column, expected)) in row.cleared.iter().zip(&cleared) {
            assert_eq!(column, expected_column);
            assert_same_polynomial(polynomial(actual), polynomial(expected));
        }
        let factors = row
            .denominator_factors
            .iter()
            .map(|(factor, multiplicity)| (polynomial(factor), *multiplicity))
            .collect::<Vec<_>>();
        assert_factors(&factors, &denominator);
        Ok(prepared)
    }

    #[test]
    fn row_factor_reuse_overlapping_powers_preserve_exact_row() -> Result<()> {
        let prepared = assert_row(
            vec![
                parse!("1/(2*(row_factor_reuse_x-1)^2*(row_factor_reuse_x+2))"),
                parse!("3/(4*(row_factor_reuse_x-1)*(row_factor_reuse_x+2)^3)"),
            ],
            parse!("4*(row_factor_reuse_x-1)^2*(row_factor_reuse_x+2)^3"),
            vec![
                (0, parse!("2*(row_factor_reuse_x+2)^2")),
                (1, parse!("3*(row_factor_reuse_x-1)")),
            ],
        )?;
        let row = &prepared.rows[0];
        for (entry, numerator, denominator) in [
            (
                &row.matrix[0],
                Atom::one(),
                parse!("2*(row_factor_reuse_x-1)^2*(row_factor_reuse_x+2)"),
            ),
            (
                &row.matrix[1],
                Atom::num(3),
                parse!("4*(row_factor_reuse_x-1)*(row_factor_reuse_x+2)^3"),
            ),
        ] {
            assert_same_polynomial(polynomial(&entry.numerator), polynomial(&numerator));
            assert_same_polynomial(polynomial(&entry.denominator), polynomial(&denominator));
        }
        for (expected, power) in [
            (parse!("row_factor_reuse_x-1"), 2),
            (parse!("row_factor_reuse_x+2"), 3),
        ] {
            assert!(
                row.denominator_factors
                    .iter()
                    .any(|(factor, multiplicity)| {
                        *multiplicity == power && polynomial(factor).to_expression() == expected
                    })
            );
        }
        Ok(())
    }

    #[test]
    fn row_factor_reuse_integer_content_is_not_max_merged() -> Result<()> {
        assert_row(
            vec![parse!("1/2"), parse!("1/4")],
            Atom::num(4),
            vec![(0, Atom::num(2)), (1, Atom::one())],
        )?;
        // The PreparedRows constant-only case can use an entry cache hit; also
        // exercise certification directly so that 2 and 4 are not multiplied
        // or merged as unrelated constant-factor keys on a cache miss.
        let constants = vec![
            Arc::new(vec![(polynomial(&Atom::num(2)), 1)]),
            Arc::new(vec![(polynomial(&Atom::num(4)), 1)]),
        ];
        let factors = reuse_row_factor_certificate(
            &polynomial(&Atom::num(4)),
            &constants,
            &RunContext::default(),
        )?
        .expect("constant residual is the complete content");
        assert_factors(&factors, &Atom::num(4));
        assert_eq!(factors.len(), 1);
        let prepared = assert_row(
            vec![
                parse!("1/(2*(row_factor_reuse_x-1))"),
                parse!("1/(3*(row_factor_reuse_x+1))"),
                parse!("1/5"),
            ],
            parse!("30*(row_factor_reuse_x-1)*(row_factor_reuse_x+1)"),
            vec![
                (0, parse!("15*(row_factor_reuse_x+1)")),
                (1, parse!("10*(row_factor_reuse_x-1)")),
                (2, parse!("6*(row_factor_reuse_x-1)*(row_factor_reuse_x+1)")),
            ],
        )?;
        assert!(
            prepared.rows[0]
                .denominator_factors
                .iter()
                .any(|(factor, multiplicity)| { *multiplicity == 1 && *factor == Atom::num(30) })
        );
        Ok(())
    }

    #[test]
    fn row_factor_reuse_signs_associates_and_variable_maps() -> Result<()> {
        let context = RunContext::default();
        let f = parse!("row_factor_reuse_x+row_factor_reuse_y");
        let negative = polynomial(&(-&f));
        let variables = Arc::new(polynomial(&f).variables().iter().rev().cloned().collect());
        let reordered: RationalPolynomial<IntegerRing, u16> = f
            .try_to_rational_polynomial(&Q, &Z, Some(variables))
            .unwrap();
        assert_ne!(reordered.numerator.variables(), polynomial(&f).variables());
        // Complete entry certificates with opposite associates and different
        // maps. Encounter order deliberately consumes the negative associate.
        let entries = vec![
            Arc::new(vec![(negative.clone(), 1), (polynomial(&Atom::num(2)), 1)]),
            Arc::new(vec![
                (reordered.numerator, 2),
                (polynomial(&Atom::num(3)), 1),
            ]),
            Arc::new(vec![(polynomial(&parse!("row_factor_reuse_x-2")), 2)]),
        ];
        let expected =
            parse!("-60*(row_factor_reuse_x+row_factor_reuse_y)^3*(row_factor_reuse_x-2)^2");
        let factors = reuse_row_factor_certificate(&polynomial(&expected), &entries, &context)?
            .expect("complete candidates");
        assert_factors(&factors, &expected);
        assert_same_polynomial(factors[0].0.clone(), negative);
        assert_eq!(factors[0].1, 3);
        assert_eq!(factors.iter().filter(|(f, _)| !f.is_constant()).count(), 2);
        assert_same_polynomial(
            factors.last().unwrap().0.clone(),
            polynomial(&Atom::num(60)),
        );
        let mut reversed = entries.clone();
        reversed.swap(0, 1);
        let factors = reuse_row_factor_certificate(&polynomial(&expected), &reversed, &context)?
            .expect("reordered associate also certifies the same denominator");
        assert_factors(&factors, &expected);
        assert_eq!(factors[0].1, 3);
        assert_same_polynomial(factors[0].0.clone(), polynomial(&f));
        assert_eq!(factors.iter().filter(|(f, _)| !f.is_constant()).count(), 2);
        assert_same_polynomial(
            factors.last().unwrap().0.clone(),
            polynomial(&Atom::num(-60)),
        );

        let expected = parse!("-42*(row_factor_reuse_x-1)^3");
        let factors = reuse_row_factor_certificate(
            &polynomial(&expected),
            &[Arc::new(vec![(
                polynomial(&parse!("row_factor_reuse_x-1")),
                1,
            )])],
            &context,
        )?
        .expect("negative content is retained");
        assert_factors(&factors, &expected);
        assert_same_polynomial(
            factors.last().unwrap().0.clone(),
            polynomial(&Atom::num(-42)),
        );
        assert_row(
            vec![
                parse!("1/(1-row_factor_reuse_x)"),
                parse!("1/(row_factor_reuse_x-1)^2"),
            ],
            parse!("(row_factor_reuse_x-1)^2"),
            vec![(0, parse!("1-row_factor_reuse_x")), (1, Atom::one())],
        )?;
        Ok(())
    }

    #[test]
    fn row_factor_reuse_zero_and_constant_entries() -> Result<()> {
        assert_row(vec![Atom::Zero; 3], Atom::one(), Vec::new())?;
        assert_row(
            vec![Atom::Zero, Atom::num(5), parse!("-3/4")],
            Atom::num(4),
            vec![(1, Atom::num(20)), (2, Atom::num(-3))],
        )?;
        assert_row(
            vec![
                Atom::Zero,
                parse!("1/(row_factor_reuse_x-1)"),
                Atom::num(-2),
            ],
            parse!("row_factor_reuse_x-1"),
            vec![(1, Atom::one()), (2, parse!("-2*(row_factor_reuse_x-1)"))],
        )?;
        Ok(())
    }

    #[test]
    fn row_factor_reuse_multivariate_formal_imaginary_parameter() -> Result<()> {
        let x = Atom::var(symbol!("row_factor_reuse_x"));
        let a = Atom::var(symbol!("row_factor_reuse_a"));
        let j = Atom::var(crate::family::imaginary_parameter());
        let f = &x + &j * &a;
        let g = &x + &a;
        let prepared = assert_row(
            vec![
                Atom::one() / (Atom::num(2) * f.clone().pow(2)),
                Atom::one() / (Atom::num(3) * &f * &g),
            ],
            Atom::num(6) * f.clone().pow(2) * &g,
            vec![(0, Atom::num(3) * &g), (1, Atom::num(2) * &f)],
        )?;
        assert!(
            prepared.rows[0]
                .denominator
                .contains_symbol(crate::family::imaginary_parameter())
        );
        Ok(())
    }

    #[test]
    fn row_factor_reuse_incomplete_candidates_fall_back_and_cache_is_preferred() -> Result<()> {
        let context = RunContext::default();
        let expected = parse!("6*(row_factor_reuse_x-1)^2*(row_factor_reuse_x+3)");
        let denominator = polynomial(&expected);
        let incomplete = vec![Arc::new(vec![(
            polynomial(&parse!("row_factor_reuse_x-1")),
            1,
        )])];
        assert!(reuse_row_factor_certificate(&denominator, &incomplete, &context)?.is_none());
        let mut cache = ahash::HashMap::default();
        let fallback = row_denominator_factors(&denominator, &incomplete, &mut cache, &context, 0)?;
        assert_factors(&fallback, &expected);
        assert_eq!(cache.len(), 1);
        let cached = row_denominator_factors(&denominator, &[], &mut cache, &context, 0)?;
        assert!(Arc::ptr_eq(&fallback, &cached));

        // A verified reuse certificate is cached as well; no common-denominator
        // factorization is needed even when it was absent from the entry cache.
        cache.clear();
        let complete = vec![
            incomplete[0].clone(),
            Arc::new(vec![(polynomial(&parse!("row_factor_reuse_x+3")), 1)]),
        ];
        let reused = row_denominator_factors(&denominator, &complete, &mut cache, &context, 0)?;
        assert_factors(&reused, &expected);
        assert!(Arc::ptr_eq(&reused, cache.get(&denominator).unwrap()));
        let cached = row_denominator_factors(&denominator, &[], &mut cache, &context, 0)?;
        assert!(Arc::ptr_eq(&reused, &cached));
        assert!(
            reuse_row_factor_certificate(&polynomial(&Atom::Zero), &complete, &context).is_err()
        );
        context.cancellation.cancel();
        assert!(row_denominator_factors(&denominator, &complete, &mut cache, &context, 0).is_err());
        Ok(())
    }
}

#[cfg(test)]
#[path = "repeated_poles_tests.rs"]
mod repeated_poles_tests;
