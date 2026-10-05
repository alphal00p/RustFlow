//! Exact reduction domains used by physical boundary reuse.
use crate::{
    DifferentialSystem, Error, Precision, Result, family::substitute, kinematics::KinematicPath,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;

type ExactPolynomial = MultivariatePolynomial<IntegerRing, u16>;

/// Validate the rational domain before canonical cancellation can erase a pole.
fn rational_expression(
    a: AtomView<'_>,
    variables: &BTreeSet<Symbol>,
    poles: &mut Vec<Atom>,
) -> Result<()> {
    match a {
        AtomView::Num(n) if matches!(n.get_coeff_view().to_owned(), Coefficient::Complex(_)) => Ok(()),
        AtomView::Var(v) if variables.contains(&v.get_symbol()) => Ok(()),
        AtomView::Add(v) => v.iter().try_for_each(|a| rational_expression(a, variables, poles)),
        AtomView::Mul(v) => v.iter().try_for_each(|a| rational_expression(a, variables, poles)),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            if let AtomView::Num(n) = exponent
                && let Coefficient::Complex(power) = n.get_coeff_view().to_owned()
                && power.im.is_zero() && power.re.is_integer()
            {
                if power.re < 0 { poles.push(base.to_owned()); }
                return rational_expression(base, variables, poles);
            }
            Err(Error::Unsupported("reduction conditions require integer rational powers".into()))
        }
        _ => Err(Error::Unsupported("reduction conditions require exact rational-complex expressions in declared physical variables and epsilon".into())),
    }
}

fn rational_parts(a: &Atom) -> Result<[Atom; 2]> {
    let rational: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(a)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::Unsupported(format!("rational reduction condition required: {e}")))?;
    let imaginary = BTreeMap::from([(
        Atom::var(crate::family::imaginary_parameter()),
        Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
    )]);
    Ok([
        rational.numerator.to_expression(),
        rational.denominator.to_expression(),
    ]
    .map(|part| substitute(&part, &imaginary).expand()))
}

/// All returned conditions are polynomials. Every original denominator is
/// retained independently before rational cancellation. Pure epsilon factors
/// remain valid generic conditions; epsilon is not set to zero here.
pub(crate) fn canonical_conditions(
    conditions: &[Atom],
    variables: &BTreeSet<Symbol>,
) -> Result<Vec<Atom>> {
    let mut candidates = conditions.to_vec();
    for condition in conditions {
        rational_expression(condition.as_view(), variables, &mut candidates)?;
    }
    let mut result = Vec::new();
    for candidate in candidates {
        for part in rational_parts(&candidate)? {
            if part.is_zero() {
                return Err(Error::InvalidInput(
                    "an identically zero reduction condition cannot define a physical domain"
                        .into(),
                ));
            }
            if matches!(part.as_view(), AtomView::Num(_)) {
                continue;
            }
            result.push(part);
        }
    }
    result.sort_by_cached_key(AtomCore::to_canonical_string);
    result.dedup();
    Ok(result)
}

/// Retain original rational denominator bases before a path expression is
/// simplified. Its value may be zero; only denominator domains are constrained.
pub(crate) fn rational_denominator_conditions(
    expressions: &[Atom],
    variables: &BTreeSet<Symbol>,
) -> Result<Vec<Atom>> {
    let mut poles = Vec::new();
    for expression in expressions {
        rational_expression(expression.as_view(), variables, &mut poles)?;
    }
    canonical_conditions(&poles, variables)
}

pub(crate) fn validate_conditions_at(
    conditions: &[Atom],
    epsilon: Symbol,
    coordinates: &BTreeMap<Symbol, Atom>,
) -> Result<()> {
    let substitutions = coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    for condition in conditions {
        let value = substitute(condition, &substitutions);
        rational_expression(value.as_view(), &BTreeSet::from([epsilon]), &mut Vec::new())?;
        if rational_parts(&value)?.iter().any(|a| a.is_zero()) {
            return Err(Error::InvalidInput(
                "physical point violates a required nonzero reduction condition".into(),
            ));
        }
    }
    Ok(())
}

/// Real path values annihilate an epsilon polynomial iff they annihilate all
/// real and imaginary epsilon coefficients. Native exact polynomial gcd finds
/// that common zero set without approximating epsilon or implementing Q(i).
fn common_real_zero_polynomial(
    polynomial: &Atom,
    epsilon: Symbol,
    parameter: Symbol,
) -> Result<Atom> {
    let mut gcd: Option<ExactPolynomial> = None;
    let variables = Arc::new(vec![PolyVariable::Symbol(parameter)]);
    for (_, coefficient) in
        crate::coefficient::exact_coefficient_list(&polynomial.expand(), &[Atom::var(epsilon)])?
    {
        let mut real = Atom::new();
        let mut imaginary = Atom::new();
        for (monomial, value) in crate::coefficient::exact_coefficient_list(
            &coefficient.expand(),
            &[Atom::var(parameter)],
        )? {
            let AtomView::Num(n) = value.as_view() else {
                return Err(Error::Unsupported(
                    "guarded physical paths require rational-complex coefficients".into(),
                ));
            };
            let Coefficient::Complex(value) = n.get_coeff_view().to_owned() else {
                return Err(Error::Unsupported(
                    "guarded physical paths require exact coefficients".into(),
                ));
            };
            real += &monomial * Atom::num(value.re);
            imaginary += &monomial * Atom::num(value.im);
        }
        for component in [real, imaginary] {
            if component.is_zero() {
                continue;
            }
            let rational: RationalPolynomial<IntegerRing, u16> = component
                .try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))
                .map_err(|e| {
                    Error::Unsupported(format!("polynomial physical path required: {e}"))
                })?;
            if !rational.denominator.is_constant() {
                return Err(Error::Unsupported(
                    "nonpolynomial condition coefficient after clearing denominators".into(),
                ));
            }
            gcd = Some(if let Some(previous) = gcd {
                previous.gcd(&rational.numerator)
            } else {
                rational.numerator
            });
            if gcd.as_ref().is_some_and(|p| p.is_constant()) {
                return Ok(Atom::one());
            }
        }
    }
    Ok(gcd.map_or_else(Atom::new, |polynomial| polynomial.to_expression()))
}

/// Admit the complete real parameter interval [0,1], including both endpoints.
/// Algebraic paths are rejected explicitly until their exact common zero set
/// can be represented by the native rational polynomial machinery.
pub(crate) fn conditions_admit_path(
    conditions: &[Atom],
    epsilon: Symbol,
    path: &KinematicPath,
    p: Precision,
    digits: u32,
) -> Result<bool> {
    path.validate()?;
    if conditions.is_empty() {
        return Ok(true);
    }
    if path.parameter == epsilon {
        return Err(Error::InvalidInput(
            "guard path parameter cannot be epsilon".into(),
        ));
    }
    let substitutions = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    let variables = BTreeSet::from([epsilon, path.parameter]);
    for condition in conditions {
        let restricted = substitute(condition, &substitutions);
        let mut candidates = vec![restricted.clone()];
        rational_expression(restricted.as_view(), &variables, &mut candidates)?;
        for candidate in candidates {
            for polynomial in rational_parts(&candidate)? {
                let common = common_real_zero_polynomial(&polynomial, epsilon, path.parameter)?;
                if common.is_zero() {
                    return Ok(false);
                }
                if common.is_one() {
                    continue;
                }
                for endpoint in [Atom::new(), Atom::one()] {
                    if substitute(
                        &common,
                        &BTreeMap::from([(Atom::var(path.parameter), endpoint)]),
                    )
                    .expand()
                    .is_zero()
                    {
                        return Ok(false);
                    }
                }
                let denominator_probe = DifferentialSystem {
                    variable: path.parameter,
                    matrix: vec![vec![Atom::one() / common]],
                };
                let compiled = denominator_probe.compile(p, &Default::default())?;
                let tolerance = p.tolerance(
                    digits
                        .checked_add(8)
                        .ok_or_else(|| Error::Limit("condition precision overflow".into()))?,
                );
                if compiled.poles.iter().any(|root| {
                    root.re >= -tolerance.clone()
                        && root.re <= p.real(1) + &tolerance
                        && root.im >= -tolerance.clone()
                        && root.im <= tolerance
                }) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

pub(crate) fn epsilon_leading_coefficient(polynomial: &Atom, epsilon: Symbol) -> Result<Atom> {
    let variable = Atom::var(epsilon);
    let mut leading: Option<(u32, Atom)> = None;
    for (monomial, coefficient) in crate::coefficient::exact_coefficient_list(
        &polynomial.expand(),
        std::slice::from_ref(&variable),
    )? {
        if coefficient.is_zero() {
            continue;
        }
        let exponent = if monomial.is_one() {
            0
        } else if monomial == variable {
            1
        } else if let AtomView::Pow(power) = monomial.as_view() {
            let (base, exponent) = power.get_base_exp();
            if base != variable.as_view() {
                return Err(Error::Unsupported(
                    "nonpolynomial epsilon denominator".into(),
                ));
            }
            exponent
                .to_string()
                .parse::<u32>()
                .map_err(|_| Error::Unsupported("noninteger epsilon denominator power".into()))?
        } else {
            return Err(Error::Unsupported(
                "nonpolynomial epsilon denominator".into(),
            ));
        };
        if leading.as_ref().is_none_or(|(old, _)| exponent < *old) {
            leading = Some((exponent, coefficient));
        }
    }
    Ok(leading.map_or_else(Atom::new, |(_, coefficient)| coefficient))
}

/// Original matrix poles and algebraic branch bases survive cancellation by a
/// stationary path coordinate. Their leading epsilon coefficient must stay
/// nonzero for one regular Laurent chart along the entire transport path.
/// This is stricter than a generic reduction guard such as s-epsilon: the
/// matrix 1/(s-epsilon) has coefficients 1/s, 1/s^2, ... and excludes s=0.
pub(crate) fn matrix_domain_conditions(
    system: &crate::kinematics::KinematicSystem,
) -> Result<Vec<Atom>> {
    fn bases(a: AtomView<'_>, out: &mut Vec<Atom>) -> Result<()> {
        match a {
            AtomView::Add(v) => v.iter().try_for_each(|a| bases(a, out)),
            AtomView::Mul(v) => v.iter().try_for_each(|a| bases(a, out)),
            AtomView::Pow(v) => {
                let (base, exponent) = v.get_base_exp();
                if let AtomView::Num(n) = exponent
                    && let Coefficient::Complex(value) = n.get_coeff_view().to_owned()
                    && value.im.is_zero()
                {
                    if value.re < 0 || !value.re.is_integer() {
                        out.push(base.to_owned());
                    }
                    return bases(base, out);
                }
                Err(Error::Unsupported(
                    "physical matrix powers require exact rational exponents".into(),
                ))
            }
            _ => Ok(()),
        }
    }
    let mut poles = Vec::new();
    for entry in system.derivatives.values().flatten().flatten() {
        bases(entry.as_view(), &mut poles)?;
    }
    let mut conditions = Vec::new();
    for pole in poles {
        for part in rational_parts(&pole)? {
            conditions.push(epsilon_leading_coefficient(&part, system.epsilon)?);
        }
    }
    Ok(conditions)
}
