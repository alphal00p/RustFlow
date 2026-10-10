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

