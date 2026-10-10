//! Exact literal coefficient extraction using Symbolica's polynomial owners.
//!
//! Symbolica's convenience `coefficient_list` currently selects statistical
//! zero testing internally. An inconclusive floating-point test can then drop
//! exact coefficients outside its numerical range. Select an exact expression
//! field explicitly and let the native polynomial owner perform the grouping.
use crate::{Error, Result};
use std::sync::Arc;
use symbolica::prelude::*;

const MAX_COEFFICIENT_EXPONENT: u32 = 100_000;

fn exponent_limit() -> Error {
    Error::Limit(format!(
        "coefficient conversion degree exceeds {MAX_COEFFICIENT_EXPONENT}"
    ))
}

/// Bound the degree that native conversion creates from a power. Positive
/// integer powers expand their base; other powers are opaque coordinates, with
/// the integer content of their exponent becoming a polynomial exponent.
fn power_degree(base: AtomView<'_>, exponent: AtomView<'_>) -> Result<(u32, u32)> {
    if let AtomView::Add(terms) = exponent {
        let mut expanded = 0_u32;
        let mut atomic = 0_u32;
        for term in terms {
            let (e, a) = power_degree(base, term)?;
            expanded = expanded.checked_add(e).ok_or_else(exponent_limit)?;
            atomic = atomic.checked_add(a).ok_or_else(exponent_limit)?;
        }
        if u64::from(expanded) + u64::from(atomic) > u64::from(MAX_COEFFICIENT_EXPONENT) {
            return Err(exponent_limit());
        }
        return Ok((expanded, atomic));
    }
    if let Ok(integer) = Integer::try_from(exponent) {
        // A negative power of a sum/product is one independent coordinate.
        // Its base is never expanded by native polynomial conversion.
        if integer < 0 && PolyVariable::try_from(base.to_owned()).is_err() {
            return Ok((0, 1));
        }
        let degree = u32::try_from(integer.abs()).map_err(|_| exponent_limit())?;
        if degree > MAX_COEFFICIENT_EXPONENT {
            return Err(exponent_limit());
        }
        return Ok(if integer < 0 {
            (0, degree)
        } else {
            (degree, 0)
        });
    }
    let coefficient = match exponent {
        AtomView::Num(n) => Some(n.get_coeff_view()),
        AtomView::Mul(m) => m.get_coefficient().and_then(|coefficient| {
            if let AtomView::Num(n) = coefficient {
                Some(n.get_coeff_view())
            } else {
                None
            }
        }),
        _ => None,
    };
    if let Some(CoefficientView::Natural(real, _, imaginary, _)) = coefficient {
        let content = symbolica::domains::integer::gcd_signed(real, imaginary);
        // Native conversion only extracts signed integer content strictly
        // within the i32 range. Larger nonintegral powers remain atomic.
        let negative = real < 0 || real == 0 && imaginary < 0;
        if content < i32::MAX as u64 || negative && content == i32::MAX as u64 {
            let degree = content.max(1) as u32;
            if degree > MAX_COEFFICIENT_EXPONENT {
                return Err(exponent_limit());
            }
            return Ok((0, degree));
        }
    }
    Ok((0, 1))
}

/// Check exponent growth before allocating a native polynomial. The native
/// fast converter accepts u32 input exponents but stores them in i32, so a
/// post-conversion check alone cannot detect positive-to-negative wrapping.
/// Bound all converted coordinates, including unrequested coefficient symbols.
fn check_conversion_degree(expression: &Atom) -> Result<()> {
    enum Plan {
        Constant,
        Atomic,
        Maximum,
        Sum,
        Power {
            base: Atom,
            expanded: u32,
            atomic: u32,
        },
    }
    let mut pending = vec![expression.clone()];
    let mut nodes = Vec::new();
    let mut error = None;
    while let Some(expression) = pending.pop() {
        expression.visitor(&mut |node| {
            if error.is_some() {
                return false;
            }
            let plan = match node {
                AtomView::Num(_) => Plan::Constant,
                AtomView::Var(_) | AtomView::Fun(_) => Plan::Atomic,
                AtomView::Add(_) => Plan::Maximum,
                AtomView::Mul(_) => Plan::Sum,
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let (expanded, atomic) = match power_degree(base, exponent) {
                        Ok(degree) => degree,
                        Err(failure) => {
                            error = Some(failure);
                            return false;
                        }
                    };
                    if expanded != 0 {
                        pending.push(base.to_owned());
                    }
                    Plan::Power {
                        base: base.to_owned(),
                        expanded,
                        atomic,
                    }
                }
            };
            let descend = matches!(plan, Plan::Maximum | Plan::Sum);
            nodes.push((node.to_owned(), plan));
            // Function arguments and nonpolynomial power bases are opaque.
            // Expanded power bases are visited separately, without entering
            // exponent expressions that native conversion keeps atomic.
            descend
        });
    }
    if let Some(error) = error {
        return Err(error);
    }
    let mut degrees = ahash::HashMap::<Atom, u64>::default();
    for (node, plan) in nodes.into_iter().rev() {
        let degree = match plan {
            Plan::Constant => 0,
            Plan::Atomic => 1,
            Plan::Maximum => node
                .children()
                .map(|child| degrees[&child.to_owned()])
                .max()
                .unwrap_or(0),
            Plan::Sum => node
                .children()
                .map(|child| degrees[&child.to_owned()])
                .sum(),
            Plan::Power {
                base,
                expanded,
                atomic,
            } => {
                u64::from(atomic)
                    + if expanded == 0 {
                        0
                    } else {
                        u64::from(expanded) * degrees[&base]
                    }
            }
        };
        if degree > u64::from(MAX_COEFFICIENT_EXPONENT) {
            return Err(exponent_limit());
        }
        degrees.insert(node, degree);
    }
    Ok(())
}

/// Collect literal powers of the requested symbols or function/power markers.
///
/// Signed integral powers are retained, so this also groups Laurent monomials.
/// A coefficient may still contain a requested symbol inside an unrequested
/// function or power: collecting `x` in `x + f(x)` returns `f(x)` as the constant
/// coefficient. Callers that require polynomial dependence must check for that
/// hidden dependence and reject negative powers themselves.
///
/// Native conversion may expand polynomial expressions in other variables.
/// If its rational-power basis change removes a requested marker, reject that
/// conversion instead of silently changing the requested literal coordinates.
pub(crate) fn exact_coefficient_list(
    expression: &Atom,
    variables: &[Atom],
) -> Result<Vec<(Atom, Atom)>> {
    if variables.is_empty() {
        return Ok(if expression.is_zero() {
            Vec::new()
        } else {
            vec![(Atom::num(1), expression.clone())]
        });
    }
    let mut requested = Vec::<PolyVariable>::with_capacity(variables.len());
    for variable in variables {
        let variable = PolyVariable::try_from(variable.clone())
            .map_err(|error| Error::InvalidInput(format!("invalid coefficient marker: {error}")))?;
        if requested.contains(&variable) {
            return Err(Error::InvalidInput(
                "coefficient markers must be distinct".into(),
            ));
        }
        requested.push(variable);
    }
    check_conversion_degree(expression)?;
    let field = AtomField {
        statistical_zero_test: false,
        cancel_check_on_division: false,
        custom_normalization: None,
    };
    let polynomial = expression
        .try_to_polynomial::<_, i32>(&field, Some(Arc::new(requested.clone())))
        .map_err(|error| {
            Error::InvalidInput(format!("exact coefficient conversion failed: {error}"))
        })?;
    let indices = requested
        .iter()
        .map(|variable| {
            polynomial
                .variables()
                .iter()
                .position(|candidate| candidate == variable)
                .ok_or_else(|| {
                    Error::Unsupported(format!(
                        "polynomial conversion changed coefficient marker {variable}"
                    ))
                })
        })
        .collect::<Result<Vec<_>>>()?;
    for exponents in polynomial.exponents_iter() {
        if exponents
            .iter()
            .any(|exponent| exponent.unsigned_abs() > MAX_COEFFICIENT_EXPONENT)
        {
            return Err(Error::Limit(format!(
                "coefficient exponent magnitude exceeds {MAX_COEFFICIENT_EXPONENT}"
            )));
        }
    }
    // Native grouping is sparse, including for signed exponents. Sort the
    // native hash-map result in the caller's marker order for reproducibility.
    let mut grouped = polynomial
        .to_multivariate_polynomial_list(&indices, true)
        .into_iter()
        .map(|(exponents, coefficient)| {
            let degrees = indices
                .iter()
                .map(|&index| exponents[index])
                .collect::<Vec<_>>();
            (degrees, coefficient.flatten(false))
        })
        .collect::<Vec<_>>();
    grouped.sort_by(|(left, _), (right, _)| left.cmp(right));
    Ok(grouped
        .into_iter()
        .map(|(degrees, coefficient)| {
            let monomial = variables
                .iter()
                .zip(degrees)
                .filter(|(_, degree)| *degree != 0)
                .fold(Atom::num(1), |monomial, (variable, degree)| {
                    monomial * variable.clone().pow(i64::from(degree))
                });
            (monomial, coefficient)
        })
        .collect())
}

pub(crate) fn factors(a: &Atom) -> Vec<(Atom, i64)> {
    let list = if let AtomView::Mul(m) = a.as_view() {
        m.iter().map(|v| v.to_owned()).collect()
    } else {
        vec![a.clone()]
    };
    list.into_iter()
        .map(|v| {
            if let AtomView::Pow(w) = v.as_view() {
                let (base, exponent) = w.get_base_exp();
                if let Ok(power) = exponent.to_string().parse::<i64>() {
                    return (base.to_owned(), power);
                }
            }
            (v, 1)
        })
        .collect()
}

pub(crate) fn powers(monomial: &Atom, variables: &[Atom]) -> Result<Vec<i16>> {
    let mut out = vec![0_i16; variables.len()];
    for (base, power) in factors(monomial) {
        if base.is_one() {
            continue;
        }
        let j = variables
            .iter()
            .position(|v| v == &base)
            .ok_or_else(|| Error::Unsupported("nonpolynomial boundary numerator".into()))?;
        out[j] =
            i16::try_from(power).map_err(|_| Error::Limit("numerator power overflow".into()))?;
    }
    Ok(out)
}

