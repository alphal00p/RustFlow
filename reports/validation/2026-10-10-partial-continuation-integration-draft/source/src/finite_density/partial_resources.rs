//! Metadata-only preflight shared by isolated partial-continuation owners.
//! Bounds inputs and conservative rational polynomial envelopes. These are not
//! allocator quotas and cannot interrupt one native CAS operation mid-call.
use crate::{Error, Result, RunContext};
use serde::Serialize;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct ExpressionLimits {
    pub bytes: usize,
    pub nodes: usize,
    pub depth: usize,
    pub exponent: u32,
    pub expanded_terms: usize,
}
impl Default for ExpressionLimits {
    fn default() -> Self {
        Self {
            bytes: 2 * 1024 * 1024,
            nodes: 65536,
            depth: 64,
            exponent: 1024,
            expanded_terms: 200000,
        }
    }
}
fn limit() -> Error {
    Error::Limit("partial continuation expression preflight".into())
}
#[derive(Clone, Copy)]
struct Shape {
    numerator: usize,
    denominator: usize,
}
fn plus(a: Shape, b: Shape) -> Option<Shape> {
    Some(Shape {
        numerator: a
            .numerator
            .checked_mul(b.denominator)?
            .checked_add(b.numerator.checked_mul(a.denominator)?)?,
        denominator: a.denominator.checked_mul(b.denominator)?,
    })
}
fn times(a: Shape, b: Shape) -> Option<Shape> {
    Some(Shape {
        numerator: a.numerator.checked_mul(b.numerator)?,
        denominator: a.denominator.checked_mul(b.denominator)?,
    })
}
fn shape(
    view: AtomView<'_>,
    depth: usize,
    nodes: &mut usize,
    limits: ExpressionLimits,
    run: &RunContext,
) -> Result<Shape> {
    run.cancellation.check()?;
    *nodes = nodes.checked_add(1).ok_or_else(limit)?;
    if *nodes > limits.nodes || depth > limits.depth {
        return Err(limit());
    }
    let one = Shape {
        numerator: 1,
        denominator: 1,
    };
    let value = match view {
        AtomView::Num(_) | AtomView::Var(_) => one,
        // HEPKit scalar-product markers are atomic polynomial variables. Walk
        // arguments only for metadata size; this preflight grants no rational
        // function/kinematic admission, which the exact owners check separately.
        AtomView::Fun(f) => {
            for arg in f.iter() {
                shape(arg, depth + 1, nodes, limits, run)?;
            }
            one
        }
        AtomView::Add(a) => a.iter().try_fold(
            Shape {
                numerator: 0,
                denominator: 1,
            },
            |s, v| plus(s, shape(v, depth + 1, nodes, limits, run)?).ok_or_else(limit),
        )?,
        AtomView::Mul(a) => a.iter().try_fold(one, |s, v| {
            times(s, shape(v, depth + 1, nodes, limits, run)?).ok_or_else(limit)
        })?,
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            let power = Rational::try_from(exponent).map_err(|_| {
                Error::Unsupported(
                    "partial continuation requires integer rational exponents".into(),
                )
            })?;
            if !power.is_integer() {
                return Err(Error::Unsupported(
                    "partial continuation requires integer rational exponents".into(),
                ));
            }
            let count = power
                .clone()
                .abs()
                .to_string()
                .parse::<u32>()
                .map_err(|_| limit())?;
            if count > limits.exponent {
                return Err(limit());
            }
            let s = shape(base, depth + 1, nodes, limits, run)?;
            let s = Shape {
                numerator: s.numerator.checked_pow(count).ok_or_else(limit)?,
                denominator: s.denominator.checked_pow(count).ok_or_else(limit)?,
            };
            if power < 0 {
                Shape {
                    numerator: s.denominator,
                    denominator: s.numerator,
                }
            } else {
                s
            }
        }
    };
    if value.numerator > limits.expanded_terms || value.denominator > limits.expanded_terms {
        return Err(limit());
    }
    Ok(value)
}
pub(crate) fn expression(a: &Atom, limits: ExpressionLimits, run: &RunContext) -> Result<usize> {
    run.cancellation.check()?;
    if limits.bytes == 0
        || limits.nodes == 0
        || limits.depth == 0
        || limits.expanded_terms == 0
        || a.as_view().get_byte_size() > limits.bytes
    {
        return Err(limit());
    }
    let mut nodes = 0;
    shape(a.as_view(), 0, &mut nodes, limits, run)?;
    Ok(nodes)
}
/// Number of possible monomials of total degree <=degree in variables variables.
/// Conservative preflight only; never used as a polynomial identity.
pub(crate) fn dense_envelope(variables: usize, degree: usize, cap: usize) -> Result<usize> {
    let n = variables.checked_add(degree).ok_or_else(limit)?;
    let k = variables.min(degree);
    let mut count = 1usize;
    for j in 1..=k {
        count = count.checked_mul(n - k + j).ok_or_else(limit)? / j;
        if count > cap {
            return Err(limit());
        }
    }
    Ok(count)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn pre_cancelled_and_expansion_growth_decline_before_cas() {
        let run = RunContext::default();
        run.cancellation.cancel();
        assert!(expression(&Atom::one(), Default::default(), &run).is_err());
        let x = Atom::var(symbol!("partial_resource_x"));
        let a = (Atom::one() + x).pow(30);
        assert!(expression(&a, Default::default(), &RunContext::default()).is_err());
        assert_eq!(dense_envelope(4, 3, 100).unwrap(), 35);
        assert!(dense_envelope(4, 3, 34).is_err());
    }
}
