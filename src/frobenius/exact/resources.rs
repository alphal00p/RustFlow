//! Metadata admission before native rational operations. This performs no
//! coefficient arithmetic, factorization, polynomial products or root finding.
use super::*;

fn degree_bounds(
    view: AtomView<'_>,
    variable: Symbol,
    context: &RunContext,
) -> Result<(usize, usize)> {
    crate::stack_resources::run(|| degree_bounds_inner(view, variable, context))
}

fn degree_bounds_inner(
    view: AtomView<'_>,
    variable: Symbol,
    context: &RunContext,
) -> Result<(usize, usize)> {
    context.cancellation.check()?;
    Ok(match view {
        AtomView::Num(_) => (0, 0),
        AtomView::Var(v) => (usize::from(v.get_symbol() == variable), 0),
        AtomView::Add(sum) => sum.iter().try_fold((0_usize, 0_usize), |(n, d), term| {
            let (a, b) = degree_bounds(term, variable, context)?;
            Ok::<_, Error>((
                scalar_count(n.checked_add(b))?.max(scalar_count(a.checked_add(d))?),
                scalar_count(d.checked_add(b))?,
            ))
        })?,
        AtomView::Mul(product) => product
            .iter()
            .try_fold((0_usize, 0_usize), |(n, d), term| {
                let (a, b) = degree_bounds(term, variable, context)?;
                Ok::<_, Error>((
                    scalar_count(n.checked_add(a))?,
                    scalar_count(d.checked_add(b))?,
                ))
            })?,
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let AtomView::Num(number) = exponent else {
                return Ok((usize::MAX, usize::MAX));
            };
            let symbolica::coefficient::Coefficient::Complex(value) =
                number.get_coeff_view().to_owned()
            else {
                return Ok((usize::MAX, usize::MAX));
            };
            let magnitude = value
                .re
                .clone()
                .abs()
                .to_string()
                .parse::<usize>()
                .map_err(|_| Error::Limit("exact rational degree exponent exceeds usize".into()))?;
            let (n, d) = degree_bounds(base, variable, context)?;
            let (n, d) = (
                scalar_count(n.checked_mul(magnitude))?,
                scalar_count(d.checked_mul(magnitude))?,
            );
            if value.re < Rational::zero() {
                (d, n)
            } else {
                (n, d)
            }
        }
        _ => (usize::MAX, usize::MAX),
    })
}

/// `factor` covers the known products/derivatives/common-denominator additions
/// in the upcoming operation. `copies` is its matrix-entry work multiplier.
/// Estimates deliberately ignore cancellations. Explicit caps are diagnostic;
/// shared RAM/time checks do not interrupt a native operation mid-call.
pub(crate) fn preflight_operation(
    expressions: &[&Atom],
    variable: Symbol,
    factor: usize,
    copies: usize,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<()> {
    limits.validate()?;
    if factor == 0 || copies == 0 {
        return Err(Error::InvalidInput(
            "exact operation multipliers must be positive".into(),
        ));
    }
    process_resources::checkpoint()?;
    let domain = ExactDomain::dense(limits, context);
    let mut degree = 0_usize;
    let mut height = 0_u64;
    let mut bytes = 0_usize;
    for expression in expressions {
        context.cancellation.check()?;
        // Validate exact syntax before the metadata-only degree traversal.
        domain.preflight_series(expression, variable, 0)?;
        let mut nodes = 0;
        height = coefficient_bits(
            height
                .checked_add(domain.series_bound(expression.as_view(), variable, 0, &mut nodes)?)
                .and_then(|v| v.checked_add(2)),
        )?;
        let (n, d) = degree_bounds(expression.as_view(), variable, context)?;
        degree = scalar_count(degree.checked_add(n).and_then(|v| v.checked_add(d)))?;
        bytes = scalar_count(bytes.checked_add(expression.as_view().get_byte_size()))?;
    }
    let predicted_degree = scalar_count(degree.checked_mul(factor))?;
    let coefficient_bits = coefficient_bits(
        height
            .checked_add(2)
            .and_then(|v| v.checked_mul(degree.checked_add(1)?.try_into().ok()?))
            .and_then(|v| v.checked_mul(factor.try_into().ok()?)),
    )?;
    let cells = scalar_count(
        predicted_degree
            .checked_add(1)
            .and_then(|v| v.checked_mul(copies)),
    )?;
    if limits
        .max_order
        .is_some_and(|limit| predicted_degree > limit)
        || limits
            .max_coefficient_bits
            .is_some_and(|limit| coefficient_bits > limit)
        || limits.max_scalar_cells.is_some_and(|limit| cells > limit)
    {
        return Err(Error::Limit(
            "exact root rational-operation degree/height/cell estimate exceeds diagnostic cap"
                .into(),
        ));
    }
    domain.check_allocation(cells)?;
    let expression_bytes = allocation_bytes(bytes, scalar_count(factor.checked_mul(copies))?)?;
    let scalar_bytes = allocation_bytes(cells, std::mem::size_of::<Gaussian>())?;
    let additional_bytes = expression_bytes
        .checked_add(scalar_bytes)
        .ok_or_else(|| Error::Limit("exact rational-operation allocation size overflow".into()))?;
    process_resources::check_allocation(additional_bytes)?;
    Ok(())
}
