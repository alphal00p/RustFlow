//! HEPKit owns rational boundary partial fractions, completion and rewriting.
//!
//! This adapter changes only scalar-product notation and preserves the original
//! family coordinates and normalization. Inputs outside the native coefficient
//! field are identified before any native decomposition; errors in an admitted
//! native calculation never fall back to a second decomposition algorithm.
use super::*;
use feynkit_graph::{IntegralFamily as NativeFamily, IntegralFamilyError};
use feynkit_kinematics::Kinematics;
use std::collections::BTreeSet;

pub(super) fn validate_coordinates(variables: &[Atom], template: &IntegralFamily) -> Result<()> {
    let loops = template.loops.len();
    let count = loops
        .checked_add(1)
        .and_then(|n| loops.checked_mul(n))
        .map(|n| n / 2)
        .and_then(|n| {
            loops
                .checked_mul(template.external.len())
                .and_then(|e| n.checked_add(e))
        })
        .ok_or_else(|| Error::Limit("boundary scalar-product count overflow".into()))?;
    if loops == 0
        || variables.len() != count
        || variables
            .iter()
            .any(|a| !matches!(a.as_view(), AtomView::Var(_)))
        || variables.iter().collect::<BTreeSet<_>>().len() != variables.len()
    {
        return Err(Error::InvalidInput(
            "boundary coordinates must be distinct symbols spanning the declared loop/external space".into(),
        ));
    }
    Ok(())
}

// Do not use a failed polynomial conversion as a fallback signal: an exponent
// overflow or native resource limit is not a coefficient-field mismatch.
fn rational_coefficient(a: AtomView<'_>) -> bool {
    match a {
        AtomView::Num(_) => Rational::try_from(a).is_ok(),
        AtomView::Var(v) => v.get_symbol() != crate::family::imaginary_parameter(),
        AtomView::Add(v) => v.iter().all(rational_coefficient),
        AtomView::Mul(v) => v.iter().all(rational_coefficient),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            Integer::try_from(exponent).is_ok() && rational_coefficient(base)
        }
        _ => false,
    }
}

pub(super) fn supports(denominators: &[Propagator]) -> bool {
    denominators.iter().all(|d| {
        std::iter::once(&d.constant)
            .chain(&d.scalar_products)
            .all(|a| rational_coefficient(a.as_view()))
    })
}

fn native_error(error: IntegralFamilyError) -> Error {
    match error {
        IntegralFamilyError::PartialFractionLimit(_) | IntegralFamilyError::PowerOverflow => {
            Error::Limit(format!("native boundary partial fractions: {error}"))
        }
        _ => Error::InvalidInput(format!("native boundary family: {error}")),
    }
}

pub(super) fn to_integrals(
    denominators: &[Propagator],
    indices: &[i16],
    numerator: &Atom,
    variables: &[Atom],
    template: &IntegralFamily,
    budget: usize,
) -> Result<Vec<IntegralTerm>> {
    if budget == 0 {
        return Err(Error::Limit(
            "native boundary partial fractions: zero state budget".into(),
        ));
    }
    let sources = boundary_sources(numerator, denominators, variables);
    let loops = (0..template.loops.len())
        .map(|i| {
            fresh_boundary_symbol(
                &format!("symbolica_amflow::boundary_native_loop_{i}"),
                &sources,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let external = (0..template.external.len())
        .map(|i| {
            fresh_boundary_symbol(
                &format!("symbolica_amflow::boundary_native_external_{i}"),
                &sources,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let dimension = fresh_boundary_symbol("symbolica_amflow::boundary_native_dimension", &sources)?;
    let kinematics = Kinematics::in_dimension(&dimension)
        .map_err(|e| Error::InvalidInput(format!("native boundary dimension: {e}")))?
        .with_momenta(loops.iter().chain(&external).cloned())
        .map_err(|e| Error::InvalidInput(format!("native boundary kinematics: {e}")))?;
    // Match exact product identities: native HEPKit interleaves LL/LE products,
    // whereas the boundary input groups all LL products before all LE products.
    let pairs = loops
        .iter()
        .enumerate()
        .flat_map(|(i, left)| loops[i..].iter().map(move |right| (left, right)))
        .chain(
            loops
                .iter()
                .flat_map(|left| external.iter().map(move |right| (left, right))),
        );
    let products = pairs
        .map(|(left, right)| {
            kinematics
                .scalar_product(left, right)
                .map_err(|e| Error::InvalidInput(format!("native boundary scalar product: {e}")))
        })
        .collect::<Result<Vec<_>>>()?;
    let to_native = variables
        .iter()
        .cloned()
        .zip(products.iter().cloned())
        .collect();
    let from_native = products
        .iter()
        .cloned()
        .zip(variables.iter().cloned())
        .collect();
    let expressions = denominators
        .iter()
        .map(|d| {
            d.scalar_products
                .iter()
                .zip(&products)
                .fold(d.constant.clone(), |sum, (coefficient, product)| {
                    sum + coefficient * product
                })
                .expand()
        })
        .collect();
    let native =
        NativeFamily::new(loops, external, expressions, &kinematics).map_err(native_error)?;
    let powers = indices.iter().map(|&n| i32::from(n)).collect::<Vec<_>>();
    let parts = native
        .partial_fraction(&powers, budget)
        .map_err(native_error)?;
    let numerator = substitute(numerator, &to_native);
    let labels = denominator_labels(&numerator, denominators, variables)?;
    let mut out = Vec::new();
    for (coefficient, powers) in parts {
        if products
            .iter()
            .any(|product| coefficient.contains(product.as_view()))
        {
            return Err(Error::Numerical(
                "native boundary partial fractions returned a loop-dependent coefficient".into(),
            ));
        }
        let physical = powers.iter().filter(|&&n| n > 0).count();
        let completed = native
            .sector(&powers)
            .and_then(|f| f.complete(&products))
            .map_err(native_error)?;
        let mut family = template.clone();
        family.propagators = completed
            .denominators()
            .iter()
            .map(|d| affine(&substitute(d, &from_native), variables))
            .collect::<Result<Vec<_>>>()?;
        family.physical_propagators = physical;
        let mut active = powers
            .iter()
            .filter(|&&n| n > 0)
            .map(|&n| {
                i16::try_from(n).map_err(|_| {
                    Error::Limit("native boundary denominator power exceeds i16".into())
                })
            })
            .collect::<Result<Vec<_>>>()?;
        active.resize(variables.len(), 0);
        let rewritten = completed
            .rewrite_numerator(&(&numerator * coefficient), &labels)
            .map_err(native_error)?;
        append_numerator_terms(&rewritten, &labels, &family, &active, &mut out)?;
    }
    Ok(out)
}
