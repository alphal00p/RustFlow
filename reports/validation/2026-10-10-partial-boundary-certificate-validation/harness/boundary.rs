use crate::{Error,Result};
use crate::coefficient::{exact_coefficient_list,powers};
use symbolica::prelude::*;
pub use symbolica_amflow::finite_density::boundary::{OccupiedBoundaryDistribution,OccupiedBoundaryLimits};
fn validate_distribution(
    d: &OccupiedBoundaryDistribution,
    limits: OccupiedBoundaryLimits,
) -> Result<()> {
    if d.upper_index < 0
        || d.lower_index < 0
        || d.shell.mass_squared < 0
        || d.shell.chemical_potential < 0
    {
        return Err(Error::InvalidInput(
            "invalid occupied boundary distribution or support".into(),
        ));
    }
    if usize::try_from(d.cut_index.max(0)).unwrap()
        + d.upper_index as usize
        + d.lower_index as usize
        > limits.max_distribution_order
    {
        return Err(Error::Limit(
            "occupied boundary distribution order limit".into(),
        ));
    }
    Ok(())
}

fn polynomial_terms(
    expression: &Atom,
    variables: &[Atom],
    budget: usize,
) -> Result<Vec<(Vec<i16>, Atom)>> {
    laurent_terms(expression, variables, &vec![false; variables.len()], budget)
}

fn laurent_terms(
    expression: &Atom,
    variables: &[Atom],
    inverse_allowed: &[bool],
    budget: usize,
) -> Result<Vec<(Vec<i16>, Atom)>> {
    if variables.len() != inverse_allowed.len() {
        return Err(Error::InvalidInput(
            "compact Laurent coordinate dimensions".into(),
        ));
    }
    let terms = exact_coefficient_list(expression, variables)?;
    if terms.len() > budget {
        return Err(Error::Limit(
            "occupied polynomial term budget exhausted".into(),
        ));
    }
    terms
        .into_iter()
        .map(|(monomial, coefficient)| {
            let degrees = powers(&monomial, variables)?;
            if degrees
                .iter()
                .zip(inverse_allowed)
                .any(|(&degree, &allowed)| degree < 0 && !allowed)
                || variables.iter().any(|variable| {
                    let AtomView::Var(v) = variable.as_view() else {
                        return true;
                    };
                    !coefficient.derivative(v.get_symbol()).is_zero()
                })
            {
                return Err(Error::Unsupported(
                    "nonpolynomial occupied soft factor requires recursive weighted reduction"
                        .into(),
                ));
            }
            Ok((degrees, coefficient))
        })
        .collect()
}
#[path="virtual_soft.rs"] mod virtual_soft;
