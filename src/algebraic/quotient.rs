use super::{SquareRoot, normalized_polynomial, rational};
use crate::{Error, Result};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};
use symbolica::domains::algebraic::AlgebraicQuotient;
use symbolica::domains::rational_polynomial::RationalPolynomialField;
use symbolica::prelude::*;
pub(super) fn recompose(terms: &BTreeMap<Vec<usize>, Atom>, roots: &[SquareRoot]) -> Atom {
    terms
        .iter()
        .fold(Atom::new(), |sum, (indices, coefficient)| {
            sum + indices.iter().fold(coefficient.clone(), |value, &i| {
                value * Atom::var(roots[i].symbol)
            })
        })
}

/// Formal inversion preserves the named generators and their independent sheets.
/// Every elimination uses Symbolica's fallible quotient-ring inversion.
pub(super) fn invert_denominator(
    denominator: &Atom,
    roots: &[SquareRoot],
    allowed: &BTreeSet<Atom>,
) -> Result<BTreeMap<Vec<usize>, Atom>> {
    let root_atoms = roots
        .iter()
        .map(|r| Atom::var(r.symbol))
        .collect::<Vec<_>>();
    let base = RationalPolynomialField::<IntegerRing, u16>::new(Z);
    let mut value = Atom::one() / denominator;
    for (index, root) in roots.iter().enumerate() {
        let fraction = rational(&value, allowed)?;
        let numerator =
            normalized_polynomial(&fraction.numerator.to_expression(), &root_atoms, roots)?;
        let denominator =
            normalized_polynomial(&fraction.denominator.to_expression(), &root_atoms, roots)?;
        if denominator.is_empty() {
            return Err(Error::InvalidInput(
                "denominator vanishes under root relations".into(),
            ));
        }
        let mut constant = Atom::new();
        let mut linear = Atom::new();
        for (indices, coefficient) in &denominator {
            let term = indices
                .iter()
                .filter(|&&i| i != index)
                .fold(coefficient.clone(), |a, &i| a * &root_atoms[i]);
            if indices.contains(&index) {
                linear += term;
            } else {
                constant += term;
            }
        }
        if linear.is_zero() {
            value = (recompose(&numerator, roots) / recompose(&denominator, roots))
                .together()
                .cancel();
            continue;
        }
        let prototype =
            MultivariatePolynomial::<_, u16>::new(&base, None, Arc::new(vec![root.symbol.into()]));
        let defining = prototype.monomial(base.one(), vec![2])
            - prototype.constant(rational(&root.radicand, allowed)?);
        let quotient = AlgebraicQuotient::new(defining);
        let element = quotient.add(
            &quotient.constant(rational(&constant, allowed)?),
            &quotient.mul(
                &quotient.constant(rational(&linear, allowed)?),
                &quotient.generator(),
            ),
        );
        let inverse = quotient.try_inv(&element).ok_or_else(|| Error::Unsupported(
            "root denominator is a nonunit in the formal quotient; a sheet-specific domain is required".into()
        ))?;
        let inverse_expression = inverse.poly().into_iter().fold(Atom::new(), |sum, term| {
            sum + term.coefficient.to_expression()
                * Atom::var(root.symbol).pow(i64::from(term.exponents[0]))
        });
        value = (recompose(&numerator, roots) * inverse_expression)
            .together()
            .cancel();
    }
    let fraction = rational(&value, allowed)?;
    let denominator =
        normalized_polynomial(&fraction.denominator.to_expression(), &root_atoms, roots)?;
    let Some(coefficient) = denominator
        .get(&Vec::new())
        .filter(|_| denominator.len() == 1)
    else {
        return Err(Error::Unsupported(
            "root denominator could not be eliminated as a formal unit".into(),
        ));
    };
    Ok(
        normalized_polynomial(&fraction.numerator.to_expression(), &root_atoms, roots)?
            .into_iter()
            .map(|(roots, value)| (roots, (value / coefficient).together().cancel()))
            .collect(),
    )
}
