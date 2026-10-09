//! Whole-amplitude convention adapters for independent occupied shells.
//!
//! Native virtual loops use d^D k/(i*pi^(D/2)). Each native occupied loop uses
//! d^D q/pi^(D/2) theta(q0) theta(mu-q0) C_n(q^2-m^2), with
//! C_n(x)=(-1)^(n-1)/(n-1)! delta^(n-1)(x). Occupations are held independent of
//! q's shell coordinate until the distribution is integrated. This measure has
//! no final-state momentum-conservation delta or total timelike channel.
//!
//! `native_measure_to_euclidean` converts measures ONLY. Euclidean quadratic
//! index signs and medium numerator Wick factors belong to target coefficients
//! and must not be applied twice, nor multiplied by another per-cut minus sign.
//! Exact routing Jacobians remain separate geometry factors. These algebraic
//! adapters do not establish common-contour or pinch admission for an amplitude.
use crate::coefficient::{exact_coefficient_list, powers};
use crate::family::substitute;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

fn count(value: usize) -> Result<i64> {
    i64::try_from(value).map_err(|_| Error::Limit("normalization loop count overflow".into()))
}

/// Convert a mixed native loop measure to dp0/(2*pi) d^(D-1)p/(2*pi)^(D-1)
/// for virtual loops and d^(D-1)q/(2*pi)^(D-1) times the occupied shell measure.
/// Every occupied loop contributes an extra 2*pi relative to a virtual loop.
/// This factor contains neither quadratic-index signs nor MSbar scale factors.
pub fn native_measure_to_euclidean(
    loops: usize,
    occupied_loops: usize,
    dimension: &Atom,
) -> Result<Atom> {
    if occupied_loops > loops {
        return Err(Error::InvalidInput(
            "more occupied loop coordinates than integration loops".into(),
        ));
    }
    let l = count(loops)?;
    let k = count(occupied_loops)?;
    let pi = Atom::var(Symbol::PI);
    Ok((Atom::num(2) * &pi).pow(k) * (Atom::num(4) * pi).pow(-Atom::num((l, 2)) * dimension))
}

/// Euclidean rho=-D for a scalar quadratic, equally for an ordinary line and a
/// required shell of positive index. Negative quadratic numerator indices use
/// the same parity. Medium-coordinate slots are excluded from this function.
/// No additional (-1)^number_of_cuts may be applied after this target phase.
pub fn quadratic_index_phase(indices: &[i16]) -> Atom {
    let odd = indices
        .iter()
        .fold(false, |odd, &n| odd ^ (n.rem_euclid(2) == 1));
    Atom::num(if odd { -1 } else { 1 })
}

/// Future-shell convention q_M=(-i*P_E0,-P_Espatial), the global reversal of
/// standard Minkowski Wick coordinates. Thus g_E=-g_M and u_E=+i*u_M; the
/// Euclidean crossed pole P_E0=+iE becomes future q_M0=E. Global loop reversal
/// preserves the virtual +i0 prescription and integration measure.
///
/// Coordinates in the returned polynomial carry the Minkowski interpretation.
/// Chemical offsets already present in the input coefficients are retained.
pub fn wick_polynomial(
    numerator: &Atom,
    scalar_products: &[Atom],
    medium_coordinates: &[Atom],
) -> Result<Atom> {
    let coordinates = scalar_products
        .iter()
        .chain(medium_coordinates)
        .cloned()
        .collect::<Vec<_>>();
    if coordinates.iter().collect::<BTreeSet<_>>().len() != coordinates.len() {
        return Err(Error::InvalidInput("repeated Wick coordinate".into()));
    }
    let symbols = coordinates
        .iter()
        .map(|a| {
            if let AtomView::Var(v) = a.as_view() {
                Ok(v.get_symbol())
            } else {
                Err(Error::InvalidInput(
                    "Wick coordinate is not a symbol".into(),
                ))
            }
        })
        .collect::<Result<Vec<_>>>()?;
    for (monomial, coefficient) in exact_coefficient_list(numerator, &coordinates)? {
        if powers(&monomial, &coordinates)?.iter().any(|&n| n < 0)
            || symbols
                .iter()
                .any(|&s| !coefficient.derivative(s).is_zero())
        {
            return Err(Error::InvalidInput(
                "Wick numerator must be polynomial in momentum coordinates".into(),
            ));
        }
    }
    let map = scalar_products
        .iter()
        .map(|g| (g.clone(), -g))
        .chain(
            medium_coordinates
                .iter()
                .map(|u| (u.clone(), Atom::i() * u)),
        )
        .collect::<BTreeMap<_, _>>();
    Ok(substitute(numerator, &map).expand())
}

/// Per-original-loop MSbar factor, independent of how many lines are occupied.
/// The physical renormalization scale must be positive; powers use that branch.
/// This common factor multiplies the assembled Euclidean amplitude once.
pub fn msbar_measure(loops: usize, epsilon: &Atom, lambda_bar: &Atom) -> Result<Atom> {
    let exponent = Atom::num(count(loops)?) * epsilon;
    let gamma = Atom::var(symbolica::transcendental::euler_gamma());
    let base = function!(Symbol::EXP, gamma) * lambda_bar.clone().pow(2)
        / (Atom::num(4) * Atom::var(Symbol::PI));
    Ok(base.pow(exponent))
}

/// Unexpanded conventional scale factor (4*pi)^(-2L)*(Lambda_bar/2)^(2L*eps).
/// This is a definition of units, with no reference integral values involved.
pub fn scale_factored_normalization(
    loops: usize,
    epsilon: &Atom,
    lambda_bar: &Atom,
) -> Result<Atom> {
    let twice = count(loops)?
        .checked_mul(2)
        .ok_or_else(|| Error::Limit("normalization exponent overflow".into()))?;
    Ok((Atom::num(4) * Atom::var(Symbol::PI)).pow(-twice)
        * (lambda_bar / Atom::num(2)).pow(Atom::num(twice) * epsilon))
}

/// Native mixed measure including MSbar and divided by the unexpanded scale
/// factor, at D=4-2*eps. Target phases/Wick factors remain separate.
pub fn native_measure_to_scale_factored_msbar(
    loops: usize,
    occupied_loops: usize,
    epsilon: &Atom,
) -> Result<Atom> {
    if occupied_loops > loops {
        return Err(Error::InvalidInput(
            "more occupied loop coordinates than integration loops".into(),
        ));
    }
    let exponent = Atom::num(count(loops)?) * epsilon;
    let gamma = Atom::var(symbolica::transcendental::euler_gamma());
    Ok(
        (Atom::num(2) * Atom::var(Symbol::PI)).pow(count(occupied_loops)?)
            * Atom::num(4).pow(exponent.clone())
            * function!(Symbol::EXP, gamma * exponent),
    )
}
