//! Analytic single-mass vacuum terminal integrals in the AMFlow normalization.
use crate::{ComplexFloat as C, Error, Precision, Result};
use symbolica::prelude::*;

/// Integral d^D l/(i pi^(D/2)) / (l^2-m^2+i0)^power.
pub fn tadpole(power: u32, mass_squared: &C, dimension: &Rational, p: Precision) -> Result<C> {
    if power == 0 {
        return Ok(p.zero());
    }
    if *mass_squared == p.zero() {
        return Ok(p.zero());
    }
    let half = dimension / &Rational::from(2);
    let argument = Rational::from(power as i64) - &half;
    let ratio = p.div(
        &p.gamma_real(&p.rational(&argument).re)?,
        &p.gamma_real(&p.real(power as i64))?,
    );
    let value = p.mul(&ratio, &p.pow(mass_squared, &p.rational(&(-argument))));
    Ok(if power.is_multiple_of(2) {
        value
    } else {
        p.neg(&value)
    })
}

/// Two-loop single-mass sunset: denominators (l1^2-m^2), l2^2,
/// (l1+l2)^2, with the corresponding positive powers [a,b,c].
/// For unit powers this is Vacuum[2,3] in AMFlow 2.0.
pub fn single_mass_sunset(
    powers: [u32; 3],
    mass_squared: &C,
    dimension: &Rational,
    p: Precision,
) -> Result<C> {
    if powers.contains(&0) || *mass_squared == p.zero() {
        return Ok(p.zero());
    }
    if powers.iter().any(|&v| v > i32::MAX as u32) {
        return Err(Error::Limit("vacuum index too large".into()));
    }
    let [a, b, c] = powers.map(|v| Rational::from(v as i64));
    let half = dimension / &Rational::from(2);
    let gamma = |v: &Rational| p.gamma_real(&p.rational(v).re);
    let numerator = p.mul(
        &p.mul(&gamma(&(&b + &c - &half))?, &gamma(&(&half - &b))?),
        &p.mul(&gamma(&(&half - &c))?, &gamma(&(&a + &b + &c - dimension))?),
    );
    let denominator = p.mul(
        &p.mul(&gamma(&a)?, &gamma(&b)?),
        &p.mul(&gamma(&c)?, &gamma(&half)?),
    );
    let value = p.mul(
        &p.div(&numerator, &denominator),
        &p.pow(mass_squared, &p.rational(&(dimension - &a - &b - &c))),
    );
    Ok(if powers.iter().map(|&v| v as u64).sum::<u64>() % 2 == 0 {
        value
    } else {
        p.neg(&value)
    })
}
