//! Exact algebraic coefficients and analytic normalizations for sampled integrals.
//!
//! Root signs are discrete data. Root magnitudes, gamma factors, and exponential
//! factors are recomputed at every working precision, before Laurent fitting.
use crate::Prescription;
use crate::algebraic::SquareRoot;
use crate::transport_cache::{RootGerm, RootSheet};
#[cfg(feature = "automatic")]
use crate::{ComplexFloat, Error, KinematicPoint, Precision, Result};
#[cfg(feature = "automatic")]
use std::collections::BTreeSet;
use symbolica::prelude::*;

/// Factors analytic and nonzero near epsilon = 0. Rational epsilon powers
/// belong in the projection rows, where their pole orders are retained.
#[derive(Clone, Debug)]
pub enum SampleNormalization {
    /// exp(coefficient * epsilon * EulerGamma).
    EulerGammaExponential(Rational),
    /// Gamma(offset + slope * epsilon)^power. Offset must be positive.
    Gamma {
        offset: Rational,
        slope: Rational,
        power: i32,
    },
    /// exp(slope * epsilon * log(base)). The prescription specifies the side
    /// of a negative real *base*, independently of the integral prescription.
    Power {
        base: Atom,
        slope: Rational,
        prescription: Prescription,
    },
}

/// A common normalization for a batch of projections and their named roots.
/// These factors multiply the sum over families, before any epsilon fit.
#[derive(Clone, Debug, Default)]
pub struct ProjectionFactors {
    pub roots: Vec<SquareRoot>,
    pub sheets: std::collections::BTreeMap<Symbol, RootSheet>,
    pub normalization: Vec<SampleNormalization>,
}

impl ProjectionFactors {
    pub fn with_roots(roots: Vec<SquareRoot>, germ: RootGerm) -> Self {
        Self {
            roots,
            sheets: germ.sheets,
            normalization: Vec::new(),
        }
    }

    #[cfg(feature = "automatic")]
    pub(crate) fn at(&self, point: &KinematicPoint) -> Result<Self> {
        let roots = self
            .roots
            .iter()
            .map(|r| SquareRoot {
                symbol: r.symbol,
                radicand: point.apply(&r.radicand).together().cancel(),
            })
            .collect::<Vec<_>>();
        let symbols = roots.iter().map(|r| r.symbol).collect::<BTreeSet<_>>();
        if symbols.len() != roots.len()
            || symbols != self.sheets.keys().copied().collect::<BTreeSet<_>>()
            || symbols.contains(&crate::family::imaginary_parameter())
        {
            return Err(Error::InvalidInput(
                "projection factors require distinct roots and exactly one sheet per root".into(),
            ));
        }
        for root in &roots {
            validate_exact_constant(&root.radicand)?;
            if root.radicand.is_zero() {
                return Err(Error::InvalidInput(
                    "projection root radicand vanishes at this point".into(),
                ));
            }
        }
        let normalization = self
            .normalization
            .iter()
            .map(|n| match n {
                SampleNormalization::Gamma { offset, .. } if offset <= &Rational::from(0) => {
                    Err(Error::InvalidInput(
                        "gamma normalization needs positive offset at epsilon = 0".into(),
                    ))
                }
                SampleNormalization::Power {
                    base,
                    slope,
                    prescription,
                } => {
                    let base = point.apply(base).together().cancel();
                    validate_exact_constant(&base)?;
                    if base.is_zero() {
                        return Err(Error::InvalidInput(
                            "normalization scale must be nonzero".into(),
                        ));
                    }
                    Ok(SampleNormalization::Power {
                        base,
                        slope: slope.clone(),
                        prescription: *prescription,
                    })
                }
                other => Ok(other.clone()),
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            roots,
            sheets: self.sheets.clone(),
            normalization,
        })
    }

    #[cfg(feature = "automatic")]
    pub(crate) fn parameters(&self, p: Precision) -> Result<ahash::HashMap<Atom, ComplexFloat>> {
        self.roots
            .iter()
            .map(|root| {
                let value = p.eval(&root.radicand, &Default::default())?;
                let root_value = crate::algebraic::principal_sqrt(p, &value)?;
                if root_value == p.zero() {
                    return Err(Error::Numerical(
                        "projection root vanished at working precision".into(),
                    ));
                }
                Ok((
                    Atom::var(root.symbol),
                    match self.sheets[&root.symbol] {
                        RootSheet::Principal => root_value,
                        RootSheet::Opposite => p.neg(&root_value),
                    },
                ))
            })
            .collect()
    }

    #[cfg(feature = "automatic")]
    pub(crate) fn evaluate(&self, epsilon: &Rational, p: Precision) -> Result<ComplexFloat> {
        let mut product = p.i(1);
        for factor in &self.normalization {
            let value = match factor {
                SampleNormalization::EulerGammaExponential(coefficient) => {
                    let euler = ComplexFloat::new(p.real(0).euler(), p.real(0));
                    p.exp(&p.mul(&p.rational(&(coefficient * epsilon)), &euler))
                }
                SampleNormalization::Gamma {
                    offset,
                    slope,
                    power,
                } => {
                    let argument = offset + &(slope * epsilon);
                    p.powi(&p.gamma_real(&p.rational(&argument).re)?, i64::from(*power))
                }
                SampleNormalization::Power {
                    base,
                    slope,
                    prescription,
                } => {
                    let base = p.eval(base, &Default::default())?;
                    let mut logarithm = p.log(&base);
                    if base.im == p.real(0) && base.re < p.real(0) {
                        let pi = p.real(0).pi();
                        logarithm.im = match prescription {
                            Prescription::PlusI0 => pi,
                            Prescription::MinusI0 => -pi,
                        };
                    }
                    p.exp(&p.mul(&p.rational(&(slope * epsilon)), &logarithm))
                }
            };
            product = p.mul(&product, &value);
            if !p.finite(&product) || product == p.zero() {
                return Err(Error::Numerical(
                    "nonfinite or vanishing projection normalization".into(),
                ));
            }
        }
        Ok(product)
    }
}

#[cfg(feature = "automatic")]
fn validate_exact_constant(value: &Atom) -> Result<()> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(value.as_view(), &mut symbols)?;
    if symbols
        .iter()
        .any(|s| *s != Atom::var(crate::family::imaginary_parameter()))
    {
        return Err(Error::InvalidInput(
            "projection factors require exact, fully specialized radicands and scales".into(),
        ));
    }
    let _: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(value)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| {
            Error::Unsupported(format!("projection scale or radicand is not rational: {e}"))
        })?;
    Ok(())
}
