//! Comparison-only independently regenerated Euclidean anchors.
//! Never import this module into a production boundary provider or cache loader.
#![allow(dead_code)]
use serde::Deserialize;
use std::collections::BTreeMap;
use symbolica_amflow::gg_hg::{HiggsJetIntegralSystem, PluginFamilyKind};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::transport_cache::{RootGerm, RootSheet};
use symbolica_amflow::{ComplexFloat, Error, KinematicPoint, Precision, Result};

#[derive(Deserialize)]
struct DecimalComplex {
    real: String,
    imaginary: String,
}
#[derive(Deserialize)]
struct Record {
    family: String,
    basis_order: Vec<usize>,
    #[serde(alias = "variables")]
    coordinates: BTreeMap<String, String>,
    #[serde(alias = "coefficients_by_epsilon")]
    coefficients: Vec<Vec<DecimalComplex>>,
    #[serde(alias = "source_absolute_errors_by_epsilon")]
    absolute_errors: Vec<Vec<String>>,
    #[serde(alias = "evidence_digits")]
    verified_digits: u32,
    root_germs: BTreeMap<String, String>,
    normalization: String,
    report_path: String,
    report_sha256: String,
    #[serde(alias = "evidence")]
    provenance: String,
    epsilon_range: Option<[i32; 2]>,
    leading_epsilon_power: Option<i32>,
    last_epsilon_power: Option<i32>,
}

pub struct EuclideanAnchor {
    pub point: KinematicPoint,
    pub root_germ: RootGerm,
    /// coefficients[epsilon_power][canonical_component], powers zero through four.
    pub coefficients: Vec<Vec<ComplexFloat>>,
    pub absolute_errors: Vec<Vec<Float>>,
    pub verified_digits: u32,
    pub provenance: String,
}

pub struct TruncatedValue {
    pub values: Vec<ComplexFloat>,
    /// Propagated source-coefficient errors ONLY. This excludes the unknown
    /// Taylor remainder and must never be used as a full finite-epsilon bound.
    pub coefficient_errors: Vec<Float>,
    pub first_unknown_power: i32,
}

impl EuclideanAnchor {
    pub fn load(input: &HiggsJetIntegralSystem, p: Precision) -> Result<Self> {
        let text = match input.basis_map().kind {
            PluginFamilyKind::Planar => {
                include_str!("../../fixtures/mg5/planar-regenerated-anchor.json")
            }
            PluginFamilyKind::Nonplanar => {
                include_str!("../../fixtures/mg5/nonplanar-regenerated-anchor.json")
            }
        };
        let data: Record = serde_json::from_str(text)
            .map_err(|e| Error::InvalidInput(format!("comparison-only ggHg anchor: {e}")))?;
        let n = input.basis_map().kind.dimension();
        let range = data
            .epsilon_range
            .or_else(|| Some([data.leading_epsilon_power?, data.last_epsilon_power?]));
        if data.family != input.basis_map().kind.id()
            || data.basis_order != (1..=n).collect::<Vec<_>>()
            || range != Some([0, 4])
            || data.coefficients.len() != 5
            || data.absolute_errors.len() != 5
            || data.coefficients.iter().any(|r| r.len() != n)
            || data.absolute_errors.iter().any(|r| r.len() != n)
            || data.verified_digits != 40
            || data.root_germs.len() != input.basis_map().roots.len()
            || data.coordinates.len() != 3
            || !data.normalization.contains("no additional gamma")
        {
            return Err(Error::InvalidInput(
                "incompatible Euclidean anchor convention".into(),
            ));
        }
        let coordinates = ["s", "t", "b"]
            .into_iter()
            .zip(input.basis_map().coordinates)
            .map(|(name, symbol)| {
                let text = data
                    .coordinates
                    .get(name)
                    .ok_or_else(|| Error::InvalidInput("missing anchor coordinate".into()))?;
                let atom = Atom::parse(text, "gg_hg_anchor_constant", Default::default())
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                Ok((Atom::var(symbol), atom))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let sheets = input
            .basis_map()
            .roots
            .iter()
            .enumerate()
            .map(|(i, root)| {
                if data
                    .root_germs
                    .get(&format!("root{}", i + 1))
                    .map(String::as_str)
                    != Some("principal")
                {
                    return Err(Error::InvalidInput("unexpected anchor root germ".into()));
                }
                Ok((root.symbol, RootSheet::Principal))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let coefficients = data
            .coefficients
            .iter()
            .map(|row| row.iter().map(|v| p.parse(&v.real, &v.imaginary)).collect())
            .collect::<Result<Vec<Vec<_>>>>()?;
        let absolute_errors = data
            .absolute_errors
            .iter()
            .map(|row| row.iter().map(|v| Ok(p.parse(v, "0")?.re)).collect())
            .collect::<Result<Vec<Vec<_>>>>()?;
        if coefficients.iter().flatten().any(|v| !p.finite(v))
            || absolute_errors
                .iter()
                .flatten()
                .any(|v| !v.is_finite() || v < &p.real(0))
        {
            return Err(Error::InvalidInput(
                "nonfinite anchor value or negative allowance".into(),
            ));
        }
        Ok(Self {
            point: KinematicPoint(coordinates),
            root_germ: RootGerm { sheets },
            coefficients,
            absolute_errors,
            verified_digits: data.verified_digits,
            provenance: format!(
                "{}; {} sha256={}; {}",
                data.family, data.report_path, data.report_sha256, data.provenance
            ),
        })
    }

    /// Evaluate only the known Taylor polynomial. No bound for O(epsilon^5)
    /// follows from the stored coefficient errors, even at a small epsilon.
    pub fn truncated(&self, epsilon: &Rational, p: Precision) -> TruncatedValue {
        let e = p.rational(epsilon);
        let mut values = vec![p.zero(); self.coefficients[0].len()];
        let mut coefficient_errors = vec![p.real(0); values.len()];
        for (power, row) in self.coefficients.iter().enumerate() {
            let weight = p.powi(&e, power as i64);
            let absolute_weight = p.norm(&weight);
            for (i, value) in row.iter().enumerate() {
                values[i] = p.add(&values[i], &p.mul(&weight, value));
                coefficient_errors[i] +=
                    absolute_weight.clone() * self.absolute_errors[power][i].clone();
            }
        }
        TruncatedValue {
            values,
            coefficient_errors,
            first_unknown_power: 5,
        }
    }
}

/// Independent all-epsilon formula for the first planar canonical component.
/// The physical integral factorizes into a massless bubble and a squared massive
/// tadpole. This is a comparison formula, never a production boundary provider.
pub fn first_planar_canonical(
    epsilon: &Rational,
    b: &Rational,
    p: Precision,
) -> Result<ComplexFloat> {
    if b >= &Rational::from(0)
        || epsilon <= &Rational::from(0)
        || epsilon >= &Rational::from((1, 2))
    {
        return Err(Error::InvalidInput(
            "analytic planar check needs b<0 and 0<epsilon<1/2".into(),
        ));
    }
    let one = Rational::from(1);
    let plus = p.gamma_real(&p.rational(&(&one + epsilon)).re)?;
    let minus = p.gamma_real(&p.rational(&(&one - epsilon)).re)?;
    let denominator = p.gamma_real(&p.rational(&(&one - &(epsilon * &Rational::from(2)))).re)?;
    let euler = ComplexFloat::new(
        Float::with_val(p.bits, rug::float::Constant::Euler),
        p.real(0),
    );
    let log_scale = p.log(&p.rational(&(-b.clone())));
    let exponential = p.exp(&p.mul(
        &p.rational(epsilon),
        &p.sub(&p.scale(&euler, 2, 1), &log_scale),
    ));
    Ok(p.mul(
        &exponential,
        &p.div(&p.mul(&p.powi(&plus, 2), &p.powi(&minus, 2)), &denominator),
    ))
}
