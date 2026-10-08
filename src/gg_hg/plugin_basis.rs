//! Exact physical-integral maps for the 48/61 canonical systems of arXiv:2112.07578.
use super::{CanonicalTerm, SourceEvidence};
use crate::algebraic::SquareRoot;
use crate::family::LinearCombination;
use crate::transport_cache::RootGerm;
use crate::{
    Error, Integral, IntegralFamily, ProjectionFactors, Propagator, Result, SampleNormalization,
};
use serde::Deserialize;
use std::collections::BTreeSet;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PluginFamilyKind {
    Planar,
    Nonplanar,
}
impl PluginFamilyKind {
    pub fn id(self) -> &'static str {
        match self {
            Self::Planar => "Planar_EW1",
            Self::Nonplanar => "NP_EW1",
        }
    }
    pub fn dimension(self) -> usize {
        match self {
            Self::Planar => 48,
            Self::Nonplanar => 61,
        }
    }
}

/// Proofs checked during mathematical-data extraction, with both source hashes.
/// Native exact inverse checking is independently available through
/// `PluginBasisMap::verify_inverse` and is exercised by the module tests.
#[derive(Clone, Debug, Deserialize)]
pub struct PluginMapEvidence {
    pub source_sha256: String,
    pub published_dlog_sha256: String,
    pub plugin_dlog_sha256: String,
    pub linear_reconstruction: bool,
    pub exact_root_reconstruction: bool,
    pub right_inverse: bool,
    pub exact_dlog_equality: bool,
}

/// Unit-boson-mass integrals and the exact change to the plugin canonical basis.
///
/// With ordinary RustFlow integrals I, the canonical vector is
/// F = exp(2 eps EulerGamma) * publication_to_canonical * I.
/// Root signs are separate caller-supplied germ data. Coordinates are
/// s=(p1+p2)^2, t=(p1-p3)^2 and b=(p1+p2-p3)^2 with massless p1,p2,p3.
#[derive(Clone, Debug)]
pub struct PluginBasisMap {
    pub kind: PluginFamilyKind,
    pub coordinates: [Symbol; 3],
    pub epsilon: Symbol,
    pub integrals: Vec<Integral>,
    pub roots: Vec<SquareRoot>,
    pub publication_to_canonical: Vec<Vec<CanonicalTerm>>,
    pub canonical_to_publication: Vec<Vec<CanonicalTerm>>,
    pub paper: SourceEvidence,
    pub evidence: PluginMapEvidence,
    routing: Vec<Routing>,
}

#[derive(Deserialize)]
struct Document {
    schema: String,
    paper: SourceEvidence,
    dimension: String,
    physical_propagators: usize,
    measure: String,
    families: Vec<Record>,
}
#[derive(Deserialize)]
struct Record {
    family: String,
    dimension: usize,
    physical_integrals: Vec<Vec<i16>>,
    roots: Vec<RootRecord>,
    publication_to_canonical: Vec<Vec<TermRecord>>,
    canonical_to_publication: Vec<Vec<TermRecord>>,
    evidence: PluginMapEvidence,
    routing: Vec<Routing>,
}
#[derive(Clone, Debug, Deserialize)]
struct Routing {
    #[serde(rename = "loop")]
    loops: [i64; 2],
    external: [i64; 3],
    mass_squared: i64,
}
#[derive(Deserialize)]
struct RootRecord {
    name: String,
    radicand: String,
}
#[derive(Deserialize)]
struct TermRecord {
    column: usize,
    coefficient: String,
}

impl PluginBasisMap {
    /// Include the extracted mathematics as well as its upstream source in
    /// boundary identities: correcting an exporter must invalidate old seeds.
    pub(crate) fn input_fingerprint() -> String {
        super::data::digest("plugin-physical-map.json")
            .expect("known fixture")
            .to_owned()
    }

    /// Parse exact mathematical data in the caller's Symbolica namespace.
    /// Names s,t,b,eps and root1..root8 therefore match a system loaded there.
    pub fn load(kind: PluginFamilyKind, namespace: &str) -> Result<Self> {
        let data: Document =
            serde_json::from_str(&super::data::get("plugin-physical-map.json")?)
                .map_err(|e| Error::InvalidInput(format!("plugin physical map: {e}")))?;
        if data.schema != "rustflow-gg-hg-plugin-physical-map-v2"
            || data.dimension != "4-2*eps"
            || data.physical_propagators != 7
            || data.measure != "exp(2*eps*EulerGamma); mV2=mu2=1"
        {
            return Err(Error::InvalidInput(
                "incompatible plugin physical map".into(),
            ));
        }
        let record = data
            .families
            .into_iter()
            .find(|r| r.family == kind.id())
            .ok_or_else(|| Error::InvalidInput("missing plugin physical map".into()))?;
        let n = kind.dimension();
        if record.dimension != n
            || record.physical_integrals.len() != n
            || record.routing.len() != 9
            || !record.evidence.linear_reconstruction
            || !record.evidence.exact_root_reconstruction
            || !record.evidence.right_inverse
            || !record.evidence.exact_dlog_equality
        {
            return Err(Error::InvalidInput(
                "uncertified plugin physical map".into(),
            ));
        }
        let parse = |text: &str| {
            Atom::parse(text, namespace.to_owned(), Default::default())
                .map_err(|e| Error::InvalidInput(e.to_string()))
        };
        let symbol = |text: &str| match parse(text)?.as_view() {
            AtomView::Var(v) => Ok(v.get_symbol()),
            _ => Err(Error::InvalidInput(
                "plugin scalar name is not a symbol".into(),
            )),
        };
        let coordinates = [symbol("s")?, symbol("t")?, symbol("b")?];
        let epsilon = symbol("eps")?;
        let mut allowed = coordinates
            .iter()
            .map(|&s| Atom::var(s))
            .collect::<BTreeSet<_>>();
        let roots = record
            .roots
            .into_iter()
            .map(|r| {
                let radicand = parse(&r.radicand)?;
                validate_scalar(&radicand, &allowed)?;
                Ok(SquareRoot {
                    symbol: symbol(&r.name)?,
                    radicand,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        allowed.extend(roots.iter().map(|r| Atom::var(r.symbol)));
        allowed.insert(Atom::var(epsilon));
        let rows = |rows: Vec<Vec<TermRecord>>| -> Result<Vec<Vec<CanonicalTerm>>> {
            if rows.len() != n {
                return Err(Error::InvalidInput("invalid plugin matrix size".into()));
            }
            rows.into_iter()
                .map(|row| {
                    let mut columns = BTreeSet::new();
                    row.into_iter()
                        .map(|t| {
                            if t.column >= n || !columns.insert(t.column) {
                                return Err(Error::InvalidInput(
                                    "invalid plugin sparse column".into(),
                                ));
                            }
                            let coefficient = parse(&t.coefficient)?;
                            validate_scalar(&coefficient, &allowed)?;
                            Ok(CanonicalTerm {
                                integral: t.column,
                                coefficient,
                            })
                        })
                        .collect()
                })
                .collect()
        };
        let integrals = record
            .physical_integrals
            .into_iter()
            .map(|p| {
                if p.len() != 9 || p[7..].iter().any(|&v| v > 0) {
                    return Err(Error::InvalidInput("invalid plugin physical powers".into()));
                }
                Ok(Integral(p))
            })
            .collect::<Result<Vec<_>>>()?;
        if integrals.iter().collect::<BTreeSet<_>>().len() != n {
            return Err(Error::InvalidInput(
                "duplicate plugin physical integral".into(),
            ));
        }
        Ok(Self {
            kind,
            coordinates,
            epsilon,
            integrals,
            roots,
            publication_to_canonical: rows(record.publication_to_canonical)?,
            canonical_to_publication: rows(record.canonical_to_publication)?,
            paper: data.paper,
            evidence: record.evidence,
            routing: record.routing,
        })
    }

    /// Native Symbolica matrix product and exact polynomial cancellation.
    /// No sampling, floating-point rank test, or root-sign assumption is used.
    pub fn verify_inverse(&self) -> Result<()> {
        let product = crate::algebra::matmul(
            &self.dense_matrix(&self.canonical_to_publication),
            &self.dense_matrix(&self.publication_to_canonical),
        );
        if product.iter().enumerate().any(|(i, row)| {
            row.iter()
                .enumerate()
                .any(|(j, v)| !(v - Atom::num(i64::from(i == j))).is_zero())
        }) {
            return Err(Error::InvalidInput(
                "plugin basis inverse identity failed".into(),
            ));
        }
        Ok(())
    }

    pub(crate) fn dense_matrix(&self, rows: &[Vec<CanonicalTerm>]) -> Vec<Vec<Atom>> {
        let n = self.kind.dimension();
        let mut matrix = vec![vec![Atom::new(); n]; n];
        for (i, row) in rows.iter().enumerate() {
            for term in row {
                matrix[i][term.integral] = term.coefficient.clone();
            }
        }
        matrix
    }

    /// Ordinary q²-m² denominators, including the signed source numerator powers.
    /// The paper uses the same denominator sign, so there is no (-1)^sum(a) factor.
    pub fn family(&self) -> Result<IntegralFamily> {
        let [s, t, b] = self.coordinates.map(Atom::var);
        let zero = Atom::new();
        let hs = &s / Atom::num(2);
        let ht = -&t / Atom::num(2);
        let hu = (&s + &t - &b) / Atom::num(2);
        let gram = vec![
            vec![zero.clone(), hs.clone(), ht.clone()],
            vec![hs, zero.clone(), hu.clone()],
            vec![ht, hu, zero],
        ];
        let propagators = self
            .routing
            .iter()
            .map(|r| Propagator::quadratic(&r.loops, &r.external, Atom::num(r.mass_squared), &gram))
            .collect::<Result<Vec<_>>>()?;
        let family = IntegralFamily {
            name: format!("gg_hg_2112_{}", self.kind.id()),
            loops: vec!["k1".into(), "k2".into()],
            external: vec!["p1".into(), "p2".into(), "p3".into()],
            external_gram: gram,
            propagators,
            physical_propagators: 7,
            epsilon: self.epsilon,
            dimension: 4,
        };
        family.validate()?;
        for integral in &self.integrals {
            family.validate_integral(integral)?;
        }
        Ok(family)
    }

    /// One exact projection per canonical master; shared targets are retained.
    pub fn canonical_rows(&self) -> Vec<LinearCombination> {
        self.publication_to_canonical
            .iter()
            .map(|row| {
                row.iter()
                    .map(|t| (self.integrals[t.integral].clone(), t.coefficient.clone()))
                    .collect()
            })
            .collect()
    }

    /// Publication measure at unit boson mass and mu²=1. Evaluated at each exact
    /// epsilon sample and precision before the native Laurent reconstruction.
    pub fn factors(&self, germ: RootGerm) -> ProjectionFactors {
        let mut factors = ProjectionFactors::with_roots(self.roots.clone(), germ);
        factors.normalization = vec![SampleNormalization::EulerGammaExponential(Rational::from(
            2,
        ))];
        factors
    }
}

fn validate_scalar(a: &Atom, allowed: &BTreeSet<Atom>) -> Result<()> {
    let mut found = BTreeSet::new();
    crate::family::scalar_symbols(a.as_view(), &mut found)?;
    if !found.is_subset(allowed) {
        return Err(Error::InvalidInput("undeclared plugin map scalar".into()));
    }
    let _: RationalPolynomial<IntegerRing, u16> = a
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::InvalidInput(format!("nonrational plugin weight: {e}")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "automatic")]
    use crate::KinematicPoint;
    #[cfg(feature = "automatic")]
    use crate::transport_cache::RootSheet;
    #[test]
    fn published_plugin_inverse_is_exact_in_native_symbolica() -> Result<()> {
        for (kind, entries) in [
            (PluginFamilyKind::Planar, 201),
            (PluginFamilyKind::Nonplanar, 453),
        ] {
            let map = PluginBasisMap::load(kind, "rustflow_plugin_map_test")?;
            assert_eq!(
                map.publication_to_canonical
                    .iter()
                    .map(Vec::len)
                    .sum::<usize>(),
                entries
            );
            map.verify_inverse()?;
        }
        Ok(())
    }
    #[test]
    #[cfg(feature = "automatic")]
    fn physical_plugin_routings_and_epsilon_normalization_are_explicit() -> Result<()> {
        for kind in [PluginFamilyKind::Planar, PluginFamilyKind::Nonplanar] {
            let map = PluginBasisMap::load(kind, "rustflow_plugin_route_test")?;
            let f = map.family()?;
            if kind == PluginFamilyKind::Planar {
                // J_1 is the product of a massless bubble (powers 2,1) and
                // a massive tadpole (power 2). Its canonical prefactor fixes
                // the denominator sign independently of matrix inversion.
                let row = &map.publication_to_canonical[0];
                assert_eq!(row.len(), 1);
                assert_eq!(row[0].integral, 0);
                assert_eq!(map.integrals[0], Integral(vec![0, 0, 0, 2, 1, 0, 2, 0, 0]));
                let expected = -Atom::var(map.coordinates[2]) * Atom::var(map.epsilon).pow(2);
                assert!((&row[0].coefficient - expected).expand().is_zero());
                // The primary J_2 rule contains two DISTINCT roots. An inverse
                // identity alone cannot detect a shared erroneous root rename.
                let row = &map.canonical_to_publication[1];
                assert_eq!(row.len(), 1);
                assert_eq!(row[0].integral, 10);
                assert_eq!(map.integrals[1], Integral(vec![0, 0, 0, 2, 1, 1, 2, 0, 0]));
                let b = Atom::var(map.coordinates[2]);
                let expected = -Atom::var(map.roots[0].symbol) * Atom::var(map.roots[1].symbol)
                    / (Atom::var(map.epsilon).pow(2) * b.clone().pow(2) * (Atom::num(4) - b));
                assert!(
                    (&row[0].coefficient - expected)
                        .together()
                        .cancel()
                        .is_zero()
                );
            }
            for i in 0..9 {
                assert_eq!(f.mass_squared(i)?, Atom::num(i64::from(i == 5 || i == 6)));
            }
            assert_eq!(f.propagators[3].constant, Atom::var(map.coordinates[0]));
            assert!(map.integrals.iter().any(|i| i.0.iter().any(|&p| p < 0)));
            let point = KinematicPoint(
                map.coordinates
                    .into_iter()
                    .map(Atom::var)
                    .zip([Atom::num(-2), Atom::num(-3), Atom::num(-7)])
                    .collect(),
            );
            let germ = RootGerm {
                sheets: map
                    .roots
                    .iter()
                    .map(|r| (r.symbol, RootSheet::Principal))
                    .collect(),
            };
            let factors = map.factors(germ).at(&point)?;
            assert_eq!(
                factors.roots.len(),
                if kind == PluginFamilyKind::Planar {
                    2
                } else {
                    8
                }
            );
            assert!(
                matches!(&factors.normalization[..],[SampleNormalization::EulerGammaExponential(v)] if v==&Rational::from(2))
            );
        }
        Ok(())
    }
}
