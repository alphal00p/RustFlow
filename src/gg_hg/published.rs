use crate::algebraic::SquareRoot;
use crate::{
    ComplexFloat, Error, Integral, IntegralFamily, Precision, Prescription, Propagator, Result,
};
use serde::Deserialize;
use std::collections::BTreeSet;
use symbolica::prelude::*;

const FAMILY_DATA: &str = include_str!("../../fixtures/gg-hg/published-physical-families.json");
const CANONICAL_DATA: &str = include_str!("../../fixtures/gg-hg/published-canonical-64.json");
const NAMESPACE: &str = "rustflow_gg_hg_published_2020";

/// Primary-source provenance of independently serialized mathematical inputs.
#[derive(Clone, Debug, Deserialize)]
pub struct SourceEvidence {
    pub url: String,
    pub sha256: String,
}

/// Exact family definitions; these names never denote plugin master indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PublishedFamilyKind {
    Planar2018,
    Planar2020,
    Nonplanar2020,
}

impl PublishedFamilyKind {
    pub fn id(self) -> &'static str {
        match self {
            Self::Planar2018 => "1810.05138v1/PL",
            Self::Planar2020 => "2007.09813v2/PL",
            Self::Nonplanar2020 => "2007.09813v2/NP",
        }
    }
}

/// The paper's multiplicative measure relative to RustFlow's ordinary measure.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PublishedMeasure {
    GammaOnePlusEpsilonInversePerLoop,
    ExpEpsilonEulerGammaTimesMassOverMuSquaredToEpsilonPerLoop,
}

/// Native ordinary integrals, with the paper's normalization kept separately.
/// All stored denominators have the q²-m² sign and use the usual +i0 convention.
#[derive(Clone, Debug)]
pub struct PublishedFamily {
    pub family: IntegralFamily,
    pub targets: Vec<Integral>,
    pub measure: PublishedMeasure,
    pub source: SourceEvidence,
    published_denominator_sign: i64,
}

#[derive(Deserialize)]
struct FamilyDocument {
    schema: String,
    physical_propagators: usize,
    dimension: i64,
    families: Vec<FamilyRecord>,
    provenance: Vec<SourceEvidence>,
}
#[derive(Deserialize)]
struct FamilyRecord {
    id: String,
    paper: String,
    routing: Vec<Routing>,
    published_denominator_sign: i64,
    measure: PublishedMeasure,
    targets: Vec<Vec<i16>>,
}
#[derive(Deserialize)]
struct Routing {
    #[serde(rename = "loop")]
    loops: [i64; 2],
    external: [i64; 3],
    mass_squared_multiple: i64,
}

impl PublishedFamily {
    /// Use s=(p1+p2)², t=(p1+p3)², u=(p2+p3)² and massless external p1,p2,p3.
    /// The Higgs virtuality is s+t+u; mv2 denotes a squared vector-boson mass.
    /// Inputs remain native exact Atoms, with no numerical specialization here.
    pub fn new(
        kind: PublishedFamilyKind,
        s: Atom,
        t: Atom,
        u: Atom,
        mv2: Atom,
        epsilon: Symbol,
    ) -> Result<Self> {
        let document: FamilyDocument = serde_json::from_str(FAMILY_DATA)
            .map_err(|e| Error::InvalidInput(format!("published family data: {e}")))?;
        if document.schema != "rustflow-published-gg-hg-families-v1"
            || document.physical_propagators != 7
            || document.dimension != 4
        {
            return Err(Error::InvalidInput(
                "incompatible published ggHg family data".into(),
            ));
        }
        let record = document
            .families
            .into_iter()
            .find(|r| r.id == kind.id())
            .ok_or_else(|| Error::InvalidInput("missing published ggHg family".into()))?;
        if record.routing.len() != 9 || ![-1, 1].contains(&record.published_denominator_sign) {
            return Err(Error::InvalidInput(
                "invalid published propagator routing".into(),
            ));
        }
        let source = document
            .provenance
            .into_iter()
            .find(|p| p.url.ends_with(&record.paper))
            .ok_or_else(|| Error::InvalidInput("missing primary-source provenance".into()))?;
        let zero = Atom::new();
        let (hs, ht, hu) = (s / Atom::num(2), t / Atom::num(2), u / Atom::num(2));
        let gram = vec![
            vec![zero.clone(), hs.clone(), ht.clone()],
            vec![hs, zero.clone(), hu.clone()],
            vec![ht, hu, zero],
        ];
        let propagators = record
            .routing
            .iter()
            .map(|r| {
                Propagator::quadratic(
                    &r.loops,
                    &r.external,
                    &mv2 * Atom::num(r.mass_squared_multiple),
                    &gram,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let family = IntegralFamily {
            name: format!(
                "published_{}_positive_denominators",
                record.id.replace(['.', '/'], "_")
            ),
            loops: vec!["k".into(), "l".into()],
            external: vec!["p1".into(), "p2".into(), "p3".into()],
            external_gram: gram,
            propagators,
            physical_propagators: 7,
            epsilon,
            dimension: 4,
        };
        let targets = record.targets.into_iter().map(Integral).collect::<Vec<_>>();
        for integral in &targets {
            family.validate_integral(integral)?;
        }
        family.validate()?;
        Ok(Self {
            family,
            targets,
            measure: record.measure,
            source,
            published_denominator_sign: record.published_denominator_sign,
        })
    }

    /// Exact denominator-sign factor; the publication's measure is additional.
    /// Signed numerator powers participate too. No noninteger power is implied.
    pub fn denominator_factor(&self, integral: &Integral) -> Result<Atom> {
        self.family.validate_integral(integral)?;
        let odd = integral.0.iter().map(|&n| i64::from(n)).sum::<i64>() % 2 != 0;
        Ok(Atom::num(if odd {
            self.published_denominator_sign
        } else {
            1
        }))
    }
}

/// Source labels are preserved verbatim. The crossed labels do not by
/// themselves certify a permutation of external momenta or invariant variables.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
pub enum CanonicalIntegralHead {
    #[serde(rename = "PL")]
    Planar,
    #[serde(rename = "NP")]
    Nonplanar,
    #[serde(rename = "PLx12")]
    PlanarX12,
    #[serde(rename = "PLx123")]
    PlanarX123,
}

#[derive(Clone, Debug)]
pub struct CanonicalIntegral {
    pub head: CanonicalIntegralHead,
    pub integral: Integral,
    pub dimension: i64,
}
#[derive(Clone, Debug)]
pub struct CanonicalTerm {
    /// Zero-based position in the published integral reference list.
    pub integral: usize,
    pub coefficient: Atom,
}

/// The 64-component published basis from arXiv:2007.09813v2.
/// Its coefficients are rational in sh,y,z,rho,eps and named square roots.
/// This is not a certified physical-to-plugin 48/61 map.
#[derive(Clone, Debug)]
pub struct PublishedCanonicalBasis {
    pub integrals: Vec<CanonicalIntegral>,
    pub rows: Vec<Vec<CanonicalTerm>>,
    pub roots: Vec<SquareRoot>,
    pub epsilon: Symbol,
    pub variables: [Symbol; 4],
    pub normalization: PublishedNormalization,
    pub source: SourceEvidence,
}

#[derive(Deserialize)]
struct CanonicalDocument {
    schema: String,
    paper: String,
    dimension: usize,
    normalization: String,
    roots: Vec<RootRecord>,
    integrals: Vec<IntegralRecord>,
    rows: Vec<Vec<TermRecord>>,
    provenance: SourceEvidence,
    extraction: Extraction,
}
#[derive(Deserialize)]
struct RootRecord {
    name: String,
    radicand: String,
}
#[derive(Deserialize)]
struct IntegralRecord {
    head: CanonicalIntegralHead,
    powers: Vec<i16>,
    dimension: i64,
}
#[derive(Deserialize)]
struct TermRecord {
    integral: usize,
    coefficient: String,
}
#[derive(Deserialize)]
struct Extraction {
    exact_linear_reconstruction: bool,
    exact_root_reconstruction: bool,
    source_sha256: String,
}

fn atom(text: &str) -> Result<Atom> {
    Atom::parse(text, NAMESPACE, Default::default()).map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(text: &str) -> Result<Symbol> {
    match atom(text)?.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(Error::InvalidInput(
            "published identifier is not a symbol".into(),
        )),
    }
}

impl PublishedCanonicalBasis {
    pub fn load() -> Result<Self> {
        let data: CanonicalDocument = serde_json::from_str(CANONICAL_DATA)
            .map_err(|e| Error::InvalidInput(format!("published canonical data: {e}")))?;
        if data.schema != "rustflow-published-gg-hg-canonical-v1"
            || data.paper != "2007.09813v2"
            || data.dimension != 64
            || data.rows.len() != 64
            || data.integrals.len() != 64
            || data.normalization != "(-mh2)^(2*eps)/Gamma(1+eps)^2"
            || !data.extraction.exact_linear_reconstruction
            || !data.extraction.exact_root_reconstruction
            || data.extraction.source_sha256 != data.provenance.sha256
        {
            return Err(Error::InvalidInput(
                "incompatible published canonical data".into(),
            ));
        }
        let epsilon = symbol("eps")?;
        let variables = [symbol("sh")?, symbol("y")?, symbol("z")?, symbol("rho")?];
        let mut allowed = variables
            .iter()
            .map(|&v| Atom::var(v))
            .collect::<BTreeSet<_>>();
        let mut roots = Vec::new();
        for root in data.roots {
            let radicand = atom(&root.radicand)?;
            let mut found = BTreeSet::new();
            crate::family::scalar_symbols(radicand.as_view(), &mut found)?;
            if !found.is_subset(&allowed) {
                return Err(Error::InvalidInput(
                    "published root depends on undeclared variables".into(),
                ));
            }
            roots.push(SquareRoot {
                symbol: symbol(&root.name)?,
                radicand,
            });
        }
        allowed.extend(roots.iter().map(|r| Atom::var(r.symbol)));
        allowed.insert(Atom::var(epsilon));
        let integrals = data
            .integrals
            .into_iter()
            .map(|r| {
                if r.dimension != 4 || r.powers.len() != 9 || r.powers[7..].iter().any(|&n| n > 0) {
                    return Err(Error::InvalidInput(
                        "invalid published canonical integral".into(),
                    ));
                }
                Ok(CanonicalIntegral {
                    head: r.head,
                    integral: Integral(r.powers),
                    dimension: r.dimension,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let rows = data
            .rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|term| {
                        if term.integral >= integrals.len() {
                            return Err(Error::InvalidInput(
                                "published canonical index out of range".into(),
                            ));
                        }
                        let coefficient = atom(&term.coefficient)?;
                        let mut found = BTreeSet::new();
                        crate::family::scalar_symbols(coefficient.as_view(), &mut found)?;
                        if !found.is_subset(&allowed) {
                            return Err(Error::InvalidInput(
                                "published coefficient uses unknown variables".into(),
                            ));
                        }
                        let _: RationalPolynomial<IntegerRing, u16> = coefficient
                            .try_to_rational_polynomial(&Q, &Z, None)
                            .map_err(|e| {
                                Error::InvalidInput(format!("nonrational canonical weight: {e}"))
                            })?;
                        Ok(CanonicalTerm {
                            integral: term.integral,
                            coefficient,
                        })
                    })
                    .collect()
            })
            .collect::<Result<Vec<Vec<_>>>>()?;
        Ok(Self {
            integrals,
            rows,
            roots,
            epsilon,
            variables,
            normalization: PublishedNormalization::HiggsScaleAndGamma2020,
            source: data.provenance,
        })
    }

    /// Until primary-source crossing semantics are certified, these references
    /// cannot be silently assigned to a physical family and AMF evaluated.
    pub fn require_resolved_crossings(&self) -> Result<()> {
        if self.integrals.iter().any(|i| {
            matches!(
                i.head,
                CanonicalIntegralHead::PlanarX12 | CanonicalIntegralHead::PlanarX123
            )
        }) {
            return Err(Error::Unsupported(
                "published PLx12/PLx123 crossing semantics require an exact routing certificate; the labels alone are not one".into()));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PublishedNormalization {
    /// Equation (22) and the ancillary common (-sh)^(2 epsilon), for two loops.
    HiggsScaleAndGamma2020,
}

impl PublishedNormalization {
    /// Evaluate the common factor before fitting in epsilon. For physical
    /// mh²>0 the +i0 prescription means log(-mh²-i0), not a principal +iπ.
    /// Only real nonzero exact Higgs virtualities are admitted by this helper.
    pub fn evaluate(
        self,
        epsilon: &Rational,
        higgs_squared: &Rational,
        prescription: Prescription,
        p: Precision,
    ) -> Result<ComplexFloat> {
        if higgs_squared.is_zero() {
            return Err(Error::InvalidInput(
                "published scale normalization requires nonzero mh²".into(),
            ));
        }
        let scale = p.rational(higgs_squared);
        let logarithm = if higgs_squared > &Rational::from(0) {
            let mut logarithm = p.log(&scale);
            let pi = p.real(0).pi();
            logarithm.im = match prescription {
                Prescription::PlusI0 => -pi,
                Prescription::MinusI0 => pi,
            };
            logarithm
        } else {
            p.log(&p.neg(&scale))
        };
        let power = p.exp(&p.mul(&p.rational(&(epsilon * &Rational::from(2))), &logarithm));
        let gamma = p.gamma_real(&p.rational(&(epsilon + &Rational::from(1))).re)?;
        let result = p.div(&power, &p.mul(&gamma, &gamma));
        if !p.finite(&result) {
            return Err(Error::Numerical(
                "nonfinite published canonical normalization".into(),
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "automatic")]
    fn family(kind: PublishedFamilyKind) -> Result<PublishedFamily> {
        PublishedFamily::new(
            kind,
            Atom::num(-2),
            Atom::num(-3),
            Atom::num(-5),
            Atom::num(7),
            symbol("eps")?,
        )
    }

    #[test]
    #[cfg(feature = "automatic")]
    fn published_routings_keep_mass_assignments_numerators_and_signs() -> Result<()> {
        let planar = family(PublishedFamilyKind::Planar2020)?;
        let nonplanar = family(PublishedFamilyKind::Nonplanar2020)?;
        let old = family(PublishedFamilyKind::Planar2018)?;
        assert_eq!(
            (
                old.targets.len(),
                planar.targets.len(),
                nonplanar.targets.len()
            ),
            (48, 45, 21)
        );
        for (f, massive) in [
            (&planar, vec![1, 6]),
            (&nonplanar, vec![3, 5, 7]),
            (&old, vec![5, 6]),
        ] {
            for i in 0..9 {
                assert_eq!(
                    f.family.mass_squared(i)?,
                    Atom::num(if massive.contains(&i) { 7 } else { 0 })
                );
            }
        }
        assert_eq!(old.targets[42], Integral(vec![1, 1, 1, -1, 1, 1, 1, 0, 0]));
        assert_eq!(old.denominator_factor(&old.targets[0])?, Atom::num(-1));
        assert_eq!(old.denominator_factor(&old.targets[10])?, Atom::one());
        assert_eq!(planar.denominator_factor(&planar.targets[0])?, Atom::one());
        // (k-p1-p2-p3)^2-m² has external constant mh²-m².
        assert_eq!(planar.family.propagators[6].constant, Atom::num(-17));
        Ok(())
    }

    #[test]
    fn published_canonical_data_is_exact_but_never_claims_the_plugin_map() -> Result<()> {
        let basis = PublishedCanonicalBasis::load()?;
        assert_eq!(
            (basis.rows.len(), basis.integrals.len(), basis.roots.len()),
            (64, 64, 7)
        );
        assert_eq!(basis.rows.iter().map(Vec::len).sum::<usize>(), 191);
        let first = &basis.rows[0][0];
        assert_eq!(first.integral, 0);
        assert!(
            (&first.coefficient - atom("(1-eps)*eps^2")?)
                .expand()
                .is_zero()
        );
        assert!(matches!(
            basis.require_resolved_crossings(),
            Err(Error::Unsupported(_))
        ));
        Ok(())
    }

    #[test]
    #[cfg(feature = "native")]
    fn finite_epsilon_normalization_preserves_the_physical_log_sheet() -> Result<()> {
        let norm = PublishedNormalization::HiggsScaleAndGamma2020;
        let p = Precision::decimal(70)?;
        let epsilon = Rational::from((1, 7));
        let positive = Rational::from(4);
        let plus = norm.evaluate(&epsilon, &positive, Prescription::PlusI0, p)?;
        let minus = norm.evaluate(&epsilon, &positive, Prescription::MinusI0, p)?;
        assert!(plus.im < p.real(0));
        assert_eq!(plus.re, minus.re);
        assert_eq!(plus.im, -minus.im);
        let euclidean = norm.evaluate(&epsilon, &Rational::from(-4), Prescription::PlusI0, p)?;
        assert!(euclidean.im.is_zero());
        assert!(p.close(
            &norm.evaluate(&Rational::from(0), &positive, Prescription::PlusI0, p)?,
            &p.i(1),
            60
        ));
        let refined = Precision::decimal(90)?;
        assert!(p.close(
            &plus,
            &refined.round(&norm.evaluate(&epsilon, &positive, Prescription::PlusI0, refined)?),
            60
        ));
        assert!(
            norm.evaluate(&epsilon, &Rational::from(0), Prescription::PlusI0, p)
                .is_err()
        );
        Ok(())
    }
}
