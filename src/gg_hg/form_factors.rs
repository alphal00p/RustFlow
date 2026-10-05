//! Generic finite form factors from exact primary expressions and certified masters.
use super::{HiggsJetIntegralSystem, PluginFamilyKind, SourceEvidence};
use crate::algebraic::SquareRoot;
use crate::transport_cache::{BoundaryIdentity, CachedBoundary, RootSheet};
use crate::{ComplexFloat, Error, KinematicPoint, Precision, Result};
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

const DATA: &str = include_str!("../../fixtures/gg-hg/generic-form-factors.json");

#[derive(Clone, Debug)]
pub struct HiggsJetFormFactorResult {
    pub values: [ComplexFloat; 4],
    pub absolute_errors: [Float; 4],
    /// Conservative relative digits conditional on the supplied master errors.
    /// An unresolved zero has no relative precision claim.
    pub verified_relative_digits: [Option<u32>; 4],
    pub provenance: String,
}

/// The four finite scalar form factors, valid at variable exact real kinematics.
/// Coefficients contain no numerical boundary data. Their physical square-root
/// sides come from the source's explicit infinitesimals, evaluated as limits.
pub struct HiggsJetFormFactors {
    pub coordinates: [Symbol; 3],
    pub source: SourceEvidence,
    rows: Vec<Vec<Weight>>,
    roots: Vec<CoefficientRoot>,
    identities: [BoundaryIdentity; 2],
    canonical_roots: [Vec<SquareRoot>; 2],
}

struct Weight {
    family: usize,
    permutation: usize,
    index: usize,
    epsilon: i32,
    coefficient: Atom,
}
struct CoefficientRoot {
    symbol: Symbol,
    radicand: Atom,
    i0_slope: Atom,
}
#[derive(Deserialize)]
struct Document {
    schema: String,
    source: SourceEvidence,
    coordinates: Vec<String>,
    normalization: String,
    crossings: Vec<Vec<String>>,
    roots: Vec<RootRecord>,
    rows: Vec<Vec<WeightRecord>>,
    evidence: Evidence,
}
#[derive(Deserialize)]
struct Evidence {
    exact_primary_equality_at_zero_i0: bool,
    exact_linear_reconstruction: bool,
    exact_root_reconstruction: bool,
}
#[derive(Deserialize)]
struct RootRecord {
    name: String,
    radicand: String,
    i0_slope: String,
}
#[derive(Deserialize)]
struct WeightRecord {
    family: String,
    permutation: usize,
    index: usize,
    epsilon: i32,
    coefficient: String,
}

fn parse(text: &str, namespace: &str) -> Result<Atom> {
    Atom::parse(text, namespace.to_owned(), Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(text: &str, namespace: &str) -> Result<Symbol> {
    match parse(text, namespace)?.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(Error::InvalidInput("form-factor scalar identifier".into())),
    }
}
fn rational(value: &Atom) -> Result<Rational> {
    let value = value.together().cancel();
    if let AtomView::Num(n) = value.as_view()
        && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
        && c.im.is_zero()
    {
        return Ok(c.re);
    }
    Err(Error::Unsupported(
        "Higgs-jet form factors require exact real rational physical inputs".into(),
    ))
}
fn exact_scalar(value: &Atom, allowed: &BTreeSet<Atom>) -> Result<()> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(value.as_view(), &mut symbols)?;
    if !symbols.is_subset(allowed) {
        return Err(Error::InvalidInput("undeclared form-factor scalar".into()));
    }
    let _: RationalPolynomial<IntegerRing, u16> = value
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::InvalidInput(format!("nonrational form-factor coefficient: {e}")))?;
    Ok(())
}

impl HiggsJetFormFactors {
    pub fn load(namespace: &str) -> Result<Self> {
        let data: Document = serde_json::from_str(DATA)
            .map_err(|e| Error::InvalidInput(format!("form-factor data: {e}")))?;
        if data.schema != "rustflow-gg-hg-generic-form-factors-v1"
            || data.coordinates != ["s", "t", "b"]
            || data.normalization != "-1/(mV2^2*(4*pi)^4)"
            || data.rows.len() != 4
            || !data.evidence.exact_primary_equality_at_zero_i0
            || !data.evidence.exact_linear_reconstruction
            || !data.evidence.exact_root_reconstruction
        {
            return Err(Error::InvalidInput(
                "uncertified generic form-factor data".into(),
            ));
        }
        let planar = HiggsJetIntegralSystem::load(PluginFamilyKind::Planar, namespace)?;
        let nonplanar = HiggsJetIntegralSystem::load(PluginFamilyKind::Nonplanar, namespace)?;
        let coordinates = planar.map.coordinates;
        let mut allowed = coordinates
            .into_iter()
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        let roots = data
            .roots
            .into_iter()
            .map(|r| {
                let radicand = parse(&r.radicand, namespace)?;
                let i0_slope = parse(&r.i0_slope, namespace)?;
                exact_scalar(&radicand, &allowed)?;
                exact_scalar(&i0_slope, &allowed)?;
                Ok(CoefficientRoot {
                    symbol: symbol(&r.name, namespace)?,
                    radicand,
                    i0_slope,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        allowed.extend(roots.iter().map(|r| Atom::var(r.symbol)));
        let source_i = Atom::var(symbol("imaginary_unit", namespace)?);
        let native_i = Atom::var(crate::family::imaginary_parameter());
        allowed.insert(native_i.clone());
        let rename = BTreeMap::from([(source_i, native_i)]);
        let rows = data
            .rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|w| {
                        let (family, dimension) = match w.family.as_str() {
                            "Planar_EW1" => (0, 48),
                            "NP_EW1" => (1, 61),
                            _ => {
                                return Err(Error::InvalidInput(
                                    "unknown form-factor family".into(),
                                ));
                            }
                        };
                        if !(1..=4).contains(&w.permutation)
                            || w.index >= dimension
                            || !(0..=4).contains(&w.epsilon)
                        {
                            return Err(Error::InvalidInput(
                                "invalid form-factor master index".into(),
                            ));
                        }
                        let coefficient =
                            crate::family::substitute(&parse(&w.coefficient, namespace)?, &rename);
                        exact_scalar(&coefficient, &allowed)?;
                        Ok(Weight {
                            family,
                            permutation: w.permutation - 1,
                            index: w.index,
                            epsilon: w.epsilon,
                            coefficient,
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let expected = Self::crossings(coordinates.map(Atom::var));
        if data.crossings.len() != 4 || data.crossings.iter().any(|c| c.len() != 3) {
            return Err(Error::InvalidInput("invalid crossing data".into()));
        }
        for (source, expected) in data.crossings.into_iter().zip(expected) {
            for (source, expected) in source.into_iter().zip(expected) {
                if !(parse(&source, namespace)? - expected).expand().is_zero() {
                    return Err(Error::InvalidInput(
                        "form-factor crossing convention mismatch".into(),
                    ));
                }
            }
        }
        Ok(Self {
            coordinates,
            source: data.source,
            rows,
            roots,
            identities: [
                planar.transport.identity().clone(),
                nonplanar.transport.identity().clone(),
            ],
            canonical_roots: [planar.map.roots, nonplanar.map.roots],
        })
    }

    /// Crossings, in the source's order: (s,u,b), (s,t,b), (u,t,b), (t,s,b).
    pub fn crossings([s, t, b]: [Atom; 3]) -> [[Atom; 3]; 4] {
        let u = (&b - &s - &t).expand();
        [
            [s.clone(), u.clone(), b.clone()],
            [s.clone(), t.clone(), b.clone()],
            [u, t.clone(), b.clone()],
            [t, s, b],
        ]
    }

    /// All physical inputs remain exact until numerical coefficient evaluation.
    /// Each topology slice has exactly four boundaries in `crossings` order.
    /// Boundary identities, coordinates, epsilon ranges and root germs are checked
    /// before contraction; the returned errors propagate every input component.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        physical_s: &Atom,
        physical_t: &Atom,
        higgs_mass_squared: &Atom,
        vector_mass_squared: &Atom,
        planar: &[CachedBoundary],
        nonplanar: &[CachedBoundary],
        p: Precision,
    ) -> Result<HiggsJetFormFactorResult> {
        let (s, t, h, mv) = (
            rational(physical_s)?,
            rational(physical_t)?,
            rational(higgs_mass_squared)?,
            rational(vector_mass_squared)?,
        );
        if h <= 0 || mv <= 0 {
            return Err(Error::InvalidInput(
                "Higgs and vector mass squares must be positive".into(),
            ));
        }
        let dimensionless = [
            Atom::num(&s / &mv),
            Atom::num(&t / &mv),
            Atom::num(&h / &mv),
        ];
        let crossing_points = Self::crossings(dimensionless.clone());
        let boundaries = [planar, nonplanar];
        let mut input_cap = u32::MAX;
        for (kind, boundaries) in boundaries.iter().enumerate() {
            if boundaries.len() != 4 {
                return Err(Error::InvalidInput(
                    "four crossed boundaries are required per topology".into(),
                ));
            }
            for (permutation, boundary) in boundaries.iter().enumerate() {
                boundary.validate()?;
                if boundary.identity.key() != self.identities[kind].key()
                    || boundary.range.leading > 0
                    || boundary.range.last < 4
                {
                    return Err(Error::InvalidInput(
                        "form-factor boundary basis or epsilon range mismatch".into(),
                    ));
                }
                let expected = self
                    .coordinates
                    .into_iter()
                    .zip(crossing_points[permutation].clone())
                    .collect::<BTreeMap<_, _>>();
                let actual = boundary.point.restart_coordinates()?;
                if actual.len() != expected.len()
                    || expected.iter().any(|(key, value)| {
                        actual
                            .get(key)
                            .is_none_or(|a| !(a - value).together().cancel().is_zero())
                    })
                {
                    return Err(Error::InvalidInput(
                        "form-factor boundary has different physical kinematics or crossing".into(),
                    ));
                }
                self.check_germ(kind, boundary, &expected)?;
                input_cap = input_cap.min(boundary.accuracy.verified_digits());
            }
        }
        let point = KinematicPoint(
            self.coordinates
                .into_iter()
                .map(Atom::var)
                .zip(dimensionless)
                .collect(),
        );
        let parameters = self.parameters(&point, p)?;
        let pi = ComplexFloat::new(Float::with_val(p.bits, rug::float::Constant::Pi), p.real(0));
        let normalization = p.div(
            &p.i(-1),
            &p.mul(
                &p.powi(&p.rational(&mv), 2),
                &p.powi(&p.scale(&pi, 4, 1), 4),
            ),
        );
        let absolute_normalization = p.norm(&normalization);
        let arithmetic_digits = ((p.bits.saturating_sub(16) as u64) * 1000 / 3322) as u32;
        let reserve = p.tolerance(arithmetic_digits.saturating_sub(12));
        let mut values = std::array::from_fn(|_| p.zero());
        let mut absolute_errors = std::array::from_fn(|_| p.real(0));
        for (row_index, row) in self.rows.iter().enumerate() {
            let mut sum = p.zero();
            let mut error = p.real(0);
            let mut magnitude = p.real(0);
            for w in row {
                let boundary = &boundaries[w.family][w.permutation];
                let order = (w.epsilon - boundary.range.leading) as usize;
                let coefficient = p.eval(&w.coefficient, &parameters)?;
                let value = p.mul(&coefficient, &boundary.coefficients[order][w.index]);
                sum = p.add(&sum, &value);
                magnitude += p.norm(&value);
                error +=
                    p.norm(&coefficient) * &boundary.accuracy.comparison_errors()[order][w.index];
            }
            values[row_index] = p.mul(&normalization, &sum);
            absolute_errors[row_index] = absolute_normalization.clone()
                * (error + reserve.clone() * (magnitude + p.real(1)));
        }
        let verified_relative_digits = std::array::from_fn(|i| {
            let norm = p.norm(&values[i]);
            if norm <= absolute_errors[i] {
                return None;
            }
            (1..=input_cap.min(arithmetic_digits))
                .rev()
                .find(|&digits| absolute_errors[i] <= p.tolerance(digits) * &norm)
        });
        Ok(HiggsJetFormFactorResult {
            values,
            absolute_errors,
            verified_relative_digits,
            provenance: format!(
                "Exact 2112.07578 form factors; source {}; explicit i0 side limits; weighted canonical boundary errors plus arithmetic reserve; physical factor -1/(mV²)^2/(4π)^4",
                self.source.sha256
            ),
        })
    }

    fn parameters(
        &self,
        point: &KinematicPoint,
        p: Precision,
    ) -> Result<ahash::HashMap<Atom, ComplexFloat>> {
        let mut parameters = point
            .0
            .iter()
            .map(|(s, v)| Ok((s.clone(), p.eval(v, &Default::default())?)))
            .collect::<Result<ahash::HashMap<_, _>>>()?;
        for root in &self.roots {
            let radicand = rational(&point.apply(&root.radicand))?;
            let slope = rational(&point.apply(&root.i0_slope))?;
            if radicand.is_zero() {
                return Err(Error::Unsupported(
                    "form-factor root is on a branch point".into(),
                ));
            }
            let principal = crate::algebraic::principal_sqrt(p, &p.rational(&radicand))?;
            let value = if radicand < 0 && slope < 0 {
                p.neg(&principal)
            } else {
                principal
            };
            parameters.insert(Atom::var(root.symbol), value);
        }
        Ok(parameters)
    }

    fn check_germ(
        &self,
        kind: usize,
        boundary: &CachedBoundary,
        point: &BTreeMap<Symbol, Atom>,
    ) -> Result<()> {
        let germ = boundary.point.root_germ().ok_or_else(|| {
            Error::InvalidInput("canonical form-factor boundary needs an explicit root germ".into())
        })?;
        let [s, t, b] = self.coordinates.map(Atom::var);
        let variables = [s.clone(), t.clone(), b.clone(), &b - &s - &t];
        let thresholds = [vec![0, 1], vec![0, 1], vec![0, 1, 4], vec![0, 1, 4]];
        let kinematics = KinematicPoint(
            point
                .iter()
                .map(|(&s, a)| (Atom::var(s), a.clone()))
                .collect(),
        );
        for root in &self.canonical_roots[kind] {
            let mut expected = RootSheet::Principal;
            for (variable, thresholds) in variables.iter().zip(&thresholds) {
                for &threshold in thresholds {
                    if (&root.radicand - (Atom::num(threshold) - variable))
                        .expand()
                        .is_zero()
                        && rational(&kinematics.apply(variable))? > threshold
                    {
                        expected = RootSheet::Opposite;
                    }
                }
            }
            if germ.sheets.get(&root.symbol) != Some(&expected) {
                return Err(Error::InvalidInput(
                    "canonical root germ disagrees with the physical form-factor convention".into(),
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generic_weights_match_independent_point_specializations() -> Result<()> {
        // The numerical file is test-only: production loads mathematical data.
        let reference: serde_json::Value =
            serde_json::from_str(include_str!("../../fixtures/gg-hg/projections.json")).unwrap();
        let namespace = "rustflow_gg_hg_generic_weights_test";
        let projector = HiggsJetFormFactors::load(namespace)?;
        assert_eq!(projector.roots.len(), 14);
        let p = Precision::decimal(100)?;
        for case in reference["form_factors"].as_array().unwrap() {
            let point = KinematicPoint(
                projector
                    .coordinates
                    .into_iter()
                    .map(Atom::var)
                    .zip(
                        case["point"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .map(|s| parse(s.as_str().unwrap(), namespace))
                            .collect::<Result<Vec<_>>>()?,
                    )
                    .collect(),
            );
            let parameters = projector.parameters(&point, p)?;
            for (row, original) in projector
                .rows
                .iter()
                .zip(case["form_factors"].as_array().unwrap())
            {
                let expected = original["weights"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|w| {
                        let family = if w["family"] == "Planar_EW1" { 0 } else { 1 };
                        (
                            (
                                family,
                                w["permutation"].as_u64().unwrap() as usize - 1,
                                w["index"].as_u64().unwrap() as usize - 1,
                                w["epsilon"].as_i64().unwrap() as i32,
                            ),
                            w,
                        )
                    })
                    .collect::<BTreeMap<_, _>>();
                assert_eq!(row.len(), expected.len());
                for weight in row {
                    let key = (
                        weight.family,
                        weight.permutation,
                        weight.index,
                        weight.epsilon,
                    );
                    let r = &expected[&key]["coefficient_value"];
                    let original = p.parse(
                        r["real"].as_str().unwrap(),
                        r["imaginary"].as_str().unwrap(),
                    )?;
                    let actual = p.eval(&weight.coefficient, &parameters)?;
                    assert!(
                        p.norm(&p.sub(&actual, &original))
                            <= p.tolerance(70) * (p.norm(&original) + p.real(1)),
                        "{} coefficient {key:?}",
                        case["mass"]
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn coefficient_roots_use_exact_side_limits_and_crossings() -> Result<()> {
        let namespace = "rustflow_gg_hg_coefficient_roots_test";
        let projector = HiggsJetFormFactors::load(namespace)?;
        let p = Precision::decimal(70)?;
        for b in [2, 5] {
            let point = KinematicPoint(
                projector
                    .coordinates
                    .into_iter()
                    .map(Atom::var)
                    .zip([Atom::num(10), Atom::num(-3), Atom::num(b)])
                    .collect(),
            );
            let parameters = projector.parameters(&point, p)?;
            let negative = projector
                .roots
                .iter()
                .filter(|root| rational(&point.apply(&root.radicand)).unwrap() < Rational::from(0))
                .collect::<Vec<_>>();
            assert!(!negative.is_empty());
            for root in negative {
                let value = &parameters[&Atom::var(root.symbol)];
                let lower = rational(&point.apply(&root.i0_slope))? < Rational::from(0);
                assert_eq!(value.im < p.real(0), lower);
            }
        }
        assert_eq!(
            HiggsJetFormFactors::crossings([Atom::num(10), Atom::num(-3), Atom::num(2)]),
            [
                [Atom::num(10), Atom::num(-5), Atom::num(2)],
                [Atom::num(10), Atom::num(-3), Atom::num(2)],
                [Atom::num(-5), Atom::num(-3), Atom::num(2)],
                [Atom::num(-3), Atom::num(10), Atom::num(2)],
            ]
        );
        Ok(())
    }
}
