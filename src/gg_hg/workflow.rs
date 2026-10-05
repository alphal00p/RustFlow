//! Native canonical boundary generation and physical transport for Higgs + jet.
//! The embedded inputs contain exact mathematics and coordinates only. Numerical
//! reference values are deliberately held in a different, unlinked fixture.
use super::{PluginBasisMap, PluginFamilyKind};
use crate::algebraic::{CanonicalAlgebraicSystem, SquareRoot};
use crate::transport_cache::{
    BoundaryAccuracy, CachedBoundary, CachedPoint, EpsilonRange, PointKind, RootGerm, RootSheet,
    ScaledDistance,
};
use crate::{Error, FlowOptions, KinematicPoint, Precision, Result, RunContext, RustFlowCache};
use serde::Deserialize;
use std::collections::BTreeMap;
use symbolica::prelude::*;

const SYSTEMS: &str = include_str!("../../fixtures/gg-hg/integral-systems.json");
const CONFIGURATIONS: &str = include_str!("../../fixtures/gg-hg/physical-configurations.json");

/// A single canonical system shared by every mass and crossing configuration.
/// All coordinates are in units mV² = 1; b is mH²/mV².
pub struct HiggsJetIntegralSystem {
    pub map: PluginBasisMap,
    pub canonical: CanonicalAlgebraicSystem,
    pub transport: crate::RustFlow<CanonicalAlgebraicSystem>,
}

#[derive(Clone, Debug)]
pub struct HiggsJetConfiguration {
    pub label: String,
    pub mass: String,
    pub family: PluginFamilyKind,
    pub permutation: usize,
    pub start: BTreeMap<Symbol, Atom>,
    pub destination: BTreeMap<Symbol, Atom>,
    pub germ: RootGerm,
}

#[derive(Deserialize)]
struct Systems {
    schema: String,
    systems: Vec<System>,
}
#[derive(Deserialize)]
struct System {
    family: String,
    dimension: usize,
    coordinate_names: Vec<String>,
    letters: Vec<String>,
    roots: Vec<Root>,
    matrix_entries: Vec<Matrix>,
}
#[derive(Deserialize)]
struct Root {
    name: String,
    radicand: String,
}
#[derive(Deserialize)]
struct Matrix {
    letter: usize,
    entries: Vec<(usize, usize, String)>,
}
#[derive(Deserialize)]
struct Configurations {
    schema: String,
    cases: Vec<Configuration>,
}
#[derive(Deserialize)]
struct Configuration {
    label: String,
    mass: String,
    family: String,
    permutation: usize,
    dimension: usize,
    start: Vec<String>,
    destination: Vec<String>,
    source_root_germs: Vec<String>,
}

fn atom(text: &str, namespace: &str) -> Result<Atom> {
    Atom::parse(text, namespace.to_owned(), Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(text: &str, namespace: &str) -> Result<Symbol> {
    match atom(text, namespace)?.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(Error::InvalidInput(
            "scientific identifier must be a symbol".into(),
        )),
    }
}

impl HiggsJetIntegralSystem {
    pub fn load(kind: PluginFamilyKind, namespace: &str) -> Result<Self> {
        let map = PluginBasisMap::load(kind, namespace)?;
        let document: Systems =
            serde_json::from_str(SYSTEMS).map_err(|e| Error::InvalidInput(e.to_string()))?;
        if document.schema != "higgs-jet-integral-systems-v1" {
            return Err(Error::InvalidInput(
                "incompatible Higgs-jet system data".into(),
            ));
        }
        let system = document
            .systems
            .into_iter()
            .find(|s| s.family == kind.id())
            .ok_or_else(|| Error::InvalidInput("missing canonical system".into()))?;
        if system.dimension != kind.dimension() || system.coordinate_names != ["s", "t", "b"] {
            return Err(Error::InvalidInput(
                "canonical coordinates or dimension disagree with basis map".into(),
            ));
        }
        let roots = system
            .roots
            .iter()
            .map(|r| {
                Ok(SquareRoot {
                    symbol: symbol(&r.name, namespace)?,
                    radicand: atom(&r.radicand, namespace)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        if roots.len() != map.roots.len()
            || roots.iter().zip(&map.roots).any(|(a, b)| {
                a.symbol != b.symbol || !(&a.radicand - &b.radicand).together().cancel().is_zero()
            })
        {
            return Err(Error::InvalidInput(
                "canonical roots disagree with certified basis map".into(),
            ));
        }
        let letters = system
            .letters
            .iter()
            .map(|s| atom(s, namespace))
            .collect::<Result<Vec<_>>>()?;
        let n = kind.dimension();
        let mut matrices = vec![vec![vec![Atom::new(); n]; n]; letters.len()];
        for matrix in system.matrix_entries {
            if matrix.letter == 0 || matrix.letter > letters.len() {
                return Err(Error::InvalidInput(
                    "canonical letter index out of bounds".into(),
                ));
            }
            for (i, j, c) in matrix.entries {
                if i == 0 || j == 0 || i > n || j > n {
                    return Err(Error::InvalidInput(
                        "canonical matrix index out of bounds".into(),
                    ));
                }
                matrices[matrix.letter - 1][i - 1][j - 1] = atom(&c, namespace)?;
            }
        }
        let canonical = CanonicalAlgebraicSystem::new(
            map.epsilon,
            &map.coordinates,
            &letters,
            &matrices,
            roots,
        )?;
        let basis = (1..=n)
            .map(|i| atom(&format!("{}_F{i}", kind.id()), namespace))
            .collect::<Result<Vec<_>>>()?;
        let transport = crate::RustFlow::new_canonical(
            canonical.clone(),
            &basis,
            &atom("exp(2*eps*EulerGamma)", namespace)?,
            crate::Prescription::PlusI0,
            &format!(
                "arXiv:2112.07578 canonical ordinary integrals; unit vector mass; physical +i0; explicit root germ; physical-map-source={}",
                map.evidence.source_sha256
            ),
        )?;
        Ok(Self {
            map,
            canonical,
            transport,
        })
    }

    /// All eight exact physical configurations for this topology, including the
    /// W/Z mass ratios at one common physical point. No numerical seeds are read.
    pub fn configurations(&self, namespace: &str) -> Result<Vec<HiggsJetConfiguration>> {
        let document: Configurations =
            serde_json::from_str(CONFIGURATIONS).map_err(|e| Error::InvalidInput(e.to_string()))?;
        if document.schema != "higgs-jet-physical-configurations-v1" {
            return Err(Error::InvalidInput(
                "incompatible physical configurations".into(),
            ));
        }
        document
            .cases
            .into_iter()
            .filter(|c| c.family == self.map.kind.id())
            .map(|c| {
                if c.dimension != self.map.kind.dimension()
                    || c.start.len() != 3
                    || c.destination.len() != 3
                    || c.source_root_germs.len() != self.map.roots.len()
                {
                    return Err(Error::InvalidInput(
                        "incompatible physical configuration dimensions".into(),
                    ));
                }
                let coordinates = |values: Vec<String>| {
                    self.map
                        .coordinates
                        .iter()
                        .copied()
                        .zip(values)
                        .map(|(s, v)| Ok((s, atom(&v, namespace)?)))
                        .collect::<Result<BTreeMap<_, _>>>()
                };
                let germ = RootGerm {
                    sheets: self
                        .map
                        .roots
                        .iter()
                        .zip(c.source_root_germs)
                        .map(|(r, g)| {
                            Ok((
                                r.symbol,
                                match g.as_str() {
                                    "principal" => RootSheet::Principal,
                                    "opposite_principal" => RootSheet::Opposite,
                                    _ => {
                                        return Err(Error::InvalidInput(
                                            "unknown physical root germ".into(),
                                        ));
                                    }
                                },
                            ))
                        })
                        .collect::<Result<_>>()?,
                };
                Ok(HiggsJetConfiguration {
                    label: c.label,
                    mass: c.mass,
                    family: self.map.kind,
                    permutation: c.permutation,
                    start: coordinates(c.start)?,
                    destination: coordinates(c.destination)?,
                    germ,
                })
            })
            .collect()
    }

    /// Generate every canonical component from ordinary integrals using the
    /// automatic boundary recursion and independent epsilon-fit refinement.
    /// Force recomputation with reuse=false; exact reduction caching remains active.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_boundary(
        &self,
        cache: &mut RustFlowCache,
        coordinates: &BTreeMap<Symbol, Atom>,
        germ: &RootGerm,
        last: i32,
        options: &FlowOptions,
        backend: &dyn crate::ReductionBackend,
        context: &RunContext,
        reuse: bool,
    ) -> Result<CachedBoundary> {
        options.validate()?;
        context.cancellation.check()?;
        if options.dimension != 4 || options.prescription != crate::Prescription::PlusI0 {
            return Err(Error::Unsupported(
                "the certified Higgs-jet canonical basis uses D=4-2 epsilon and +i0".into(),
            ));
        }
        let point = CachedPoint::Exact(coordinates.clone()).with_root_germ(germ.clone())?;
        self.transport.identity().validate_point(&point)?;
        let range = EpsilonRange::new(0, last)?;
        let point_key = point.key()?;
        if reuse {
            for boundary in cache.entries() {
                if boundary.identity.key() == self.transport.identity().key()
                    && boundary.point.key()? == point_key
                    && boundary.range.leading == 0
                    && boundary.range.last >= last
                    && boundary.accuracy.verified_digits() >= options.digits
                {
                    return Ok(boundary.clone());
                }
            }
        }
        let mut options = options.clone();
        options.reuse_samples &= reuse;
        let kinematics = KinematicPoint(
            coordinates
                .iter()
                .map(|(&s, a)| (Atom::var(s), a.clone()))
                .collect(),
        );
        let fits = crate::solve_integral_projections_normalized(
            &[(self.map.family()?, self.map.canonical_rows())],
            &kinematics,
            last,
            &options,
            backend,
            context,
            &self.map.factors(germ.clone()),
        )?;
        let bits = fits
            .iter()
            .map(|f| f.working_bits)
            .min()
            .ok_or_else(|| Error::Numerical("empty canonical fit".into()))?;
        let digits = fits
            .iter()
            .map(|f| {
                f.verified_digits
                    .ok_or_else(|| Error::Accuracy("unverified canonical fit".into()))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .min()
            .unwrap();
        let p = Precision { bits };
        for fit in &fits {
            if fit.coefficients.iter().any(|(&power, c)| {
                power < 0 && p.norm(c) + &fit.comparison_errors[&power] > p.tolerance(digits)
            }) {
                return Err(Error::Accuracy(
                    "canonical negative-epsilon poles did not cancel".into(),
                ));
            }
        }
        let coefficients = (0..=last)
            .map(|power| {
                fits.iter()
                    .map(|f| p.round(&f.coefficients[&power]))
                    .collect()
            })
            .collect();
        let errors = (0..=last)
            .map(|power| {
                fits.iter()
                    .map(|f| {
                        let magnitude = p.norm(&f.coefficients[&power]);
                        let scale = if magnitude > p.real(1) {
                            magnitude
                        } else {
                            p.real(1)
                        };
                        f.comparison_errors[&power].clone()
                            + p.tolerance(options.digits + options.guard_digits / 2) * scale
                    })
                    .collect()
            })
            .collect();
        let provenance = format!(
            "Native automatic auxiliary-mass boundary; {}; certified exact map {}; independent epsilon samples and precision/order refinement; {} sample sets; no numerical reference seed",
            self.map.kind.id(),
            self.map.evidence.source_sha256,
            fits[0].refinements + 1
        );
        let boundary = CachedBoundary {
            identity: self.transport.identity().clone(),
            point,
            kind: PointKind::Physical,
            range,
            coefficients,
            accuracy: BoundaryAccuracy::supplied(digits, bits, errors, &provenance)?,
        };
        cache.insert(boundary.clone())?;
        Ok(boundary)
    }

    /// Straight regular transport within the same real invariant-sign chamber.
    /// Root-domain and sheet checks remain those of the shared continuation engine.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: &RootGerm,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<crate::physical_transport::PhysicalResult> {
        let policy = ScaledDistance {
            scales: BTreeMap::new(),
            admissible: |source: &CachedBoundary, target: &CachedPoint| {
                Ok(self.chamber(&source.point.restart_coordinates()?)?
                    == self.chamber(&target.restart_coordinates()?)?)
            },
        };
        self.transport.evaluate_to(
            cache,
            destination,
            Some(germ),
            EpsilonRange::new(0, last)?,
            options,
            context,
            &policy,
        )
    }

    fn chamber(&self, point: &BTreeMap<Symbol, Atom>) -> Result<Vec<std::cmp::Ordering>> {
        let [s, t, b] = self.map.coordinates.map(|s| {
            point
                .get(&s)
                .cloned()
                .ok_or_else(|| Error::InvalidInput("missing physical coordinate".into()))
        });
        let (s, t, b) = (s?, t?, b?);
        [s.clone(),t.clone(),&b-&s-&t,b.clone(),Atom::num(4)-b].into_iter().map(|a| {
            let a = a.together().cancel();
            if let AtomView::Num(n) = a.as_view()
                && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
                && c.im.is_zero() && !c.re.is_zero()
            { return Ok(c.re.cmp(&Rational::from(0))); }
            Err(Error::Unsupported("the Higgs-jet straight-path policy requires nonzero exact real invariants away from b=4".into()))
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_sixteen_configurations_share_one_exact_physical_point() -> Result<()> {
        let namespace = "higgs_jet_configuration_test";
        let s = atom("7173070292440521/111284741846000", namespace)?;
        let t = atom("-12058167788971/339319588980", namespace)?;
        let u = Atom::num(1) - &s - &t;
        let expected = [(&s, &u), (&s, &t), (&u, &t), (&t, &s)];
        let mut count = 0;
        for kind in [PluginFamilyKind::Planar, PluginFamilyKind::Nonplanar] {
            let system = HiggsJetIntegralSystem::load(kind, namespace)?;
            let configurations = system.configurations(namespace)?;
            assert_eq!(configurations.len(), 8);
            for c in configurations {
                let [ss, tt, b] = system.map.coordinates.map(|s| c.destination[&s].clone());
                let (se, te) = expected[c.permutation - 1];
                assert!((ss / &b - se).together().cancel().is_zero());
                assert!((tt / b - te).together().cancel().is_zero());
                assert_eq!(c.germ.sheets.len(), system.map.roots.len());
                assert_eq!(system.chamber(&c.start)?, system.chamber(&c.destination)?);
                system
                    .transport
                    .identity()
                    .validate_point(&CachedPoint::Exact(c.start).with_root_germ(c.germ)?)?;
                count += 1;
            }
        }
        assert_eq!(count, 16);
        Ok(())
    }
}
