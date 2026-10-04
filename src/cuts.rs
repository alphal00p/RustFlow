//! Typed reverse-unitarity cuts and positive-energy phase-space metadata.
//!
//! The squared denominator does not determine the energy orientation of an
//! on-shell particle. Positive-energy cuts retain their exact momentum routing;
//! algebraic distribution cuts make no claim about a physical energy channel.
//! This module configures existing RustRed facilities instead of introducing
//! another cut reduction or graph-enumeration algorithm.
use crate::{Error, Integral, IntegralFamily, KinematicPoint, Result};
use rustred::sector::{CutConstraint, Pattern, Restrictions};
use rustred::solver::SectorConfig;
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LoopPrescription {
    PlusI0,
    MinusI0,
    Insensitive,
}

/// The physical momentum `sum loops[i] l_i + sum external[j] p_j`.
/// A positive-energy cut uses this oriented momentum, not its negative.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MomentumRouting {
    pub loops: Vec<Rational>,
    pub external: Vec<Rational>,
}

/// An external total momentum explicitly declared future timelike by the caller.
/// Scalar invariants alone cannot distinguish it from the past-directed vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FutureTimelikeChannel {
    pub external: Vec<Rational>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CutDefinition {
    /// `delta_+(q^2-m^2)` with the ordinary Lorentz-invariant phase-space measure.
    /// A positive power n uses (-1)^(n-1)/(n-1)! times the (n-1)th
    /// derivative of the on-shell delta; differentiating with respect to the
    /// squared mass raises its power with coefficient n.
    PositiveEnergy { momentum: MomentumRouting },
    /// Algebraic discontinuity of a denominator, without a positive-energy
    /// condition. Suitable for cut IBPs; not automatically a phase-space volume.
    Distribution,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CutLine {
    pub propagator: usize,
    pub definition: CutDefinition,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CutMetadata {
    lines: BTreeMap<usize, CutDefinition>,
    loop_prescriptions: Vec<LoopPrescription>,
}
impl CutMetadata {
    pub fn new(
        lines: impl IntoIterator<Item = CutLine>,
        loop_prescriptions: Vec<LoopPrescription>,
    ) -> Result<Self> {
        let mut by_slot = BTreeMap::new();
        for line in lines {
            if by_slot.insert(line.propagator, line.definition).is_some() {
                return Err(Error::InvalidInput("duplicate cut propagator slot".into()));
            }
        }
        if by_slot.is_empty() {
            return Err(Error::InvalidInput(
                "a cut family requires at least one cut".into(),
            ));
        }
        Ok(Self {
            lines: by_slot,
            loop_prescriptions,
        })
    }
    pub fn lines(&self) -> &BTreeMap<usize, CutDefinition> {
        &self.lines
    }
    pub fn loop_prescriptions(&self) -> &[LoopPrescription] {
        &self.loop_prescriptions
    }
    pub fn is_cut(&self, slot: usize) -> bool {
        self.lines.contains_key(&slot)
    }
}

/// An ordinary algebraic family together with its integration measure.
/// There is deliberately no `Deref<IntegralFamily>`: passing a cut problem to
/// an uncut evaluator must not silently discard its measure. Validation checks
/// algebraic routing and restrictions; it does not certify physical sheets for
/// general mixed real/virtual families or arbitrary denominator rescalings.
#[derive(Clone, Debug)]
pub struct CutFamily {
    family: IntegralFamily,
    cuts: CutMetadata,
}
impl CutFamily {
    pub fn new(family: IntegralFamily, cuts: CutMetadata) -> Result<Self> {
        let result = Self { family, cuts };
        result.validate()?;
        Ok(result)
    }
    pub fn family(&self) -> &IntegralFamily {
        &self.family
    }
    pub fn cuts(&self) -> &CutMetadata {
        &self.cuts
    }
    pub fn into_parts(self) -> (IntegralFamily, CutMetadata) {
        (self.family, self.cuts)
    }

    pub fn validate(&self) -> Result<()> {
        self.family.validate()?;
        if self.cuts.loop_prescriptions.len() != self.family.loops.len() {
            return Err(Error::InvalidInput(
                "one prescription is required per loop momentum".into(),
            ));
        }
        for (&slot, definition) in &self.cuts.lines {
            if slot >= self.family.physical_propagators {
                return Err(Error::InvalidInput(
                    "cuts must refer to physical propagator slots".into(),
                ));
            }
            if let CutDefinition::PositiveEnergy { momentum } = definition {
                self.validate_mass_shell(slot, momentum)?;
            }
        }
        for slot in 0..self.family.physical_propagators {
            if !self.cuts.is_cut(slot) {
                self.propagator_loop_prescription(slot)?;
            }
        }
        Ok(())
    }

    fn validate_mass_shell(&self, slot: usize, momentum: &MomentumRouting) -> Result<()> {
        let loops = self.family.loops.len();
        let external = self.family.external.len();
        if momentum.loops.len() != loops || momentum.external.len() != external {
            return Err(Error::InvalidInput(
                "cut momentum routing dimensions".into(),
            ));
        }
        if momentum.loops.iter().all(Rational::is_zero) {
            return Err(Error::Unsupported(
                "a positive-energy cut must contain an integration momentum".into(),
            ));
        }
        let mut expected = Vec::new();
        for (i, a) in momentum.loops.iter().enumerate() {
            for (j, b) in momentum.loops.iter().enumerate().skip(i) {
                expected.push(Atom::num(a.clone() * b) * Atom::num(if i == j { 1 } else { 2 }));
            }
        }
        for a in &momentum.loops {
            for b in &momentum.external {
                expected.push(Atom::num(a.clone() * b) * Atom::num(2));
            }
        }
        if expected
            .iter()
            .zip(&self.family.propagators[slot].scalar_products)
            .any(|(a, b)| !(a - b).together().cancel().is_zero())
        {
            return Err(Error::InvalidInput(
                "positive-energy cut routing does not match a normalized quadratic denominator"
                    .into(),
            ));
        }
        Ok(())
    }

    /// Intrinsic squared mass using the retained, normalized physical momentum.
    /// Its reality and threshold support belong to physical-point validation.
    pub fn cut_mass_squared(&self, slot: usize) -> Result<Atom> {
        let Some(CutDefinition::PositiveEnergy { momentum }) = self.cuts.lines.get(&slot) else {
            return Err(Error::Unsupported(
                "distribution cut has no positive-energy mass shell".into(),
            ));
        };
        let mut mass = -self.family.propagators[slot].constant.clone();
        for (i, a) in momentum.external.iter().enumerate() {
            for (j, b) in momentum.external.iter().enumerate() {
                mass += Atom::num(a.clone() * b) * &self.family.external_gram[i][j];
            }
        }
        Ok(mass.together().cancel())
    }

    pub fn at(&self, point: &KinematicPoint) -> Result<Self> {
        Self::new(self.family.at(point), self.cuts.clone())
    }

    /// A cut with nonpositive integer power vanishes by its distributional
    /// definition. This is separate from a RustRed scaleless-sector certificate.
    pub fn is_cut_zero(&self, integral: &Integral) -> Result<bool> {
        self.family.validate_integral(integral)?;
        Ok(self.cuts.lines.keys().any(|&slot| integral.0[slot] <= 0))
    }

    pub fn validate_auxiliary_mask(&self, shifted: &[bool]) -> Result<()> {
        if shifted.len() != self.family.propagators.len() {
            return Err(Error::InvalidInput("auxiliary mass mask dimensions".into()));
        }
        if self.cuts.lines.keys().any(|&slot| shifted[slot]) {
            return Err(Error::InvalidInput(
                "auxiliary mass cannot shift a cut denominator".into(),
            ));
        }
        Ok(())
    }

    /// Required-active constraints for RustRed's zero analysis and symmetry
    /// filtering. An excluded sector is not thereby a scaleless certificate.
    pub fn native_restrictions(&self) -> Result<Restrictions> {
        let arity = self.family.propagators.len();
        let cuts = CutConstraint::try_from_positions(arity, self.cuts.lines.keys().copied())
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        let pattern = Pattern::any(arity).map_err(|e| Error::InvalidInput(e.to_string()))?;
        Restrictions::try_new(cuts, pattern).map_err(|e| Error::InvalidInput(e.to_string()))
    }

    /// Configure the existing native solver for ordinary, possibly raised cuts.
    /// Removed cuts require a separate authenticated source transformation and
    /// are deliberately not inferred here. The caller's resource limits and
    /// runtime arity dispatch remain unchanged.
    pub fn configure_native<const N: usize>(
        &self,
        mut config: SectorConfig<N>,
    ) -> Result<SectorConfig<N>> {
        if N != self.family.propagators.len() {
            return Err(Error::InvalidInput("native cut configuration arity".into()));
        }
        if config.integral_order.is_some() {
            return Err(Error::Unsupported("RustRed programmable ordering does not support cuts; use its native cut-aware ordering".into()));
        }
        if config.removed_deltas.iter().any(|&x| x) {
            return Err(Error::Unsupported(
                "removed cuts require an explicitly prepared native source system".into(),
            ));
        }
        let deltas = std::array::from_fn(|slot| self.cuts.is_cut(slot));
        if config.deltas.iter().any(|&x| x) && config.deltas != deltas {
            return Err(Error::InvalidInput(
                "native cut configuration conflicts with family metadata".into(),
            ));
        }
        config.deltas = deltas;
        Ok(config)
    }

    /// Read the per-loop causal convention on which this propagator depends.
    /// Opposite nonzero conventions in one denominator are inconsistent.
    ///
    /// This precedes denominator normalization: for `D=c(q²-m²)`, the physical
    /// `i0` multiplying D also contains sign(c). General linear denominators
    /// need their own orientation. This method is therefore not an inferred
    /// per-denominator prescription for arbitrary rescaled input. Numerical
    /// mixed-cut adapters must validate and retain those normalizations.
    pub fn propagator_loop_prescription(&self, slot: usize) -> Result<LoopPrescription> {
        let propagator = self
            .family
            .propagators
            .get(slot)
            .ok_or_else(|| Error::InvalidInput("propagator slot out of range".into()))?;
        let loops = self.family.loops.len();
        let mut active = vec![false; loops];
        let mut coordinate = 0;
        for i in 0..loops {
            for j in i..loops {
                if !propagator.scalar_products[coordinate].is_zero() {
                    active[i] = true;
                    active[j] = true;
                }
                coordinate += 1;
            }
        }
        for slot in &mut active {
            for _ in &self.family.external {
                if !propagator.scalar_products[coordinate].is_zero() {
                    *slot = true;
                }
                coordinate += 1;
            }
        }
        let mut result = LoopPrescription::Insensitive;
        for (active, prescription) in active.iter().zip(&self.cuts.loop_prescriptions) {
            if !*active || *prescription == LoopPrescription::Insensitive {
                continue;
            }
            if result != LoopPrescription::Insensitive && result != *prescription {
                return Err(Error::InvalidInput(
                    "denominator mixes opposite loop prescriptions".into(),
                ));
            }
            result = *prescription;
        }
        Ok(result)
    }

    /// Value identity includes cuts and orientation; it cannot collide with an
    /// uncut family or the opposite positive-energy channel.
    pub fn fingerprint(&self) -> Result<String> {
        let algebraic = self.family.convert()?.family;
        let content = format!(
            "cut-family-v1:{}:{}:{:?}",
            algebraic.fingerprint(),
            self.family.physical_propagators,
            self.cuts
        );
        Ok(blake3::hash(content.as_bytes()).to_hex().to_string())
    }
}

/// Internal view used to reuse the ordinary authenticated reduction-cache codec.
/// Its identity contains the measure, and every reduction reattaches that measure
/// explicitly before calling the underlying cut-aware backend.
pub(crate) struct CutBackendView<'a> {
    pub backend: &'a dyn crate::ReductionBackend,
    pub cuts: &'a CutMetadata,
}
impl crate::ReductionBackend for CutBackendView<'_> {
    fn identity(&self) -> String {
        let measure = blake3::hash(format!("cut-measure-v1:{:?}", self.cuts).as_bytes());
        format!("{}:cut-measure-v1:{measure}", self.backend.identity())
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &crate::RunContext,
    ) -> Result<crate::reduction::Reduction> {
        let family = CutFamily::new(family.clone(), self.cuts.clone())?;
        self.backend.reduce_cut(&family, targets, context)
    }
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &crate::RunContext,
    ) -> Result<crate::reduction::Reduction> {
        let family = CutFamily::new(family.clone(), self.cuts.clone())?;
        self.backend
            .reduce_cut_at_epsilon(&family, targets, epsilon, context)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Propagator;
    use rustred::solver::{CoordinateCase, SearchOptions, SectorSolver, SourceSystem};

    fn phase_bubble() -> CutFamily {
        let gram = vec![vec![Atom::num(25)]];
        let family = IntegralFamily {
            name: "two_body_phase_space".into(),
            loops: vec!["k".into()],
            external: vec!["p".into()],
            external_gram: gram.clone(),
            propagators: vec![
                Propagator::quadratic(&[1], &[0], Atom::num(1), &gram).unwrap(),
                Propagator::quadratic(&[-1], &[1], Atom::num(4), &gram).unwrap(),
            ],
            physical_propagators: 2,
            epsilon: symbol!("cuts_test::eps"),
            dimension: 4,
        };
        CutFamily::new(
            family,
            CutMetadata::new(
                [
                    CutLine {
                        propagator: 0,
                        definition: CutDefinition::PositiveEnergy {
                            momentum: MomentumRouting {
                                loops: vec![1.into()],
                                external: vec![0.into()],
                            },
                        },
                    },
                    CutLine {
                        propagator: 1,
                        definition: CutDefinition::PositiveEnergy {
                            momentum: MomentumRouting {
                                loops: vec![(-1).into()],
                                external: vec![1.into()],
                            },
                        },
                    },
                ],
                vec![LoopPrescription::Insensitive],
            )
            .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn cut_constraints_preserve_raised_powers_and_distinguish_measure_identity() {
        let family = phase_bubble();
        assert_eq!(family.cut_mass_squared(0).unwrap(), Atom::num(1));
        assert_eq!(family.cut_mass_squared(1).unwrap(), Atom::num(4));
        assert!(!family.is_cut_zero(&Integral(vec![2, 1])).unwrap());
        assert!(family.is_cut_zero(&Integral(vec![0, 1])).unwrap());
        assert!(family.is_cut_zero(&Integral(vec![1, -3])).unwrap());
        assert!(family.validate_auxiliary_mask(&[false, false]).is_ok());
        assert!(family.validate_auxiliary_mask(&[true, false]).is_err());
        let configuration = family
            .configure_native(SectorConfig::<2>::default())
            .unwrap();
        assert_eq!(configuration.deltas, [true, true]);
        assert_eq!(configuration.removed_deltas, [false, false]);
        let restrictions = family.native_restrictions().unwrap();
        assert!(
            restrictions
                .exclusion(&rustred::sector::Mask::try_new([false, true]).unwrap())
                .unwrap()
                .is_some()
        );
        let mut reversed = family.cuts.clone();
        for definition in reversed.lines.values_mut() {
            if let CutDefinition::PositiveEnergy { momentum } = definition {
                for coefficient in momentum.loops.iter_mut().chain(&mut momentum.external) {
                    *coefficient = -coefficient.clone();
                }
            }
        }
        let reversed = CutFamily::new(family.family.clone(), reversed).unwrap();
        assert_ne!(
            family.fingerprint().unwrap(),
            reversed.fingerprint().unwrap()
        );
        let mut distribution = family.cuts.clone();
        distribution.lines.insert(0, CutDefinition::Distribution);
        let distribution = CutFamily::new(family.family.clone(), distribution).unwrap();
        assert_ne!(
            family.fingerprint().unwrap(),
            distribution.fingerprint().unwrap()
        );
        assert!(distribution.cut_mass_squared(0).is_err());
    }

    #[test]
    fn malformed_cut_metadata_is_rejected_before_reduction() {
        let family = phase_bubble();
        let line = CutLine {
            propagator: 0,
            definition: CutDefinition::Distribution,
        };
        assert!(
            CutMetadata::new([line.clone(), line], vec![LoopPrescription::Insensitive]).is_err()
        );
        let mut cuts = family.cuts.clone();
        cuts.loop_prescriptions.clear();
        assert!(CutFamily::new(family.family.clone(), cuts).is_err());
        let mut cuts = family.cuts.clone();
        cuts.lines.insert(2, CutDefinition::Distribution);
        assert!(CutFamily::new(family.family.clone(), cuts).is_err());
        let mut cuts = family.cuts.clone();
        cuts.lines.insert(
            0,
            CutDefinition::PositiveEnergy {
                momentum: MomentumRouting {
                    loops: vec![2.into()],
                    external: vec![0.into()],
                },
            },
        );
        assert!(CutFamily::new(family.family.clone(), cuts).is_err());
        let config = SectorConfig::<2> {
            deltas: [true, false],
            ..Default::default()
        };
        assert!(family.configure_native(config).is_err());
        let config = SectorConfig::<2> {
            deltas: [true, true],
            removed_deltas: [true, false],
            ..Default::default()
        };
        assert!(family.configure_native(config).is_err());
    }

    #[test]
    fn native_cut_ibps_match_mass_derivatives_of_two_body_phase_space() {
        let cut = phase_bubble();
        let converted = cut
            .family
            .convert_at_epsilon(Some(&Rational::from(0)))
            .unwrap();
        let sources = SourceSystem::<2>::from_family(&converted.family).unwrap();
        let solver = SectorSolver::new(
            &sources,
            [true, true],
            cut.configure_native(SectorConfig::default()).unwrap(),
        )
        .unwrap();
        let mut reduction = crate::reduction::Reduction {
            residuals: vec![Integral(vec![1, 1])],
            ..Default::default()
        };
        for powers in [[2, 1], [1, 2]] {
            let candidate = solver
                .solve_case(
                    CoordinateCase::new(powers.map(Some)).unwrap(),
                    SearchOptions {
                        max_depth: Some(3),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(
                candidate.target,
                rustred::solver::Integral::numeric(powers).unwrap()
            );
            let mut terms = BTreeMap::new();
            for term in candidate.rhs {
                let powers = term
                    .integral
                    .powers()
                    .iter()
                    .map(|p| {
                        assert!(!p.is_symbolic());
                        assert!(p.value() > 0);
                        p.value()
                    })
                    .collect();
                let coefficient = crate::family::substitute(
                    &term.coefficient.to_expression(),
                    &converted.reverse,
                );
                terms.insert(Integral(powers), coefficient);
            }
            reduction.rules.insert(Integral(powers.to_vec()), terms);
        }
        // For s=25,m1²=1,m2²=4, lambda=384. At D=4 the raised
        // cuts are d/dm_i² of phase space, namely (m_i²-s-m_j²)/lambda.
        for (powers, expected) in [
            (vec![2, 1], Atom::num((-7, 96))),
            (vec![1, 2], Atom::num((-11, 192))),
        ] {
            let result = reduction.expand(&Integral(powers)).unwrap();
            assert_eq!(result.len(), 1);
            assert!(
                (&result[&Integral(vec![1, 1])] - expected)
                    .together()
                    .cancel()
                    .is_zero()
            );
        }
    }
}
