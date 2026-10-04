//! Positive-energy two-body phase space using HEPKit normalization and native cut IBPs.
//!
//! This is deliberately restricted to a genuine all-cut one-loop/two-line
//! component. It is not a whole-graph imaginary-part prescription for arbitrary
//! cuts or a replacement for mixed real/virtual AMF boundary construction.
use crate::cuts::{CutDefinition, CutFamily, FutureTimelikeChannel, LoopPrescription};
use crate::family::substitute;
use crate::{
    ComplexFloat as C, Error, FlowOptions, Integral, KinematicPoint, LaurentExpansion, Precision,
    Progress, ReductionBackend, Result, RunContext,
};
use feynkit_kinematics::Kinematics;
use rayon::prelude::*;
use std::collections::BTreeMap;
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;

fn real_rational(value: &Atom, description: &str) -> Result<Rational> {
    let normalized = value.together().cancel();
    if let AtomView::Num(number) = normalized.as_view()
        && let Coefficient::Complex(value) = number.get_coeff_view().to_owned()
        && value.im.is_zero()
    {
        return Ok(value.re);
    }
    Err(Error::Unsupported(format!(
        "{description} must specialize to an exact real rational value"
    )))
}

#[derive(Clone, Debug)]
struct TwoBodyGeometry {
    invariant: Rational,
    kallen: Rational,
    loop_scale: Rational,
    four_dimensional_volume: Atom,
    supported: bool,
}
impl TwoBodyGeometry {
    fn new(family: &CutFamily, channel: &FutureTimelikeChannel) -> Result<Self> {
        let ordinary = family.family();
        if ordinary.loops.len() != 1
            || ordinary.external.len() != 1
            || ordinary.physical_propagators != 2
            || ordinary.propagators.len() != 2
            || family.cuts().lines().len() != 2
        {
            return Err(Error::Unsupported("two-body terminal requires one loop, one external channel, and exactly two physical cut denominators".into()));
        }
        if family.cuts().loop_prescriptions() != [LoopPrescription::Insensitive] {
            return Err(Error::Unsupported(
                "a pure phase-space terminal requires a prescription-insensitive integration loop"
                    .into(),
            ));
        }
        if channel.external.len() != 1 || channel.external[0].is_zero() {
            return Err(Error::InvalidInput(
                "future channel must be a nonzero external momentum of this family".into(),
            ));
        }
        let mut momenta = Vec::new();
        for definition in family.cuts().lines().values() {
            let CutDefinition::PositiveEnergy { momentum } = definition else {
                return Err(Error::Unsupported(
                    "algebraic distribution cuts do not define positive-energy phase space".into(),
                ));
            };
            momenta.push(momentum);
        }
        if !(&momenta[0].loops[0] + &momenta[1].loops[0]).is_zero() {
            return Err(Error::Unsupported(
                "cut momenta do not conserve an external total momentum".into(),
            ));
        }
        let total = &momenta[0].external[0] + &momenta[1].external[0];
        let forward = total == channel.external[0];
        if !forward && total != -channel.external[0].clone() {
            return Err(Error::InvalidInput(
                "cut total momentum is not the specified oriented channel".into(),
            ));
        }
        let gram = real_rational(&ordinary.external_gram[0][0], "channel Gram invariant")?;
        let invariant = channel.external[0].clone() * &channel.external[0] * gram;
        if invariant <= 0 {
            return Err(Error::InvalidInput(
                "the declared future channel is not timelike".into(),
            ));
        }
        let first = real_rational(&family.cut_mass_squared(0)?, "first cut mass squared")?;
        let second = real_rational(&family.cut_mass_squared(1)?, "second cut mass squared")?;
        if first < 0 || second < 0 {
            return Err(Error::Unsupported(
                "physical phase-space masses must be nonnegative".into(),
            ));
        }
        let gap = &(&invariant - &first) - &second;
        let kallen = gap.clone() * &gap - Rational::from(4) * &first * &second;
        // This exact test distinguishes the physical threshold from the
        // pseudothreshold; lambda>0 alone would admit negative-energy solutions.
        let above = gap >= 0 && kallen > 0;
        if forward && gap >= 0 && kallen.is_zero() {
            return Err(Error::Unsupported("the two-body threshold is a singular point for raised cuts; approach it with an endpoint expansion".into()));
        }
        let coefficient = momenta[0].loops[0].clone();
        let loop_scale = if coefficient < 0 {
            -coefficient
        } else {
            coefficient
        };
        let supported = forward && above;
        let four_dimensional_volume = if supported {
            let first_momentum = Atom::var(symbol!("symbolica_amflow::phase_space_first"));
            let second_momentum = Atom::var(symbol!("symbolica_amflow::phase_space_second"));
            let native_error = |error: feynkit_kinematics::SymbolicKinematicsError| {
                Error::InvalidInput(format!("native two-body phase space: {error}"))
            };
            let kinematics = Kinematics::new()
                .with_mass_squared(&first_momentum, Atom::num(first))
                .map_err(native_error)?
                .with_mass_squared(&second_momentum, Atom::num(second))
                .map_err(native_error)?
                .with_scalar_product(
                    &first_momentum,
                    &second_momentum,
                    Atom::num(gap / Rational::from(2)),
                )
                .map_err(native_error)?;
            // Native HEPKit owns dPhi_2/dOmega at D=4. Integrate its rest-frame
            // solid angle exactly; no flux or identical-particle factor belongs
            // in this cut integral. The D-dependent ratio is applied below.
            Atom::num(4)
                * Atom::var(Symbol::PI)
                * kinematics
                    .two_body_phase_space(&first_momentum, &second_momentum)
                    .map_err(native_error)?
        } else {
            Atom::new()
        };
        Ok(Self {
            invariant,
            kallen,
            loop_scale,
            four_dimensional_volume,
            supported,
        })
    }

    fn value(&self, dimension: &Rational, p: Precision) -> Result<C> {
        if !self.supported {
            return Ok(p.zero());
        }
        let half = Rational::from((1, 2));
        let pi = C::new(p.real(1).pi(), p.real(0));
        // Continue the native four-dimensional angular measure to D. This
        // ratio is one at D=4 and retains the dimensional phase-space measure,
        // while the last factor accounts for q = c*k + external routing.
        let factors = [
            p.pow(
                &p.i(2),
                &p.rational(&(Rational::from(8) - Rational::from(2) * dimension)),
            ),
            p.pow(&pi, &p.rational(&((Rational::from(4) - dimension) * &half))),
            p.pow(
                &p.rational(&self.kallen),
                &p.rational(&((dimension - &Rational::from(4)) * &half)),
            ),
            p.pow(
                &p.rational(&self.invariant),
                &p.rational(&((Rational::from(4) - dimension) * &half)),
            ),
            p.pow(
                &p.rational(&self.loop_scale),
                &p.rational(&(-dimension.clone())),
            ),
        ];
        let numerator = factors.iter().fold(
            p.eval(&self.four_dimensional_volume, &Default::default())?,
            |product, factor| p.mul(&product, factor),
        );
        let gamma_four = p.gamma_real(&p.rational(&Rational::from((3, 2))).re)?;
        let gamma_dimension =
            p.gamma_real(&p.rational(&((dimension - &Rational::from(1)) * half)).re)?;
        let value = p.div(&p.mul(&numerator, &gamma_four), &gamma_dimension);
        if !p.finite(&value) {
            return Err(Error::Numerical(
                "nonfinite two-body phase-space terminal".into(),
            ));
        }
        Ok(value)
    }
}

pub struct PreparedTwoBodyPhaseSpace {
    family: CutFamily,
    coefficients: Vec<Atom>,
    nonzero_conditions: Vec<Atom>,
    geometry: TwoBodyGeometry,
}
impl PreparedTwoBodyPhaseSpace {
    pub fn new(
        family: &CutFamily,
        channel: &FutureTimelikeChannel,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        context.cancellation.check()?;
        if targets.is_empty() {
            return Err(Error::InvalidInput("empty phase-space target list".into()));
        }
        let (mut ordinary, cuts) = family.at(point)?.into_parts();
        ordinary.dimension = options.dimension;
        let family = CutFamily::new(ordinary, cuts)?;
        for target in targets {
            family.family().validate_integral(target)?;
        }
        let geometry = TwoBodyGeometry::new(&family, channel)?;
        if !geometry.supported {
            return Ok(Self {
                family,
                coefficients: vec![Atom::new(); targets.len()],
                nonzero_conditions: vec![],
                geometry,
            });
        }
        let reduction = backend.reduce_cut(&family, targets, context)?;
        let unit = Integral(vec![1, 1]);
        let coefficients = targets
            .iter()
            .map(|target| {
                let expanded = reduction.expand(target)?;
                if expanded.keys().any(|integral| *integral != unit) {
                    return Err(Error::IncompleteReduction(
                        "cut reduction did not close onto the unit two-body phase-space volume"
                            .into(),
                    ));
                }
                Ok(expanded.get(&unit).cloned().unwrap_or_default())
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            family,
            coefficients,
            nonzero_conditions: reduction.nonzero_conditions,
            geometry,
        })
    }

    pub fn family(&self) -> &CutFamily {
        &self.family
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<C>> {
        options.validate()?;
        context.cancellation.check()?;
        if options.dimension != self.family.family().dimension {
            return Err(Error::InvalidInput(
                "prepared phase-space dimension differs from evaluation options".into(),
            ));
        }
        if epsilon.is_zero() {
            return Err(Error::InvalidInput(
                "epsilon samples must be nonzero".into(),
            ));
        }
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let rules = BTreeMap::from([(
            Atom::var(self.family.family().epsilon),
            Atom::num(epsilon.clone()),
        )]);
        for condition in &self.nonzero_conditions {
            let condition = substitute(condition, &rules).together().cancel();
            if condition.is_zero() {
                return Err(Error::Numerical(
                    "cut reduction nonzero condition vanishes at the epsilon sample".into(),
                ));
            }
            let value = p.eval(&condition, &Default::default())?;
            if !p.finite(&value) || value == p.zero() {
                return Err(Error::Numerical(
                    "unresolved cut reduction nonzero condition".into(),
                ));
            }
        }
        let dimension =
            Rational::from(self.family.family().dimension) - Rational::from(2) * epsilon;
        let volume = self.geometry.value(&dimension, p)?;
        self.coefficients
            .iter()
            .map(|coefficient| {
                let coefficient = substitute(coefficient, &rules).together().cancel();
                let value = p.mul(&p.eval(&coefficient, &Default::default())?, &volume);
                if !p.finite(&value) {
                    return Err(Error::Numerical(
                        "nonfinite reduced phase-space sample".into(),
                    ));
                }
                Ok(value)
            })
            .collect()
    }

    pub fn evaluate_samples(
        &self,
        samples: &[Rational],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<Vec<C>>> {
        options.validate()?;
        let evaluate = |(index, epsilon): (usize, &Rational)| {
            context.emit(Progress::Sample {
                index,
                total: samples.len(),
            })?;
            self.evaluate(epsilon, options, context)
        };
        if options.workers == 1 {
            samples.iter().enumerate().map(evaluate).collect()
        } else {
            rayon::ThreadPoolBuilder::new()
                .num_threads(options.workers)
                .build()
                .map_err(|e| Error::InvalidInput(e.to_string()))?
                .install(|| samples.par_iter().enumerate().map(evaluate).collect())
        }
    }

    /// All off-threshold two-body volume derivatives are analytic at epsilon=0.
    /// Reuse the library's independent-grid/precision fit, with leading power 0.
    pub fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        crate::engine::fit_samples_refined_leading(
            self.coefficients.len(),
            0,
            last,
            options,
            |samples, refined| self.evaluate_samples(samples, refined, context),
        )
    }
}

#[allow(clippy::too_many_arguments)] // Ordinary solve arguments plus the required oriented physical channel.
pub fn solve_two_body_phase_space(
    family: &CutFamily,
    channel: &FutureTimelikeChannel,
    targets: &[Integral],
    point: &KinematicPoint,
    last: i32,
    options: &FlowOptions,
    backend: &dyn ReductionBackend,
    context: &RunContext,
) -> Result<Vec<LaurentExpansion>> {
    PreparedTwoBodyPhaseSpace::new(family, channel, targets, point, backend, options, context)?
        .solve(last, options, context)
}
