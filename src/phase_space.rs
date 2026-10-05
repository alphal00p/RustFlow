//! Positive-energy phase-space terminals using HEPKit normalization and native cut IBPs.
//!
//! Supported terminals are genuine all-cut two-body components and massless
//! N-body components. This is not a whole-graph imaginary-part prescription for arbitrary
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

fn native_two_body_volume(
    first: &Rational,
    second: &Rational,
    invariant: &Rational,
) -> Result<Atom> {
    let first_momentum = Atom::var(symbol!("symbolica_amflow::phase_space_first"));
    let second_momentum = Atom::var(symbol!("symbolica_amflow::phase_space_second"));
    let native_error = |error: feynkit_kinematics::SymbolicKinematicsError| {
        Error::InvalidInput(format!("native two-body phase space: {error}"))
    };
    let kinematics = Kinematics::new()
        .with_mass_squared(&first_momentum, Atom::num(first.clone()))
        .map_err(native_error)?
        .with_mass_squared(&second_momentum, Atom::num(second.clone()))
        .map_err(native_error)?
        .with_scalar_product(
            &first_momentum,
            &second_momentum,
            Atom::num((invariant - first - second) / Rational::from(2)),
        )
        .map_err(native_error)?;
    // Native HEPKit owns dPhi_2/dOmega at D=4. Integrate its rest-frame
    // solid angle exactly; no flux or identical-particle factor belongs
    // in this cut integral. The D-dependent ratio is applied below.
    Ok(Atom::num(4)
        * Atom::var(Symbol::PI)
        * kinematics
            .two_body_phase_space(&first_momentum, &second_momentum)
            .map_err(native_error)?)
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
            native_two_body_volume(&first, &second, &invariant)?
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

/// All L+1 cut momenta form a complete massless final state. External shifts
/// are affine translations; only the exact loop-routing determinant changes
/// the integration measure. No graph enumeration or tensor reduction occurs.
#[derive(Clone, Debug)]
struct MasslessGeometry {
    particles: usize,
    invariant: Rational,
    determinant: Rational,
    supported: bool,
    two_body_base: TwoBodyGeometry,
}
impl MasslessGeometry {
    fn new(family: &CutFamily, channel: &FutureTimelikeChannel) -> Result<Self> {
        let ordinary = family.family();
        let loops = ordinary.loops.len();
        let particles = loops
            .checked_add(1)
            .ok_or_else(|| Error::Limit("phase-space multiplicity overflow".into()))?;
        if loops == 0
            || ordinary.physical_propagators != particles
            || family.cuts().lines().len() != particles
        {
            return Err(Error::Unsupported(
                "massless N-body terminal requires exactly L+1 physical denominators, all cut"
                    .into(),
            ));
        }
        if family
            .cuts()
            .loop_prescriptions()
            .iter()
            .any(|p| *p != LoopPrescription::Insensitive)
        {
            return Err(Error::Unsupported(
                "a pure phase-space terminal requires prescription-insensitive integration loops"
                    .into(),
            ));
        }
        if channel.external.len() != ordinary.external.len()
            || channel.external.iter().all(Rational::is_zero)
        {
            return Err(Error::InvalidInput(
                "future channel must be a nonzero external momentum of this family".into(),
            ));
        }
        let gram = ordinary
            .external_gram
            .iter()
            .map(|row| {
                row.iter()
                    .map(|value| real_rational(value, "physical channel Gram matrix"))
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut invariant = Rational::zero();
        for (i, a) in channel.external.iter().enumerate() {
            for (j, b) in channel.external.iter().enumerate() {
                invariant += a * b * &gram[i][j];
            }
        }
        if invariant <= 0 {
            return Err(Error::InvalidInput(
                "the declared future channel is not timelike".into(),
            ));
        }
        let mut momenta = Vec::with_capacity(particles);
        for (&slot, definition) in family.cuts().lines() {
            let CutDefinition::PositiveEnergy { momentum } = definition else {
                return Err(Error::Unsupported(
                    "algebraic distribution cuts do not define positive-energy phase space".into(),
                ));
            };
            if !family.cut_mass_squared(slot)?.is_zero() {
                return Err(Error::Unsupported(
                    "massless N-body terminal requires every cut mass squared to vanish exactly"
                        .into(),
                ));
            }
            momenta.push(momentum);
        }
        for index in 0..loops {
            if !momenta
                .iter()
                .fold(Rational::zero(), |sum, q| sum + &q.loops[index])
                .is_zero()
            {
                return Err(Error::Unsupported(
                    "cut momenta do not conserve an external total momentum".into(),
                ));
            }
        }
        let total = (0..ordinary.external.len())
            .map(|index| {
                momenta
                    .iter()
                    .fold(Rational::zero(), |sum, q| sum + &q.external[index])
            })
            .collect::<Vec<_>>();
        let forward = total == channel.external;
        if !forward
            && total
                != channel
                    .external
                    .iter()
                    .map(|x| -x.clone())
                    .collect::<Vec<_>>()
        {
            return Err(Error::InvalidInput(
                "cut total momentum is not the specified oriented channel".into(),
            ));
        }
        // All omitted-row minors agree up to sign because sum(q_i)=P. The
        // first L cut momenta therefore suffice to certify a complete routing.
        let routing = momenta[..loops].iter().map(|q| q.loops.clone()).collect();
        let determinant = symbolica::tensors::matrix::Matrix::from_nested_vec(routing, Q)
            .map_err(|e| Error::InvalidInput(format!("phase-space routing matrix: {e}")))?
            .det()
            .map_err(|e| Error::InvalidInput(format!("phase-space routing determinant: {e}")))?
            .abs();
        if determinant.is_zero() {
            return Err(Error::Unsupported(
                "cut momenta do not span the independent integration loops".into(),
            ));
        }
        let two_body_base = TwoBodyGeometry {
            invariant: Rational::one(),
            kallen: Rational::one(),
            loop_scale: Rational::one(),
            four_dimensional_volume: native_two_body_volume(
                &Rational::zero(),
                &Rational::zero(),
                &Rational::one(),
            )?,
            supported: true,
        };
        Ok(Self {
            particles,
            invariant,
            determinant,
            supported: forward,
            two_body_base,
        })
    }

    fn value(&self, dimension: &Rational, p: Precision) -> Result<C> {
        if !self.supported {
            return Ok(p.zero());
        }
        let n = i64::try_from(self.particles).map_err(|_| {
            Error::Limit("phase-space multiplicity exceeds signed arithmetic".into())
        })?;
        let a = (dimension - &Rational::from(2)) / Rational::from(2);
        let gamma = |argument: Rational| p.gamma_real(&p.rational(&argument).re);
        let pi = C::new(p.real(1).pi(), p.real(0));
        let base = self.two_body_base.value(dimension, p)?;
        // For Re(a)>0 this is repeated dt/(2pi) phase-space factorization,
        // with beta((n-2)a,2a) at each step. The gamma quotient supplies its
        // meromorphic dimensional continuation away from exceptional samples.
        // For N=2 the gamma quotient is exactly one, including removable
        // exceptional points such as D=2. Do not separately evaluate its
        // divergent gamma factors instead of the finite native two-body value.
        let mut value = base.clone();
        if n > 2 {
            value = p.div(&p.powi(&base, n - 1), &p.powi(&p.scale(&pi, 2, 1), n - 2));
            value = p.mul(&value, &gamma(a.clone())?);
            value = p.mul(&value, &p.powi(&gamma(Rational::from(2) * &a)?, n - 1));
            value = p.div(&value, &gamma(Rational::from(n - 1) * &a)?);
            value = p.div(&value, &gamma(Rational::from(n) * &a)?);
        }
        value = p.mul(
            &value,
            &p.pow(
                &p.rational(&self.invariant),
                &p.rational(&(Rational::from(n - 1) * a - Rational::one())),
            ),
        );
        value = p.mul(
            &value,
            &p.pow(
                &p.rational(&self.determinant),
                &p.rational(&(-dimension.clone())),
            ),
        );
        if !p.finite(&value) {
            return Err(Error::Numerical(
                "nonfinite massless phase-space terminal".into(),
            ));
        }
        Ok(value)
    }
}

#[derive(Clone, Copy)]
enum TerminalKind {
    TwoBody,
    Massless,
}
enum Geometry {
    TwoBody(TwoBodyGeometry),
    Massless(MasslessGeometry),
}
impl Geometry {
    fn supported(&self) -> bool {
        match self {
            Self::TwoBody(g) => g.supported,
            Self::Massless(g) => g.supported,
        }
    }
    fn value(&self, dimension: &Rational, p: Precision) -> Result<C> {
        match self {
            Self::TwoBody(g) => g.value(dimension, p),
            Self::Massless(g) => g.value(dimension, p),
        }
    }
}

/// Coefficient poles supplement the generic -2L dimensional pole allowance.
/// Exact specialized rational weights are required; no sampled pole inference.
fn reduction_pole_order(coefficients: &[Atom], epsilon: Symbol) -> Result<i32> {
    let mut poles = 0_i64;
    for coefficient in coefficients {
        let coefficient: RationalPolynomial<IntegerRing, u16> = coefficient
            .try_to_rational_polynomial(&Q, &Z, None)
            .map_err(|e| {
                Error::Unsupported(format!(
                    "rational phase-space reduction weight required: {e}"
                ))
            })?;
        if coefficient.numerator.is_zero() {
            continue;
        }
        if coefficient
            .numerator
            .variables()
            .iter()
            .any(|v| *v != PolyVariable::Symbol(epsilon))
        {
            return Err(Error::Unsupported(
                "phase-space reduction weights must specialize to rational functions of epsilon"
                    .into(),
            ));
        }
        let order = |poly: &MultivariatePolynomial<IntegerRing, u16>| {
            poly.variables()
                .iter()
                .position(|v| *v == PolyVariable::Symbol(epsilon))
                .map_or(0, |index| i64::from(poly.degree_bounds(index).0))
        };
        poles = poles.max(order(&coefficient.denominator) - order(&coefficient.numerator));
    }
    i32::try_from(poles)
        .map_err(|_| Error::Limit("phase-space reduction pole order overflow".into()))
}

struct PreparedPhaseSpaceTerminal {
    family: CutFamily,
    coefficients: Vec<Atom>,
    nonzero_conditions: Vec<Atom>,
    geometry: Geometry,
    leading: i32,
}
impl PreparedPhaseSpaceTerminal {
    #[allow(clippy::too_many_arguments)]
    fn new(
        kind: TerminalKind,
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
        let geometry = match kind {
            TerminalKind::TwoBody => Geometry::TwoBody(TwoBodyGeometry::new(&family, channel)?),
            TerminalKind::Massless => Geometry::Massless(MasslessGeometry::new(&family, channel)?),
        };
        let mut coefficients = vec![Atom::new(); targets.len()];
        let mut nonzero_conditions = Vec::new();
        let unit = Integral(
            (0..family.family().propagators.len())
                .map(|slot| i16::from(family.cuts().is_cut(slot)))
                .collect(),
        );
        if geometry.supported() {
            let mut pending = Vec::new();
            for (index, target) in targets.iter().enumerate() {
                if family.is_cut_zero(target)? {
                    continue;
                }
                if matches!(kind, TerminalKind::Massless)
                    && target
                        .0
                        .iter()
                        .enumerate()
                        .any(|(slot, &power)| !family.cuts().is_cut(slot) && power > 0)
                {
                    return Err(Error::Unsupported("massless phase-space terminal does not admit additional uncut denominators".into()));
                }
                if matches!(kind, TerminalKind::Massless) && *target == unit {
                    coefficients[index] = Atom::one();
                } else {
                    pending.push((index, target.clone()));
                }
            }
            if !pending.is_empty() {
                let reduction = backend.reduce_cut(
                    &family,
                    &pending
                        .iter()
                        .map(|(_, target)| target.clone())
                        .collect::<Vec<_>>(),
                    context,
                )?;
                for (index, target) in pending {
                    let expanded = reduction.expand(&target)?;
                    if expanded.keys().any(|integral| *integral != unit) {
                        return Err(Error::IncompleteReduction(
                            "cut reduction did not close onto the proved unit phase-space volume"
                                .into(),
                        ));
                    }
                    coefficients[index] = expanded.get(&unit).cloned().unwrap_or_default();
                }
                nonzero_conditions = reduction.nonzero_conditions;
            }
        }
        let leading = match kind {
            TerminalKind::TwoBody => 0,
            TerminalKind::Massless => {
                let loops = i32::try_from(family.family().loops.len())
                    .map_err(|_| Error::Limit("phase-space loop count overflow".into()))?;
                let extra_poles = reduction_pole_order(&coefficients, family.family().epsilon)?;
                loops
                    .checked_mul(-2)
                    .and_then(|bound| bound.checked_sub(extra_poles))
                    .ok_or_else(|| Error::Limit("phase-space Laurent bound overflow".into()))?
            }
        };
        Ok(Self {
            family,
            coefficients,
            nonzero_conditions,
            geometry,
            leading,
        })
    }
    fn evaluate(
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

    fn evaluate_samples(
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

    // Both terminal classes use the same independent-grid/precision fit.
    fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        crate::engine::fit_samples_refined_leading(
            self.coefficients.len(),
            self.leading,
            last,
            options,
            |samples, refined| self.evaluate_samples(samples, refined, context),
        )
    }
}

/// Genuine all-cut two-body phase space, including unequal nonnegative masses.
/// The established API and analytic-at-epsilon-zero fitting convention are kept.
pub struct PreparedTwoBodyPhaseSpace {
    terminal: PreparedPhaseSpaceTerminal,
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
        Ok(Self {
            terminal: PreparedPhaseSpaceTerminal::new(
                TerminalKind::TwoBody,
                family,
                channel,
                targets,
                point,
                backend,
                options,
                context,
            )?,
        })
    }
    pub fn family(&self) -> &CutFamily {
        &self.terminal.family
    }
    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<C>> {
        self.terminal.evaluate(epsilon, options, context)
    }
    pub fn evaluate_samples(
        &self,
        samples: &[Rational],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<Vec<C>>> {
        self.terminal.evaluate_samples(samples, options, context)
    }
    pub fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        self.terminal.solve(last, options, context)
    }
}

/// Massless N-body positive-energy phase space, N=L+1>=2. Cut routing and the
/// future channel are exact. Raised cuts and polynomial numerators are admitted
/// only when the caller's cut-aware backend reduces them onto the unit volume.
/// Mixed cuts, additional uncut denominators and massive N>2 states are excluded.
pub struct PreparedMasslessPhaseSpace {
    terminal: PreparedPhaseSpaceTerminal,
}
impl PreparedMasslessPhaseSpace {
    pub fn new(
        family: &CutFamily,
        channel: &FutureTimelikeChannel,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        Ok(Self {
            terminal: PreparedPhaseSpaceTerminal::new(
                TerminalKind::Massless,
                family,
                channel,
                targets,
                point,
                backend,
                options,
                context,
            )?,
        })
    }
    pub fn family(&self) -> &CutFamily {
        &self.terminal.family
    }
    /// Conservative generic loop pole bound, supplemented by exact poles of
    /// reduction weights. This is not inferred from a unit-volume finite value.
    pub fn leading_power(&self) -> i32 {
        self.terminal.leading
    }
    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<C>> {
        self.terminal.evaluate(epsilon, options, context)
    }
    pub fn evaluate_samples(
        &self,
        samples: &[Rational],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<Vec<C>>> {
        self.terminal.evaluate_samples(samples, options, context)
    }
    pub fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        self.terminal.solve(last, options, context)
    }
}

#[allow(clippy::too_many_arguments)]
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

#[allow(clippy::too_many_arguments)]
pub fn solve_massless_phase_space(
    family: &CutFamily,
    channel: &FutureTimelikeChannel,
    targets: &[Integral],
    point: &KinematicPoint,
    last: i32,
    options: &FlowOptions,
    backend: &dyn ReductionBackend,
    context: &RunContext,
) -> Result<Vec<LaurentExpansion>> {
    PreparedMasslessPhaseSpace::new(family, channel, targets, point, backend, options, context)?
        .solve(last, options, context)
}
