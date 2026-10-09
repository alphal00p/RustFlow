//! Shared numerical connection transport for ordinary and distribution-weighted AMF.
//! Physical admission and source replay belong to the calling measure owner.
use super::*;
use crate::frobenius::FrobeniusBasis;

#[derive(Default)]
pub(crate) struct ConnectionTransport {
    pub(super) preparations: std::sync::Arc<super::frobenius_cache::FrobeniusPreparations>,
}

pub(crate) struct ConnectionRequest<'a> {
    pub system: &'a DifferentialSystem,
    pub reduced: &'a ReducedSystem,
    pub epsilon: Symbol,
    pub loops: usize,
    pub positive_mass_contour: bool,
    /// For occupied domains proved analytic only for Re(eta)>=0, every straight
    /// transport segment must stay in that convex domain.
    pub restrict_to_right_half_plane: bool,
    pub start_scale: u32,
}

impl ConnectionTransport {
    pub(crate) fn evaluate(
        &self,
        request: ConnectionRequest<'_>,
        eps: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        lift_exponents: impl Fn(&mut FrobeniusBasis) -> Result<()>,
        boundary: impl FnOnce(&FrobeniusBasis, Precision) -> Result<Vec<ComplexFloat>>,
    ) -> Result<Vec<ComplexFloat>> {
        options.validate()?;
        if request.start_scale < 4 {
            return Err(Error::InvalidInput(
                "auxiliary start scale must be at least four".into(),
            ));
        }
        context.cancellation.check()?;
        let digits = options
            .digits
            .checked_add(options.guard_digits)
            .ok_or_else(|| Error::InvalidInput("precision overflow".into()))?;
        let p = Precision::decimal(digits)?;
        let parameters = ahash::HashMap::from_iter([(Atom::var(request.epsilon), p.rational(eps))]);
        // Preserve original rational denominator domains before cancellation,
        // then specialize epsilon exactly. Floating substitution can miss an
        // identically zero eta-dependent condition at a rational sample.
        let epsilon_rule = BTreeMap::from([(Atom::var(request.epsilon), Atom::num(eps.clone()))]);
        let guards = crate::physical_conditions::canonical_conditions(
            &request.reduced.nonzero_conditions,
            &std::collections::BTreeSet::from([
                request.system.variable,
                request.epsilon,
                crate::family::imaginary_parameter(),
            ]),
        )?
        .iter()
        .map(|condition| {
            let specialized = crate::family::substitute(condition, &epsilon_rule)
                .together()
                .cancel();
            if specialized.is_zero() {
                return Err(Error::Reduction(format!(
                    "a reduction nonzero condition vanishes identically at epsilon={eps}"
                )));
            }
            Ok(specialized)
        })
        .collect::<Result<Vec<_>>>()?;
        if request.reduced.basis.is_empty() {
            return Ok(vec![p.zero(); request.reduced.targets.len()]);
        }
        let stage = |name: &str| {
            context.emit(Progress::Stage {
                name: format!(
                    "{name} ({} basis integrals, {} loops)",
                    request.reduced.basis.len(),
                    request.loops
                ),
            })
        };
        stage("compiling the auxiliary-mass differential equation")?;
        let mut compiled = request.system.compile(p, &parameters)?;
        compiled.exclude_polynomials(&guards)?;
        stage("constructing the expansion at auxiliary-mass infinity")?;
        let mut infinity = self
            .preparations
            .get(&request.system, true, context)?
            .evaluate(p, &parameters, options.series_order, context)?;
        lift_exponents(&mut infinity)?;
        ensure_generic_indicial(&infinity, &parameters)?;
        stage("generating and matching native asymptotic boundary regions")?;
        let constants = boundary(&infinity, p)?;
        let maximum = compiled
            .poles
            .iter()
            .map(|v| p.norm(v))
            .fold(p.real(1), |a, b| if a > b { a } else { b });
        let direction = if options.prescription == Prescription::PlusI0 {
            -1
        } else {
            1
        };
        let positive_mass_contour = request.positive_mass_contour;
        let start_scale = i64::from(request.start_scale);
        let start = if positive_mass_contour {
            ComplexFloat::new(maximum * start_scale, p.real(0))
        } else {
            ComplexFloat::new(p.real(0), maximum * (start_scale * direction))
        };
        let at_start = infinity.evaluate(&p.div(&p.i(1), &start), &parameters)?;
        let initial = at_start
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&constants)
                    .fold(p.zero(), |s, (v, c)| p.add(&s, &p.mul(v, c)))
            })
            .collect();
        let minimum = compiled
            .poles
            .iter()
            .map(|v| p.norm(v))
            .filter(|v| *v > p.real(0))
            .fold(p.real(1), |a, b| if a < b { a } else { b });
        let end = if positive_mass_contour {
            ComplexFloat::new(minimum / 8, p.real(0))
        } else {
            ComplexFloat::new(p.real(0), (minimum / 8) * direction)
        };
        let path = compiled.plan_path(&start, &end, -direction)?;
        if request.restrict_to_right_half_plane
            && path
                .iter()
                .chain([&start, &end])
                .any(|point| point.re < p.real(0))
        {
            return Err(Error::Unsupported(
                "planned auxiliary path leaves the certified Re(eta)>=0 occupied continuation domain".into(),
            ));
        }
        stage("transporting the auxiliary-mass solution")?;
        let transported = compiled.transport(
            &BoundaryData {
                point: start,
                values: initial,
            },
            &path,
            options,
            context,
        )?;
        stage("constructing the physical endpoint expansion")?;
        let mut endpoint = self
            .preparations
            .get(&request.system, false, context)?
            .evaluate(p, &parameters, options.series_order, context)?;
        lift_exponents(&mut endpoint)?;
        ensure_generic_indicial(&endpoint, &parameters)?;
        let constants = endpoint.match_values(&end, &transported.values, &parameters)?;
        stage("extracting the dimensionally regulated physical limit")?;
        // Reduction coefficients are multiplied into endpoint series before
        // selecting the physical constant; poles in these coefficients matter.
        request
            .reduced
            .targets
            .iter()
            .map(|target| {
                let weights = request
                    .reduced
                    .basis
                    .iter()
                    .map(|i| target.get(i).cloned().unwrap_or_default())
                    .collect::<Vec<_>>();
                project_limit(
                    &endpoint,
                    &constants,
                    &weights,
                    request.epsilon,
                    request.system.variable,
                    &parameters,
                )
            })
            .collect()
    }
}
