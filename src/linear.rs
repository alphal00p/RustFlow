//! Linear denominators through a quadratic deformation and its exact endpoint.
//!
//! A denominator `D = b + (u.l)(v.p)` is replaced by
//! `D(x) = D + x (u.l)^2`. The branch `u` is normalized to have its first
//! nonzero component equal to one. This normalization changes the auxiliary
//! coordinate, but the denominator at `x=0` is exactly the original one.
//! Positive real `x` preserves the prescription of the finite poles. The
//! additional hard momentum regions have dimension-dependent powers of `x`;
//! the endpoint projector removes those sectors before extracting the finite
//! constant. It never substitutes a small numerical value for the linear limit.
//!
//! The outer system retains symbolic epsilon. Its reference values are obtained
//! automatically from an ordinary quadratic AMF problem at an exact positive
//! rational point inside the endpoint disk of convergence.

use crate::family::{imaginary_parameter, substitute};
use crate::reduction::{ReducedSystem, ReductionBackend};
use crate::{
    ComplexFloat, DifferentialSystem, Error, FlowOptions, Integral, IntegralFamily, KinematicPoint,
    Precision, RecursionMode, Result, RunContext,
};
use std::collections::BTreeMap;
use symbolica::domains::float::Complex;
use symbolica::prelude::*;

/// One physical linear slot and its normalized loop-momentum branch.
#[derive(Clone, Debug)]
pub struct LinearLine {
    pub propagator: usize,
    pub loop_coefficients: Vec<Atom>,
}

/// Exact family deformation. Irreducible numerator slots are unchanged.
#[derive(Clone, Debug)]
pub struct LinearDeformation {
    pub variable: Symbol,
    pub family: IntegralFamily,
    pub lines: Vec<LinearLine>,
    /// Original and deformed denominator-map conditions. These are checked at
    /// the regular reference point, not imposed at the singular endpoint.
    pub nonzero_conditions: Vec<Atom>,
}

fn input_conditions(family: &IntegralFamily) -> Result<Vec<Atom>> {
    let converted = family.convert()?;
    let imaginary = BTreeMap::from([(
        Atom::var(imaginary_parameter()),
        Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
    )]);
    converted
        .family
        .domain()
        .conditions()
        .map(|condition| {
            let value = substitute(
                &substitute(&condition.polynomial().to_expression(), &converted.reverse),
                &imaginary,
            )
            .together()
            .cancel();
            if value.is_zero() {
                Err(Error::InvalidInput(
                    "vanishing linear-family domain condition".into(),
                ))
            } else {
                Ok(value)
            }
        })
        .collect()
}

fn real_rational(a: &Atom) -> bool {
    matches!(a.as_view(), AtomView::Num(n)
        if matches!(n.get_coeff_view().to_owned(),
            symbolica::coefficient::Coefficient::Complex(c) if c.im.is_zero()))
}

impl LinearDeformation {
    /// Return `None` if every physical slot already has a quadratic loop term.
    ///
    /// A linear denominator must factor into one real rational loop branch
    /// times an external vector. General sums of independent loop/external
    /// bilinears and complex loop routings require a different contour proof.
    pub fn new(family: &IntegralFamily, variable: Symbol) -> Result<Option<Self>> {
        let mut nonzero_conditions = input_conditions(family)?;
        if variable == family.epsilon {
            return Err(Error::InvalidInput(
                "linear deformation variable equals epsilon".into(),
            ));
        }
        if family
            .propagators
            .iter()
            .flat_map(|d| std::iter::once(&d.constant).chain(&d.scalar_products))
            .chain(family.external_gram.iter().flatten())
            .any(|a| !a.derivative(variable).is_zero())
        {
            return Err(Error::InvalidInput(
                "linear deformation variable is already in the family".into(),
            ));
        }
        let loops = family.loops.len();
        let external = family.external.len();
        let offset = loops * (loops + 1) / 2;
        let mut deformed = family.clone();
        let mut lines = Vec::new();
        for (index, original) in family.propagators[..family.physical_propagators]
            .iter()
            .enumerate()
        {
            if original.scalar_products[..offset]
                .iter()
                .any(|a| !a.together().cancel().is_zero())
            {
                continue;
            }
            let matrix = &original.scalar_products[offset..];
            let pivot = matrix
                .iter()
                .position(|a| !a.together().cancel().is_zero())
                .ok_or_else(|| {
                    Error::Unsupported(
                        "constant physical denominator is not a linear momentum propagator".into(),
                    )
                })?;
            let row = pivot / external;
            let column = pivot % external;
            let branch = (0..loops)
                .map(|i| {
                    (&matrix[i * external + column] / &matrix[pivot])
                        .together()
                        .cancel()
                })
                .collect::<Vec<_>>();
            for (i, coefficient) in branch.iter().enumerate() {
                if !real_rational(coefficient) {
                    return Err(Error::Unsupported(
                        "linear deformation requires a real rational loop branch".into(),
                    ));
                }
                for a in 0..external {
                    if !(&matrix[i * external + a] - coefficient * &matrix[row * external + a])
                        .together()
                        .cancel()
                        .is_zero()
                    {
                        return Err(Error::Unsupported("linear denominator couples independent loop branches to external vectors".into()));
                    }
                }
            }
            let mut slot = 0;
            for (i, left) in branch.iter().enumerate() {
                for (j, right) in branch.iter().enumerate().skip(i) {
                    deformed.propagators[index].scalar_products[slot] = (Atom::var(variable)
                        * left
                        * right
                        * Atom::num(if i == j { 1 } else { 2 }))
                    .together()
                    .cancel();
                    slot += 1;
                }
            }
            lines.push(LinearLine {
                propagator: index,
                loop_coefficients: branch,
            });
        }
        if lines.is_empty() {
            return Ok(None);
        }
        nonzero_conditions.extend(input_conditions(&deformed)?);
        nonzero_conditions.sort();
        nonzero_conditions.dedup();
        Ok(Some(Self {
            variable,
            family: deformed,
            lines,
            nonzero_conditions,
        }))
    }
}

/// A symbolic-epsilon differential system in the linear deformation parameter.
/// This independent interface performs automatic quadratic reference solves;
/// ordinary [`crate::PreparedFlow`] does not dispatch to it implicitly.
pub struct PreparedLinearFlow {
    pub deformation: LinearDeformation,
    pub reduced: ReducedSystem,
    pub system: DifferentialSystem,
}

/// Numerical linear limits plus the exact regular point used for matching.
#[derive(Clone, Debug)]
pub struct LinearEvaluation {
    pub values: Vec<ComplexFloat>,
    pub reference_point: Rational,
    pub reference_basis_size: usize,
}

impl PreparedLinearFlow {
    pub fn new(
        family: &IntegralFamily,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        context.cancellation.check()?;
        if options.recursion != RecursionMode::Amf {
            return Err(Error::Unsupported(
                "linear quadratic-deformation flow currently uses AMF recursion".into(),
            ));
        }
        if targets.is_empty() {
            return Err(Error::InvalidInput("empty linear-flow target list".into()));
        }
        let mut family = family.at(point);
        family.dimension = options.dimension;
        for target in targets {
            family.validate_integral(target)?;
        }
        let variable = symbol!("symbolica_amflow::linear_x");
        let deformation = LinearDeformation::new(&family, variable)?.ok_or_else(|| {
            Error::InvalidInput("family has no physical linear propagators".into())
        })?;
        let session = crate::reduction::ReductionSession::new(backend);
        let mut reduced = if options.skip_reduction {
            crate::reduction::parameter_differential_system_skip_initial(
                &session,
                &deformation.family,
                targets,
                variable,
                12,
                context,
            )?
        } else {
            crate::reduction::parameter_differential_system(
                &session,
                &deformation.family,
                targets,
                variable,
                12,
                context,
            )?
        };
        reduced
            .nonzero_conditions
            .extend(deformation.nonzero_conditions.iter().cloned());
        reduced.nonzero_conditions.sort();
        reduced.nonzero_conditions.dedup();
        if options.refine_basis {
            crate::refine::refine_basis(&mut reduced, variable, family.epsilon, 32)?;
        }
        let system = DifferentialSystem {
            variable,
            matrix: reduced.matrix.clone(),
        };
        Ok(Self {
            deformation,
            reduced,
            system,
        })
    }

    /// Evaluate from automatically generated quadratic-family boundary data.
    /// Call again with increased precision and order to assess achieved accuracy.
    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        backend: &dyn ReductionBackend,
        context: &RunContext,
    ) -> Result<LinearEvaluation> {
        options.validate()?;
        context.cancellation.check()?;
        if epsilon.is_zero() {
            return Err(Error::InvalidInput(
                "linear-flow epsilon sample must be nonzero".into(),
            ));
        }
        if options.recursion != RecursionMode::Amf
            || options.dimension != self.deformation.family.dimension
        {
            return Err(Error::InvalidInput(
                "linear-flow recursion or dimension differs from preparation".into(),
            ));
        }
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        if self.reduced.basis.is_empty() {
            return Ok(LinearEvaluation {
                values: vec![p.zero(); self.reduced.targets.len()],
                reference_point: Rational::from(1),
                reference_basis_size: 0,
            });
        }
        let eps_symbol = self.deformation.family.epsilon;
        let parameters = ahash::HashMap::from_iter([(Atom::var(eps_symbol), p.rational(epsilon))]);
        let exact_epsilon = BTreeMap::from([(Atom::var(eps_symbol), Atom::num(epsilon.clone()))]);
        let conditions = self
            .reduced
            .nonzero_conditions
            .iter()
            .map(|a| substitute(a, &exact_epsilon).together().cancel())
            .collect::<Vec<_>>();
        if conditions.iter().any(Atom::is_zero) {
            return Err(Error::Reduction(
                "linear-flow nonzero condition vanishes at this epsilon sample".into(),
            ));
        }
        let compiled = self.system.compile(p, &parameters)?;
        let radius = compiled
            .poles
            .iter()
            .map(|z| p.norm(z))
            .filter(|r| *r > p.real(0))
            .fold(
                p.real(1),
                |left, right| if left < right { left } else { right },
            );
        let mut reference_point = Rational::from(1);
        let mut reference = None;
        for _ in 0..options.max_steps {
            context.cancellation.check()?;
            reference_point = &reference_point / &Rational::from(2);
            if p.rational(&reference_point).re >= radius.clone() / 8 {
                continue;
            }
            let point = KinematicPoint(BTreeMap::from([(
                Atom::var(self.system.variable),
                Atom::num(reference_point.clone()),
            )]));
            if conditions
                .iter()
                .any(|a| point.apply(a).together().cancel().is_zero())
            {
                continue;
            }
            // Check every remaining scalar condition numerically before any
            // reference IBP solve, including external Gram/input domains.
            for condition in &conditions {
                if p.eval(&point.apply(condition), &parameters)? == p.zero() {
                    return Err(Error::Reduction(
                        "vanishing linear reference-point condition".into(),
                    ));
                }
            }
            let family = self.deformation.family.at(&point);
            family.validate()?;
            reference = Some(family);
            break;
        }
        let family = reference
            .ok_or_else(|| Error::Limit("linear reference-point search exhausted".into()))?;
        let endpoint = self
            .system
            .frobenius(p, &parameters, options.series_order)?;
        crate::engine::ensure_generic_indicial(&endpoint, &parameters)?;
        // Keeping epsilon symbolic in the inner reference solve also supports
        // rational regulator values for which sampled exponents are ambiguous.
        let flow = crate::PreparedFlow::new(
            &family,
            &self.reduced.basis,
            &KinematicPoint::default(),
            backend,
            options,
            context,
        )?;
        let boundary = crate::recursive::RecursiveBoundary::new(backend, options, context);
        let values = flow.evaluate(epsilon, options, &boundary, context)?;
        let constants =
            endpoint.match_values(&p.rational(&reference_point), &values, &parameters)?;
        let values = self
            .reduced
            .targets
            .iter()
            .map(|target| {
                let weights = self
                    .reduced
                    .basis
                    .iter()
                    .map(|integral| target.get(integral).cloned().unwrap_or_default())
                    .collect::<Vec<_>>();
                crate::engine::project_limit(
                    &endpoint,
                    &constants,
                    &weights,
                    eps_symbol,
                    self.system.variable,
                    &parameters,
                )
            })
            .collect::<Result<_>>()?;
        Ok(LinearEvaluation {
            values,
            reference_point,
            reference_basis_size: flow.reduced.basis.len(),
        })
    }
}
