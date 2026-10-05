//! Physical differential equations and AMF-generated seeds for one exact family.
use crate::{
    kinematic_derivative::KinematicDerivative,
    kinematics::KinematicSystem,
    reduction::LinearCombination,
    transport_cache::{BoundaryAccuracy, CachedBoundary, CachedPoint, EpsilonRange, PointKind},
    *,
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

/// A common derivative-closed physical basis, its exact target reductions and
/// the AMF family that can generate fresh numerical starting values.
///
/// Fixed parameters must be specialized before construction. Every remaining
/// physical scalar variable needs a derivative, including mass parameters.
/// Boundary values use the library's ordinary `d^D l/(i*pi^(D/2))` convention.
pub struct PreparedPhysicalFamily {
    family: IntegralFamily,
    basis: Vec<Integral>,
    targets: Vec<LinearCombination>,
    conditions: Vec<Atom>,
    variables: BTreeSet<Symbol>,
    prescription: Prescription,
    transformations: Vec<reduction::BasisTransformation>,
    refinement: Option<refine::RefinementReport>,
    flow: RustFlow,
    sheared: Option<EpsilonShearedFlow>,
}

impl PreparedPhysicalFamily {
    pub fn new(
        family: &IntegralFamily,
        targets: &[Integral],
        variables: &[Symbol],
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        branch_domain: &str,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        if targets.is_empty() || variables.is_empty() {
            return Err(Error::InvalidInput(
                "physical preparation requires targets and variables".into(),
            ));
        }
        let mut family = family.clone();
        family.dimension = options.dimension;
        family.validate()?;
        let variable_set = variables.iter().copied().collect::<BTreeSet<_>>();
        if variable_set.len() != variables.len() || variable_set.contains(&family.epsilon) {
            return Err(Error::InvalidInput(
                "physical variables must be distinct and exclude epsilon".into(),
            ));
        }
        // A parameter can cancel out of the final differential matrix while
        // still changing an integral normalization. Check the family itself.
        let mut allowed = variable_set.clone();
        allowed.insert(family.epsilon);
        let mut symbols = BTreeSet::new();
        for atom in family.external_gram.iter().flatten().chain(
            family
                .propagators
                .iter()
                .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
        ) {
            crate::family::scalar_symbols(atom.as_view(), &mut symbols)?;
        }
        symbols.remove(&Atom::var(crate::family::imaginary_parameter()));
        if symbols
            .iter()
            .any(|s| !matches!(s.as_view(), AtomView::Var(v) if allowed.contains(&v.get_symbol())))
        {
            return Err(Error::InvalidInput("every remaining family parameter must have a physical derivative or be specialized before preparation".into()));
        }
        let derivatives = variables
            .iter()
            .map(|&s| Ok((s, KinematicDerivative::new(&family, s)?)))
            .collect::<Result<BTreeMap<_, _>>>()?;
        // The shared closure builder sees every derivative target separately.
        // Its summed matrix is unused: each coordinate matrix below is formed
        // from the same returned reductions in the same ordered basis.
        let derivative = |integral: &Integral| {
            let mut terms = Vec::new();
            for derivative in derivatives.values() {
                terms.extend(derivative.integral(integral)?);
            }
            Ok(terms)
        };
        let mut closed = if options.skip_reduction {
            crate::reduction::build_retained_system(
                backend, &family, targets, 32, context, derivative,
            )?
        } else {
            crate::reduction::build_differential_system(
                backend, &family, targets, 32, context, derivative,
            )?
        };
        // Native rational reduction represents i as a formal coefficient.
        // Decode it before forming the physical connection and its domain.
        let imaginary = BTreeMap::from([(
            Atom::var(crate::family::imaginary_parameter()),
            Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
        )]);
        for coefficient in closed
            .candidates
            .values_mut()
            .flat_map(|terms| terms.values_mut())
            .chain(
                closed
                    .targets
                    .iter_mut()
                    .flat_map(|terms| terms.values_mut()),
            )
        {
            *coefficient = crate::family::substitute(coefficient, &imaginary)
                .together()
                .cancel();
        }
        for condition in &mut closed.nonzero_conditions {
            *condition = crate::family::substitute(condition, &imaginary);
        }
        if closed.basis.is_empty() {
            return Err(Error::Unsupported(
                "all targets vanish; use solve_integrals for the exact zero family".into(),
            ));
        }
        let mut matrices = BTreeMap::new();
        let mut derivative_conditions = Vec::new();
        for (&variable, derivative) in &derivatives {
            let mut matrix = vec![vec![Atom::new(); closed.basis.len()]; closed.basis.len()];
            for (row, integral) in closed.basis.iter().enumerate() {
                for (integral, coefficient) in derivative.integral(integral)? {
                    let expansion = closed.candidates.get(&integral).ok_or_else(|| {
                        Error::IncompleteReduction(
                            "physical derivative target absent after common closure".into(),
                        )
                    })?;
                    for (master, weight) in expansion {
                        let column = closed.basis.binary_search(master).map_err(|_| {
                            Error::IncompleteReduction(
                                "physical derivative leaves the common basis".into(),
                            )
                        })?;
                        matrix[row][column] += &coefficient * weight;
                    }
                }
                for entry in &mut matrix[row] {
                    *entry = entry.together().cancel();
                }
            }
            matrices.insert(variable, matrix);
            derivative_conditions.extend_from_slice(derivative.nonzero_conditions());
        }
        let refinement = if options.refine_basis {
            let first_variable = *derivatives.keys().next().unwrap();
            closed.matrix = matrices[&first_variable].clone();
            let report =
                crate::refine::refine_basis(&mut closed, first_variable, family.epsilon, 32)?;
            // A basis change depends on every coordinate. Apply the same
            // ordered transformations with each coordinate's own derivative.
            for (&variable, matrix) in &mut matrices {
                if variable == first_variable {
                    *matrix = closed.matrix.clone();
                } else {
                    for transformation in &closed.transformations {
                        *matrix = DifferentialSystem {
                            variable,
                            matrix: std::mem::take(matrix),
                        }
                        .change_basis(&transformation.matrix)?
                        .matrix;
                    }
                }
            }
            Some(report)
        } else {
            None
        };
        let mut conditions = closed.nonzero_conditions;
        conditions.extend(derivative_conditions);
        let conditions = crate::physical_conditions::canonical_conditions(&conditions, &allowed)?;
        let system = KinematicSystem {
            epsilon: family.epsilon,
            derivatives: matrices,
        };
        // The native family fingerprint anonymizes coefficient names. Keep
        // all original names, slot order and guards in physical master labels.
        let fingerprint = physical_family_fingerprint(&family, &variable_set, &conditions)?;
        let head = symbol!(&format!("symbolica_amflow::physical_family_{fingerprint}"));
        let labels = closed
            .basis
            .iter()
            .map(|integral| {
                head.call_args(integral.0.iter().map(|&power| Atom::num(i64::from(power))))
            })
            .collect::<Vec<_>>();
        let flow = RustFlow::with_conditions(
            system,
            &labels,
            &Atom::one(),
            options.prescription,
            branch_domain,
            &conditions,
        )?;
        Ok(Self {
            family,
            basis: closed.basis,
            targets: closed.targets,
            conditions,
            variables: derivatives.into_keys().collect(),
            prescription: options.prescription,
            transformations: closed.transformations,
            refinement,
            flow,
            sheared: None,
        })
    }
    /// Rescale the common integral basis by exact epsilon powers when one
    /// diagonal transformation regularizes every physical partial. Original
    /// integral labels/reductions remain available; `flow()` exposes the
    /// explicitly rescaled cache basis, and seed/projection methods map it.
    pub fn with_epsilon_shearing(mut self, context: &RunContext) -> Result<Self> {
        if self.sheared.is_none() {
            self.sheared = Some(self.flow.regularize_epsilon(context)?);
        }
        Ok(self)
    }
    pub fn epsilon_shearing(&self) -> Option<&EpsilonShearing> {
        self.sheared.as_ref().map(EpsilonShearedFlow::shearing)
    }
    pub fn family(&self) -> &IntegralFamily {
        &self.family
    }
    /// Original integral basis. With epsilon shearing enabled, cached columns
    /// are epsilon^-weight times these integrals; see `epsilon_shearing()`.
    pub fn basis(&self) -> &[Integral] {
        &self.basis
    }
    /// Exact ordered basis changes, with previous = matrix * current.
    pub fn basis_transformations(&self) -> &[reduction::BasisTransformation] {
        &self.transformations
    }
    pub fn basis_refinement(&self) -> Option<&refine::RefinementReport> {
        self.refinement.as_ref()
    }
    pub fn target_reductions(&self) -> &[LinearCombination] {
        &self.targets
    }
    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.conditions
    }
    pub fn flow(&self) -> &RustFlow {
        self.sheared
            .as_ref()
            .map_or(&self.flow, EpsilonShearedFlow::flow)
    }
    pub(crate) fn cache_target_weight(&self, column: usize, weight: &Atom) -> Atom {
        if let Some(shearing) = self.epsilon_shearing() {
            weight * Atom::var(self.family.epsilon).pow(i64::from(shearing.weights()[column]))
        } else {
            weight.clone()
        }
    }

    /// Compute an automatic AMF or FT boundary and append it to the growing bank.
    /// Extra requested seed digits reserve accuracy for subsequent transport;
    /// independent fitting still determines the stored verified accuracy.
    #[allow(clippy::too_many_arguments)] // Automatic solve inputs plus explicit seed accuracy reserve.
    pub fn seed_cache(
        &self,
        cache: &mut RustFlowCache,
        point: &BTreeMap<Symbol, Atom>,
        last: i32,
        extra_seed_digits: u32,
        options: &FlowOptions,
        backend: &dyn ReductionBackend,
        context: &RunContext,
    ) -> Result<CachedBoundary> {
        options.validate()?;
        context.cancellation.check()?;
        if !point.keys().eq(self.variables.iter())
            || options.dimension != self.family.dimension
            || options.prescription != self.prescription
        {
            return Err(Error::InvalidInput("physical seed coordinates, dimension or prescription do not match the prepared family".into()));
        }
        let cached_point = CachedPoint::Exact(point.clone());
        self.flow().identity().validate_point(&cached_point)?;
        let leading = self.ordinary_master_leading()?;
        let range = EpsilonRange::new(leading, last)?;
        let point = KinematicPoint(
            point
                .iter()
                .map(|(&s, a)| (Atom::var(s), a.clone()))
                .collect(),
        );
        let mut options = options.clone();
        options.digits = options
            .digits
            .checked_add(extra_seed_digits)
            .ok_or_else(|| Error::Limit("physical seed precision overflow".into()))?;
        options.validate()?;
        let expansions = solve_integrals(
            &self.family,
            &self.basis,
            &point,
            last,
            &options,
            backend,
            context,
        )?;
        if expansions.len() != self.basis.len() {
            return Err(Error::Numerical(
                "automatic seed returned the wrong number of master values".into(),
            ));
        }
        let verified = expansions
            .iter()
            .map(|e| {
                e.verified_digits.ok_or_else(|| {
                    Error::Accuracy("automatic seed lacks independent fit verification".into())
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .min()
            .unwrap();
        let bits = expansions.iter().map(|e| e.working_bits).min().unwrap();
        let mut coefficients = Vec::new();
        let mut errors = Vec::new();
        for power in leading..=last {
            coefficients.push(
                expansions
                    .iter()
                    .map(|e| {
                        e.coefficients.get(&power).cloned().ok_or_else(|| {
                            Error::Numerical(format!(
                                "automatic seed lacks epsilon coefficient {power}"
                            ))
                        })
                    })
                    .collect::<Result<_>>()?,
            );
            errors.push(
                expansions
                    .iter()
                    .map(|e| {
                        e.comparison_errors.get(&power).cloned().ok_or_else(|| {
                            Error::Accuracy(format!(
                                "automatic seed lacks comparison error at epsilon power {power}"
                            ))
                        })
                    })
                    .collect::<Result<_>>()?,
            );
        }
        let method = match options.recursion {
            RecursionMode::Amf => "AMF",
            RecursionMode::Ft => "FT",
        };
        let provenance = format!(
            "Automatic {method} master values (including quadratic deformation for linear families when applicable), checked by independent epsilon grids and precision/order refinement"
        );
        let boundary = CachedBoundary {
            // Automatic AMF/FT always evaluates actual integrals. Only the
            // exact, checked conversion below assigns a rescaled identity.
            identity: self.flow.identity().clone(),
            point: cached_point,
            kind: PointKind::Physical,
            range,
            coefficients,
            accuracy: BoundaryAccuracy::supplied(verified, bits, errors, &provenance)?,
        };
        let boundary = if let Some(sheared) = &self.sheared {
            sheared.to_sheared_boundary(&boundary, last)?
        } else {
            boundary
        };
        cache.insert(boundary.clone())?;
        Ok(boundary)
    }
}

pub(crate) fn physical_family_fingerprint(
    family: &IntegralFamily,
    variables: &BTreeSet<Symbol>,
    conditions: &[Atom],
) -> Result<String> {
    let mut variables = variables
        .iter()
        .map(|&s| Atom::var(s).to_canonical_string())
        .collect::<Vec<_>>();
    variables.sort();
    let metadata = (
        "physical-family-v1",
        &family.name,
        &family.loops,
        &family.external,
        family
            .external_gram
            .iter()
            .map(|row| {
                row.iter()
                    .map(AtomCore::to_canonical_string)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>(),
        family
            .propagators
            .iter()
            .map(|p| {
                (
                    p.constant.to_canonical_string(),
                    p.scalar_products
                        .iter()
                        .map(AtomCore::to_canonical_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>(),
        family.physical_propagators,
        family.dimension,
        Atom::var(family.epsilon).to_canonical_string(),
        variables,
        conditions
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>(),
    );
    Ok(
        blake3::hash(&serde_json::to_vec(&metadata).map_err(|e| Error::Cache(e.to_string()))?)
            .to_hex()
            .to_string(),
    )
}
