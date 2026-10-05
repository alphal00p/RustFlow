//! Exact caller-declared auxiliary-mass systems, without an additional IBP search.
use super::*;
use crate::reduction::LinearCombination;
use std::collections::BTreeSet;

/// A closed ordinary-integral system declared by its provider.
///
/// The matrix and target maps are already specialized to the kinematic point
/// passed to [`PreparedFlow::from_supplied`]. They may depend only on `variable`
/// and the family's epsilon symbol. This declaration is a trust boundary:
/// structural validation does not independently prove the supplied IBP identities.
#[derive(Clone, Debug)]
pub struct SuppliedAuxiliarySystem {
    pub variable: Symbol,
    pub reduced: ReducedSystem,
    /// One slot per propagator, including false slots for irreducible numerators.
    /// A true slot declares the deformation `D_i - variable`.
    pub deformation_mask: Vec<bool>,
    /// The source and certification of the supplied exact connection and maps.
    pub provenance: String,
}

pub(super) struct SuppliedSeal {
    pub(super) identity: String,
    pub(super) provenance: String,
    targets: Vec<Integral>,
}

impl PreparedFlow {
    /// Prepare an explicitly declared closed system in the ordinary integral basis.
    ///
    /// Only the undeformed family is specialized with `point`. The supplied
    /// connection and maps must already use that same point, dimension, routing,
    /// normalization, and deformation. Their mathematical correctness is the
    /// provider's responsibility; this constructor checks the structural contract.
    pub fn from_supplied(
        family: &IntegralFamily,
        targets: &[Integral],
        point: &KinematicPoint,
        supplied: SuppliedAuxiliarySystem,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        context.cancellation.check()?;
        options.validate()?;
        if options.recursion != RecursionMode::Amf {
            return Err(Error::Unsupported(
                "supplied auxiliary systems require AMF recursion".into(),
            ));
        }
        if options.refine_basis || !supplied.reduced.transformations.is_empty() {
            return Err(Error::Unsupported(
                "supplied auxiliary systems require explicit ordinary integral bases without further basis refinement".into(),
            ));
        }
        if family.dimension != options.dimension {
            return Err(Error::InvalidInput(
                "supplied system dimension differs from evaluation dimension".into(),
            ));
        }
        if supplied.provenance.trim().is_empty() || targets.is_empty() {
            return Err(Error::InvalidInput(
                "supplied systems require targets and source provenance".into(),
            ));
        }
        let family = family.at(point);
        family.validate()?;
        crate::boundary::validate_deformation_mask(&family, &supplied.deformation_mask)?;
        let variable = supplied.variable;
        if variable == family.epsilon || variable == crate::family::imaginary_parameter() {
            return Err(Error::InvalidInput(
                "auxiliary variable collides with epsilon or the reserved imaginary unit".into(),
            ));
        }
        for value in family.external_gram.iter().flatten().chain(
            family
                .propagators
                .iter()
                .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
        ) {
            validate_coefficient(value, &[family.epsilon])?;
        }
        let reduced = supplied.reduced;
        let basis = reduced.basis.iter().cloned().collect::<BTreeSet<_>>();
        if basis.len() != reduced.basis.len() {
            return Err(Error::InvalidInput(
                "duplicate supplied basis integral".into(),
            ));
        }
        for integral in targets
            .iter()
            .chain(&reduced.basis)
            .chain(reduced.candidates.keys())
        {
            family.validate_integral(integral)?;
        }
        if reduced.targets.len() != targets.len()
            || reduced.matrix.len() != basis.len()
            || reduced.matrix.iter().any(|row| row.len() != basis.len())
        {
            return Err(Error::InvalidInput(
                "supplied basis, matrix, or target-map dimensions differ".into(),
            ));
        }
        for row in reduced.targets.iter().chain(reduced.candidates.values()) {
            if row.keys().any(|integral| !basis.contains(integral)) {
                return Err(Error::IncompleteReduction(
                    "supplied map contains an integral outside its declared closed basis".into(),
                ));
            }
        }
        for coefficient in reduced
            .matrix
            .iter()
            .flatten()
            .chain(reduced.targets.iter().flat_map(|row| row.values()))
            .chain(reduced.candidates.values().flat_map(|row| row.values()))
            .chain(&reduced.nonzero_conditions)
        {
            validate_coefficient(coefficient, &[variable, family.epsilon])?;
        }
        if reduced
            .nonzero_conditions
            .iter()
            .any(|condition| condition.together().cancel().is_zero())
        {
            return Err(Error::Reduction(
                "a supplied nonzero condition is identically zero".into(),
            ));
        }
        let principal_mass_constraints =
            crate::vacuum_peel::principal_mass_constraints(&family, targets)?;
        crate::vacuum_peel::validate_principal_mass_constraints(
            &principal_mass_constraints,
            options.prescription,
        )?;
        let system = DifferentialSystem {
            variable,
            matrix: reduced.matrix.clone(),
        };
        let blocks = if basis.is_empty() {
            vec![]
        } else {
            system.blocks()?.iter().map(Vec::len).collect()
        };
        let mut flow = Self {
            family,
            reduced,
            system,
            epsilon_sample: None,
            deformation_mask: supplied.deformation_mask,
            supplied: None,
            principal_mass_constraints,
            basis_refinement: None,
        };
        let identity = declaration_identity(&flow, targets, &supplied.provenance)?;
        flow.supplied = Some(SuppliedSeal {
            identity,
            provenance: supplied.provenance,
            targets: targets.to_vec(),
        });
        context.emit(Progress::Prepared {
            basis_size: basis.len(),
            blocks,
        })?;
        Ok(flow)
    }

    /// Content identity for caller-declared systems, retained in sample caches.
    pub(crate) fn source_identity(&self) -> Option<&str> {
        self.supplied.as_ref().map(|seal| seal.identity.as_str())
    }

    /// Declared exact source, when preparation bypassed a new reduction search.
    pub fn supplied_provenance(&self) -> Option<&str> {
        self.supplied.as_ref().map(|seal| seal.provenance.as_str())
    }

    /// The actual fixed auxiliary deformation used for preparation and boundaries.
    pub fn deformation_mask(&self) -> &[bool] {
        &self.deformation_mask
    }

    pub(crate) fn validate_supplied_contract(&self, options: &FlowOptions) -> Result<()> {
        if let Some(seal) = &self.supplied {
            if options.dimension != self.family.dimension {
                return Err(Error::InvalidInput(
                    "evaluation dimension differs from supplied connection dimension".into(),
                ));
            }
            if self.basis_refinement.is_some()
                || declaration_identity(self, &seal.targets, &seal.provenance)? != seal.identity
            {
                return Err(Error::InvalidInput("a prepared supplied-system declaration was modified; prepare a new system to retain a valid source identity".into()));
            }
        }
        Ok(())
    }
}

fn validate_coefficient(value: &Atom, allowed: &[Symbol]) -> Result<()> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(value.as_view(), &mut symbols)?;
    if symbols.iter().any(|symbol| {
        *symbol != Atom::var(crate::family::imaginary_parameter())
            && !allowed.iter().any(|&allowed| *symbol == Atom::var(allowed))
    }) {
        return Err(Error::InvalidInput(format!(
            "supplied system has an undeclared scalar parameter in {value}"
        )));
    }
    let rational: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(value)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|error| {
            Error::Unsupported(format!(
                "supplied coefficients must be exact rational functions: {error}"
            ))
        })?;
    // Symbolica can introduce a formal variable for a nonrational subtree.
    // That facility is useful for algebra, but does not certify a rational ODE.
    if rational.numerator.variables().iter().any(|variable| {
        !matches!(variable, PolyVariable::Symbol(symbol)
            if *symbol == crate::family::imaginary_parameter() || allowed.contains(symbol))
    }) {
        return Err(Error::Unsupported(
            "supplied ordinary connections must be rational in eta and epsilon; algebraic generators must be eliminated first".into(),
        ));
    }
    Ok(())
}

fn declaration_identity(
    flow: &PreparedFlow,
    targets: &[Integral],
    provenance: &str,
) -> Result<String> {
    let expression = AtomCore::to_canonical_string;
    let matrix = |rows: &[Vec<Atom>]| {
        rows.iter()
            .map(|row| row.iter().map(expression).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    };
    let combination = |row: &LinearCombination| {
        row.iter()
            .map(|(i, a)| (i.0.clone(), expression(a)))
            .collect::<Vec<_>>()
    };
    let metadata = serde_json::json!({
        "version": "supplied-auxiliary-system-v1",
        "family": crate::physical_family::physical_family_fingerprint(&flow.family, &BTreeSet::new(), &[])?,
        "requested_targets": targets.iter().map(|i| &i.0).collect::<Vec<_>>(),
        "variable": Atom::var(flow.system.variable).to_canonical_string(),
        "system_matrix": matrix(&flow.system.matrix),
        "reduced_matrix": matrix(&flow.reduced.matrix),
        "basis": flow.reduced.basis.iter().map(|i| &i.0).collect::<Vec<_>>(),
        "target_maps": flow.reduced.targets.iter().map(combination).collect::<Vec<_>>(),
        "candidate_maps": flow.reduced.candidates.iter().map(|(i, row)| (&i.0, combination(row))).collect::<Vec<_>>(),
        "conditions": flow.reduced.nonzero_conditions.iter().map(expression).collect::<Vec<_>>(),
        "transformations": flow.reduced.transformations.iter().map(|t| (
            t.previous.iter().map(|i| &i.0).collect::<Vec<_>>(),
            t.current.iter().map(|i| &i.0).collect::<Vec<_>>(), matrix(&t.matrix),
        )).collect::<Vec<_>>(),
        "deformation": flow.deformation_mask,
        "epsilon_sample": flow.epsilon_sample.as_ref().map(ToString::to_string),
        "provenance": provenance,
    });
    Ok(blake3::hash(
        &serde_json::to_vec(&metadata).map_err(|error| Error::Cache(error.to_string()))?,
    )
    .to_hex()
    .to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn supplied_identity_covers_connection_maps_mask_guards_and_origin() {
        let epsilon = symbol!("supplied_identity::epsilon");
        let variable = symbol!("supplied_identity::eta");
        let family = IntegralFamily {
            name: "fingerprint tadpoles".into(),
            loops: vec!["l".into()],
            external: vec!["p".into()],
            external_gram: vec![vec![Atom::num(-1)]],
            propagators: vec![
                Propagator {
                    constant: Atom::num(-2),
                    scalar_products: vec![Atom::num(1), Atom::num(0)],
                },
                Propagator {
                    constant: Atom::num(-4),
                    scalar_products: vec![Atom::num(1), Atom::num(2)],
                },
            ],
            physical_propagators: 2,
            epsilon,
            dimension: 4,
        };
        let integral = Integral(vec![1, 0]);
        let input = SuppliedAuxiliarySystem {
            variable,
            reduced: ReducedSystem {
                basis: vec![integral.clone()],
                matrix: vec![vec![
                    (Atom::num(1) - Atom::var(epsilon)) / (Atom::num(2) + Atom::var(variable)),
                ]],
                targets: vec![BTreeMap::from([(integral.clone(), Atom::num(1))])],
                candidates: BTreeMap::new(),
                nonzero_conditions: vec![],
                transformations: vec![],
            },
            deformation_mask: vec![true, false],
            provenance: "declared analytic connection v1".into(),
        };
        let identity = |input| {
            PreparedFlow::from_supplied(
                &family,
                std::slice::from_ref(&integral),
                &KinematicPoint::default(),
                input,
                &FlowOptions::default(),
                &RunContext::default(),
            )
            .unwrap()
            .source_identity()
            .unwrap()
            .to_owned()
        };
        let original = identity(input.clone());
        assert_eq!(identity(input.clone()), original);
        let mut variants = vec![input.clone(); 5];
        variants[0].reduced.matrix[0][0] *= 2;
        variants[1].reduced.targets[0].insert(integral.clone(), Atom::num(2));
        variants[2].deformation_mask[1] = true;
        variants[3]
            .reduced
            .nonzero_conditions
            .push(Atom::var(epsilon));
        variants[4].provenance.push_str(" revised");
        let identities = variants
            .into_iter()
            .map(identity)
            .chain([original])
            .collect::<BTreeSet<_>>();
        assert_eq!(identities.len(), 6);
    }
}
