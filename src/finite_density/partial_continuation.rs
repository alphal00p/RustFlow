//! Private continuation bridge for one virtual loop and two occupied massless shells.
//! No public admission and no source/region/transport dispatch live here. Exact
//! native closure and complete symbolic Frobenius projection remain consumers.
use super::super::guarded::GuardedContext;
use super::super::source_class::{
    SourceImageBinding, SourceImageBudget, SourceImageClassCertificate, certify_source_image_class,
};
use super::endpoint::{PartialEndpointAudit, PartialEndpointBudget, PartialEndpointEvidence};
use super::*;
use crate::{Integral, reduction::ReducedSystem};

const CONTINUATION_VERSION: &str = "partial-continuation-h1-k2-v1";
const CONSUMER_CONTRACT: &str = "original-raw-domains;universal-source-image-class;explicit-final-labels;native-replayed-weighted-closure;all-hard-ordinary-plus-certified-virtual-soft;actual-convergent-regular-singular-frobenius;combine-rational-targets-before-shared-project-limit;right-half-plane-eta";

#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct PartialContinuationBudget {
    pub endpoint: PartialEndpointBudget,
    pub source_image_budget: SourceImageBudget,
    pub raw_expressions: usize,
    pub reduced_expressions: usize,
    pub target_text_bytes: usize,
}
impl Default for PartialContinuationBudget {
    fn default() -> Self {
        Self {
            endpoint: Default::default(),
            source_image_budget: Default::default(),
            raw_expressions: 65536,
            reduced_expressions: 262144,
            target_text_bytes: 65536,
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub(crate) struct PartialContinuationAudit {
    pub continuation_identity: String,
    pub endpoint: PartialEndpointAudit,
    pub combined_target_signature: Vec<String>,
    pub reduced_system_signature: String,
    pub nonzero_conditions: Vec<String>,
    #[serde(skip)]
    retained_condition_atoms: Vec<Atom>,
    pub enumerated_raw_derivative_labels: usize,
    pub native_weighted_reconstruction_required: bool,
    pub actual_symbolic_frobenius_required: bool,
    pub combined_target_projection_required: bool,
    pub universal_source_image_certificate_required: bool,
    pub boundary_region_certificates_required: bool,
    budget: PartialContinuationBudget,
}
impl PartialContinuationAudit {
    pub(crate) fn validate_reduced_system(
        &self,
        reduced: &ReducedSystem,
        run: &RunContext,
    ) -> Result<()> {
        run.cancellation.check()?;
        preflight_reduced(reduced, self.budget, run)?;
        if self
            .retained_condition_atoms
            .iter()
            .any(|c| !reduced.nonzero_conditions.contains(c))
        {
            return Err(invalid(
                "partial consumer omitted retained coefficient domains",
            ));
        }
        if self.reduced_system_signature != reduced_signature(reduced, run)? {
            return Err(invalid("partial reduced consumer binding changed"));
        }
        Ok(())
    }
}
/// Sealed proof bundle. Raw input domains are captured before any support-zero
/// classification. The narrow constructor cannot grant singleton Ward authority.
#[derive(Clone, Debug)]
pub(crate) struct PartialContinuation {
    origin: PartialOriginCapability,
    endpoint: PartialEndpointEvidence,
    budget: PartialContinuationBudget,
    source_identity: String,
    origin_identity: String,
    endpoint_identity: String,
    origin_loops: Vec<usize>,
    raw_nonzero_conditions: Vec<Atom>,
    epsilon: Symbol,
}
fn preflight_expressions<'a>(
    expressions: impl IntoIterator<Item = &'a Atom>,
    count_cap: usize,
    budget: PartialOriginBudget,
    run: &RunContext,
) -> Result<()> {
    let mut count = 0usize;
    let mut nodes = 0usize;
    let mut bytes = 0usize;
    for a in expressions {
        run.cancellation.check()?;
        count = count
            .checked_add(1)
            .ok_or_else(|| limit("expression count overflow"))?;
        if count > count_cap {
            return Err(limit("expression count budget"));
        }
        bytes = bytes
            .checked_add(a.as_view().get_byte_size())
            .ok_or_else(|| limit("aggregate expression bytes overflow"))?;
        if bytes > budget.aggregate_expression_bytes {
            return Err(limit("aggregate expression bytes budget"));
        }
        nodes = nodes
            .checked_add(resources::expression(a, budget.expressions, run)?)
            .ok_or_else(|| limit("aggregate expression overflow"))?;
        if nodes > budget.operations {
            return Err(limit("aggregate expression preflight"));
        }
    }
    Ok(())
}
/// Same ordering/domain algorithm as the reviewed singleton owner: validate raw
/// numerators and converted coefficients before cut removal or physical masses.
fn original_conditions(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
    b: PartialContinuationBudget,
    epsilon: Symbol,
    run: &RunContext,
) -> Result<Vec<Atom>> {
    run.cancellation.check()?;
    if b.raw_expressions == 0 || b.target_text_bytes == 0 {
        return Err(limit("zero raw-target budget"));
    }
    let mut expressions = Vec::new();
    let mut bytes = 0usize;
    for target in &input.input().targets {
        run.cancellation.check()?;
        bytes = bytes
            .checked_add(target.numerator.len())
            .ok_or_else(|| limit("raw target text overflow"))?;
        if bytes > b.target_text_bytes || expressions.len() >= b.raw_expressions {
            return Err(limit("raw target preflight"));
        }
        expressions.push(
            Atom::parse(&target.numerator, "rustflow_density", Default::default())
                .map_err(|e| invalid(&e.to_string()))?,
        );
    }
    for a in input.targets().iter().flat_map(|t| t.values()) {
        if expressions.len() >= b.raw_expressions {
            return Err(limit("raw converted target budget"));
        }
        expressions.push(a.clone());
    }
    preflight_expressions(&expressions, b.raw_expressions, b.endpoint.origin, run)?;
    let vars = input
        .basis()
        .coordinates()
        .iter()
        .flat_map(|a| a.get_all_symbols(true))
        .chain(input.independent_masses().iter().copied())
        .chain([epsilon])
        .collect::<BTreeSet<_>>();
    let masses = input
        .independent_masses()
        .iter()
        .zip(input.physical_masses())
        .map(|(&m, v)| (Atom::var(m), v.clone()))
        .collect::<BTreeMap<_, _>>();
    let raw = crate::physical_conditions::rational_denominator_conditions(&expressions, &vars)?;
    run.cancellation.check()?;
    let mut conditions = raw
        .iter()
        .map(|c| crate::family::substitute(c, &masses))
        .collect::<Vec<_>>();
    let converted = family
        .targets()
        .iter()
        .flat_map(|t| t.values())
        .collect::<Vec<_>>();
    preflight_expressions(
        converted.iter().copied(),
        b.raw_expressions,
        b.endpoint.origin,
        run,
    )?;
    let scalar_vars = BTreeSet::from([epsilon]);
    conditions.extend(crate::physical_conditions::rational_denominator_conditions(
        &converted.into_iter().cloned().collect::<Vec<_>>(),
        &scalar_vars,
    )?);
    preflight_expressions(&conditions, b.raw_expressions, b.endpoint.origin, run)?;
    let conditions = crate::physical_conditions::canonical_conditions(&conditions, &scalar_vars)?;
    run.cancellation.check()?;
    Ok(conditions)
}
fn combination_signature(
    combination: &crate::family::LinearCombination,
    run: &RunContext,
) -> Result<Vec<String>> {
    combination
        .iter()
        .map(|(i, c)| {
            run.cancellation.check()?;
            Ok(format!("{:?}:{}", i.0, c.to_canonical_string()))
        })
        .collect()
}
fn preflight_reduced(
    reduced: &ReducedSystem,
    budget: PartialContinuationBudget,
    run: &RunContext,
) -> Result<()> {
    run.cancellation.check()?;
    if !reduced.transformations.is_empty() {
        return Err(unsupported(
            "partial continuation requires an integral basis",
        ));
    }
    let mut count = 0usize;
    let mut axes = 0usize;
    for label in reduced
        .basis
        .iter()
        .chain(reduced.targets.iter().flat_map(|t| t.keys()))
        .chain(reduced.candidates.keys())
        .chain(reduced.candidates.values().flat_map(|t| t.keys()))
    {
        run.cancellation.check()?;
        count = count
            .checked_add(1)
            .ok_or_else(|| limit("reduced label count overflow"))?;
        axes = axes
            .checked_add(label.0.len())
            .ok_or_else(|| limit("reduced axes overflow"))?;
        if count > budget.endpoint.labels || axes > budget.endpoint.total_label_axes {
            return Err(limit("reduced label preflight"));
        }
    }
    preflight_expressions(
        reduced
            .matrix
            .iter()
            .flatten()
            .chain(reduced.targets.iter().flat_map(|t| t.values()))
            .chain(reduced.candidates.values().flat_map(|t| t.values()))
            .chain(&reduced.nonzero_conditions),
        budget.reduced_expressions,
        budget.endpoint.origin,
        run,
    )
}
fn reduced_signature(reduced: &ReducedSystem, run: &RunContext) -> Result<String> {
    run.cancellation.check()?;
    let targets = reduced
        .targets
        .iter()
        .map(|t| combination_signature(t, run))
        .collect::<Result<Vec<_>>>()?;
    let candidates = reduced
        .candidates
        .iter()
        .map(|(i, t)| Ok((i.0.clone(), combination_signature(t, run)?)))
        .collect::<Result<Vec<_>>>()?;
    let transforms=reduced.transformations.iter().map(|t|serde_json::json!({"previous":t.previous.iter().map(|i|&i.0).collect::<Vec<_>>(),"current":t.current.iter().map(|i|&i.0).collect::<Vec<_>>(),"matrix":t.matrix.iter().map(|r|r.iter().map(Atom::to_canonical_string).collect::<Vec<_>>()).collect::<Vec<_>>()})).collect::<Vec<_>>();
    let body = serde_json::json!({"basis":reduced.basis.iter().map(|i|&i.0).collect::<Vec<_>>(),"matrix":reduced.matrix.iter().map(|r|r.iter().map(Atom::to_canonical_string).collect::<Vec<_>>()).collect::<Vec<_>>(),"targets":targets,"conditions":reduced.nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),"candidates":candidates,"transformations":transforms});
    let bytes = serde_json::to_vec(&body).map_err(|e| Error::Cache(e.to_string()))?;
    run.cancellation.check()?;
    Ok(blake3::hash(&bytes).to_string())
}
impl PartialContinuation {
    pub(crate) fn new(
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        budget: PartialContinuationBudget,
        run: &RunContext,
    ) -> Result<Self> {
        run.cancellation.check()?;
        if family.shells().len() != 2
            || family.loops().checked_sub(family.shells().len()) != Some(1)
        {
            return Err(unsupported(
                "partial occupied continuation currently admits one virtual loop and two occupied shells",
            ));
        }
        if input.input().loops > budget.endpoint.origin.loops
            || input.input().edges.len() > budget.endpoint.origin.physical_rows
            || input.input().targets.len() > budget.raw_expressions
        {
            return Err(limit("input structure preflight"));
        }
        preflight_family(family, budget.endpoint.origin, run)?;
        let epsilon = symbol!("rustflow_occupied::epsilon");
        // This is intentionally before construction/classification of any zero support.
        let raw_nonzero_conditions = original_conditions(input, family, budget, epsilon, run)?;
        let origin = PartialOriginCapability::new_with_context(
            input,
            family,
            shifted,
            options,
            budget.endpoint.origin,
            run,
        )?;
        let endpoint = PartialEndpointEvidence::new_with_context(
            input,
            family,
            shifted,
            options,
            &origin,
            budget.endpoint,
            run,
        )?;
        let origin_identity = origin.identity()?;
        let endpoint_identity = endpoint.identity()?;
        let identity_body = serde_json::json!({"version":CONTINUATION_VERSION,"contract":CONSUMER_CONTRACT,"origin":origin_identity,"endpoint":endpoint_identity,"source_options":options,"raw_conditions":raw_nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),"budget":budget});
        let source_identity = format!(
            "{CONTINUATION_VERSION}:{}",
            blake3::hash(
                &serde_json::to_vec(&identity_body).map_err(|e| Error::Cache(e.to_string()))?
            )
        );
        run.cancellation.check()?;
        Ok(Self {
            origin,
            endpoint,
            budget,
            source_identity,
            origin_identity,
            endpoint_identity,
            origin_loops: family.shells().iter().map(|s| s.loop_index).collect(),
            raw_nonzero_conditions,
            epsilon,
        })
    }
    pub(crate) fn source_identity(&self) -> &str {
        &self.source_identity
    }
    pub(crate) fn origin_identity(&self) -> &str {
        &self.origin_identity
    }
    pub(crate) fn source_options(&self) -> WeightedSourceOptions {
        self.origin.source_options
    }
    pub(crate) fn certified_origin_loops(&self) -> &[usize] {
        &self.origin_loops
    }
    pub(crate) fn raw_nonzero_conditions(&self) -> &[Atom] {
        &self.raw_nonzero_conditions
    }
    pub(crate) fn validate_family_binding(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        run: &RunContext,
    ) -> Result<()> {
        preflight_family(family, self.budget.endpoint.origin, run)?;
        if self.origin.family_signature != family_signature(family)
            || self.origin.shifted != shifted
            || self.origin.source_options != options
            || self.origin.version != VERSION
            || self.origin.prescription != PRESCRIPTION
        {
            return Err(invalid("partial bundle family/placement/options binding"));
        }
        run.cancellation.check()?;
        Ok(())
    }
    fn source_binding(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        run: &RunContext,
    ) -> Result<SourceImageBinding> {
        self.validate_family_binding(family, shifted, options, run)?;
        Ok(SourceImageBinding {
            input_identity: self.origin.input_identity.clone(),
            family_signature: self.origin.family_signature.clone(),
            origin_identity: self.origin_identity.clone(),
            shifted_slots: shifted.to_vec(),
            source_options: options,
            expected_roles: family.roles().to_vec(),
            physical_arity: family.factors().len(),
            physical_slots: family.physical_slots(),
            input_slots: family.input_slots(),
        })
    }
    pub(crate) fn certify_sources<const N: usize>(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        context: &GuardedContext<N>,
        admitted: &IndexDomain<N>,
        run: &RunContext,
    ) -> Result<SourceImageClassCertificate> {
        let binding = self.source_binding(family, shifted, options, run)?;
        certify_source_image_class(
            context,
            admitted,
            &binding,
            self.budget.source_image_budget,
            run,
        )
    }
    pub(crate) fn validate_sources<const N: usize>(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        context: &GuardedContext<N>,
        admitted: &IndexDomain<N>,
        certificate: &SourceImageClassCertificate,
        run: &RunContext,
    ) -> Result<()> {
        let binding = self.source_binding(family, shifted, options, run)?;
        certificate.validate_context(context, admitted, &binding, run)
    }
    pub(crate) fn free_virtual_zero_domains<const N: usize>(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        run: &RunContext,
    ) -> Result<Vec<IndexDomain<N>>> {
        self.validate_family_binding(family, shifted, options, run)?;
        // The origin API fixes the leading shell-count entries to its lower-contact
        // boxes. Source preparation owns those; append only independent virtual zeros.
        let boxes = self.origin.zero_domains_with_context(
            family,
            shifted,
            options,
            self.budget.endpoint.origin,
            run,
        )?;
        Ok(boxes.into_iter().skip(family.shells().len()).collect())
    }
    pub(crate) fn audit_labels(
        &self,
        family: &OccupiedCutFamily,
        labels: &[Integral],
        run: &RunContext,
    ) -> Result<PartialEndpointAudit> {
        self.validate_family_binding(
            family,
            &self.origin.shifted,
            self.origin.source_options,
            run,
        )?;
        if labels.len() > self.budget.endpoint.labels {
            return Err(limit("bundle finite label budget"));
        }
        let labels = labels
            .iter()
            .map(|i| i.0.iter().map(|n| i64::from(*n)).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        self.endpoint
            .audit_labels_with_context(family, &labels, self.budget.endpoint, run)
    }
    /// Bind every raw domain to the actual object consumed by transport, then
    /// hash that exact object. Work completes on a bounded clone before commit.
    pub(crate) fn audit_and_bind_reduced_system(
        &self,
        family: &OccupiedCutFamily,
        reduced: &mut ReducedSystem,
        eta: Symbol,
        epsilon: Symbol,
        run: &RunContext,
    ) -> Result<PartialContinuationAudit> {
        let mut audit = self.audit_reduced_system(family, reduced, eta, epsilon, run)?;
        let mut bound = reduced.clone();
        bound.nonzero_conditions = audit.retained_condition_atoms.clone();
        let signature = reduced_signature(&bound, run)?;
        run.cancellation.check()?;
        reduced.nonzero_conditions = bound.nonzero_conditions;
        audit.reduced_system_signature = signature;
        audit.endpoint.original_coefficient_conditions_checked = true;
        Ok(audit)
    }
    pub(crate) fn audit_reduced_system(
        &self,
        family: &OccupiedCutFamily,
        reduced: &ReducedSystem,
        eta: Symbol,
        epsilon: Symbol,
        run: &RunContext,
    ) -> Result<PartialContinuationAudit> {
        self.validate_family_binding(
            family,
            &self.origin.shifted,
            self.origin.source_options,
            run,
        )?;
        if epsilon != self.epsilon
            || eta == epsilon
            || reduced.targets.len() != family.targets().len()
            || reduced.matrix.len() != reduced.basis.len()
            || reduced
                .matrix
                .iter()
                .any(|r| r.len() != reduced.basis.len())
        {
            return Err(invalid("partial reduced-system dimensions/parameters"));
        }
        if !reduced.transformations.is_empty() {
            return Err(unsupported(
                "partial continuation requires the native integral basis without a later basis transform",
            ));
        }
        preflight_reduced(reduced, self.budget, run)?;
        let expressions = reduced
            .matrix
            .iter()
            .flatten()
            .chain(reduced.targets.iter().flat_map(|t| t.values()))
            .chain(reduced.candidates.values().flat_map(|t| t.values()))
            .chain(&reduced.nonzero_conditions)
            .collect::<Vec<_>>();
        preflight_expressions(
            expressions.iter().copied(),
            self.budget.reduced_expressions,
            self.budget.endpoint.origin,
            run,
        )?;
        // Capture reduced target/matrix denominator domains BEFORE cancellation, too.
        let vars = BTreeSet::from([epsilon, eta]);
        let mut conditions = self.raw_nonzero_conditions.clone();
        conditions.extend(reduced.nonzero_conditions.iter().cloned());
        conditions.extend(crate::physical_conditions::rational_denominator_conditions(
            &expressions.into_iter().cloned().collect::<Vec<_>>(),
            &vars,
        )?);
        preflight_expressions(
            &conditions,
            self.budget.reduced_expressions,
            self.budget.endpoint.origin,
            run,
        )?;
        let nonzero = crate::physical_conditions::canonical_conditions(&conditions, &vars)?;
        run.cancellation.check()?;
        let mut labels = BTreeSet::new();
        let input_labels = family
            .targets()
            .iter()
            .flat_map(|t| t.keys())
            .chain(&reduced.basis)
            .chain(reduced.targets.iter().flat_map(|t| t.keys()))
            .chain(reduced.candidates.keys())
            .chain(reduced.candidates.values().flat_map(|t| t.keys()));
        let mut count = 0usize;
        for label in input_labels {
            run.cancellation.check()?;
            count = count
                .checked_add(1)
                .ok_or_else(|| limit("consumer label count overflow"))?;
            if count > self.budget.endpoint.labels {
                return Err(limit("consumer label count budget"));
            }
            labels.insert(label.clone());
        }
        let mut derivative_count = 0usize;
        for label in &reduced.basis {
            if label.0.len() != family.factors().len() {
                return Err(invalid("physical basis arity"));
            }
            for &slot in &self.origin.shifted {
                if label.0[slot] == 0 {
                    continue;
                }
                let mut raised = label.clone();
                raised.0[slot] = raised.0[slot]
                    .checked_add(1)
                    .ok_or_else(|| limit("raw derivative index overflow"))?;
                derivative_count = derivative_count
                    .checked_add(1)
                    .ok_or_else(|| limit("raw derivative count overflow"))?;
                if count
                    .checked_add(derivative_count)
                    .is_none_or(|n| n > self.budget.endpoint.labels)
                {
                    return Err(limit("raw derivative label budget"));
                }
                labels.insert(raised);
            }
        }
        let endpoint = self.audit_labels(family, &labels.into_iter().collect::<Vec<_>>(), run)?;
        let combined_target_signature = reduced
            .targets
            .iter()
            .enumerate()
            .map(|(i, t)| {
                Ok(format!(
                    "target={i}:{}",
                    serde_json::to_string(&combination_signature(t, run)?)
                        .map_err(|e| Error::Cache(e.to_string()))?
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(PartialContinuationAudit {
            continuation_identity: self.source_identity.clone(),
            endpoint,
            combined_target_signature,
            reduced_system_signature: reduced_signature(reduced, run)?,
            nonzero_conditions: nonzero.iter().map(Atom::to_canonical_string).collect(),
            retained_condition_atoms: nonzero,
            enumerated_raw_derivative_labels: derivative_count,
            native_weighted_reconstruction_required: true,
            actual_symbolic_frobenius_required: true,
            combined_target_projection_required: true,
            universal_source_image_certificate_required: true,
            boundary_region_certificates_required: true,
            budget: self.budget,
        })
    }
    pub(crate) fn report(&self) -> serde_json::Value {
        serde_json::json!({"version":CONTINUATION_VERSION,"source_identity":self.source_identity,"origin_identity":self.origin_identity,"endpoint_identity":self.endpoint_identity,"input_identity":self.origin.input_identity,"family_signature":self.origin.family_signature,"shifted_slots":self.origin.shifted,"source_options":self.origin.source_options,"raw_nonzero_conditions":self.raw_nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),"consumer_contract":CONSUMER_CONTRACT,"native_closure":false,"transport_performed":false,"raw_ward_authority":false,"infinity_boundary_authority":false,"endpoint_requirements":self.endpoint.consumer_requirements(),"limits":self.budget,"cooperative_cancellation":true,"internal_cas_work_bound":false})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::{DensityInput, DensityTarget};
    fn input() -> DensityInput {
        serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/finite_density/massless_three_loop_chain.json"
        )))
        .unwrap()
    }
    fn options() -> WeightedSourceOptions {
        WeightedSourceOptions {
            free_virtual_zero_sectors: true,
            ..Default::default()
        }
    }
    fn construct(
        raw: DensityInput,
        shifted: &[usize],
    ) -> (PreparedDensityInput, OccupiedCutFamily, PartialContinuation) {
        let input = raw.prepare().unwrap();
        let family = input
            .occupied_cut(&[0, 3], 1024)
            .unwrap()
            .at_physical_masses();
        let evidence = PartialContinuation::new(
            &input,
            &family,
            shifted,
            options(),
            Default::default(),
            &RunContext::default(),
        )
        .unwrap();
        (input, family, evidence)
    }
    fn reduced(f: &OccupiedCutFamily) -> ReducedSystem {
        let basis = f
            .targets()
            .iter()
            .flat_map(|t| t.keys())
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        ReducedSystem {
            matrix: vec![vec![Atom::zero(); basis.len()]; basis.len()],
            basis,
            targets: f.targets().to_vec(),
            nonzero_conditions: vec![],
            candidates: BTreeMap::new(),
            transformations: vec![],
        }
    }
    #[test]
    fn narrow_actual_masks_bind_distinct_authorities_and_no_lower_zero_duplicates() {
        for shift in [vec![1, 2], vec![2, 4]] {
            let (_, f, p) = construct(input(), &shift);
            let run = RunContext::default();
            assert_ne!(p.source_identity(), p.origin_identity());
            assert_eq!(p.certified_origin_loops().len(), 2);
            let actual = p
                .free_virtual_zero_domains::<16>(&f, &shift, options(), &run)
                .unwrap();
            let complete = p
                .origin
                .zero_domains_with_context::<16>(
                    &f,
                    &shift,
                    options(),
                    p.budget.endpoint.origin,
                    &run,
                )
                .unwrap();
            assert_eq!(actual.len() + f.shells().len(), complete.len());
            for domain in actual {
                for s in f.shells() {
                    assert!(domain.bounds()[s.lower_slot].contains(0));
                }
            }
            assert_eq!(p.report()["raw_ward_authority"], false);
            assert_eq!(p.report()["infinity_boundary_authority"], false);
            let wrong = if shift == [1, 2] {
                vec![2, 4]
            } else {
                vec![1, 2]
            };
            assert!(
                p.validate_family_binding(&f, &wrong, options(), &run)
                    .is_err()
            );
            assert!(
                p.validate_family_binding(&f, &shift, WeightedSourceOptions::default(), &run)
                    .is_err()
            );
        }
        let i = input().prepare().unwrap();
        let f = i.occupied_cut(&[3], 1024).unwrap().at_physical_masses();
        assert!(
            PartialContinuation::new(
                &i,
                &f,
                &[0],
                options(),
                Default::default(),
                &RunContext::default()
            )
            .is_err()
        );
    }
    #[test]
    fn raw_removed_cut_poles_and_symbols_are_checked_before_zero_support() {
        let mut raw = input();
        raw.targets = vec![DensityTarget {
            powers: vec![0, 1, 1, 1, 1],
            numerator: "(rustflow_occupied::epsilon^2-1)/(rustflow_occupied::epsilon-1)".into(),
        }];
        let (_, f, p) = construct(raw.clone(), &[1, 2]);
        assert!(f.targets()[0].is_empty());
        assert!(!p.raw_nonzero_conditions().is_empty());
        let eps = symbol!("rustflow_occupied::epsilon");
        assert!(
            crate::physical_conditions::validate_conditions_at(
                p.raw_nonzero_conditions(),
                eps,
                &BTreeMap::from([(eps, Atom::one())])
            )
            .is_err()
        );
        for bad in [
            "rustflow_occupied::eta",
            "1/rustflow_occupied::eta",
            "not_assigned",
        ] {
            raw.targets[0].numerator = bad.into();
            let i = raw.prepare().unwrap();
            let f = i.occupied_cut(&[0, 3], 1024).unwrap().at_physical_masses();
            assert!(
                PartialContinuation::new(
                    &i,
                    &f,
                    &[1, 2],
                    options(),
                    Default::default(),
                    &RunContext::default()
                )
                .is_err()
            );
        }
    }
    #[test]
    fn final_audit_includes_raw_derivatives_candidates_and_combined_pole_weights() {
        let (_, f, p) = construct(input(), &[1, 2]);
        let mut r = reduced(&f);
        let eta = symbol!("rustflow_occupied::eta");
        let eps = symbol!("rustflow_occupied::epsilon");
        let original = r.basis[0].clone();
        let mut candidate = original.clone();
        candidate.0[f.shells()[0].upper_slot] = 2;
        r.candidates.insert(
            candidate.clone(),
            BTreeMap::from([(original.clone(), Atom::var(eta).pow(-2))]),
        );
        // This fixture exercises auditing, not a claimed native differential system.
        r.targets[0].insert(original, Atom::var(eta).pow(-1));
        let unbound = p
            .audit_reduced_system(&f, &r, eta, eps, &RunContext::default())
            .unwrap();
        assert!(
            unbound
                .validate_reduced_system(&r, &RunContext::default())
                .is_err()
        );
        let a = p
            .audit_and_bind_reduced_system(&f, &mut r, eta, eps, &RunContext::default())
            .unwrap();
        assert!(
            r.nonzero_conditions
                .iter()
                .any(|c| c.contains(Atom::var(eta).as_view()))
        );
        assert!(a.enumerated_raw_derivative_labels > 0);
        assert!(a.nonzero_conditions.iter().any(|s| s.contains("eta")));
        assert!(a.endpoint.labels.iter().any(|b| {
            b.label
                == candidate
                    .0
                    .iter()
                    .map(|n| i64::from(*n))
                    .collect::<Vec<_>>()
        }));
        assert!(
            a.actual_symbolic_frobenius_required
                && a.native_weighted_reconstruction_required
                && a.combined_target_projection_required
        );
        a.validate_reduced_system(&r, &RunContext::default())
            .unwrap();
        r.targets[0].clear();
        assert!(
            a.validate_reduced_system(&r, &RunContext::default())
                .is_err()
        );
    }
    #[test]
    fn invalid_roles_and_tail_reject_before_required_cut_zero() {
        let (_, f, p) = construct(input(), &[1, 2]);
        let mut label = f.targets()[0].keys().next().unwrap().clone();
        let s = &f.shells()[0];
        label.0[s.physical_slot] = 0;
        let mut inverse = label.clone();
        inverse.0[f.physical_slots()] = 1;
        assert!(
            p.audit_labels(&f, &[inverse], &RunContext::default())
                .is_err()
        );
        let mut occupation = label.clone();
        occupation.0[s.lower_slot] = -1;
        assert!(
            p.audit_labels(&f, &[occupation], &RunContext::default())
                .is_err()
        );
        let mut tail = label;
        tail.0.push(1);
        assert!(p.audit_labels(&f, &[tail], &RunContext::default()).is_err());
    }
    #[test]
    fn pre_cancelled_construction_and_consumers_publish_nothing() {
        let (i, f, p) = construct(input(), &[1, 2]);
        let run = RunContext::default();
        run.cancellation.cancel();
        assert!(
            PartialContinuation::new(&i, &f, &[1, 2], options(), Default::default(), &run).is_err()
        );
        assert!(
            p.validate_family_binding(&f, &[1, 2], options(), &run)
                .is_err()
        );
        assert!(
            p.free_virtual_zero_domains::<16>(&f, &[1, 2], options(), &run)
                .is_err()
        );
        assert!(p.audit_labels(&f, &[], &run).is_err());
        assert!(
            p.audit_reduced_system(
                &f,
                &reduced(&f),
                symbol!("rustflow_occupied::eta"),
                symbol!("rustflow_occupied::epsilon"),
                &run
            )
            .is_err()
        );
    }
    #[test]
    fn labels_expressions_and_raw_target_budgets_decline_without_partial_evidence() {
        let (i, f, mut p) = construct(input(), &[1, 2]);
        let mut b = PartialContinuationBudget::default();
        b.raw_expressions = 1;
        assert!(
            PartialContinuation::new(&i, &f, &[1, 2], options(), b, &RunContext::default())
                .is_err()
        );
        p.budget.endpoint.labels = 0;
        assert!(p.audit_labels(&f, &[], &RunContext::default()).is_err());
        let (_, f, mut p) = construct(input(), &[1, 2]);
        p.budget.reduced_expressions = 0;
        assert!(
            p.audit_reduced_system(
                &f,
                &reduced(&f),
                symbol!("rustflow_occupied::eta"),
                symbol!("rustflow_occupied::epsilon"),
                &RunContext::default()
            )
            .is_err()
        );
    }
    #[test]
    fn typed_source_wrapper_binds_real_rows_and_rejects_wrong_origin_or_domain() {
        use super::super::super::guarded::GuardedMeasureIdentity;
        for shifted in [vec![1, 2], vec![2, 4]] {
            let i = input().prepare().unwrap();
            let f = i.occupied_cut(&[0, 3], 1024).unwrap().at_physical_masses();
            let options = WeightedSourceOptions::default();
            let run = RunContext::default();
            let p = PartialContinuation::new(&i, &f, &shifted, options, Default::default(), &run)
                .unwrap();
            let id = GuardedMeasureIdentity {
                measure: p.source_identity().to_owned(),
                support: "polynomial class; no source theorem inferred by checker".into(),
                orientation: "existing occupied routing".into(),
                normalization: "native C/H".into(),
                branch: PRESCRIPTION.into(),
                deformation: format!("fixed shells; shifted={shifted:?}"),
            };
            let sources = f
                .guarded_sources::<16>(
                    symbol!("rustflow_occupied::epsilon"),
                    4,
                    symbol!("rustflow_occupied::eta"),
                    &shifted,
                    4096,
                    vec![],
                    id,
                )
                .unwrap();
            let cert = p
                .certify_sources(
                    &f,
                    &shifted,
                    options,
                    &sources.context,
                    sources.deformation.admitted_domain(),
                    &run,
                )
                .unwrap();
            p.validate_sources(
                &f,
                &shifted,
                options,
                &sources.context,
                sources.deformation.admitted_domain(),
                &cert,
                &run,
            )
            .unwrap();
            let wrong = if shifted == [1, 2] {
                vec![2, 4]
            } else {
                vec![1, 2]
            };
            assert!(
                p.validate_sources(
                    &f,
                    &wrong,
                    options,
                    &sources.context,
                    sources.deformation.admitted_domain(),
                    &cert,
                    &run
                )
                .is_err()
            );
            assert!(
                p.certify_sources(
                    &f,
                    &shifted,
                    options,
                    &sources.context,
                    &IndexDomain::for_roles(sources.context.sources().roles()),
                    &run
                )
                .is_err()
            );
            assert_eq!(cert.report()["individual_application_trace"], false);
        }
    }
}
