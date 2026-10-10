//! Isolated finite-label endpoint bound owner, nested under the frozen origin owner.
//! No values, source zeros, mode projection, boundary, or production dispatch.
//! RunContext cancellation and internal CAS work limits are not implemented here.
use super::*;

const ENDPOINT_VERSION: &str = "partial-endpoint-rank-two-finite-JN-holder-v1";
const ENDPOINT_ORDER: &str = "same-sign-independent-line-masses;future-real-spatial-shells;right-half-plane-eta;fixed-positive-T-pole-free-regulator-removal;joint-eta-regulator-kernel-bound;real-Fermi-limit;then-meromorphic-D";

#[derive(Clone, Copy, Debug, Serialize)]
pub struct PartialEndpointBudget {
    pub origin: PartialOriginBudget,
    /// Per audit call, before deduplication. No partial audit on exhaustion.
    pub labels: usize,
    pub total_label_axes: usize,
    pub index_magnitude: u64,
    pub degree_sum: u64,
    pub dimension_integer_bound: u64,
}
impl Default for PartialEndpointBudget {
    fn default() -> Self {
        Self {
            origin: Default::default(),
            labels: 65536,
            total_label_axes: 1048576,
            index_magnitude: 65536,
            degree_sum: 1048576,
            dimension_integer_bound: 16777216,
        }
    }
}
fn endpoint_limit(message: &str) -> Error {
    Error::Limit(format!("partial endpoint: {message}"))
}
fn endpoint_invalid(message: &str) -> Error {
    Error::InvalidInput(format!("partial endpoint: {message}"))
}
fn integer(value: u64) -> Result<Rational> {
    Ok(Rational::from(i64::try_from(value).map_err(|_| {
        endpoint_limit("rational integer overflow")
    })?))
}
#[derive(Clone, Debug, Serialize)]
pub struct StrictEndpointMargin {
    pub stratum: String,
    pub value: String,
}
fn positive_margin(stratum: String, value: Rational) -> Result<StrictEndpointMargin> {
    if value <= 0 {
        return Err(endpoint_invalid(&format!(
            "nonpositive witness at {stratum}"
        )));
    }
    Ok(StrictEndpointMargin {
        stratum,
        value: value.to_string(),
    })
}
#[derive(Clone, Debug, Serialize)]
pub struct PartialEndpointLabelAudit {
    pub label: Vec<i64>,
    pub origin_audit: FiniteLabelOriginAudit,
    pub positive_virtual_degree: u64,
    pub positive_pure_transfer_degree: u64,
    pub polynomial_degree: u64,
    pub compact_jet_orders: Vec<u64>,
    pub total_jet_order: u64,
    pub dimension_witness: String,
    pub taylor_terms: u64,
    pub sigma: String,
    pub transfer_holder_power: String,
    pub strict_margins: Vec<StrictEndpointMargin>,
    pub physical_endpoint_class: String,
}
#[derive(Clone, Debug, Serialize)]
pub struct PartialEndpointAudit {
    pub endpoint_identity: String,
    pub labels: Vec<PartialEndpointLabelAudit>,
    pub combined_rational_target_checked: bool,
    pub original_coefficient_conditions_checked: bool,
    pub native_mode_authority: bool,
}
#[derive(Clone, Debug, Serialize)]
struct LabelData {
    label: Vec<i64>,
    origin: FiniteLabelOriginAudit,
    p: u64,
    p0: u64,
    r: u64,
    jets: Vec<u64>,
    total_jets: u64,
    cuts: Vec<u64>,
    uppers: Vec<u64>,
    lowers: Vec<u64>,
    physical_class: String,
}
/// Privately populated proof/bound data. No public deserialization or constructor
/// from positivity strings. Every construction replays the actual native origin.
#[derive(Clone, Debug, Serialize)]
pub struct PartialEndpointEvidence {
    endpoint_version: &'static str,
    endpoint_order: &'static str,
    origin_identity: String,
    origin: PartialOriginCapability,
    original_target_label_bounds: Vec<PartialEndpointLabelAudit>,
    native_mode_authority: bool,
}
impl PartialEndpointEvidence {
    pub fn new(
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        origin: &PartialOriginCapability,
        budget: PartialEndpointBudget,
    ) -> Result<Self> {
        if budget.labels == 0
            || budget.total_label_axes == 0
            || budget.index_magnitude == 0
            || budget.degree_sum == 0
            || budget.dimension_integer_bound == 0
        {
            return Err(endpoint_limit("zero endpoint budget"));
        }
        origin.validate_binding(input, family, shifted, options)?;
        // Complete native/full-channel/chart replay is mandatory, even if all
        // initial labels vanish. No initial-target-only or positivity shortcut.
        origin.replay(input, family, budget.origin)?;
        let mut evidence = Self {
            endpoint_version: ENDPOINT_VERSION,
            endpoint_order: ENDPOINT_ORDER,
            origin_identity: origin.identity()?,
            origin: origin.clone(),
            original_target_label_bounds: Vec::new(),
            native_mode_authority: false,
        };
        let mut labels = Vec::new();
        for target in family.targets() {
            for label in target.keys() {
                if labels.len() >= budget.labels {
                    return Err(endpoint_limit("original target label budget"));
                }
                labels.push(label.0.iter().map(|x| i64::from(*x)).collect());
            }
        }
        evidence.original_target_label_bounds = evidence.audit_bounds(family, &labels, budget)?;
        Ok(evidence)
    }
    pub fn identity(&self) -> Result<String> {
        let bytes = serde_json::to_vec(self).map_err(|e| Error::Cache(e.to_string()))?;
        Ok(format!("{ENDPOINT_VERSION}:{}", blake3::hash(&bytes)))
    }
    pub fn origin_identity(&self) -> &str {
        &self.origin_identity
    }
    pub fn validate_binding(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
    ) -> Result<()> {
        self.origin
            .validate_binding(input, family, shifted, options)?;
        if self.endpoint_version != ENDPOINT_VERSION
            || self.endpoint_order != ENDPOINT_ORDER
            || self.origin_identity != self.origin.identity()?
            || self.native_mode_authority
        {
            return Err(endpoint_invalid(
                "endpoint/origin theorem identity mismatch",
            ));
        }
        Ok(())
    }
    pub fn replay(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        budget: PartialEndpointBudget,
    ) -> Result<()> {
        let rebuilt = Self::new(
            input,
            family,
            &self.origin.shifted,
            self.origin.source_options,
            &self.origin,
            budget,
        )?;
        if rebuilt.identity()? != self.identity()? {
            return Err(endpoint_invalid("endpoint replay mismatch"));
        }
        Ok(())
    }
    /// This deliberately does not authorize a transport projector. The actual
    /// ODE owner must check its regular-singular representation and joint target.
    pub fn consumer_requirements(&self) -> serde_json::Value {
        serde_json::json!({"endpoint_theorem":ENDPOINT_VERSION,"origin_identity":self.origin_identity,
            "audit_every_target_basis_candidate_and_source_application_label":true,
            "retain_raw_target_conditions_before_cancellation_or_zero_classification":true,
            "check_original_conditions_at_each_epsilon_sample":true,
            "audit_combined_rational_reconstruction_with_all_eta_poles":true,
            "require_actual_regular_singular_frobenius_representation_and_sufficient_series":true,
            "retain_native_joint_target_divergence_checks":true,
            "right_half_plane_eta_only":true,
            "native_mode_authority":false,"boundary_authority":false,"source_zero_authority":false,
            "singleton_raw_ward_authority":false,"production_admission":false,
            "cancellation_supported":false,"internal_cas_work_bound":false})
    }
    pub fn audit_labels(
        &self,
        family: &OccupiedCutFamily,
        labels: &[Vec<i64>],
        budget: PartialEndpointBudget,
    ) -> Result<PartialEndpointAudit> {
        let bounds = self.audit_bounds(family, labels, budget)?;
        Ok(PartialEndpointAudit {
            endpoint_identity: self.identity()?,
            labels: bounds,
            combined_rational_target_checked: false,
            original_coefficient_conditions_checked: false,
            native_mode_authority: false,
        })
    }
    fn audit_bounds(
        &self,
        family: &OccupiedCutFamily,
        labels: &[Vec<i64>],
        budget: PartialEndpointBudget,
    ) -> Result<Vec<PartialEndpointLabelAudit>> {
        if self.origin.family_signature != family_signature(family) {
            return Err(endpoint_invalid("finite label family binding"));
        }
        if budget.labels == 0
            || budget.total_label_axes == 0
            || budget.index_magnitude == 0
            || budget.degree_sum == 0
            || budget.dimension_integer_bound == 0
        {
            return Err(endpoint_limit("zero endpoint audit budget"));
        }
        if labels.len() > budget.labels {
            return Err(endpoint_limit("finite label count"));
        }
        let mut axes = 0usize;
        let mut result = Vec::new();
        for label in labels {
            axes = axes
                .checked_add(label.len())
                .ok_or_else(|| endpoint_limit("label axes overflow"))?;
            if axes > budget.total_label_axes
                || label
                    .iter()
                    .any(|x| x.unsigned_abs() > budget.index_magnitude)
            {
                return Err(endpoint_limit("label axes/index magnitude"));
            }
            let data = self.label_data(family, label, budget)?;
            let mut bound = mul(2, add(add(add(data.p, data.r)?, data.total_jets)?, 3)?)?;
            bound = bound.max(mul(2, add(add(data.p0, data.total_jets)?, 2)?)?);
            bound = bound.max(add(self.origin.compact_loops as u64, 2)?);
            for i in 0..data.jets.len() {
                bound = bound.max(add(
                    add(add(3, mul(2, data.jets[i])?)?, data.p0)?,
                    data.total_jets,
                )?);
                let lower = add(
                    add(
                        add(mul(2, data.cuts[i])?, data.lowers[i])?,
                        mul(2, data.uppers[i].saturating_sub(1))?,
                    )?,
                    add(add(data.p0, data.total_jets)?, 1)?,
                )?;
                bound = bound.max(lower);
            }
            if bound > budget.dimension_integer_bound {
                return Err(endpoint_limit("dimension witness budget"));
            }
            let numerator = add(mul(5, bound)?, 1)?;
            let d = &integer(numerator)? / &Rational::from(5);
            let n = if self.origin.virtual_loops == 2 {
                add(add(bound / 2, data.r)?, 1)?
            } else {
                0
            };
            result.push(self.bounds_at_dimension(data, &d, n)?);
        }
        Ok(result)
    }
    fn label_data(
        &self,
        family: &OccupiedCutFamily,
        label: &[i64],
        budget: PartialEndpointBudget,
    ) -> Result<LabelData> {
        // Origin validates all H/completion/storage roles before any zero branch.
        let origin = self.origin.audit_label(family, label)?;
        let mut p0 = 0;
        for &slot in &self.origin.pure_compact_slots {
            p0 = add(p0, label[slot].max(0) as u64)?;
        }
        let mut jets = Vec::new();
        let mut cuts = Vec::new();
        let mut uppers = Vec::new();
        let mut lowers = Vec::new();
        let mut total = 0;
        for shell in family.shells() {
            let n = label[shell.physical_slot].max(1) as u64;
            let s = label[shell.upper_slot] as u64;
            let ell = label[shell.lower_slot] as u64;
            let j = add(add(n - 1, s.saturating_sub(1))?, ell.saturating_sub(1))?;
            total = add(total, j)?;
            cuts.push(n);
            uppers.push(s);
            lowers.push(ell);
            jets.push(j);
        }
        let p = origin.positive_virtual_degree;
        let r = origin.polynomial_degree;
        if add(add(add(p, p0)?, r)?, total)? > budget.degree_sum {
            return Err(endpoint_limit("combined finite degree budget"));
        }
        let support = self
            .origin
            .complete_supports
            .iter()
            .find(|s| s.active == origin.active_virtual_slots)
            .ok_or_else(|| endpoint_invalid("missing origin support"))?;
        let physical_class = if origin.cut_vanishing {
            "required-cut-zero"
        } else if lowers.iter().any(|n| *n > 0) {
            "continued-lower-contact-zero"
        } else {
            match support.class {
                SupportClass::FreeVirtual { .. } => "ordinary-free-polynomial-zero",
                SupportClass::RegulatedVacuumZero => "regulated-vacuum-zero",
                SupportClass::FiniteEtaJets {
                    kinematics_vanish: true,
                    ..
                } => "vanishing-physical-endpoint-not-finite-eta-zero",
                SupportClass::FiniteEtaJets {
                    kinematics_vanish: false,
                    ..
                } => "continued-hard-limit-with-positive-holder-remainder",
            }
        }
        .to_owned();
        Ok(LabelData {
            label: label.to_vec(),
            origin,
            p,
            p0,
            r,
            jets,
            total_jets: total,
            cuts,
            uppers,
            lowers,
            physical_class,
        })
    }
    fn bounds_at_dimension(
        &self,
        data: LabelData,
        d: &Rational,
        n: u64,
    ) -> Result<PartialEndpointLabelAudit> {
        let sigma = Rational::from((1, 2));
        let q = &integer(add(data.p0, data.total_jets)?)? + &sigma;
        let mut margins = Vec::new();
        let virtual_power = &integer(self.origin.virtual_loops as u64)? * d / &Rational::from(2)
            - integer(add(add(data.p, data.total_jets)?, n)?)?;
        margins.push(positive_margin(
            "virtual finite J+N Holder remainder".into(),
            &virtual_power - &sigma,
        )?);
        if self.origin.virtual_loops == 2 {
            margins.push(positive_margin(
                "UV Taylor remainder integrability".into(),
                integer(n)? - d / &Rational::from(2) - integer(data.r)?,
            )?);
            margins.push(positive_margin(
                "finite Taylor upper-order window".into(),
                d / &Rational::from(2) + integer(add(data.r, 2)?)? - integer(n)?,
            )?);
        } else if n != 0 {
            return Err(endpoint_invalid(
                "one-loop chart has no D-growing Taylor subtraction",
            ));
        }
        margins.push(positive_margin(
            "pair-angle Holder integrability".into(),
            d / &Rational::from(2) - Rational::one() - &q,
        )?);
        margins.push(positive_margin(
            "positive compact Gram-angle measure".into(),
            d - &integer(add(self.origin.compact_loops as u64, 2)?)?,
        )?);
        for i in 0..data.jets.len() {
            margins.push(positive_margin(
                format!("compact radial difference loop {i}"),
                d - &integer(add(2, mul(2, data.jets[i])?)?)? - &q - &sigma,
            )?);
            let lower = add(
                add(mul(2, data.cuts[i])?, data.lowers[i])?,
                mul(2, data.uppers[i].saturating_sub(1))?,
            )?;
            margins.push(positive_margin(
                format!("compact lower-contact jets loop {i}"),
                d - &integer(lower)? - &q - &sigma,
            )?);
        }
        Ok(PartialEndpointLabelAudit {
            label: data.label,
            origin_audit: data.origin,
            positive_virtual_degree: data.p,
            positive_pure_transfer_degree: data.p0,
            polynomial_degree: data.r,
            compact_jet_orders: data.jets,
            total_jet_order: data.total_jets,
            dimension_witness: d.to_string(),
            taylor_terms: n,
            sigma: sigma.to_string(),
            transfer_holder_power: q.to_string(),
            strict_margins: margins,
            physical_endpoint_class: data.physical_class,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn input() -> PreparedDensityInput {
        serde_json::from_str::<crate::finite_density::DensityInput>(include_str!(
            "/common/dev/rustflow_fermi/examples/finite_density/massless_three_loop_chain.json"
        ))
        .unwrap()
        .prepare()
        .unwrap()
    }
    fn options() -> WeightedSourceOptions {
        WeightedSourceOptions {
            free_virtual_zero_sectors: true,
            ..Default::default()
        }
    }
    fn construct(
        cuts: &[usize],
        shifted: &[usize],
    ) -> (
        PreparedDensityInput,
        OccupiedCutFamily,
        PartialOriginCapability,
        PartialEndpointEvidence,
    ) {
        let i = input();
        let f = i.occupied_cut(cuts, 1024).unwrap().at_physical_masses();
        let o =
            PartialOriginCapability::new(&i, &f, shifted, options(), Default::default()).unwrap();
        let e = PartialEndpointEvidence::new(&i, &f, shifted, options(), &o, Default::default())
            .unwrap();
        (i, f, o, e)
    }
    fn label(f: &OccupiedCutFamily) -> Vec<i64> {
        let mut v = vec![0; 16];
        for x in &mut v[..f.physical_slots()] {
            *x = 1;
        }
        v
    }
    #[test]
    fn actual_singleton_and_double_original_labels_replay_and_keep_authorities_separate() {
        for (cuts, shifted) in [
            (vec![3], vec![0]),
            (vec![0, 3], vec![1, 2]),
            (vec![0, 3], vec![2, 4]),
        ] {
            let (i, f, o, e) = construct(&cuts, &shifted);
            e.replay(&i, &f, Default::default()).unwrap();
            assert_ne!(e.identity().unwrap(), o.identity().unwrap());
            assert!(!e.original_target_label_bounds.is_empty());
            assert_eq!(e.consumer_requirements()["native_mode_authority"], false);
            assert_eq!(e.consumer_requirements()["boundary_authority"], false);
            assert!(
                e.original_target_label_bounds
                    .iter()
                    .all(|a| !a.strict_margins.is_empty())
            );
            if let Ok(dir) = std::env::var("PARTIAL_ENDPOINT_REPORT") {
                std::fs::create_dir_all(&dir).unwrap();
                std::fs::write(
                    std::path::Path::new(&dir)
                        .join(format!("cuts-{cuts:?}-shift-{shifted:?}.json")),
                    serde_json::to_vec_pretty(
                        &serde_json::json!({"evidence":e,"requirements":e.consumer_requirements()}),
                    )
                    .unwrap(),
                )
                .unwrap();
            }
        }
    }
    #[test]
    fn every_finite_cut_upper_lower_and_transfer_jet_has_positive_recorded_margins() {
        let (_, f, _, e) = construct(&[0, 3], &[1, 2]);
        let mut labels = Vec::new();
        for n in [1, 3, 9] {
            for s in [0, 1, 7] {
                for ell in [0, 1, 5] {
                    let mut a = label(&f);
                    for shell in f.shells() {
                        a[shell.physical_slot] = n;
                        a[shell.upper_slot] = s;
                        a[shell.lower_slot] = ell;
                    }
                    a[2] = 4;
                    a[f.physical_slots()] = -11;
                    labels.push(a);
                }
            }
        }
        let report = e.audit_labels(&f, &labels, Default::default()).unwrap();
        assert_eq!(report.labels.len(), 27);
        for bound in &report.labels {
            assert_eq!(bound.positive_pure_transfer_degree, 4);
            assert_eq!(bound.compact_jet_orders.len(), 2);
            assert_eq!(bound.sigma, "1/2");
            assert!(bound.polynomial_degree >= 11);
            assert!(
                !report.native_mode_authority
                    && !report.combined_rational_target_checked
                    && !report.original_coefficient_conditions_checked
            );
            for m in &bound.strict_margins {
                assert!(
                    rat(&Atom::parse(&m.value, "endpoint_test", Default::default()).unwrap())
                        .unwrap()
                        > 0
                );
            }
        }
        if let Ok(dir) = std::env::var("PARTIAL_ENDPOINT_REPORT") {
            std::fs::write(
                std::path::Path::new(&dir).join("finite-label-audits.json"),
                serde_json::to_vec_pretty(&report).unwrap(),
            )
            .unwrap();
        }
    }
    #[test]
    fn wrong_origin_family_placement_and_options_reject() {
        let (i, f, o, e) = construct(&[3], &[0]);
        let other = i.occupied_cut(&[0], 1024).unwrap().at_physical_masses();
        assert!(
            PartialEndpointEvidence::new(&i, &other, &[0], options(), &o, Default::default())
                .is_err()
        );
        assert!(
            PartialEndpointEvidence::new(&i, &f, &[0, 1], options(), &o, Default::default())
                .is_err()
        );
        let opt = WeightedSourceOptions::default();
        assert!(PartialEndpointEvidence::new(&i, &f, &[0], opt, &o, Default::default()).is_err());
        let mut foreign = o.clone();
        foreign.input_identity = "foreign".into();
        assert!(
            PartialEndpointEvidence::new(&i, &f, &[0], options(), &foreign, Default::default())
                .is_err()
        );
        assert!(e.validate_binding(&i, &other, &[0], options()).is_err());
        assert!(e.audit_labels(&other, &[], Default::default()).is_err());
    }
    #[test]
    fn missing_charts_full_off_null_or_support_data_never_authorize_endpoint() {
        let (i, f, o, _) = construct(&[3], &[0]);
        let mut missing = o.clone();
        missing.complete_supports.pop();
        assert!(
            PartialEndpointEvidence::new(&i, &f, &[0], options(), &missing, Default::default())
                .is_err()
        );
        let mut missing = o.clone();
        missing
            .complete_supports
            .iter_mut()
            .find(|s| !s.charts.is_empty())
            .unwrap()
            .charts
            .clear();
        assert!(
            PartialEndpointEvidence::new(&i, &f, &[0], options(), &missing, Default::default())
                .is_err()
        );
        let mut missing = o.clone();
        missing
            .complete_supports
            .iter_mut()
            .find_map(|s| s.full.as_mut())
            .unwrap()
            .diagonal
            .clear();
        assert!(
            PartialEndpointEvidence::new(&i, &f, &[0], options(), &missing, Default::default())
                .is_err()
        );
    }
    #[test]
    fn role_inverse_and_storage_fail_even_for_cut_zero() {
        let (_, f, _, e) = construct(&[3], &[0]);
        let mut a = label(&f);
        a[3] = 0;
        a[f.shells()[0].lower_slot] = -1;
        assert!(
            e.audit_labels(&f, &[a.clone()], Default::default())
                .is_err()
        );
        a[f.shells()[0].lower_slot] = 0;
        a[f.physical_slots()] = 1;
        assert!(
            e.audit_labels(&f, &[a.clone()], Default::default())
                .is_err()
        );
        a[f.physical_slots()] = 0;
        a[15] = 1;
        assert!(e.audit_labels(&f, &[a], Default::default()).is_err());
    }
    #[test]
    fn endpoint_budgets_and_arithmetic_are_explicit_no_partial_result() {
        let (i, f, o, e) = construct(&[3], &[0]);
        let b = PartialEndpointBudget::default();
        assert!(matches!(
            PartialEndpointEvidence::new(
                &i,
                &f,
                &[0],
                options(),
                &o,
                PartialEndpointBudget {
                    origin: PartialOriginBudget {
                        supports: 1,
                        ..b.origin
                    },
                    ..b
                }
            ),
            Err(Error::Limit(_))
        ));
        let a = label(&f);
        for bad in [
            PartialEndpointBudget { labels: 0, ..b },
            PartialEndpointBudget {
                total_label_axes: 1,
                ..b
            },
            PartialEndpointBudget {
                index_magnitude: 0,
                ..b
            },
            PartialEndpointBudget { degree_sum: 0, ..b },
            PartialEndpointBudget {
                dimension_integer_bound: 1,
                ..b
            },
        ] {
            assert!(matches!(
                e.audit_labels(&f, &[a.clone()], bad),
                Err(Error::Limit(_))
            ));
        }
        let mut huge = a;
        huge[0] = i64::MAX;
        assert!(matches!(
            e.audit_labels(
                &f,
                &[huge],
                PartialEndpointBudget {
                    index_magnitude: u64::MAX,
                    ..b
                }
            ),
            Err(Error::Limit(_))
        ));
        assert!(add(u64::MAX, 1).is_err());
    }
    #[test]
    fn nonpositive_or_tampered_witness_is_rejected() {
        let (i, f, _, e) = construct(&[3], &[0]);
        let data = e.label_data(&f, &label(&f), Default::default()).unwrap();
        assert!(matches!(
            e.bounds_at_dimension(data, &Rational::one(), 0),
            Err(Error::InvalidInput(_))
        ));
        let mut altered = e.clone();
        altered.original_target_label_bounds[0].dimension_witness = "1".into();
        assert!(altered.replay(&i, &f, Default::default()).is_err());
        let mut altered = e.clone();
        altered.native_mode_authority = true;
        assert!(altered.validate_binding(&i, &f, &[0], options()).is_err());
    }
    #[test]
    fn shifted_v_zero_endpoint_is_distinct_from_finite_eta_zero() {
        let i = input();
        let f = i.occupied_cut(&[3], 1024).unwrap().at_physical_masses();
        let o =
            PartialOriginCapability::new(&i, &f, &[0, 1], options(), Default::default()).unwrap();
        let e = PartialEndpointEvidence::new(&i, &f, &[0, 1], options(), &o, Default::default())
            .unwrap();
        let mut a = label(&f);
        a[0] = 0;
        let audit = e
            .audit_labels(&f, &[a.clone()], Default::default())
            .unwrap();
        assert_eq!(
            audit.labels[0].physical_endpoint_class,
            "vanishing-physical-endpoint-not-finite-eta-zero"
        );
        assert_eq!(
            audit.labels[0].origin_audit.classification,
            "finite-positive-eta-jets"
        );
        let a: [i64; 16] = a.try_into().unwrap();
        assert!(
            !o.zero_domains::<16>(&f, &[0, 1], options())
                .unwrap()
                .iter()
                .any(|d| d.contains(&a))
        );
    }
}
