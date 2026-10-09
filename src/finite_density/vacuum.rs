//! The zero-cut contribution uses the ordinary native AMF owner.
use super::PreparedDensityInput;
use super::normalization::native_measure_to_euclidean;
use crate::family::{Integral, IntegralFamily, LinearCombination, substitute};
use crate::recursive::{RecursiveBoundary, RecursiveTerminalPolicy};
use crate::{
    ComplexFloat, Error, FlowOptions, KinematicPoint, Precision, PreparedFlow, Result, RunContext,
    RustRedBackend,
};
use std::collections::BTreeMap;
use symbolica::prelude::*;

/// Native exact proof for one undeformed, unoccupied vacuum sector.
/// The certificate is valid only on its retained nonzero coefficient domain.
#[derive(Debug)]
pub struct VacuumZeroCertificate {
    integral: Integral,
    native: rustred::sector::zero::Certificate,
    nonzero_conditions: Vec<Atom>,
}

impl VacuumZeroCertificate {
    pub fn integral(&self) -> &Integral {
        &self.integral
    }
    pub fn native(&self) -> &rustred::sector::zero::Certificate {
        &self.native
    }
    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.nonzero_conditions
    }

    /// Exact native proof metadata for the undeformed zero-cut contribution.
    /// Kernel entries are decimal strings so arbitrarily large integers remain exact.
    pub fn report(&self) -> serde_json::Value {
        serde_json::json!({
            "scope": "undeformed zero-cut vacuum; native ProvedZero certificate",
            "integral": self.integral.0,
            "family_fingerprint": self.native.family_fingerprint(),
            "raw_sector": self.native.raw_sector().active_bits(),
            "effective_sector": self.native.effective_sector().active_bits(),
            "active_parameter_order": self.native.active_parameter_order(),
            "primitive_kernel": self.native.primitive_kernel().iter().map(ToString::to_string).collect::<Vec<_>>(),
            "rank": self.native.rank(),
            "exponent_row_count": self.native.exponent_row_count(),
            "nonzero_conditions": self.nonzero_conditions.iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>(),
        })
    }
}

fn undeformed_zero_certificates(
    family: &IntegralFamily,
    requested: &[Integral],
    context: &RunContext,
) -> Result<Vec<VacuumZeroCertificate>> {
    use rustred::sector::{
        Mask,
        zero::{Analyzer, Decision},
    };
    if requested.is_empty() {
        return Ok(Vec::new());
    }
    context.cancellation.check()?;
    let converted = family.convert()?;
    let analyzer = Analyzer::try_unrestricted(&converted.family)
        .map_err(|error| Error::Reduction(error.to_string()))?;
    let mut proofs = Vec::new();
    for integral in requested {
        context.cancellation.check()?;
        family.validate_integral(integral)?;
        let sector = Mask::try_new(integral.0.iter().map(|&power| power > 0))
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        if let Decision::ProvedZero(native) = analyzer
            .analyze(&sector)
            .map_err(|error| Error::Reduction(error.to_string()))?
        {
            let nonzero_conditions = native
                .domain()
                .conditions()
                .iter()
                .map(|condition| {
                    substitute(&condition.polynomial().to_expression(), &converted.reverse)
                })
                .collect();
            proofs.push(VacuumZeroCertificate {
                integral: integral.clone(),
                native,
                nonzero_conditions,
            });
        }
    }
    Ok(proofs)
}

pub struct PreparedDensityVacuum {
    flow: Option<PreparedFlow>,
    loops: usize,
    dimension: i64,
    targets: usize,
    target_conditions: Vec<Atom>,
    zero_certificates: Vec<VacuumZeroCertificate>,
}

impl PreparedDensityVacuum {
    pub fn prepare(
        input: &PreparedDensityInput,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        let geometry = input.continued_cut(&[], 1)?.at_physical_masses();
        let family =
            geometry.region_family(symbol!("rustflow_occupied::epsilon"), options.dimension)?;
        let mut target_conditions = crate::physical_conditions::rational_denominator_conditions(
            &geometry
                .targets()
                .iter()
                .flat_map(|target| target.values().cloned())
                .collect::<Vec<_>>(),
            &std::collections::BTreeSet::from([family.epsilon]),
        )?;
        let requested = geometry
            .targets()
            .iter()
            .flat_map(|t| t.keys().cloned())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        // Analyze the actual zero-cut physical family BEFORE introducing eta.
        // Occupied sectors have different support and never pass this path.
        // Keep raw target poles even when every term is proved zero.
        let zero_certificates = undeformed_zero_certificates(&family, &requested, context)?;
        let zeros = zero_certificates
            .iter()
            .map(|proof| proof.integral.clone())
            .collect::<std::collections::BTreeSet<_>>();
        for proof in &zero_certificates {
            target_conditions.extend(proof.nonzero_conditions.clone());
        }
        target_conditions = crate::physical_conditions::canonical_conditions(
            &target_conditions,
            &std::collections::BTreeSet::from([family.epsilon]),
        )?;
        let requested = requested
            .into_iter()
            .filter(|integral| !zeros.contains(integral))
            .collect::<Vec<_>>();
        let flow = if requested.is_empty() {
            None
        } else {
            let backend = RustRedBackend {
                bubble_subloops: false,
                ..Default::default()
            };
            let mut flow = PreparedFlow::new(
                &family,
                &requested,
                &KinematicPoint::default(),
                &backend,
                options,
                context,
            )?;
            let indexed = requested
                .iter()
                .cloned()
                .zip(flow.reduced.targets.clone())
                .collect::<BTreeMap<_, _>>();
            // Compose the original fixed polynomial coefficients with the
            // native reductions BEFORE rational endpoint projection. This
            // retains cancellation between singular reduction weights.
            flow.reduced.targets = geometry
                .targets()
                .iter()
                .map(|target| {
                    let mut out = LinearCombination::new();
                    for (integral, coefficient) in target {
                        if zeros.contains(integral) {
                            continue;
                        }
                        for (basis, value) in &indexed[integral] {
                            *out.entry(basis.clone()).or_default() += coefficient * value;
                        }
                    }
                    for value in out.values_mut() {
                        *value = value.together().cancel();
                    }
                    out.retain(|_, c| !c.is_zero());
                    out
                })
                .collect();
            flow.reduced
                .nonzero_conditions
                .extend(target_conditions.clone());
            Some(flow)
        };
        Ok(Self {
            flow,
            loops: geometry.loops(),
            dimension: options.dimension,
            targets: geometry.targets().len(),
            target_conditions,
            zero_certificates,
        })
    }

    /// Inspect proof-bearing native zero sectors retained before deformation.
    pub fn zero_certificates(&self) -> &[VacuumZeroCertificate] {
        &self.zero_certificates
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        options.validate()?;
        if options.dimension != self.dimension {
            return Err(Error::InvalidInput(
                "zero-cut evaluation dimension changed after preparation".into(),
            ));
        }
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("vacuum precision overflow".into()))?,
        )?;
        crate::physical_conditions::validate_conditions_at(
            &self.target_conditions,
            symbol!("rustflow_occupied::epsilon"),
            &BTreeMap::from([(
                symbol!("rustflow_occupied::epsilon"),
                Atom::num(epsilon.clone()),
            )]),
        )?;
        let Some(flow) = &self.flow else {
            return Ok(vec![p.zero(); self.targets]);
        };
        let backend = RustRedBackend {
            bubble_subloops: false,
            ..Default::default()
        };
        let boundary = RecursiveBoundary::new(&backend, options, context)
            .with_terminal_policy(RecursiveTerminalPolicy::TadpolesOnly);
        let native = flow.evaluate(epsilon, options, &boundary, context)?;
        let dimension = Atom::num(self.dimension) - Atom::num(2) * Atom::num(epsilon.clone());
        let measure = p.eval(
            &native_measure_to_euclidean(self.loops, 0, &dimension)?,
            &ahash::HashMap::default(),
        )?;
        Ok(native.iter().map(|v| p.mul(v, &measure)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    #[test]
    fn zero_vacuum_projection_retains_original_coefficient_poles() {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_one_loop_tadpole.json"
        ))
        .unwrap();
        input.targets.truncate(1);
        input.targets[0].numerator = "u1/(rustflow_occupied::epsilon-1/2)".into();
        let input = input.prepare().unwrap();
        let options = FlowOptions::default();
        let context = RunContext::default();
        let prepared = PreparedDensityVacuum::prepare(&input, &options, &context).unwrap();
        assert!(
            prepared
                .evaluate(&Rational::from((1, 2)), &options, &context)
                .is_err()
        );
        assert!(
            prepared
                .evaluate(&Rational::from((1, 3)), &options, &context)
                .unwrap()[0]
                .is_zero()
        );
    }
    #[test]
    fn massless_four_loop_vacuum_is_proved_zero_before_amf_and_keeps_poles() {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap();
        input.targets.truncate(1);
        input.targets[0].numerator = "1/(rustflow_occupied::epsilon-1/2)".into();
        let input = input.prepare().unwrap();
        let options = FlowOptions::default();
        let context = RunContext::default();
        let prepared = PreparedDensityVacuum::prepare(&input, &options, &context).unwrap();
        assert!(prepared.flow.is_none());
        assert!(!prepared.zero_certificates().is_empty());
        assert!(
            prepared
                .zero_certificates()
                .iter()
                .all(|p| !p.native.family_fingerprint().is_empty())
        );
        assert!(
            prepared
                .evaluate(&Rational::from((1, 2)), &options, &context)
                .is_err()
        );
        assert!(
            prepared
                .evaluate(&Rational::from((1, 3)), &options, &context)
                .unwrap()[0]
                .is_zero()
        );
    }

    #[test]
    fn mixed_vacuum_requests_retain_nonzero_flow_and_original_target_order() {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_one_loop_tadpole.json"
        ))
        .unwrap();
        input.targets[0].powers[0] = 0;
        let options = FlowOptions::default();
        let context = RunContext::default();
        let prepared =
            PreparedDensityVacuum::prepare(&input.prepare().unwrap(), &options, &context).unwrap();
        let flow = prepared.flow.as_ref().unwrap();
        assert_eq!(prepared.zero_certificates().len(), 1);
        assert_eq!(flow.reduced.targets.len(), 2);
        assert!(flow.reduced.targets[0].is_empty());
        assert!(!flow.reduced.targets[1].is_empty());
    }

    #[test]
    fn native_zero_certificate_retains_symbolic_family_domain() {
        let scale = symbol!("vacuum_zero_test::scale");
        let epsilon = symbol!("vacuum_zero_test::epsilon");
        let family = IntegralFamily {
            name: "scaled massless vacuum".into(),
            loops: vec!["k".into()],
            external: vec![],
            external_gram: vec![],
            propagators: vec![crate::Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::one() / Atom::var(scale)],
            }],
            physical_propagators: 1,
            epsilon,
            dimension: 4,
        };
        let proofs =
            undeformed_zero_certificates(&family, &[Integral(vec![1])], &RunContext::default())
                .unwrap();
        assert_eq!(proofs.len(), 1);
        let conditions = &proofs[0].nonzero_conditions;
        assert!(!conditions.is_empty());
        assert!(
            crate::physical_conditions::validate_conditions_at(
                conditions,
                epsilon,
                &BTreeMap::from([(scale, Atom::new())]),
            )
            .is_err()
        );
        assert!(
            crate::physical_conditions::validate_conditions_at(
                conditions,
                epsilon,
                &BTreeMap::from([(scale, Atom::one())]),
            )
            .is_ok()
        );
    }
}
