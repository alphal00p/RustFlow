//! Physical singleton zeros from the uniform massless finite-jet theorem.
//!
//! The proof uses the all-uncut positive-eta homotopy, but evaluates no period
//! and supplies neither finite-eta zero rules nor a closed connection. Only
//! original, eta-independent physical targets can enter this owner.
use super::PreparedDensityInput;
use super::flow::{OccupiedFlowEvaluation, validate_options};
use super::flow_boundary::OccupiedBoundaryProvenance;
use super::geometry::OccupiedCutFamily;
use super::guarded::IndexRole;
use super::massless_endpoint::{MasslessFlowEvidence, MasslessLabelAudit};
use super::preparation::WeightedSourceOptions;
use crate::{Error, FlowOptions, Precision, Result, RunContext};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

const CONSTRUCTION: &str = "certified_massless_singleton_physical_zero";
const PROOF: &str = "original-eta-independent-singleton-endpoint-zero-v1";

/// A sealed zero for the original targets on their retained coefficient domain.
/// Other cuts, massive lines and inverse completions cannot obtain this result.
pub struct PreparedSingletonZero {
    family: OccupiedCutFamily,
    evidence: MasslessFlowEvidence,
    audit: MasslessLabelAudit,
    nonzero_conditions: Vec<Atom>,
    dimension: i64,
    epsilon: Symbol,
}

impl PreparedSingletonZero {
    /// Try the proved singleton class. Unsupported geometry continues through
    /// the usual occupied owner; malformed inputs and resource errors propagate.
    pub fn prepare(
        input: &PreparedDensityInput,
        cuts: &[usize],
        options: &FlowOptions,
        source_options: WeightedSourceOptions,
        context: &RunContext,
    ) -> Result<Option<Self>> {
        validate_options(options)?;
        context.cancellation.check()?;
        if cuts.len() != 1
            || input.input().loops < 2
            || input.physical_masses().iter().any(|mass| !mass.is_zero())
            || source_options.positive_compact_energy_powers
        {
            return Ok(None);
        }
        let family = input.occupied_cut(cuts, 65536)?.at_physical_masses();
        if input.input().targets.len() > 65536
            || input.targets().iter().map(|t| t.len()).sum::<usize>() > 65536
            || family.targets().iter().map(|t| t.len()).sum::<usize>() > 65536
        {
            return Err(Error::Limit(
                "singleton zero target-label budget exceeded".into(),
            ));
        }
        let epsilon = symbol!("rustflow_occupied::epsilon");
        // Capture raw poles before a zero classification, including terms
        // removed by cutting or by physical mass assignment. Original raised
        // shells retain C_n times H_0: their on-shell action includes support
        // derivatives, so adding explicit surface terms here would double count.
        let mut expressions = input
            .targets()
            .iter()
            .flat_map(|target| target.values().cloned())
            .collect::<Vec<_>>();
        for target in &input.input().targets {
            context.cancellation.check()?;
            expressions.push(
                Atom::parse(&target.numerator, "rustflow_density", Default::default())
                    .map_err(|error| Error::InvalidInput(error.to_string()))?,
            );
        }
        let variables = input
            .basis()
            .coordinates()
            .iter()
            .flat_map(|a| a.get_all_symbols(true))
            .chain(input.independent_masses().iter().copied())
            .chain([epsilon])
            .collect::<BTreeSet<_>>();
        let assignments = input
            .independent_masses()
            .iter()
            .zip(input.physical_masses())
            .map(|(&mass, value)| (Atom::var(mass), value.clone()))
            .collect::<BTreeMap<_, _>>();
        let mut conditions =
            crate::physical_conditions::rational_denominator_conditions(&expressions, &variables)?
                .iter()
                .map(|condition| crate::family::substitute(condition, &assignments))
                .collect::<Vec<_>>();
        let coefficients = family
            .targets()
            .iter()
            .flat_map(|target| target.values().cloned())
            .collect::<Vec<_>>();
        // Only epsilon is an unassigned scalar. In particular, no eta-pole
        // weight or externally supplied reduced combination can enter here.
        let variables = BTreeSet::from([epsilon]);
        conditions.extend(crate::physical_conditions::rational_denominator_conditions(
            &coefficients,
            &variables,
        )?);
        let nonzero_conditions =
            crate::physical_conditions::canonical_conditions(&conditions, &variables)?;
        context.cancellation.check()?;
        let shifted = (0..family.physical_slots())
            .filter(|&slot| family.roles()[slot] == IndexRole::Ordinary)
            .collect::<Vec<_>>();
        let evidence = match MasslessFlowEvidence::new(input, &family, &shifted, source_options) {
            Ok(evidence) => evidence,
            Err(Error::Unsupported(_)) => return Ok(None),
            Err(error) => return Err(error),
        };
        let audit = evidence.singleton_physical_zero_audit(&family)?;
        context.cancellation.check()?;
        Ok(Some(Self {
            family,
            evidence,
            audit,
            nonzero_conditions,
            dimension: options.dimension,
            epsilon,
        }))
    }

    /// The homotopy and finite-jet audit are proof data, not a solved flow.
    pub fn report(&self) -> serde_json::Value {
        serde_json::json!({
            "construction": CONSTRUCTION,
            "proof_version": PROOF,
            "scope": "original eta-independent physical singleton targets only",
            "native_closure": false,
            "transport_performed": false,
            "input_identity": self.evidence.input_identity(),
            "family_signature": super::massless_endpoint::bound_family_signature(&self.family),
            "proof_homotopy_shifted_slots": self.evidence.shifted_slots(),
            "proof_homotopy_identity": self.evidence.source_identity(),
            "normalization": "zero in the full Euclidean convention; original Wick target coefficients and routing determinant retained in family binding",
            "label_audit": self.audit,
            "nonzero_conditions": self.nonzero_conditions.iter()
                .map(Atom::to_canonical_string).collect::<Vec<_>>(),
        })
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<OccupiedFlowEvaluation> {
        validate_options(options)?;
        context.cancellation.check()?;
        if epsilon.is_zero() || options.dimension != self.dimension {
            return Err(Error::InvalidInput(
                "singleton zero requires nonzero epsilon and the prepared dimension".into(),
            ));
        }
        self.evidence
            .validate_family(&self.family, self.evidence.shifted_slots())?;
        crate::physical_conditions::validate_conditions_at(
            &self.nonzero_conditions,
            self.epsilon,
            &BTreeMap::from([(self.epsilon, Atom::num(epsilon.clone()))]),
        )?;
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("singleton zero precision overflow".into()))?,
        )?;
        Ok(OccupiedFlowEvaluation {
            construction: CONSTRUCTION,
            values: vec![p.zero(); self.family.targets().len()],
            boundary: OccupiedBoundaryProvenance {
                massless_origin: Some(self.evidence.origin_identity()),
                ..Default::default()
            },
            nonzero_conditions: self.nonzero_conditions.clone(),
            contour_admission: format!(
                "{PROOF}; {}; ordered fixed-T/positive-eta regulator removal, real thermal limit, high-D eta endpoint, then meromorphic D continuation; no finite-eta zero or native closure assertion",
                self.evidence.source_identity()
            ),
            shifted_slots: vec![],
            basis_size: 0,
            physical_arity: self.family.factors().len(),
            native_storage_capacity: None,
            source_options: None,
            massless_endpoint: Some(self.audit.clone()),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::{DensityInput, DensityTarget};

    fn input() -> DensityInput {
        serde_json::from_str(include_str!(
            "../../examples/finite_density/massless_three_loop_chain.json"
        ))
        .unwrap()
    }
    fn prepare(input: &DensityInput, cut: usize) -> PreparedSingletonZero {
        PreparedSingletonZero::prepare(
            &input.prepare().unwrap(),
            &[cut],
            &FlowOptions::default(),
            WeightedSourceOptions::default(),
            &RunContext::default(),
        )
        .unwrap()
        .unwrap()
    }

    #[test]
    fn original_raised_medium_targets_retain_independent_energy_jets_and_open_witnesses() {
        let mut input = input();
        input.targets.push(DensityTarget {
            powers: vec![1, 1, 1, 3, 1],
            numerator: "(g1_2+u1*u2)^2".into(),
        });
        for cut in [0, 3] {
            let prepared = prepare(&input, cut);
            let shell = &prepared.family.shells()[0];
            assert!(
                prepared
                    .audit
                    .labels
                    .iter()
                    .any(|label| label.indices[shell.physical_slot] > 1)
            );
            assert!(prepared.audit.labels.iter().any(|label| {
                label.indices[prepared.family.physical_slots()..prepared.family.input_slots()]
                    .iter()
                    .any(|n| *n < 0)
            }));
            for label in &prepared.audit.labels {
                // Geometry retains the independent-energy distributions. The
                // C_n action produces its Fermi surface after energy pairing;
                // the germ bound must already count all n-1 mass/energy jets.
                assert_eq!(label.indices[shell.upper_slot], 0);
                assert_eq!(label.indices[shell.lower_slot], 0);
                assert_eq!(
                    label.mass_and_upper_jet_orders,
                    [i32::from(label.indices[shell.physical_slot] - 1)]
                );
                let d = Rational::try_from(
                    Atom::parse(
                        &label.high_dimension_witness,
                        "singleton_zero_test",
                        Default::default(),
                    )
                    .unwrap()
                    .as_view(),
                )
                .unwrap();
                assert!(
                    label
                        .inequalities
                        .iter()
                        .all(|bound| bound.degree.at(&d) > 0)
                );
            }
            let result = prepared
                .evaluate(
                    &Rational::from((1, 7)),
                    &FlowOptions::default(),
                    &RunContext::default(),
                )
                .unwrap();
            assert_eq!(result.values.len(), input.targets.len());
            assert!(
                result
                    .values
                    .iter()
                    .all(|value| value.re.is_zero() && value.im.is_zero())
            );
            assert_eq!(result.construction, CONSTRUCTION);
            assert!(result.source_options.is_none() && result.native_storage_capacity.is_none());
            assert_eq!(prepared.report()["native_closure"], false);
            assert_eq!(prepared.report()["transport_performed"], false);
        }
    }

    #[test]
    fn coefficient_poles_survive_cancellation_and_a_removed_cut_target() {
        let mut input = input();
        input.targets = vec![
            DensityTarget {
                powers: vec![0, 1, 1, 1, 1],
                numerator: "(rustflow_occupied::epsilon^2-1)/(rustflow_occupied::epsilon-1)".into(),
            },
            DensityTarget {
                powers: vec![1, 1, 1, 1, 1],
                numerator: "(rustflow_occupied::epsilon^2-1)/(rustflow_occupied::epsilon-1)".into(),
            },
        ];
        let prepared = prepare(&input, 0);
        assert!(prepared.family.targets()[0].is_empty());
        assert!(!prepared.family.targets()[1].is_empty());
        assert!(!prepared.nonzero_conditions.is_empty());
        assert!(
            matches!(prepared.evaluate(&Rational::one(), &FlowOptions::default(), &RunContext::default()),
            Err(Error::InvalidInput(message)) if message.contains("nonzero reduction condition"))
        );
        prepared
            .evaluate(
                &Rational::from((1, 2)),
                &FlowOptions::default(),
                &RunContext::default(),
            )
            .unwrap();
    }

    #[test]
    fn unrelated_contours_and_inverse_completion_options_do_not_receive_zero_authority() {
        let mut input = input();
        let opts = FlowOptions::default();
        let context = RunContext::default();
        let prepared = input.prepare().unwrap();
        assert!(
            PreparedSingletonZero::prepare(&prepared, &[0, 3], &opts, Default::default(), &context)
                .unwrap()
                .is_none()
        );
        assert!(
            PreparedSingletonZero::prepare(
                &prepared,
                &[0],
                &opts,
                WeightedSourceOptions {
                    positive_compact_energy_powers: true,
                    ..Default::default()
                },
                &context
            )
            .unwrap()
            .is_none()
        );
        input.edges[1].mass_squared = "1".into();
        assert!(
            PreparedSingletonZero::prepare(
                &input.prepare().unwrap(),
                &[0],
                &opts,
                Default::default(),
                &context
            )
            .unwrap()
            .is_none()
        );
        // A typed rank-one-block permit is not a singleton zero certificate.
        let family = prepared
            .occupied_cut(&[0, 3], 16)
            .unwrap()
            .at_physical_masses();
        let block =
            MasslessFlowEvidence::new(&prepared, &family, &[1, 2, 4], Default::default()).unwrap();
        assert!(block.singleton_physical_zero_audit(&family).is_err());
        let mut external: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut external.edges {
            edge.mass_squared = "0".into();
        }
        external.edges[1].routing = vec!["1".into(), "0".into()];
        external.edges[2].vertices = [0, 0];
        external.edges[2].routing = vec!["0".into(), "1".into()];
        external.loop_charges[1][0] = 0;
        assert!(
            PreparedSingletonZero::prepare(
                &external.prepare().unwrap(),
                &[0],
                &opts,
                Default::default(),
                &context
            )
            .unwrap()
            .is_none()
        );
    }

    #[test]
    fn original_auxiliary_or_unassigned_coefficients_reject_before_cut_erasure() {
        for numerator in [
            "rustflow_occupied::eta",
            "1/rustflow_occupied::eta",
            "unassigned_coefficient",
        ] {
            let mut input = input();
            input.targets = vec![DensityTarget {
                powers: vec![0, 1, 1, 1, 1],
                numerator: numerator.into(),
            }];
            assert!(matches!(
                PreparedSingletonZero::prepare(
                    &input.prepare().unwrap(),
                    &[0],
                    &FlowOptions::default(),
                    Default::default(),
                    &RunContext::default()
                ),
                Err(Error::Unsupported(_))
            ));
        }
    }

    #[test]
    fn zero_evaluation_preserves_cancellation_and_prepared_dimension_contract() {
        let definition = input();
        let prepared = prepare(&definition, 3);
        let context = RunContext::default();
        context.cancellation.cancel();
        assert!(matches!(
            PreparedSingletonZero::prepare(
                &definition.prepare().unwrap(),
                &[3],
                &FlowOptions::default(),
                Default::default(),
                &context
            ),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            prepared.evaluate(&Rational::from((1, 7)), &FlowOptions::default(), &context),
            Err(Error::Cancelled)
        ));
        let options = FlowOptions {
            dimension: 6,
            ..Default::default()
        };
        assert!(
            prepared
                .evaluate(&Rational::from((1, 7)), &options, &RunContext::default())
                .is_err()
        );
        assert!(
            prepared
                .evaluate(
                    &Rational::zero(),
                    &FlowOptions::default(),
                    &RunContext::default()
                )
                .is_err()
        );
    }
}
