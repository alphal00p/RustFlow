//! Shared JSON request handling for the native CLI and embedding Python host.
//!
//! Answers are serialized only after the complete native operation succeeds.
//! A failed occupied sector is never replaced by a partial amplitude.
use super::assembly::{DensityEvaluation, PreparedDensityFlow};
use super::guarded::GuardedDiscoveryOptions;
use super::reduction::WeightedClosureOptions;
use super::{DensityInput, PreparedDensityInput};
use crate::{ComplexFloat, Error, FlowOptions, LaurentExpansion, MassMode, Result, RunContext};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub fn default_options(input: &DensityInput) -> FlowOptions {
    FlowOptions {
        digits: input.digits,
        mass_mode: MassMode::All,
        ..Default::default()
    }
}

pub fn default_closure_options() -> WeightedClosureOptions {
    WeightedClosureOptions {
        max_rounds: 12,
        discovery: GuardedDiscoveryOptions {
            max_depth: 3,
            max_domains: 8192,
            ..Default::default()
        },
        ..Default::default()
    }
}

/// Parse a nonzero real rational before any potentially expensive preparation.
/// Approximate decimal coefficients and symbolic expressions are rejected.
pub fn parse_epsilon(text: &str) -> Result<Rational> {
    let atom = Atom::parse(text, "rustflow_density", Default::default())
        .map_err(|e| Error::InvalidInput(format!("epsilon: {e}")))?;
    if let AtomView::Num(number) = atom.as_view()
        && let symbolica::coefficient::Coefficient::Complex(value) =
            number.get_coeff_view().to_owned()
        && value.im.is_zero()
    {
        if value.re.is_zero() {
            return Err(Error::InvalidInput(
                "epsilon sample must be nonzero; omit epsilon for Laurent reconstruction".into(),
            ));
        }
        return Ok(value.re);
    }
    Err(Error::InvalidInput(
        "epsilon sample must be an exact real rational such as 4/5".into(),
    ))
}

fn complex_json(value: &ComplexFloat) -> Value {
    json!({"real":value.re.as_raw().to_string(),"imaginary":value.im.as_raw().to_string()})
}

fn header(input: &PreparedDensityInput, options: &FlowOptions, operation: &str) -> Value {
    json!({
        "schema_version":1,
        "operation":operation,
        "name":input.input().name,
        "input_identity":input.identity(),
        "normalization":"unscaled Euclidean amplitude",
        "full_amplitude":true,
        "dimension_base":options.dimension,
        "requested_digits":options.digits,
        "guard_digits":options.guard_digits,
        "series_order":options.series_order,
    })
}

fn sample_report(
    input: &PreparedDensityInput,
    epsilon: &Rational,
    options: &FlowOptions,
    start_scale: u32,
    result: DensityEvaluation,
) -> Value {
    let mut output = header(input, options, "finite-density-sample");
    output["epsilon"] = epsilon.to_string().into();
    output["occupied_start_scale"] = start_scale.into();
    // A single regulator evaluation has no independent accuracy certificate.
    output["verified_digits"] = Value::Null;
    output["values"] = json!(result.values.iter().map(complex_json).collect::<Vec<_>>());
    output["targets"] = json!(input.input().targets);
    output["contributions"] = json!(
        result
            .contributions
            .iter()
            .map(|(cuts, values)| json!({
                "cut_slots":cuts,
                "values":values.iter().map(complex_json).collect::<Vec<_>>(),
            }))
            .collect::<Vec<_>>()
    );
    output["empty_support"] = json!(
        result
            .empty_support
            .iter()
            .map(|(cuts, proof)| json!({
                "cut_slots":cuts,"proof":proof,
            }))
            .collect::<Vec<_>>()
    );
    output["occupied_reports"] = json!(result.occupied_reports.iter().map(|(cuts, report)| json!({
        "cut_slots":cuts,
        "construction":report.construction,
        "source_options":report.source_options.as_ref().map(|source| json!({
            "policy":source.policy.as_str(),
            "positive_compact_energy_powers":source.positive_compact_energy_powers,
        })),
        "basis_size":report.basis_size,
        "shifted_slots":report.shifted_slots,
        "contour_admission":report.contour_admission,
        "nonzero_conditions":report.nonzero_conditions.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
        "boundary":{
            "integrated_coefficients":report.boundary.integrated_coefficients,
            "integrated_products":report.boundary.integrated_products,
            "scaleless_noncompact_products":report.boundary.scaleless_noncompact_products,
            "vanishing_required_cut_coefficients":report.boundary.vanishing_required_cut_coefficients,
            "regions":report.boundary.regions.iter().map(|region| json!({
                "hard":region.hard,
                "transformation":region.transformation.iter().map(|row| row.iter().map(Atom::to_canonical_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "jacobian_determinant":region.jacobian_determinant.to_canonical_string(),
                "exponents":region.exponents.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
                "half_orders":region.half_orders,
            })).collect::<Vec<_>>(),
        },
    })).collect::<Vec<_>>());
    output
}

fn laurent_report(
    input: &PreparedDensityInput,
    options: &FlowOptions,
    expansions: &[LaurentExpansion],
) -> Value {
    let mut output = header(input, options, "finite-density");
    output["laurent_orders"] = json!(input.input().laurent_orders);
    output["targets"] = json!(input.input().targets);
    output["expansions"] = json!(expansions.iter().map(|result| json!({
        "verified_digits":result.verified_digits,
        "working_bits":result.working_bits,
        "samples":result.samples,
        "validation_samples":result.validation_samples,
        "refinements":result.refinements,
        "coefficients":result.coefficients.iter().map(|(k,v)| (k.to_string(),complex_json(v))).collect::<BTreeMap<_,_>>(),
        "comparison_errors":result.comparison_errors.iter().map(|(k,v)| (k.to_string(),v.as_raw().to_string())).collect::<BTreeMap<_,_>>(),
    })).collect::<Vec<_>>());
    output
}

/// Evaluate all connected cut contributions through native AMF, or reconstruct
/// their summed Laurent expansion when `epsilon` is absent. This entry point
/// inherits the assembly owner's numerical admission and compiled arity limits.
pub fn evaluate_request(
    input: &DensityInput,
    epsilon: Option<&str>,
    options: &FlowOptions,
    closure: WeightedClosureOptions,
    context: &RunContext,
    start_scale: u32,
) -> Result<Value> {
    let epsilon = epsilon.map(parse_epsilon).transpose()?;
    if start_scale < 4 {
        return Err(Error::InvalidInput(
            "occupied start scale must be at least four".into(),
        ));
    }
    options.validate()?;
    let prepared = PreparedDensityFlow::prepare(input, options, closure, context)?;
    let mut report = if let Some(epsilon) = epsilon {
        let result = prepared.evaluate(&epsilon, options, context, start_scale)?;
        sample_report(prepared.input(), &epsilon, options, start_scale, result)
    } else {
        let result = prepared.solve_with_start_scale(options, context, start_scale)?;
        let mut report = laurent_report(prepared.input(), options, &result);
        report["occupied_start_scale"] = start_scale.into();
        report
    };
    report["closure"] = json!(
        prepared
            .closure_diagnostics()
            .iter()
            .map(|(cuts, diagnostics)| json!({
                "cut_slots":cuts,
                "rounds":diagnostics.rounds,
                "requested":diagnostics.requested,
                "native_frontier_requests":diagnostics.native_frontier_requests,
                "provisional_sizes":diagnostics.provisional_sizes,
                "native_rules":diagnostics.native_rules,
                "native_rule_applications":diagnostics.native_rule_applications,
                "uncovered_discovery_domains":diagnostics.uncovered_discovery_domains,
                "guard_refinement_passes":diagnostics.guard_refinement_passes,
                "guard_refinement_added_domains":diagnostics.guard_refinement_added_domains,
                "guard_refinement_budget_exhausted":diagnostics.guard_refinement_budget_exhausted,
                "guard_refinements":diagnostics.guard_refinements,
            }))
            .collect::<Vec<_>>()
    );
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Precision;

    fn input() -> DensityInput {
        serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap()
    }

    #[test]
    fn epsilon_is_exact_and_rejected_before_graph_preparation() {
        assert_eq!(parse_epsilon("4/5").unwrap(), Rational::from((4, 5)));
        assert_eq!(parse_epsilon("-1/7").unwrap(), Rational::from((-1, 7)));
        for bad in ["0", "0.8", "eps", "𝑖", "1/0"] {
            assert!(parse_epsilon(bad).is_err(), "accepted {bad}");
        }
        let mut invalid = input();
        invalid.loops = 0;
        let error = evaluate_request(
            &invalid,
            Some("0"),
            &default_options(&invalid),
            default_closure_options(),
            &RunContext::default(),
            8,
        )
        .unwrap_err();
        assert!(error.to_string().contains("epsilon sample must be nonzero"));
    }

    #[test]
    fn sample_json_preserves_target_order_and_native_decimal_precision() {
        let input = input().prepare().unwrap();
        let p = Precision::decimal(70).unwrap();
        let value = p.add(&p.rational(&Rational::from((1, 7))), &p.complex(0, 2));
        let result = DensityEvaluation {
            contributions: vec![
                (vec![], vec![p.zero(), p.i(3)]),
                (vec![0], vec![value.clone(), p.zero()]),
            ],
            values: vec![value.clone(), p.i(3)],
            occupied_reports: vec![],
            empty_support: vec![(vec![1], "strictly disjoint occupied support".into())],
        };
        let report = sample_report(
            &input,
            &Rational::from((4, 5)),
            &default_options(input.input()),
            12,
            result,
        );
        assert_eq!(report["values"][0]["real"], value.re.as_raw().to_string());
        assert_eq!(
            report["values"][0]["imaginary"],
            value.im.as_raw().to_string()
        );
        assert_eq!(report["values"][1], complex_json(&p.i(3)));
        assert!(report["values"][0]["real"].as_str().unwrap().len() > 60);
        assert_eq!(report["contributions"][1]["cut_slots"], json!([0]));
        assert_eq!(report["targets"][1]["powers"], json!([2, 1, 1]));
        assert!(report["verified_digits"].is_null());
        assert_eq!(report["empty_support"][0]["cut_slots"], json!([1]));
        assert_eq!(report["normalization"], "unscaled Euclidean amplitude");
    }

    #[test]
    fn laurent_json_retains_accuracy_evidence_and_signed_powers() {
        let input = input().prepare().unwrap();
        let p = Precision::decimal(60).unwrap();
        let expansion = LaurentExpansion {
            coefficients: BTreeMap::from([(-2, p.i(3)), (0, p.i(-5))]),
            verified_digits: Some(12),
            working_bits: p.bits,
            samples: 10,
            validation_samples: 14,
            refinements: 2,
            comparison_errors: BTreeMap::from([(-2, p.real(0))]),
        };
        let report = laurent_report(&input, &default_options(input.input()), &[expansion]);
        assert_eq!(
            report["expansions"][0]["coefficients"]["-2"],
            complex_json(&p.i(3))
        );
        assert_eq!(report["expansions"][0]["verified_digits"], 12);
        assert_eq!(report["expansions"][0]["validation_samples"], 14);
        assert!(report["expansions"][0]["comparison_errors"]["-2"].is_string());
    }
}
