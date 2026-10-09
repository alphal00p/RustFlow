//! Draft for the next test batch: public assembly with two independent charges.
use ahash::HashMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput, assembly::PreparedDensityFlow, reduction::WeightedClosureOptions,
};
use symbolica_amflow::*;

#[test]
fn independent_chemical_potentials_assemble_scalar_and_medium_numerators() {
    let input: DensityInput = serde_json::from_value(serde_json::json!({
        "name": "two_independent_compact_cycles",
        "loops": 2,
        "vertices": 1,
        "edges": [
            {"vertices": [0,0], "routing": ["1","0"], "mass_squared": "1/4", "charges": [1,0]},
            {"vertices": [0,0], "routing": ["0","1"], "mass_squared": "1/16", "charges": [0,1]}
        ],
        "loop_charges": [[1,0],[0,1]],
        "chemical_potentials": ["1","3/4"],
        "targets": [
            {"powers": [1,1], "numerator": "1"},
            {"powers": [1,1], "numerator": "u1*u2"}
        ],
        "numerator_convention": "shifted_euclidean",
        "laurent_orders": [-2,0],
        "digits": 18
    }))
    .unwrap();
    let options = FlowOptions {
        digits: 18,
        guard_digits: 24,
        series_order: 60,
        ..Default::default()
    };
    let context = RunContext::default();
    let prepared = PreparedDensityFlow::prepare(
        &input,
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    // The two singleton cuts have a genuine uncut virtual tadpole and use the
    // native source-generated weighted connection; only the double cut is a
    // direct compact polynomial terminal.
    assert_eq!(prepared.closure_diagnostics().len(), 2);
    let report = prepared
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    let exact = |expression: &str| {
        p.eval(
            &Atom::parse(expression, "multiple_mu_test", Default::default()).unwrap(),
            &HashMap::default(),
        )
        .unwrap()
    };
    assert_eq!(report.contributions.len(), 4);
    assert_eq!(report.occupied_reports.len(), 3);
    assert!(report.empty_support.is_empty());
    assert_eq!(
        report
            .occupied_reports
            .iter()
            .filter(|(_, r)| r.construction == "weighted_amf")
            .count(),
        2
    );
    assert_eq!(
        report
            .occupied_reports
            .iter()
            .filter(|(_, r)| r.construction == "compact_polynomial_moments")
            .count(),
        1
    );
    for (cuts, values) in &report.contributions {
        let scalar = match cuts.as_slice() {
            [] => "1/(128*pi^2)",
            [0] => "1/(128*pi^2)",
            [1] => "1/(64*pi^2)",
            [0, 1] => "1/(64*pi^2)",
            _ => panic!("unexpected cut set {cuts:?}"),
        };
        assert!(
            p.close(&values[0], &exact(scalar), 14),
            "scalar cut {cuts:?}"
        );
        if cuts.len() == 2 {
            assert!(p.close(&values[1], &exact("-3/(512*pi^2)"), 14));
        } else {
            assert!(
                p.norm(&values[1]) < p.tolerance(14),
                "odd virtual moment cut {cuts:?}"
            );
        }
    }
    assert!(p.close(&report.values[0], &exact("3/(64*pi^2)"), 14));
    assert!(p.close(&report.values[1], &exact("-3/(512*pi^2)"), 14));
}
