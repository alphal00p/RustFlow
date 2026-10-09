//! Public amplitude assembly at compact support thresholds, without a flow
//! connection for the occupied polynomial terminal.
use ahash::HashMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput, assembly::PreparedDensityFlow, reduction::WeightedClosureOptions,
    terminal::PreparedOccupiedTerminal,
};
use symbolica_amflow::*;

fn input(mu: &str, powers: &[i16]) -> DensityInput {
    serde_json::from_value(serde_json::json!({
        "name": "one_loop_compact_threshold",
        "loops": 1,
        "vertices": 1,
        "edges": [{"vertices": [0,0], "routing": ["1"], "mass_squared": "1/4", "charges": [1]}],
        "loop_charges": [[1]],
        "chemical_potentials": [mu],
        "targets": powers.iter().map(|power| serde_json::json!({"powers": [power], "numerator": "1"})).collect::<Vec<_>>(),
        "numerator_convention": "shifted_euclidean",
        "laurent_orders": [-1,0],
        "digits": 18
    })).unwrap()
}

fn options() -> FlowOptions {
    FlowOptions {
        digits: 18,
        guard_digits: 24,
        series_order: 60,
        ..Default::default()
    }
}

fn expected(expression: &str, p: Precision) -> ComplexFloat {
    p.eval(
        &Atom::parse(expression, "density_terminal_test", Default::default()).unwrap(),
        &HashMap::default(),
    )
    .unwrap()
}

#[test]
fn public_assembly_above_threshold_keeps_vacuum_and_raised_surface() {
    let options = options();
    let context = RunContext::default();
    let flow = PreparedDensityFlow::prepare(
        &input("1", &[1, 2]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    assert!(flow.closure_diagnostics().is_empty());
    let report = flow
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    assert_eq!(report.contributions.len(), 2);
    assert_eq!(report.occupied_reports.len(), 1);
    let occupied = &report.occupied_reports[0].1;
    assert_eq!(occupied.construction, "compact_polynomial_moments");
    assert!(occupied.shifted_slots.is_empty());
    assert_eq!(occupied.basis_size, 0);
    // In D=3, I1^vac=-m/(4pi), I1^occ=-(mu-m)/(4pi).
    // Holding the original numerator fixed gives I2=-dI1/d(m²): the
    // occupied raised surface exactly cancels the vacuum raised integral.
    for (actual, reference) in report.contributions[0]
        .1
        .iter()
        .zip(["-1/(8*pi)", "1/(4*pi)"])
    {
        assert!(p.close(actual, &expected(reference, p), 14));
    }
    for (actual, reference) in occupied.values.iter().zip(["-1/(8*pi)", "-1/(4*pi)"]) {
        assert!(p.close(actual, &expected(reference, p), 14));
    }
    assert!(p.close(&report.values[0], &expected("-1/(4*pi)", p), 14));
    assert!(p.norm(&report.values[1]) < p.tolerance(14));
}

#[test]
fn public_assembly_below_threshold_retains_massive_vacuum() {
    let options = options();
    let context = RunContext::default();
    let flow = PreparedDensityFlow::prepare(
        &input("1/4", &[1]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    let report = flow
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    assert_eq!(report.empty_support.len(), 1);
    assert!(report.occupied_reports.is_empty());
    assert_eq!(report.contributions[1].1[0], p.zero());
    assert!(p.close(&report.values[0], &expected("-1/(8*pi)", p), 14));
}

#[test]
fn empty_support_and_compact_terminals_keep_the_common_amf_contract() {
    let context = RunContext::default();
    for mu in ["1/4", "1"] {
        let input = input(mu, &[1]);
        let baseline = options();
        let prepared = PreparedDensityFlow::prepare(
            &input,
            &baseline,
            WeightedClosureOptions::default(),
            &context,
        )
        .unwrap();
        let inadmissible = [
            FlowOptions {
                prescription: Prescription::MinusI0,
                ..baseline.clone()
            },
            FlowOptions {
                recursion: RecursionMode::Ft,
                ..baseline.clone()
            },
            FlowOptions {
                mass_mode: MassMode::Propagator,
                ..baseline.clone()
            },
        ];
        for options in inadmissible {
            assert!(matches!(
                PreparedDensityFlow::prepare(
                    &input,
                    &options,
                    WeightedClosureOptions::default(),
                    &context,
                ),
                Err(Error::Unsupported(_))
            ));
            assert!(matches!(
                prepared.evaluate(&Rational::from((1, 2)), &options, &context, 8,),
                Err(Error::Unsupported(_))
            ));
        }
    }
}

#[test]
fn public_assembly_threshold_uses_thermal_finite_jump_and_diagnoses_singular_case() {
    let options = options();
    let context = RunContext::default();
    let flow = PreparedDensityFlow::prepare(
        &input("1/2", &[1]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    let report = flow
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    assert!(report.empty_support.is_empty());
    assert_eq!(report.occupied_reports[0].1.values[0], p.zero());
    assert!(p.close(&report.values[0], &expected("-1/(8*pi)", p), 14));
    let raised = PreparedDensityFlow::prepare(
        &input("1/2", &[2]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    let raised_report = raised
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    // The T->0+ Fermi regulator at mu=m has n_F(0)=1/2. It gives half the
    // upper-side occupied jump while preserving the nonzero vacuum term.
    assert!(p.close(
        &raised_report.occupied_reports[0].1.values[0],
        &expected("-1/(8*pi)", p),
        14
    ));
    assert!(p.close(&raised_report.values[0], &expected("1/(8*pi)", p), 14));
    let singular = PreparedDensityFlow::prepare(
        &input("1/2", &[3]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    assert!(
        matches!(singular.evaluate(&Rational::from((1, 2)), &options, &context, 8),
        Err(Error::Numerical(message)) if message.contains("raised occupied threshold") && message.contains("gap exponent -1"))
    );
}

#[test]
fn terminal_declines_virtual_poles_and_zeroes_unrestricted_polynomials() {
    let mut input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    let options = options();
    assert!(
        PreparedOccupiedTerminal::prepare(&input.prepare().unwrap(), &[0], &options)
            .unwrap()
            .is_none()
    );
    input.targets.truncate(1);
    input.targets[0].powers = vec![1, 0, 0];
    input.targets[0].numerator = "1".into();
    let terminal = PreparedOccupiedTerminal::prepare(&input.prepare().unwrap(), &[0], &options)
        .unwrap()
        .unwrap();
    let report = terminal
        .evaluate(&Rational::from((1, 2)), &options, &RunContext::default())
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    assert_eq!(report.values, vec![p.zero()]);
    assert!(report.boundary.scaleless_noncompact_products > 0);
}
