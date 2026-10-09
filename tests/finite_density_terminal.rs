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

fn massless_input(mu: &str, targets: &[(i16, &str)]) -> DensityInput {
    let mut density = input(mu, &[1]);
    density.edges[0].mass_squared = "0".into();
    density.targets = targets
        .iter()
        .map(|(power, numerator)| {
            serde_json::from_value(serde_json::json!({
                "powers": [power], "numerator": numerator,
            }))
            .unwrap()
        })
        .collect();
    density
}

#[test]
fn public_massless_terminal_assembly_matches_fixed_numerator_mass_derivatives() {
    let options = options();
    let context = RunContext::default();
    let density = massless_input(
        "1",
        &[
            (1, "1"),
            (2, "1"),
            (1, "u1"),
            (2, "u1^2"),
            (2, "g1_1-u1^2"),
            (2, "g1_1"),
        ],
    );
    let flow = PreparedDensityFlow::prepare(
        &density,
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    assert!(flow.closure_diagnostics().is_empty());
    // D=5 (d=4) is a convergent compact domain for these derivatives.
    // Independently, the regulated thermal T->0 simple scalar is
    // -1/(16*pi²)*integral_sqrt(a)^mu (E²-a)dE. Differentiate its elementary
    // antiderivative at fixed original numerator before a->0. Inserting
    // original P0=iE and |q|²=E²-a yields the other references below.
    let report = flow
        .evaluate(&Rational::from((-1, 2)), &options, &context, 8)
        .unwrap();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    assert_eq!(report.contributions.len(), 2);
    assert!(
        report.contributions[0]
            .1
            .iter()
            .all(|value| *value == p.zero())
    );
    assert!(
        report.occupied_reports[0]
            .1
            .contour_admission
            .contains("joint-high-D-massless-polynomial-origin-v1")
    );
    for (actual, expression) in report.values.iter().zip([
        "-1/(48*pi^2)",
        "-1/(16*pi^2)",
        "-1i/(64*pi^2)",
        "1/(48*pi^2)",
        "-1/(24*pi^2)",
        "-1/(48*pi^2)",
    ]) {
        assert!(
            p.close(actual, &expected(expression, p), 16),
            "{expression}"
        );
    }
    // Multiplication is applied to the original polynomial, before cut
    // derivatives: P²/(P²)^2 equals 1/P² after the complete cut decomposition.
    assert!(p.close(&report.values[0], &report.values[5], 16));
}

#[test]
fn public_massless_terminal_keeps_removable_and_true_dimensional_poles_distinct() {
    let options = options();
    let context = RunContext::default();
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    let removable = PreparedDensityFlow::prepare(
        &massless_input("1", &[(2, "u1")]),
        &options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    let report = removable
        .evaluate(&Rational::from((1, 2)), &options, &context, 8)
        .unwrap();
    // At d=2 the independent massive thermal primitive is
    // -i/(8*pi)*(mu²-a), so -d/da gives -i/(8*pi) at every 0<a<mu².
    // The continued massless radial ratio is the removable (d/2-1)/(d-2).
    assert!(p.close(&report.values[0], &expected("-1i/(8*pi)", p), 16));
    let pole_options = FlowOptions {
        dimension: 5,
        ..options.clone()
    };
    let true_pole = PreparedDensityFlow::prepare(
        &massless_input("1", &[(2, "1")]),
        &pole_options,
        WeightedClosureOptions::default(),
        &context,
    )
    .unwrap();
    let Err(assembly_error) =
        true_pole.evaluate(&Rational::from((1, 2)), &pole_options, &context, 8)
    else {
        panic!("full amplitude accepted an uncancelled dimensional pole");
    };
    eprintln!("full massless D=4 pole diagnostic: {assembly_error}");
    // The ordinary vacuum owner runs first and may reach its auxiliary
    // tadpole Gamma(2-D/2) pole or a retained native reduction condition.
    assert!(
        matches!(&assembly_error,
            Error::Numerical(message) if message.contains("uncancelled dimensional pole"))
            || matches!(&assembly_error,
                Error::Numerical(message)
                    if message == "gamma evaluation reached a pole or exceeded the numeric backend's range")
            || matches!(&assembly_error,
                Error::InvalidInput(message)
                    if message == "physical point violates a required nonzero reduction condition"),
        "unexpected full-amplitude pole error: {assembly_error}"
    );
    let terminal = PreparedOccupiedTerminal::prepare(
        &massless_input("1", &[(2, "1")]).prepare().unwrap(),
        &[0],
        &pole_options,
    )
    .unwrap()
    .unwrap();
    assert!(matches!(
        terminal.evaluate(&Rational::from((1, 2)), &pole_options, &context),
        Err(Error::Numerical(message)) if message.contains("uncancelled dimensional pole")
    ));
    assert!(matches!(
        PreparedDensityFlow::prepare(
            &massless_input("0", &[(1, "1")]), &options,
            WeightedClosureOptions::default(), &context,
        ), Err(Error::InvalidInput(message)) if message.contains("zero chemical potential")
    ));
}

#[test]
fn massless_terminal_admission_does_not_admit_flowing_virtual_poles() {
    let mut density: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    density.edges[0].mass_squared = "0".into();
    density.edges[1].mass_squared = "0".into();
    let options = options();
    assert!(
        PreparedOccupiedTerminal::prepare(&density.prepare().unwrap(), &[0], &options,)
            .unwrap()
            .is_none()
    );
    assert!(matches!(PreparedDensityFlow::prepare(
        &density, &options, WeightedClosureOptions::default(), &RunContext::default(),
    ), Err(Error::Unsupported(message)) if message.contains("strictly positive shell masses")));
}
