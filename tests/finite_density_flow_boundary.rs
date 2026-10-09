use ahash::HashMap;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::boundary::{
    IntegratedOccupiedBoundary, OccupiedBoundaryDistribution, OccupiedBoundaryLimits,
};
use symbolica_amflow::finite_density::compact::CompactShell;
use symbolica_amflow::finite_density::flow_boundary::OccupiedFlowBoundary;
use symbolica_amflow::finite_density::geometry::OccupiedCutFamily;
use symbolica_amflow::frobenius::{FrobeniusBasis, FrobeniusColumn};
use symbolica_amflow::*;

fn family(cuts: &[usize]) -> OccupiedCutFamily {
    let input: DensityInput = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/finite_density/massive_two_loop_sunset.json"
    )))
    .unwrap();
    input
        .prepare()
        .unwrap()
        .occupied_cut(cuts, 16)
        .unwrap()
        .at_physical_masses()
}

fn shift(family: &OccupiedCutFamily) -> Vec<bool> {
    (0..family.factors().len())
        .map(|slot| {
            slot < family.physical_slots()
                && !family
                    .shells()
                    .iter()
                    .any(|shell| shell.physical_slot == slot)
        })
        .collect()
}

fn reference(expression: &str, p: Precision) -> ComplexFloat {
    p.eval(&parse!(expression), &HashMap::default()).unwrap()
}

#[test]
fn occupied_matching_plans_and_integrates_subleading_coefficients() {
    let family = family(&[0, 1]);
    let original = family.targets()[0].keys().next().unwrap().clone();
    let mut pinched = original.clone();
    pinched.0[2] = 0;
    let basis = [pinched, original];
    let p = Precision::decimal(60).unwrap();
    // This deliberately simple test-only basis isolates the matching owner:
    // [J,I] = c0*[1,-z] + c1*[0,z^2] + O(z^3).
    // Its second constant cannot be fixed from the leading region terms.
    // The actual boundary coefficients still come from the physical sunset
    // integrand, graph routing, compact moments and spatial tensor projection.
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![
            FrobeniusColumn {
                exponent: Atom::zero(),
                coefficients: vec![
                    vec![vec![p.i(1), p.zero()]],
                    vec![vec![p.zero(), p.i(-1)]],
                    vec![vec![p.zero(), p.zero()]],
                ],
            },
            FrobeniusColumn {
                exponent: Atom::num(2),
                coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
            },
        ],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let eps = symbol!("occupied_flow_boundary_eps");
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        2,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let report = boundary
        .constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    assert!(p.close(&report.constants[0], &reference("1/(4*pi)", p), 40));
    assert!(p.close(&report.constants[1], &reference("13/(32*pi)", p), 40));
    assert_eq!(report.provenance.regions.len(), 1);
    assert_eq!(
        report.provenance.regions[0].exponents,
        vec![Atom::zero(), Atom::one()]
    );
    assert_eq!(
        report.provenance.regions[0].half_orders,
        vec![Some(0), Some(2)]
    );
    assert_eq!(report.provenance.integrated_coefficients, 4);
    assert_eq!(report.provenance.scaleless_noncompact_products, 0);

    let shallow = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        1,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    assert!(matches!(
        shallow.constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::IncompleteReduction(_))
    ));
    let mut short_series = solutions.clone();
    short_series.columns[0].coefficients.pop();
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &short_series,
        ),
        Err(Error::IncompleteReduction(_))
    ));
}

#[test]
fn mixed_boundary_keeps_symbolic_epsilon_power_and_deformation_roles() {
    let family = family(&[0]);
    let basis = vec![family.targets()[0].keys().next().unwrap().clone()];
    let p = Precision::decimal(60).unwrap();
    let eps = symbol!("occupied_mixed_flow_boundary_eps");
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![FrobeniusColumn {
            exponent: Atom::var(eps),
            coefficients: vec![vec![vec![p.i(1)]]],
        }],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        2,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let shifted = shift(&family);
    let report = boundary
        .constants(
            &family,
            &basis,
            &shifted,
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    assert!(p.close(&report.constants[0], &p.scale(&p.i(1), 1, 2), 40));
    assert!(report.provenance.regions.iter().any(|region| {
        region.exponents[0] == Atom::var(eps) && region.half_orders[0] == Some(0)
    }));
    assert!(
        report
            .provenance
            .regions
            .iter()
            .all(|region| !region.hard[0])
    );
    let mut moving_shell = shifted.clone();
    moving_shell[0] = true;
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &moving_shell,
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::Unsupported(_))
    ));
    let mut moving_occupation = shifted;
    moving_occupation[family.shells()[0].upper_slot] = true;
    assert!(matches!(
        boundary.constants(
            &family,
            &basis,
            &moving_occupation,
            &Rational::from((1, 2)),
            p,
            &solutions,
        ),
        Err(Error::Unsupported(_))
    ));
}
fn multihard_family(cuts: &[usize], powers: [i16; 4]) -> OccupiedCutFamily {
    // Four oriented parallel edges obey exact momentum and charge conservation:
    // q0=(q0-q1-q2)+q1+q2. Pinching the second denominator leaves independent
    // virtual loops, while the input graph and routing remain connected.
    let input: DensityInput = serde_json::from_value(serde_json::json!({
        "name": "occupied_multihard_boundary_regression",
        "loops": 3,
        "vertices": 2,
        "edges": [
            { "vertices": [0,1], "routing": ["1","0","0"], "mass_squared": "1/4", "charges": [1] },
            { "vertices": [1,0], "routing": ["1","-1","-1"], "mass_squared": "16", "charges": [1] },
            { "vertices": [1,0], "routing": ["0","1","0"], "mass_squared": "1/4", "charges": [0] },
            { "vertices": [1,0], "routing": ["0","0","1"], "mass_squared": "1/4", "charges": [0] }
        ],
        "loop_charges": [[1],[0],[0]],
        "chemical_potentials": ["1"],
        "targets": [{ "powers": powers, "numerator": "1" }],
        "numerator_convention": "shifted_euclidean",
        "laurent_orders": [-3,0],
        "digits": 12
    }))
    .unwrap();
    input
        .prepare()
        .unwrap()
        .occupied_cut(cuts, 64)
        .unwrap()
        .at_physical_masses()
}

#[test]
fn two_virtual_hard_loops_use_recursive_gaussian_boundaries() {
    let family = multihard_family(&[0], [1, 0, 2, 2]);
    let basis = vec![family.targets()[0].keys().next().unwrap().clone()];
    let eps = symbol!("multihard_occupied_eps");
    let p = Precision::decimal(60).unwrap();
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![FrobeniusColumn {
            exponent: Atom::num(2) * Atom::var(eps),
            coefficients: vec![vec![vec![p.i(1)]]],
        }],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let preparations = Arc::new(AtomicUsize::new(0));
    let observed = preparations.clone();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Prepared { .. }) {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        })),
        ..Default::default()
    };
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        0,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let report = boundary
        .constants(
            &family,
            &basis,
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    // D=3. Each hard power-two Gaussian has value sqrt(pi), and the raw
    // compact C1 is (mu-m)/sqrt(pi)=1/(2sqrt(pi)). No bubble/sunset formula.
    let expected = p.eval(&parse!("pi^(1/2)/2"), &HashMap::default()).unwrap();
    assert!(p.close(&report.constants[0], &expected, 20));
    assert!(preparations.load(Ordering::Relaxed) > 0);
    assert!(report.provenance.regions.iter().any(|region| {
        region.hard.iter().filter(|&&hard| hard).count() == 2
            && region.half_orders[0] == Some(0)
            && (region.exponents[0].clone() - Atom::num(2) * Atom::var(eps))
                .expand()
                .is_zero()
    }));
}

#[test]
fn three_compact_loops_need_only_polynomial_soft_moments() {
    // One nonbranching charged triangle and a neutral spanning tree. Cutting
    // all charged edges leaves a connected complement with three compact loops.
    let input: DensityInput = serde_json::from_value(serde_json::json!({
        "name": "three_compact_triangle_boundary",
        "loops": 3,
        "vertices": 3,
        "edges": [
            {"vertices": [0,1], "routing": ["1","0","0"], "mass_squared": "1/4", "charges": [1]},
            {"vertices": [1,2], "routing": ["0","1","0"], "mass_squared": "1/4", "charges": [1]},
            {"vertices": [2,0], "routing": ["0","0","1"], "mass_squared": "1/4", "charges": [1]},
            {"vertices": [1,0], "routing": ["1","0","-1"], "mass_squared": "16", "charges": [0]},
            {"vertices": [2,1], "routing": ["0","1","-1"], "mass_squared": "16", "charges": [0]}
        ],
        "loop_charges": [[1],[1],[1]],
        "chemical_potentials": ["1"],
        "targets": [{"powers": [1,1,1,1,1], "numerator": "1"}],
        "numerator_convention": "shifted_euclidean",
        "laurent_orders": [-3,0],
        "digits": 12
    }))
    .unwrap();
    let family = input
        .prepare()
        .unwrap()
        .occupied_cut(&[0, 1, 2], 64)
        .unwrap()
        .at_physical_masses();
    let eps = symbol!("three_compact_boundary_eps");
    let ordinary = family.region_family(eps, 4).unwrap();
    let mut target = family.targets()[0].keys().next().unwrap().clone();
    target.0.truncate(family.input_slots());
    for shell in family.shells() {
        target.0[shell.physical_slot] = 0;
    }
    let region = regions::LoopRegion {
        transformation: (0..3)
            .map(|i| (0..3).map(|j| Atom::num(i64::from(i == j))).collect())
            .collect(),
        hard: vec![false; 3],
        hard_branches: vec![],
        jacobian_determinant: Atom::one(),
    };
    let shifted = shift(&family);
    let expansion = regions::expand_region(
        &ordinary,
        &target,
        &shifted[..family.input_slots()],
        &region,
        2,
    )
    .unwrap();
    assert_eq!(expansion.eta_power, Atom::num(-2));
    let distributions = family
        .shells()
        .iter()
        .map(|shell| OccupiedBoundaryDistribution {
            source_loop_index: shell.loop_index,
            shell: CompactShell {
                mass_squared: Rational::from((1, 4)),
                chemical_potential: Rational::one(),
            },
            cut_index: 1,
            upper_index: 0,
            lower_index: 0,
        })
        .collect::<Vec<_>>();
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let p = Precision::decimal(60).unwrap();
    for (index, expected) in ["1/(8*pi^(3/2))", "0", "-133/(32*pi^(3/2))"]
        .into_iter()
        .enumerate()
    {
        let projected = integrand::projected_factor_region(
            &expansion.coefficients[index],
            &expansion.coordinates,
            &ordinary,
            &region.hard,
            1000,
        )
        .unwrap();
        let report = boundary
            .evaluate_projected_with_provenance(
                &projected,
                &distributions,
                &Rational::from((1, 2)),
                &HashMap::default(),
                p,
            )
            .unwrap();
        let expected = p.eval(&parse!(expected), &HashMap::default()).unwrap();
        assert!(p.close(&report.value, &expected, 40));
        assert_eq!(report.scaleless_noncompact_products, 0);
    }
}

#[test]
fn inverse_completion_boundary_uses_exact_positive_energy_certificate() {
    let virtual_family = family(&[0]);
    let family = family(&[0, 1]);
    let certificate = (family.physical_slots()..family.input_slots())
        .find_map(|slot| family.compact_energy_completion(slot))
        .unwrap();
    let mut integral = family.targets()[0].keys().next().unwrap().clone();
    integral.0[certificate.slot] = 1;
    let eps = symbol!("inverse_energy_flow_boundary_eps");
    let p = Precision::decimal(60).unwrap();
    let solutions = FrobeniusBasis {
        precision: p,
        columns: vec![FrobeniusColumn {
            exponent: Atom::one(),
            coefficients: vec![vec![vec![p.i(1)]]],
        }],
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        2,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    assert!(matches!(
        boundary.constants(
            &family,
            &[integral.clone()],
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions
        ),
        Err(Error::InvalidInput(_))
    ));
    let enabled = boundary.with_positive_compact_energy_powers(true);
    let report = enabled
        .constants(
            &family,
            &[integral.clone()],
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &solutions,
        )
        .unwrap();
    let expected = p.div(
        &reference("-log(2)/(2*pi)", p),
        &p.rational(&certificate.coefficient),
    );
    assert!(p.close(&report.constants[0], &expected, 40));
    assert_eq!(report.provenance.regions[0].exponents, vec![Atom::one()]);

    // An inverse compact energy also survives non-leading matching. Its
    // radial integral weights D=q1^2+q2^2-2q1.q2-1 before shell substitution.
    let mut pinched = integral.clone();
    pinched.0[2] = 0;
    let subleading_solutions = FrobeniusBasis {
        precision: p,
        columns: vec![
            FrobeniusColumn {
                exponent: Atom::zero(),
                coefficients: vec![
                    vec![vec![p.i(1), p.zero()]],
                    vec![vec![p.zero(), p.i(-1)]],
                    vec![vec![p.zero(), p.zero()]],
                ],
            },
            FrobeniusColumn {
                exponent: Atom::num(2),
                coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
            },
        ],
    };
    let subleading = enabled
        .constants(
            &family,
            &[pinched, integral],
            &shift(&family),
            &Rational::from((1, 2)),
            p,
            &subleading_solutions,
        )
        .unwrap();
    for (actual, expression) in subleading
        .constants
        .iter()
        .zip(["log(2)/(2*pi)", "(log(2)/4+3/8)/pi"])
    {
        let expected = p.div(
            &reference(expression, p),
            &p.rational(&certificate.coefficient),
        );
        assert!(p.close(actual, &expected, 40));
    }
    assert_eq!(
        subleading.provenance.regions[0].half_orders,
        vec![Some(0), Some(2)]
    );

    let non_energy = (virtual_family.physical_slots()..virtual_family.input_slots())
        .find(|&slot| virtual_family.compact_energy_completion(slot).is_none())
        .unwrap();
    let mut integral = virtual_family.targets()[0].keys().next().unwrap().clone();
    integral.0[non_energy] = 1;
    assert!(matches!(
        enabled.constants(
            &virtual_family,
            &[integral],
            &shift(&virtual_family),
            &Rational::from((1, 2)),
            p,
            &solutions
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn sealed_massless_origin_integrates_hard_factors_upper_contacts_and_zero_lower_jets() {
    use symbolica_amflow::finite_density::massless_endpoint::MasslessFlowEvidence;
    use symbolica_amflow::finite_density::preparation::WeightedSourceOptions;
    let mut input: DensityInput = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/finite_density/massive_two_loop_sunset.json"
    )))
    .unwrap();
    for edge in &mut input.edges {
        edge.mass_squared = "0".into();
    }
    let input = input.prepare().unwrap();
    let family = input.occupied_cut(&[0], 16).unwrap().at_physical_masses();
    let proof =
        MasslessFlowEvidence::new(&input, &family, &[1, 2], WeightedSourceOptions::default())
            .unwrap();
    let shell = &family.shells()[0];
    let base = family.targets()[0].keys().next().unwrap().clone();
    let mut raised = base.clone();
    raised.0[shell.physical_slot] = 2;
    let mut upper = base.clone();
    upper.0[shell.upper_slot] = 1;
    let mut lower = base;
    lower.0[shell.lower_slot] = 1;
    let basis = [raised, upper, lower];
    let eps = symbol!("massless_flow_boundary_eps");
    let p = Precision::decimal(60).unwrap();
    // The artificial identity basis isolates boundary matching. Each actual
    // leading region is eta^(-eps), with its coefficient integrated from the
    // original graph; no endpoint or differential-equation value is supplied.
    let solutions = FrobeniusBasis {
        precision: p,
        columns: (0..3)
            .map(|column| FrobeniusColumn {
                exponent: Atom::var(eps),
                coefficients: vec![vec![
                    (0..3).map(|row| p.i(i64::from(row == column))).collect(),
                ]],
            })
            .collect(),
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        eps,
        0,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let epsilon = Rational::from((4, 5));
    assert!(matches!(
        boundary.constants(&family, &basis, &shift(&family), &epsilon, p, &solutions),
        Err(Error::Unsupported(_))
    ));
    let boundary = boundary.with_massless_evidence(&proof).unwrap();
    let result = boundary
        .constants(&family, &basis, &shift(&family), &epsilon, p, &solutions)
        .unwrap();
    // At D=12/5 the C2 seed is outside its bare origin convergence range.
    // These are the continued one-shell moments times the ordinary power-two
    // Gaussian seed Gamma(4/5), not a two-loop integral formula.
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let normalization = p.mul(
        &p.pow(&pi, &p.rational(&Rational::from((1, 2)))),
        &p.gamma_real(&p.rational(&Rational::from((7, 10))).re)
            .unwrap(),
    );
    let gaussian = p
        .gamma_real(&p.rational(&Rational::from((4, 5))).re)
        .unwrap();
    let upper_expected = p.div(&gaussian, &normalization);
    let raised_expected = p.scale(&upper_expected, -3, 16);
    assert!(p.close(&result.constants[0], &raised_expected, 40));
    assert!(p.close(&result.constants[1], &upper_expected, 40));
    assert!(p.close(&result.constants[2], &p.zero(), 40));
    assert_eq!(
        result.provenance.massless_origin,
        Some(proof.origin_identity())
    );
    assert!(
        result
            .provenance
            .regions
            .iter()
            .any(|r| r.hard.iter().any(|&h| h))
    );
}
