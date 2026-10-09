//! Independent schema-3 input identities and native admission checks.
//! The supplied values below are parsed only as validation data. No numerical
//! integral predictions, boundary values, or reference comparisons are claimed.
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;
use symbolica_amflow::Integral;
use symbolica_amflow::finite_density::{DensityInput, InversePropagatorBasis, NumeratorConvention};

const ORACLE: &str = include_str!("../fixtures/finite_density/oracle_results_supplied.json");
const EXAMPLES: [&str; 4] = [
    include_str!("../examples/finite_density/chain_of_three_parallel_pairs.json"),
    include_str!("../examples/finite_density/five_vertex_eight_edge.json"),
    include_str!("../examples/finite_density/triangular_prism.json"),
    include_str!("../examples/finite_density/massive_two_loop_sunset.json"),
];

fn parse(value: &Value) -> Atom {
    Atom::parse(
        value.as_str().expect("mathematical expression string"),
        "rustflow_density",
        Default::default(),
    )
    .expect("native Symbolica expression")
}
fn powers(value: &Value) -> Vec<i16> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|n| i16::try_from(n.as_i64().unwrap()).unwrap())
        .collect()
}
fn equal(left: &Atom, right: &Atom, context: &str) {
    assert!(
        (left - right).together().cancel().is_zero(),
        "{context}: {left} != {right}"
    );
}
fn small() -> DensityInput {
    serde_json::from_str(EXAMPLES[3]).unwrap()
}

#[test]
fn all_108_supplied_target_expansions_use_the_same_inverse_propagator_slots() {
    let oracle: Value = serde_json::from_str(ORACLE).unwrap();
    assert_eq!(oracle["schema_version"], 3);
    assert_eq!(oracle["four_loop_oracles"].as_array().unwrap().len(), 95);
    assert_eq!(oracle["lower_loop_oracles"].as_array().unwrap().len(), 13);
    let mut bases = BTreeMap::new();
    for (name, value) in oracle["bases"].as_object().unwrap() {
        let map = &value["inverse_propagator_representation"];
        let coordinates = map["scalar_product_symbols"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse)
            .collect();
        let slots = map["inverse_propagator_polynomials"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse)
            .collect();
        let labels = map["inverse_propagator_symbols"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse)
            .collect();
        let basis = InversePropagatorBasis::new(coordinates, slots, labels).unwrap();
        for (coordinate, expected) in map["scalar_products_in_inverse_propagators"]
            .as_object()
            .unwrap()
        {
            equal(
                &basis.inverse_map()
                    [&Atom::parse(coordinate, "rustflow_density", Default::default()).unwrap()],
                &parse(expected),
                name,
            );
        }
        // Independently rebuild every forward quadratic from its routing.
        for (slot, routing) in value["positive_index_routings"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            let mut quadratic = parse(&map["mass_squared"][slot]);
            for (coordinate, pair) in basis
                .coordinates()
                .iter()
                .zip(map["scalar_product_pairs"].as_array().unwrap())
            {
                let a = pair[0].as_u64().unwrap() as usize - 1;
                let b = pair[1].as_u64().unwrap() as usize - 1;
                let coefficient = routing[a].as_i64().unwrap()
                    * routing[b].as_i64().unwrap()
                    * if a == b { 1 } else { 2 };
                quadratic += Atom::num(coefficient) * coordinate;
            }
            equal(&basis.slots()[slot], &quadratic, "forward routing square");
        }
        for mass in map["mass_squared"].as_array().unwrap() {
            let _ = parse(mass);
        }
        let _ = parse(&map["reference_normalization"]);
        bases.insert(name.clone(), basis);
    }
    let records = oracle["four_loop_oracles"]
        .as_array()
        .unwrap()
        .iter()
        .chain(oracle["lower_loop_oracles"].as_array().unwrap());
    let mut checked = 0;
    for record in records {
        let id = record["id"].as_str().unwrap();
        let basis = &bases[record["basis"].as_str().unwrap()];
        let definition = &record["finite_density_target"];
        let denominator_powers = powers(&definition["denominator_powers"]);
        assert!(denominator_powers.iter().all(|&n| n >= 0));
        let numerator = parse(&definition["numerator_polynomial"]);
        let legacy = powers(&record["powers"]);
        assert_eq!(
            denominator_powers,
            legacy.iter().map(|&n| n.max(0)).collect::<Vec<_>>()
        );
        let dot_pairs =
            &oracle["bases"][record["basis"].as_str().unwrap()]["negative_index_dot_pairs"];
        let mut legacy_numerator = Atom::num(1);
        for (slot, &power) in legacy.iter().enumerate() {
            if power < 0 {
                let pair = &dot_pairs[slot];
                let g = Atom::parse(
                    &format!("g{}_{}", pair[0], pair[1]),
                    "rustflow_density",
                    Default::default(),
                )
                .unwrap();
                legacy_numerator *= g.pow(-i64::from(power));
            }
        }
        equal(
            &numerator,
            &legacy_numerator,
            "legacy dotted numerator definition",
        );
        let actual = basis.convert(&numerator, &denominator_powers).unwrap();
        let mut expected = BTreeMap::<Integral, Atom>::new();
        for term in definition["indexed_terms"].as_array().unwrap() {
            *expected
                .entry(Integral(powers(&term["powers"])))
                .or_default() += parse(&term["coefficient"]);
        }
        let indices = actual
            .keys()
            .chain(expected.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        for index in indices {
            equal(
                actual.get(&index).unwrap_or(&Atom::new()),
                expected.get(&index).unwrap_or(&Atom::new()),
                id,
            );
        }
        // Also verify the complete reconstructed integrand, so target powers
        // cannot accidentally be interpreted with the legacy dotted convention.
        let direct = denominator_powers
            .iter()
            .zip(basis.slots())
            .fold(numerator, |a, (&n, rho)| a * rho.clone().pow(-i64::from(n)));
        equal(&basis.reconstruct(&actual).unwrap(), &direct, id);
        let _ = parse(&definition["dimension"]);
        for chemical in definition["chemical_potentials"].as_array().unwrap() {
            let _ = parse(chemical);
        }
        checked += 1;
    }
    assert_eq!(checked, 108);
}

#[test]
fn complete_reference_payload_is_preserved_and_parses_only_in_validation() {
    let oracle: Value = serde_json::from_str(ORACLE).unwrap();
    assert_eq!(oracle["constants"].as_object().unwrap().len(), 9);
    assert!(oracle["constants"].get("a4").is_none());
    for value in oracle["constants"].as_object().unwrap().values() {
        let _ = parse(value);
    }
    let _ = parse(&oracle["conventions"]["measure_factor"]);
    let _ = parse(&oracle["conventions"]["normalization"]);
    let mut coefficients = 0;
    let mut uncertainties = 0;
    let mut determinant_tags = 0;
    for record in oracle["four_loop_oracles"]
        .as_array()
        .unwrap()
        .iter()
        .chain(oracle["lower_loop_oracles"].as_array().unwrap())
    {
        for (order, value) in record["coefficients"].as_object().unwrap() {
            let _: i32 = order.parse().unwrap();
            let _ = parse(value);
            coefficients += 1;
        }
        if let Some(values) = record["uncertainties"].as_object() {
            for (order, value) in values {
                let _: i32 = order.parse().unwrap();
                let absolute = parse(&value["absolute"]);
                let rational = Rational::try_from(absolute.as_view()).unwrap();
                assert!(rational >= Rational::from(0));
                uncertainties += 1;
                if value.get("tag").is_some() {
                    assert_eq!(value["tag"], "det");
                    determinant_tags += 1;
                }
            }
        }
    }
    assert!(coefficients > 200 && uncertainties > 0 && determinant_tags > 0);
    // Production examples expose only input identities. Recursively forbid the
    // reference-bearing keys, including nested per-record payloads.
    let definitions: Value = serde_json::from_str(include_str!(
        "../examples/finite_density/oracle_definitions.json"
    ))
    .unwrap();
    fn no_answers(value: &Value) {
        match value {
            Value::Object(map) => {
                for (key, child) in map {
                    assert!(
                        ![
                            "constants",
                            "coefficients",
                            "uncertainties",
                            "det",
                            "reference_value"
                        ]
                        .contains(&key.as_str())
                    );
                    no_answers(child);
                }
            }
            Value::Array(array) => array.iter().for_each(no_answers),
            _ => {}
        }
    }
    no_answers(&definitions);
    assert_eq!(definitions["targets"].as_array().unwrap().len(), 108);
}

#[test]
fn user_defined_graphs_prepare_and_retain_all_connected_occupied_cuts() {
    let manifest: Value =
        serde_json::from_str(include_str!("../docs/finite-density-acceptance.json")).unwrap();
    for (i, source) in EXAMPLES.iter().enumerate() {
        let input: DensityInput = serde_json::from_str(source).unwrap();
        let prepared = input.prepare().unwrap();
        assert_eq!(prepared.basis().physical_slots(), input.edges.len());
        assert_eq!(prepared.targets().len(), 2);
        assert_eq!(prepared.independent_masses().len(), input.edges.len());
        let actual = prepared.cut_decomposition(1024).unwrap();
        let actual_slots = actual
            .iter()
            .map(|cut| {
                cut["cut_slots"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_u64().unwrap() as usize)
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>();
        let certificate = if i < 3 {
            &manifest["four_loop_families"][i]["graph_certificate"]
        } else {
            &manifest["small_graph"]["graph_certificate"]
        };
        let mut expected_slots = certificate["connected_admissible_occupied_cutsets"]
            .as_array()
            .unwrap()
            .iter()
            .map(|cut| {
                cut["slots"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_u64().unwrap() as usize)
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>();
        expected_slots.insert(vec![]); // Complete amplitude includes the vacuum.
        assert_eq!(actual_slots, expected_slots, "{}", input.name);
        for cut in actual {
            let count = cut["cut_slots"].as_array().unwrap().len();
            assert_eq!(
                cut["cut_sign"].as_i64().unwrap(),
                if count % 2 == 0 { 1 } else { -1 }
            );
            let determinant = parse(&cut["determinant"]);
            assert!(!determinant.is_zero());
        }
        assert!(prepared.cut_decomposition(1).is_err());
        let report = prepared.preparation_report(1024).unwrap();
        assert_eq!(report["status"], "algebraic_preparation_only");
        assert_eq!(report["numerical_evaluation"]["available"], false);
    }
}

#[test]
fn full_routing_rank_does_not_replace_incidence_or_chemical_evidence() {
    let mut wrong_incidence = small();
    wrong_incidence.edges[0].vertices.reverse(); // Rank of R remains two.
    let error = wrong_incidence.prepare().unwrap_err().to_string();
    assert!(error.contains("incidence"), "{error}");

    let mut wrong_charges = small();
    wrong_charges.loop_charges = vec![vec![1], vec![0]];
    let error = wrong_charges.prepare().unwrap_err().to_string();
    assert!(error.contains("chemical charge"), "{error}");

    let mut branched = small();
    branched.loop_charges = vec![vec![2], vec![1]];
    branched.edges[0].charges = vec![2];
    branched.edges[1].charges = vec![1];
    branched.edges[2].charges = vec![1];
    let error = branched.prepare().unwrap_err().to_string();
    assert!(error.contains("nonbranching"), "{error}");
}

#[test]
fn zero_power_preserves_chemical_assignments_and_medium_offsets() {
    let mut input = small();
    input.targets.truncate(1);
    input.targets[0].powers = vec![1, 0, 1];
    input.targets[0].numerator = "g1_2+u1*u2".into();
    input.numerator_convention = NumeratorConvention::UnshiftedEuclidean;
    let prepared = input.prepare().unwrap();
    assert_eq!(prepared.input().edges[1].charges, vec![1]);
    assert_eq!(
        prepared.loop_chemical_potentials(),
        &[Atom::num(1), Atom::num(1)]
    );
    let expected_numerator = parse!("g1_2+u1*u2-2", default_namespace = "rustflow_density")
        - Atom::num(2) * Atom::i() * parse!("u1+u2", default_namespace = "rustflow_density");
    let slots = prepared.basis().slots();
    let expected = expected_numerator / (&slots[0] * &slots[2]);
    equal(
        &prepared
            .basis()
            .reconstruct(&prepared.targets()[0])
            .unwrap(),
        &expected,
        "zero-power charged edge retains shifted-coordinate offsets",
    );
    assert!(
        prepared.targets()[0]
            .keys()
            .all(|index| index.0.len() == slots.len())
    );
}

#[test]
fn independent_physical_mass_derivatives_preserve_original_polynomial() {
    let prepared = small().prepare().unwrap();
    let basis = prepared.basis();
    let original = &prepared.targets()[1];
    let [first, second, ..] = prepared.independent_masses() else {
        panic!("two independent line masses")
    };
    assert_ne!(first, second);
    let derivative_first = basis.mass_derivative(original, *first).unwrap();
    let derivative_second = basis.mass_derivative(original, *second).unwrap();
    let mixed = basis.mass_derivative(&derivative_first, *second).unwrap();
    let slots = basis.slots();
    let numerator = parse!("g1_2+u1*u2", default_namespace = "rustflow_density");
    let expected_first =
        -Atom::num(2) * &numerator / (slots[0].clone().pow(3) * &slots[1] * &slots[2]);
    let expected_second =
        -&numerator / (slots[0].clone().pow(2) * slots[1].clone().pow(2) * &slots[2]);
    let expected_mixed =
        Atom::num(2) * numerator / (slots[0].clone().pow(3) * slots[1].clone().pow(2) * &slots[2]);
    equal(
        &basis.reconstruct(&derivative_first).unwrap(),
        &expected_first,
        "first physical mass",
    );
    equal(
        &basis.reconstruct(&derivative_second).unwrap(),
        &expected_second,
        "second physical mass",
    );
    let actual_mixed = basis.reconstruct(&mixed).unwrap();
    equal(&actual_mixed, &expected_mixed, "mixed independent masses");
    let specialize = |a: Atom| {
        prepared
            .independent_masses()
            .iter()
            .zip(prepared.physical_masses())
            .fold(a, |value, (&mass, physical)| {
                value.replace(Atom::var(mass)).with(physical.clone())
            })
    };
    equal(
        &specialize(actual_mixed),
        &specialize(expected_mixed),
        "equal-mass specialization after derivatives",
    );
    // Differentiating a common mass before identifying which line was raised
    // produces the sum of line derivatives, not either individual result.
    assert!(
        !(basis.reconstruct(&derivative_first).unwrap()
            - basis.reconstruct(&derivative_second).unwrap())
        .together()
        .cancel()
        .is_zero()
    );
}

#[test]
fn independent_chemical_species_retain_separate_cycle_assignments() {
    let mut input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/two_independent_chemical_cycles.json"
    ))
    .unwrap();
    let first = input.prepare().unwrap();
    assert_eq!(
        first.loop_chemical_potentials(),
        &[Atom::num(1), Atom::num((3, 2)), Atom::new()]
    );
    let cuts = first.cut_decomposition(32).unwrap();
    assert!(
        cuts.iter()
            .any(|cut| cut["cut_slots"] == serde_json::json!([0, 2]))
    );
    // All four charged lines cannot be removed: the neutral edges alone leave
    // disconnected components. Charge independence does not override the graph.
    assert!(
        !cuts
            .iter()
            .any(|cut| cut["cut_slots"] == serde_json::json!([0, 1, 2, 3]))
    );
    input.chemical_potentials = vec!["5/4".into(), "7/4".into()];
    let second = input.prepare().unwrap();
    assert_eq!(
        second.loop_chemical_potentials(),
        &[Atom::num((5, 4)), Atom::num((7, 4)), Atom::new()]
    );
    assert_ne!(first.identity(), second.identity());
    assert_eq!(first.input().edges[4].charges, vec![0, 0]);
    assert_eq!(second.input().edges[5].charges, vec![0, 0]);
}

#[test]
fn physical_connectivity_is_undirected_under_consistent_edge_reorientation() {
    let input = small();
    let original = input.prepare().unwrap();
    let mut reversed = input;
    // Reverse exactly one edge together with its physical momentum and charge.
    // Its squared denominator and all undirected connected complements agree.
    reversed.edges[0].vertices.reverse();
    for coefficient in &mut reversed.edges[0].routing {
        *coefficient = format!("-({coefficient})");
    }
    for charge in &mut reversed.edges[0].charges {
        *charge = -*charge;
    }
    let oriented = reversed.prepare().unwrap();
    let cuts = |prepared: &symbolica_amflow::finite_density::PreparedDensityInput| {
        prepared
            .cut_decomposition(32)
            .unwrap()
            .into_iter()
            .map(|entry| entry["cut_slots"].to_string())
            .collect::<BTreeSet<_>>()
    };
    assert_eq!(cuts(&original), cuts(&oriented));
    for (left, right) in original.targets().iter().zip(oriented.targets()) {
        equal(
            &original.basis().reconstruct(left).unwrap(),
            &oriented.basis().reconstruct(right).unwrap(),
            "consistent edge reorientation preserves physical quadratic",
        );
    }
}

#[test]
fn compatible_rational_routing_retains_a_nonunit_exact_jacobian() {
    let mut input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/two_independent_chemical_cycles.json"
    ))
    .unwrap();
    // Rescale the neutral loop coordinate. Its zero species charges allow this
    // rational transformation without fractional loop-charge assignments.
    for edge in &mut input.edges {
        edge.routing[2] = format!("({})/2", edge.routing[2]);
    }
    let prepared = input.prepare().unwrap();
    let cuts = prepared.cut_decomposition(32).unwrap();
    let cut = cuts
        .iter()
        .find(|cut| cut["cut_slots"] == serde_json::json!([0, 1, 2]))
        .unwrap();
    let determinant = parse(&cut["determinant"]);
    equal(
        &determinant.clone().pow(2),
        &Atom::num((1, 4)),
        "nonunit routing determinant",
    );
    assert_eq!(
        cut["jacobian_convention"],
        "abs(determinant)^(-D) before energy-shell integration"
    );
}
