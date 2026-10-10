use rustred::solver::guarded::{
    GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain,
};
use rustred::solver::{Integral, Power, Term};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::guarded::GuardedMeasureIdentity;
use symbolica_amflow::finite_density::preparation::{WeightedSourceOptions, WeightedSourcePolicy};

pub struct Prepared {
    pub sources: Arc<GuardedSourceSystem<16>>,
    pub targets: Vec<BTreeMap<[i64; 16], Atom>>,
    pub shifted: Vec<usize>,
    pub admitted: IndexDomain<16>,
}

pub fn prepare(args: &[String], out: &Path) -> Prepared {
    let mode = &args[3];
    let proof_sha = &args[4];
    let (cuts, shifted, zero_forced): (Vec<usize>, Vec<usize>, Vec<Vec<usize>>) =
        match mode.as_str() {
            "single-slot0" => (vec![3], vec![0], vec![vec![0], vec![1, 4]]),
            "double-slots12" => (vec![0, 3], vec![1, 2], vec![vec![1]]),
            "double-slots24" => (vec![0, 3], vec![2, 4], vec![vec![4]]),
            _ => panic!("only separately proved partial placements admitted to diagnostic"),
        };
    let raw: serde_json::Value = serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();
    // This is a fixture-bound validation certificate, never a graph-name admission.
    assert_eq!(raw["loops"], 3);
    assert_eq!(raw["edges"].as_array().unwrap().len(), 5);
    let expected = [
        ["1", "0", "0"],
        ["0", "1", "0"],
        ["0", "0", "1"],
        ["1", "0", "-1"],
        ["0", "1", "-1"],
    ];
    for (edge, row) in raw["edges"].as_array().unwrap().iter().zip(expected) {
        assert_eq!(edge["routing"], serde_json::json!(row));
        assert_eq!(edge["mass_squared"], "0");
    }
    let input: DensityInput = serde_json::from_value(raw).unwrap();
    let input = input.prepare().unwrap();
    let family = input
        .occupied_cut(&cuts, 4096)
        .unwrap()
        .at_physical_masses();
    assert_eq!(family.physical_slots(), 5);
    assert_eq!(family.input_slots(), 9);
    assert_eq!(family.loops(), 3);
    assert!(
        family
            .shells()
            .iter()
            .all(|s| s.mass_squared.is_zero() && s.chemical_potential > Rational::zero())
    );
    let eta = symbol!("partial_placement_closure::eta");
    let epsilon = symbol!("partial_placement_closure::epsilon");
    let identity=GuardedMeasureIdentity{
        measure:format!("standalone-partial-placement-support-v1;input={};cuts={cuts:?};shifted={shifted:?}",input.identity()),
        support:"actual positive-support classification; polynomial completions; all valid occupation orders; required cuts positive; storage tails zero".into(),
        orientation:"exact future occupied geometry; ordinary virtual translation/rank witnesses in support certificate".into(),
        normalization:"native homogeneous identities, no value or endpoint shortcut".into(),
        branch:format!("support proof SHA256={proof_sha}; fixed eta/T regulator-first; labelwise high-D finite jets; ordinary virtual vacuum zero only after uniform line-regulator proof"),
        deformation:format!("fixed shells/H; native physical factors minus eta in slots{shifted:?}; origin masks cover every actual positive virtual support"),
    };
    let prepared = family
        .guarded_sources_with_options::<16>(
            epsilon,
            4,
            eta,
            &shifted,
            4096,
            vec![],
            identity,
            WeightedSourceOptions {
                policy: WeightedSourcePolicy::PolynomialClosure,
                positive_compact_energy_powers: false,
                free_virtual_zero_sectors: false,
            },
        )
        .unwrap();
    // The frozen RustFlow library links the old native crate while discovery
    // uses the isolated verified-composition build. Convert only immutable
    // source values, preserving every index, polynomial, guard and condition;
    // no rule/proof object crosses this crate-version boundary.
    let convert_domain = |d: &symbolica_amflow::finite_density::guarded::IndexDomain<16>| {
        IndexDomain::new(std::array::from_fn(|i| {
            IndexBounds::new(d.bounds()[i].lower(), d.bounds()[i].upper()).unwrap()
        }))
        .unwrap()
    };
    let admitted = convert_domain(prepared.deformation.admitted_domain());
    let original = prepared.context.sources();
    assert!(original.zero_domains().is_empty());
    let mut positive_cuts = *admitted.bounds();
    for &cut in &cuts {
        positive_cuts[cut] = IndexBounds::new(Some(1), None).unwrap()
    }
    let mut zeros = Vec::new();
    let mut descriptions = Vec::new();
    for shell in family.shells() {
        let mut b = positive_cuts;
        b[shell.lower_slot] = IndexBounds::new(Some(1), None).unwrap();
        zeros.push(IndexDomain::new(b).unwrap());
        descriptions.push(format!(
            "all positive virtual supports certified: lower-origin loop{}; finite eta only",
            shell.loop_index
        ))
    }
    for forced in &zero_forced {
        let mut b = positive_cuts;
        for &slot in forced {
            assert!(!cuts.contains(&slot));
            b[slot] = IndexBounds::new(None, Some(0)).unwrap()
        }
        zeros.push(IndexDomain::new(b).unwrap());
        descriptions.push(format!("separated ordinary virtual zero: forced nonpositive physical slots{forced:?}; exact support certificate"))
    }
    let rows = original
        .sources()
        .iter()
        .zip(original.native_sources().rows())
        .map(|(s, r)| {
            GuardedSource::new(
                s.id.clone(),
                r.iter()
                    .map(|t| Term {
                        integral: Integral::new(std::array::from_fn(|i| {
                            let p = t.integral.powers()[i];
                            Power::new(p.is_symbolic(), p.value()).unwrap()
                        })),
                        coefficient: t.coefficient.clone(),
                    })
                    .collect(),
                convert_domain(&s.domain),
            )
            .with_nonzero_conditions(s.nonzero_conditions.clone())
        })
        .collect();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            original.measure_id(),
            original.roles().map(|r| match r {
                symbolica_amflow::finite_density::guarded::IndexRole::Ordinary => {
                    rustred::solver::guarded::IndexRole::Ordinary
                }
                symbolica_amflow::finite_density::guarded::IndexRole::RequiredCut => {
                    rustred::solver::guarded::IndexRole::RequiredCut
                }
                symbolica_amflow::finite_density::guarded::IndexRole::Occupation => {
                    rustred::solver::guarded::IndexRole::Occupation
                }
            }),
            *original.native_sources().index_variables(),
            rows,
        )
        .unwrap()
        .with_zero_domains(zeros.clone())
        .unwrap(),
    );
    let bounds = |d: &IndexDomain<16>| {
        d.bounds()
            .iter()
            .map(|b| (b.lower(), b.upper()))
            .collect::<Vec<_>>()
    };
    // Exact finite boundary examples validate both positive and negative source guards.
    let mut bulk = [0; 16];
    bulk[..5].fill(1);
    assert!(!sources.is_zero(&bulk));
    for shell in family.shells() {
        let mut p = bulk;
        p[shell.lower_slot] = 1;
        assert!(sources.is_zero(&p));
        for &slot in &shifted {
            let mut q = p;
            q[slot] = 0;
            assert!(sources.is_zero(&q))
        }
    }
    for forced in &zero_forced {
        let mut p = bulk;
        for &slot in forced {
            p[slot] = 0
        }
        assert!(sources.is_zero(&p));
        for slot in 5..9 {
            let mut q = p;
            q[slot] = 1;
            assert!(!sources.is_zero(&q))
        }
        for slot in prepared.deformation.physical_arity()..16 {
            let mut q = p;
            q[slot] = 1;
            assert!(!sources.is_zero(&q))
        }
    }
    let measure = family.deformed_measure::<16>(eta, &shifted).unwrap();
    // Ensure selected deformation derivative agrees with the later native driver.
    for target in &prepared.targets {
        for point in target.keys() {
            let d = prepared.deformation.derivative(*point).unwrap();
            assert_eq!(d, super::derivative(*point, &shifted))
        }
    }
    let metadata = serde_json::json!({"scope":"Fixture-bound mathematical support certificate for a fresh native source-only closure pilot. No production flow or boundary admission.","input_identity":input.identity(),"cuts":cuts,"shifted":shifted,"proof_binding_sha256":proof_sha,"source_identity":sources.measure_id(),"source_count":sources.sources().len(),"raw_ward_omitted":true,
        "zero_domains":zeros.iter().map(bounds).collect::<Vec<_>>(),"zero_proofs":descriptions,"admitted_domain":bounds(&admitted),"physical_arity":prepared.deformation.physical_arity(),
        "factor_basis":measure.factors().iter().map(ToString::to_string).collect::<Vec<_>>(),"coordinates":family.coordinates().iter().map(ToString::to_string).collect::<Vec<_>>(),"inverse_routing":family.inverse_routing().iter().map(|r|r.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "targets":prepared.targets.iter().map(super::display).collect::<Vec<_>>(),"original_sources":sources.sources().iter().zip(sources.native_sources().rows()).map(|(s,row)|serde_json::json!({"id":s.id,"domain":bounds(&s.domain),"conditions":s.nonzero_conditions.iter().map(|c|c.to_expression().to_string()).collect::<Vec<_>>(),"terms":row.iter().map(|t|serde_json::json!({"indices":t.integral.powers().iter().map(|p|(p.is_symbolic(),p.value())).collect::<Vec<_>>(),"coefficient":t.coefficient.to_expression().to_string()})).collect::<Vec<_>>()})).collect::<Vec<_>>()});
    std::fs::write(
        out.join("source-system.json"),
        serde_json::to_vec_pretty(&metadata).unwrap(),
    )
    .unwrap();
    let empty = GuardedProgram::new(sources.clone(), vec![], []).unwrap();
    std::fs::write(
        out.join("fresh-source-program.bin"),
        empty.encode_native(Default::default()).unwrap(),
    )
    .unwrap();
    Prepared {
        sources,
        targets: prepared.targets,
        shifted,
        admitted,
    }
}
