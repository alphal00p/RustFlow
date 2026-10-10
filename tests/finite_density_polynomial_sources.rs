//! Compare the public polynomial factory with the frozen source-only 62-row
//! experiment. Historical rules and numerical answers are never imported.
#![allow(dead_code)]

use bincode::Decode;
use rustred::persistence::{
    BinaryIoLimits, CoefficientId, DecodedCoefficientTable, SectionTag, inspect_program,
};
use rustred::solver::guarded::GuardedProgram;
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{
    DensityInput,
    guarded::{GuardedMeasureIdentity, IndexRole},
    massless_endpoint::MasslessFlowEvidence,
    preparation::{WeightedSourceOptions, WeightedSourcePolicy},
};

type Domain = Vec<(Option<i64>, Option<i64>)>;
type Label = Vec<(bool, i16)>;
#[derive(Decode)]
struct StoredTerm {
    integral: Label,
    coefficient: usize,
}
#[derive(Decode)]
struct StoredSource {
    id: String,
    domain: Domain,
    conditions: Vec<usize>,
    terms: Vec<StoredTerm>,
}
#[derive(Decode)]
struct StoredDependency {
    row: usize,
    integral: Label,
    shifts: Vec<i16>,
}
#[derive(Decode)]
struct DiscardedRule {
    fixed: Vec<Option<i16>>,
    target: Label,
    rhs: Vec<StoredTerm>,
    sources: Vec<StoredDependency>,
    domain: Domain,
    discovery_domain: Domain,
    conditions: Vec<usize>,
    sector: Vec<bool>,
    permutation: Option<Vec<usize>>,
}
#[derive(Decode)]
struct StoredProgram {
    schema: String,
    measure: String,
    roles: Vec<u8>,
    indices: Vec<usize>,
    sources: Vec<StoredSource>,
    zero_domains: Vec<Domain>,
    rules: Vec<DiscardedRule>,
    terminals: Vec<Vec<i64>>,
}

fn substitute(a: &Atom, replacements: &BTreeMap<Atom, Atom>) -> Atom {
    a.replace_map(|view, _, out| {
        if !matches!(view, AtomView::Var(_)) {
            return;
        }
        if let Some(value) = replacements.get(&view.to_owned()) {
            **out = value.clone();
        }
    })
}

fn rational_ratio(a: &Atom, b: &Atom) -> Option<Rational> {
    let ratio = (a / b).together().cancel();
    Rational::try_from(ratio.as_view())
        .ok()
        .filter(|r| !r.is_zero())
}

fn equivalent_conditions(a: &[Atom], b: &[Atom]) {
    for condition in a.iter().chain(b) {
        assert!(!condition.clone().together().cancel().is_zero());
    }
    // Preserve every actual condition in the report. Only nonzero rational
    // rescaling, repetition and verified nonzero rational constants may differ.
    for (left, right) in [(a, b), (b, a)] {
        for condition in left {
            if let Ok(value) = Rational::try_from(condition.as_view()) {
                assert!(!value.is_zero());
                continue;
            }
            assert!(
                right
                    .iter()
                    .any(|other| rational_ratio(condition, other).is_some()),
                "unmatched retained condition {condition}"
            );
        }
    }
}

fn identity(input: &str) -> GuardedMeasureIdentity {
    GuardedMeasureIdentity {
        measure: format!("public polynomial source regression; input={input}; cut=[0]"),
        support: "original polynomial completion domain and sealed origin".into(),
        orientation: "actual future occupied routing".into(),
        normalization: "homogeneous native sources; no period values".into(),
        branch: "fixed-positive-eta source theorem bound by sealed evidence".into(),
        deformation: "all uncut physical factors D-eta; distributions fixed".into(),
    }
}

#[test]
fn public_polynomial_factory_matches_frozen_complete_source_corpus() {
    const CORPUS: &[u8] = include_bytes!(
        "../reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe/strongest62-boost-ward/program.bin"
    );
    let limits = BinaryIoLimits::default();
    let envelope = inspect_program(CORPUS, limits).unwrap();
    let payload = envelope.section(SectionTag::PROGRAM).unwrap();
    let (stored, consumed): (StoredProgram, usize) =
        bincode::decode_from_slice(payload, bincode::config::standard()).unwrap();
    assert_eq!(consumed, payload.len());
    assert_eq!(stored.schema, "rustred.guarded-source-program.v2");
    assert_eq!(stored.sources.len(), 62);
    // Decode the explicit record shape only to reach source data. No stored
    // candidate, rule, terminal or expected integral value is passed to RustRed.
    let discarded_rule_count = stored.rules.len();
    assert_eq!(discarded_rule_count, 19);
    assert!(stored.terminals.is_empty());
    let table = DecodedCoefficientTable::import_generated_normalized(
        envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),
        envelope.section(SectionTag::COEFFICIENTS).unwrap(),
        limits,
    )
    .unwrap();
    let coefficient = |id| {
        let c = table
            .coefficient(CoefficientId::try_from_index(id).unwrap())
            .unwrap();
        assert!(c.denominator.is_one());
        c.numerator.clone()
    };
    let old_variables = coefficient(stored.sources[0].terms[0].coefficient)
        .variables()
        .clone();
    let symbols = old_variables
        .iter()
        .map(|v| match v {
            PolyVariable::Symbol(s) => *s,
            _ => panic!("source corpus has a non-symbol coefficient variable"),
        })
        .collect::<Vec<_>>();
    assert_eq!(stored.indices.len(), 12);
    let parameters = symbols
        .iter()
        .enumerate()
        .filter(|(i, _)| !stored.indices.contains(i))
        .map(|(_, s)| *s)
        .collect::<Vec<_>>();
    assert_eq!(parameters.len(), 2);
    let parameter = |name: &str| {
        let values = parameters
            .iter()
            .copied()
            .filter(|s| Atom::var(*s).to_string() == name)
            .collect::<Vec<_>>();
        assert_eq!(values.len(), 1);
        values[0]
    };
    let (epsilon, eta) = (parameter("epsilon"), parameter("eta"));
    let definition: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massless_three_loop_chain.json"
    ))
    .unwrap();
    let input = definition.prepare().unwrap();
    assert!(stored.measure.contains(input.identity()));
    let family = input.occupied_cut(&[0], 4096).unwrap().at_physical_masses();
    assert_eq!(family.loops(), 3);
    assert_eq!(family.factors().len(), 11);
    let shifted = (0..family.physical_slots())
        .filter(|&slot| family.roles()[slot] == IndexRole::Ordinary)
        .collect::<Vec<_>>();
    assert_eq!(shifted, [1, 2, 3, 4]);
    let options = WeightedSourceOptions {
        policy: WeightedSourcePolicy::PolynomialClosure,
        positive_compact_energy_powers: false,
        free_virtual_zero_sectors: true,
    };
    let origin = MasslessFlowEvidence::new(&input, &family, &shifted, options).unwrap();
    let prepared = family
        .guarded_sources_with_massless_origin::<12>(
            epsilon,
            4,
            eta,
            &shifted,
            4096,
            vec![],
            identity(input.identity()),
            options,
            &origin,
        )
        .unwrap();
    let context = &prepared.context;
    let native = context.sources();
    assert_eq!(context.physical_arity(), 11);
    assert_eq!(native.sources().len(), 62);
    assert_ne!(native.measure_id(), stored.measure);
    assert!(native.measure_id().contains("polynomial-closure-v1"));
    let actual_roles = native
        .roles()
        .iter()
        .map(|r| match r {
            IndexRole::Ordinary => 0,
            IndexRole::RequiredCut => 1,
            IndexRole::Occupation => 2,
        })
        .collect::<Vec<_>>();
    assert_eq!(actual_roles, stored.roles);
    let bounds = |d: &rustred::solver::guarded::IndexDomain<12>| {
        d.bounds()
            .iter()
            .map(|b| (b.lower(), b.upper()))
            .collect::<Domain>()
    };
    assert_eq!(
        native.zero_domains().iter().map(bounds).collect::<Vec<_>>(),
        stored.zero_domains
    );
    let fresh_variables = native.native_sources().coefficient_variables();
    let replacements = stored
        .indices
        .iter()
        .zip(native.native_sources().index_variables())
        .map(|(&old, &new)| {
            let PolyVariable::Symbol(new) = &fresh_variables[new] else {
                panic!()
            };
            (Atom::var(symbols[old]), Atom::var(*new))
        })
        .collect::<BTreeMap<_, _>>();
    let mut report_rows = Vec::new();
    for (ordinal, ((saved, actual), row)) in stored
        .sources
        .iter()
        .zip(native.sources())
        .zip(native.native_sources().rows())
        .enumerate()
    {
        match ordinal {
            0..=51 => assert_eq!(saved.id, actual.id),
            52..=55 => {
                assert_eq!(
                    saved.id,
                    format!(
                        "shifted-polynomial-temporal-U/lorentz/0/3/occupation-case-{}",
                        ordinal - 52
                    )
                );
                assert_eq!(
                    actual.id,
                    format!("polynomial-closure-v1/temporal-U/0/6/case-{}", ordinal - 52)
                );
            }
            56..=59 => {
                assert_eq!(
                    saved.id,
                    format!(
                        "shifted-polynomial-global-boost/global-boost/occupied-0/occupation-case-{}",
                        ordinal - 56
                    )
                );
                assert_eq!(
                    actual.id,
                    format!(
                        "polynomial-closure-v1/common-boost/0/6/case-{}",
                        ordinal - 56
                    )
                );
            }
            60..=61 => {
                assert_eq!(
                    saved.id,
                    format!(
                        "polynomial-singleton-raw-normalized-boost/{}",
                        if ordinal == 60 {
                            "bulk"
                        } else {
                            "positive-upper"
                        }
                    )
                );
                assert_eq!(
                    actual.id,
                    format!(
                        "raw-polynomial-singleton-Ward-v1/upper-{}",
                        if ordinal == 60 { "zero" } else { "positive" }
                    )
                );
            }
            _ => unreachable!(),
        }
        assert_eq!(
            bounds(&actual.domain),
            saved.domain,
            "source domain {ordinal}"
        );
        let expected = saved
            .terms
            .iter()
            .map(|t| {
                (
                    t.integral.clone(),
                    substitute(&coefficient(t.coefficient).to_expression(), &replacements),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let actual_row = row
            .iter()
            .map(|t| {
                (
                    t.integral
                        .powers()
                        .iter()
                        .map(|p| (p.is_symbolic(), p.value()))
                        .collect::<Label>(),
                    t.coefficient.to_expression(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        assert_eq!(expected.len(), saved.terms.len());
        assert_eq!(actual_row.len(), row.len());
        assert_eq!(
            actual_row.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>(),
            "source images {ordinal}"
        );
        let first = expected.keys().next().unwrap();
        let scale = rational_ratio(&actual_row[first], &expected[first])
            .unwrap_or_else(|| panic!("nonconstant row scaling at {ordinal}"));
        if ordinal < 52 {
            assert_eq!(scale, Rational::one(), "legacy row rescaled {ordinal}");
        }
        for (label, want) in &expected {
            assert!(
                (&actual_row[label] - Atom::num(scale.clone()) * want)
                    .expand()
                    .together()
                    .cancel()
                    .is_zero(),
                "coefficient mismatch in source {ordinal}, {label:?}"
            );
        }
        let expected_conditions = saved
            .conditions
            .iter()
            .map(|&id| substitute(&coefficient(id).to_expression(), &replacements))
            .collect::<Vec<_>>();
        let actual_conditions = actual
            .nonzero_conditions
            .iter()
            .map(|c| c.to_expression())
            .collect::<Vec<_>>();
        equivalent_conditions(&actual_conditions, &expected_conditions);
        report_rows.push(serde_json::json!({
            "ordinal":ordinal,"frozen_id":saved.id,"factory_id":actual.id,
            "factory_over_frozen_nonzero_rational_scale":scale.to_string(),"domain":saved.domain,
            "frozen_conditions":expected_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "factory_conditions":actual_conditions.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "terms":actual_row.iter().map(|(label,c)|serde_json::json!({"label":label,"coefficient":c.to_canonical_string()})).collect::<Vec<_>>(),
        }));
    }
    let empty = GuardedProgram::new(native.clone(), vec![], []).unwrap();
    let encoded = empty.encode_native(limits).unwrap();
    let decoded = GuardedProgram::decode_generated(&encoded, native.clone(), limits).unwrap();
    assert!(decoded.rules().is_empty());
    assert_eq!(
        decoded.sources().native_sources().rows(),
        native.native_sources().rows()
    );
    let mut invalid = [0_i64; 12];
    invalid[0] = 1;
    invalid[6] = 1;
    assert!(!native.is_zero(&invalid));
    assert!(prepared.deformation.derivative(invalid).is_err());
    if let Some(path) = std::env::var_os("RUSTFLOW_POLYNOMIAL_SOURCE_REPORT") {
        let path = std::path::PathBuf::from(path);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        std::fs::write(path, serde_json::to_vec_pretty(&serde_json::json!({
            "scope":"Public factory equals complete frozen source corpus after explicit index renaming and nonzero rational row normalization. No historical rules, periods or closure assertions.",
            "input_identity":input.identity(),"source_corpus_blake3":blake3::hash(CORPUS).to_hex().to_string(),
            "source_corpus_bytes":CORPUS.len(),"historical_rules_discarded":discarded_rule_count,
            "historical_rules_imported":false,"source_count":62,"physical_arity":11,"capacity":12,
            "roles":actual_roles,"zero_domains":stored.zero_domains,"options":options,
            "factory_measure_id":native.measure_id(),"original_measure_id":stored.measure,
            "all_original_legacy_rows_preserved":true,"positive_completion_rejected":true,
            "source_context_roundtrip":true,"rows":report_rows,
        })).unwrap()).unwrap();
    }
}
