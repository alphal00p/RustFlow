use super::tests::{restricted_sample, sample};
use super::*;
use crate::algebra::CoefficientContext;
use crate::solver::Term;
use crate::solver::guarded::{GuardedSource, IndexBounds, IndexDomain};

fn equivalent(left: &GuardedProgram<1>, right: &GuardedProgram<1>) {
    assert_eq!(
        left.encode_native(Default::default()).unwrap(),
        right.encode_native(Default::default()).unwrap()
    );
    for point in [[0], [1], [2], [3], [7]] {
        let a = left.reduce(point, Default::default()).unwrap();
        let b = right.reduce(point, Default::default()).unwrap();
        assert_eq!(a.terms, b.terms);
        assert_eq!(a.nonzero_conditions, b.nonzero_conditions);
        assert_eq!(a.rule_applications, b.rule_applications);
        assert_eq!(format!("{:?}", a.unresolved), format!("{:?}", b.unresolved));
    }
}

#[test]
fn verified_union_preserves_precedence_conditions_and_native_bytes() {
    let first = || restricted_sample(IndexBounds::new(Some(1), None).unwrap(), 1);
    let second = || restricted_sample(IndexBounds::new(Some(2), None).unwrap(), 2);
    let a = first();
    let b = second();
    assert!(!Arc::ptr_eq(a.sources(), b.sources()));
    let verified = a.union_verified(b, [[0]], 2).unwrap();
    let replayed = first().union_replayed(second(), [[0]], 2).unwrap();
    equivalent(&verified, &replayed);
    assert_eq!(
        verified.apply(&[3]).unwrap().status,
        GuardedApplicationStatus::Applied { rule: 0 }
    );
    assert_eq!(verified.terminals(), &BTreeSet::from([[0]]));
    let loaded = GuardedProgram::decode_generated(
        &verified.encode_native(Default::default()).unwrap(),
        verified.sources().clone(),
        Default::default(),
    )
    .unwrap();
    equivalent(&verified, &loaded);
}

#[test]
fn verified_terminal_rebind_replaces_stale_terminals_and_preserves_conditions() {
    for terminals in [vec![], vec![[1]], vec![[0]]] {
        let verified = sample("verified-rebind")
            .with_terminals_verified(terminals.clone(), 1)
            .unwrap();
        let replayed = sample("verified-rebind")
            .with_terminals_replayed(terminals, 1)
            .unwrap();
        equivalent(&verified, &replayed);
    }
    let result = sample("verified-no-stop")
        .with_terminals_verified([], 1)
        .unwrap()
        .reduce([1], Default::default())
        .unwrap();
    assert_eq!(result.unresolved.len(), 1);
    assert_eq!(result.unresolved[0].integral, [0]);
    assert!(!result.nonzero_conditions.is_empty());
}

#[test]
fn verified_composition_counts_duplicates_and_checks_rebind_budget() {
    assert!(
        sample("verified-budget")
            .union_verified(sample("verified-budget"), [[0]], 1)
            .is_err()
    );
    let p = sample("verified-budget")
        .union_verified(sample("verified-budget"), [[0]], 2)
        .unwrap();
    assert_eq!(p.rules().len(), 2);
    assert!(p.with_terminals_verified([[0]], 1).is_err());
}

#[test]
fn verified_union_accepts_empty_inputs_without_imposing_an_order() {
    for empty_first in [false, true] {
        let program = sample("verified-empty");
        let empty = GuardedProgram::new(program.sources().clone(), vec![], []).unwrap();
        let verified = if empty_first {
            empty.union_verified(program, [[0]], 1)
        } else {
            program.union_verified(empty, [[0]], 1)
        }
        .unwrap();
        equivalent(&verified, &sample("verified-empty"));
    }
    let source_owner = sample("verified-empty-both");
    let empty = || GuardedProgram::new(source_owner.sources().clone(), vec![], []).unwrap();
    let both = empty().union_verified(empty(), [[1]], 0).unwrap();
    assert!(both.rules().is_empty());
    assert_eq!(
        both.apply(&[1]).unwrap().status,
        GuardedApplicationStatus::Terminal
    );
}

fn cut_program() -> GuardedProgram<2> {
    let context = CoefficientContext::new(["verified_cut_n", "verified_surface_n"]);
    let roles = [IndexRole::RequiredCut, IndexRole::Occupation];
    let sources = GuardedSourceSystem::new(
        "verified-terminals",
        roles,
        [0, 1],
        vec![GuardedSource::new(
            "cut-row",
            vec![Term {
                integral: Integral::symbolic([0, 0]).unwrap(),
                coefficient: context.one().numerator,
            }],
            IndexDomain::for_roles(&roles),
        )],
    )
    .unwrap()
    .with_zero_domains(vec![
        IndexDomain::new([IndexBounds::fixed(2), IndexBounds::fixed(0)]).unwrap(),
    ])
    .unwrap();
    GuardedProgram::new(Arc::new(sources), vec![], []).unwrap()
}

#[test]
fn verified_apis_reject_invalid_and_declared_zero_terminals() {
    for terminal in [[0, 0], [-1, 0], [1, -1], [2, 0]] {
        assert!(
            cut_program()
                .with_terminals_verified([terminal], 0)
                .is_err()
        );
        assert!(
            cut_program()
                .union_verified(cut_program(), [terminal], 0)
                .is_err()
        );
    }
    let p = cut_program().with_terminals_verified([[1, 0]], 0).unwrap();
    assert_eq!(
        p.apply(&[1, 0]).unwrap().status,
        GuardedApplicationStatus::Terminal
    );
}
