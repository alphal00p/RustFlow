use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{
    Error, Integral,
    reduction::{LinearCombination, Reduction},
};
fn expand_many(
    reduction: &Reduction,
    targets: &[Integral],
) -> symbolica_amflow::Result<BTreeMap<Integral, LinearCombination>> {
    reduction.expand_many(targets)
}
fn integral(i: i16) -> Integral {
    Integral(vec![i])
}
fn row(terms: &[(i16, i64)]) -> LinearCombination {
    terms
        .iter()
        .map(|&(i, c)| (integral(i), Atom::num(c)))
        .collect()
}
fn diamond(depth: i16) -> Reduction {
    let terminal = integral(2 * depth + 1);
    let mut reduction = Reduction {
        residuals: vec![terminal.clone()],
        ..Default::default()
    };
    for level in 0..depth {
        let children = if level + 1 == depth {
            BTreeMap::from([(terminal.clone(), Atom::num(2))])
        } else {
            row(&[(2 * level + 3, 1), (2 * level + 4, 1)])
        };
        for side in [1, 2] {
            reduction
                .rules
                .insert(integral(2 * level + side), children.clone());
        }
    }
    reduction
}
#[test]
fn long_chain_does_not_use_call_stack() {
    let mut table = Reduction {
        residuals: vec![integral(30000)],
        ..Default::default()
    };
    for i in 1..30000 {
        table.rules.insert(integral(i), row(&[(i + 1, 1)]));
    }
    let result = expand_many(&table, &[integral(1), integral(15000)]).unwrap();
    for root in [integral(1), integral(15000)] {
        assert_eq!(result[&root], row(&[(30000, 1)]));
    }
}
#[test]
fn shared_diamonds_and_multiple_roots_are_memoized() {
    let shallow = diamond(8);
    let targets = [integral(1), integral(2), integral(7), integral(1)];
    let result = expand_many(&shallow, &targets).unwrap();
    assert_eq!(result.len(), 3);
    for target in &targets {
        assert_eq!(result[target], shallow.expand(target).unwrap());
    }
    let deep = diamond(100);
    let result = expand_many(&deep, &[integral(1), integral(2)]).unwrap();
    assert_eq!(result[&integral(1)][&integral(201)], Atom::num(2).pow(100));
    assert_eq!(result[&integral(1)], result[&integral(2)]);
}
#[test]
fn only_reachable_cycles_are_rejected_even_with_zero_coefficients() {
    let mut table = Reduction {
        rules: BTreeMap::from([(integral(1), row(&[(2, 1)])), (integral(3), row(&[(3, 1)]))]),
        residuals: vec![integral(2)],
        ..Default::default()
    };
    assert_eq!(
        expand_many(&table, &[integral(1)]).unwrap()[&integral(1)],
        table.expand(&integral(1)).unwrap()
    );
    assert!(expand_many(&table, &[]).unwrap().is_empty());
    table
        .rules
        .get_mut(&integral(1))
        .unwrap()
        .insert(integral(3), Atom::new());
    assert!(matches!(
        expand_many(&table, &[integral(1)]),
        Err(Error::IncompleteReduction(_))
    ));
    table.rules.remove(&integral(3)); // undeclared leaf is also checked through coefficient zero
    assert!(matches!(
        expand_many(&table, &[integral(1)]),
        Err(Error::IncompleteReduction(_))
    ));
}
#[test]
fn zero_rules_residual_priority_nontriangular_rows_and_mutations_match() {
    let parent = Integral(vec![1, 0]);
    let child = Integral(vec![1, 1]);
    let mut table = Reduction {
        rules: BTreeMap::from([(
            parent.clone(),
            BTreeMap::from([(child.clone(), parse!("(audit::x+1)/(audit::x-2)"))]),
        )]),
        residuals: vec![child.clone()],
        ..Default::default()
    };
    let targets = [parent.clone(), child.clone()];
    let result = expand_many(&table, &targets).unwrap();
    for root in &targets {
        assert_eq!(result[root], table.expand(root).unwrap());
    }
    table.rules.insert(child.clone(), BTreeMap::new()); // rule takes priority over residual declaration
    assert!(
        expand_many(&table, &targets)
            .unwrap()
            .values()
            .all(BTreeMap::is_empty)
    );
    table.rules.remove(&child);
    table
        .rules
        .get_mut(&parent)
        .unwrap()
        .insert(child.clone(), Atom::num(7));
    assert_eq!(
        expand_many(&table, &targets).unwrap()[&parent][&child],
        Atom::num(7)
    );
}
#[test]
fn exact_cancellation_and_uncovered_leaves_match_existing_table_semantics() {
    let mut table = Reduction {
        rules: BTreeMap::from([
            (integral(1), row(&[(2, 1), (3, -1)])),
            (integral(2), row(&[(4, 1)])),
            (integral(3), row(&[(4, 1)])),
        ]),
        residuals: vec![integral(4)],
        ..Default::default()
    };
    assert!(expand_many(&table, &[integral(1)]).unwrap()[&integral(1)].is_empty());
    table.residuals.clear();
    assert!(matches!(
        expand_many(&table, &[integral(1)]),
        Err(Error::IncompleteReduction(_))
    ));
}

#[test]
fn native_batched_tables_preserve_physical_differential_closure_and_guards() {
    use symbolica_amflow::{
        FlowOptions, IntegralFamily, PreparedPhysicalFamily, Propagator, ReductionBackend, Result,
        RunContext, TableBackend,
    };
    struct BatchedTable(Reduction);
    impl ReductionBackend for BatchedTable {
        fn identity(&self) -> String {
            "native batched supplied table".into()
        }
        fn reduce(
            &self,
            _: &IntegralFamily,
            targets: &[Integral],
            context: &RunContext,
        ) -> Result<Reduction> {
            context.cancellation.check()?;
            Ok(Reduction {
                rules: expand_many(&self.0, targets)?,
                residuals: self.0.residuals.clone(),
                nonzero_conditions: self.0.nonzero_conditions.clone(),
            })
        }
    }
    let (s, epsilon) = symbol!("expansion_closure::s", "expansion_closure::epsilon");
    let family = IntegralFamily {
        name: "tadpole supplied chained IBPs".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -Atom::var(s),
            scalar_products: vec![Atom::one()],
        }],
        physical_propagators: 1,
        epsilon,
        dimension: 4,
    };
    let table = Reduction {
        rules: BTreeMap::from([
            (
                integral(2),
                BTreeMap::from([(
                    integral(3),
                    parse!("-2*expansion_closure::s/expansion_closure::epsilon"),
                )]),
            ),
            (
                integral(3),
                BTreeMap::from([(
                    integral(1),
                    parse!(
                        "expansion_closure::epsilon*(expansion_closure::epsilon-1)/(2*expansion_closure::s^2)"
                    ),
                )]),
            ),
        ]),
        residuals: vec![integral(1)],
        nonzero_conditions: vec![Atom::var(s), Atom::var(epsilon)],
    };
    let baseline = TableBackend {
        name: "chained exact IBPs".into(),
        reduction: table.clone(),
    };
    let batched = BatchedTable(table);
    let options = FlowOptions::default();
    let context = RunContext::default();
    let old = PreparedPhysicalFamily::new(
        &family,
        &[integral(2)],
        &[s],
        &baseline,
        &options,
        "regular tadpole",
        &context,
    )
    .unwrap();
    let new = PreparedPhysicalFamily::new(
        &family,
        &[integral(2)],
        &[s],
        &batched,
        &options,
        "regular tadpole",
        &context,
    )
    .unwrap();
    assert_eq!(new.basis(), old.basis());
    assert_eq!(new.target_reductions(), old.target_reductions());
    assert_eq!(new.nonzero_conditions(), old.nonzero_conditions());
    assert_eq!(
        new.flow().system().derivatives,
        old.flow().system().derivatives
    );
    assert_eq!(new.flow().identity().key(), old.flow().identity().key());
    assert!(
        (&new.flow().system().derivatives[&s][0][0]
            - parse!("(1-expansion_closure::epsilon)/expansion_closure::s"))
        .together()
        .cancel()
        .is_zero()
    );
}

#[test]
fn flat_rows_cancel_exactly_but_still_authenticate_zero_rhs() {
    let x = parse!("flat_table::x");
    let weight = (&x * &x - Atom::one()) / (&x - Atom::one());
    let mut table = Reduction {
        rules: BTreeMap::from([(integral(1), BTreeMap::from([(integral(2), weight)]))]),
        residuals: vec![integral(2)],
        ..Default::default()
    };
    assert_eq!(
        expand_many(&table, &[integral(1), integral(2)]).unwrap()[&integral(1)],
        table.expand(&integral(1)).unwrap()
    );
    table
        .rules
        .get_mut(&integral(1))
        .unwrap()
        .insert(integral(3), Atom::new());
    assert!(matches!(
        expand_many(&table, &[integral(1)]),
        Err(Error::IncompleteReduction(_))
    ));
}
