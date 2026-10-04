use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn native_dependency_components_preserve_block_order_and_pivot_ties() {
    // Once component1 is emitted, 0 becomes ready and precedes independent2.
    // A FIFO topological queue would change the established public ordering.
    let matrix = vec![vec![0, 1, 0], vec![0, 0, 0], vec![0, 0, 1]];
    let system = DifferentialSystem {
        variable: symbol!("x"),
        matrix: matrix
            .into_iter()
            .map(|r| r.into_iter().map(Atom::num).collect())
            .collect(),
    };
    assert_eq!(system.blocks().unwrap(), vec![vec![1], vec![0], vec![2]]);
    let matrix = vec![
        vec![1, 0, 1, 0],
        vec![0, 1, 0, 0],
        vec![1, 1, 1, 0],
        vec![0, 0, 1, 1],
    ];
    let system = DifferentialSystem {
        variable: symbol!("x"),
        matrix: matrix
            .into_iter()
            .map(|r| r.into_iter().map(Atom::num).collect())
            .collect(),
    };
    assert_eq!(system.blocks().unwrap(), vec![vec![1], vec![0, 2], vec![3]]);
}
