//! Shared bounded integer potentials for diagonal power transformations.
use crate::{Error, Result, RunContext};
use linnet::half_edge::{HedgeGraph, builder::HedgeGraphBuilder, involution::HedgePair};

/// Generalizes the constraint relaxation formerly inside diagonal_fuchsian_form.
/// `orders[i][j] + shifts[j] - shifts[i] >= -allowed_pole`.
/// Zero entries are absent edges; Linnet owns all graph storage/incidence.
pub(crate) fn integer_potentials(
    orders: &[Vec<Option<i64>>],
    allowed_pole: i64,
    max_span: i64,
    context: Option<&RunContext>,
) -> Result<Vec<i64>> {
    let n = orders.len();
    if n == 0 || orders.iter().any(|row| row.len() != n) || max_span < 0 {
        return Err(Error::InvalidInput(
            "integer shearing dimensions or limit".into(),
        ));
    }
    let mut builder = HedgeGraphBuilder::<i64, usize>::new();
    let nodes = (0..n).map(|i| builder.add_node(i)).collect::<Vec<_>>();
    for (i, row) in orders.iter().enumerate() {
        for (j, order) in row.iter().enumerate() {
            if let Some(order) = order {
                let bound = order
                    .checked_add(allowed_pole)
                    .ok_or_else(|| Error::Limit("integer shearing bound overflow".into()))?;
                builder.add_edge(nodes[j], nodes[i], bound, true);
            }
        }
    }
    let graph: HedgeGraph<i64, usize> = builder.build();
    let mut shifts = vec![0_i64; n];
    let mut predecessor = vec![None; n];
    for iteration in 0..n {
        if let Some(context) = context {
            context.cancellation.check()?;
        }
        let mut changed = None;
        for (pair, _, edge) in graph.iter_edges() {
            let HedgePair::Paired { source, sink } = pair else {
                return Err(Error::Numerical(
                    "internal shearing graph has an open edge".into(),
                ));
            };
            let from = graph.node_id(source).0;
            let to = graph.node_id(sink).0;
            let bound = shifts[from]
                .checked_add(*edge.data)
                .ok_or_else(|| Error::Limit("integer shearing potential overflow".into()))?;
            if shifts[to] > bound {
                shifts[to] = bound;
                predecessor[to] = Some(from);
                changed = Some(to);
            }
        }
        let Some(mut at) = changed else {
            break;
        };
        if iteration == n - 1 {
            // The nth relaxation witnesses a negative cycle; trace only the
            // native graph's retained predecessor indices, not a second graph.
            for _ in 0..n {
                at = predecessor[at]
                    .ok_or_else(|| Error::Numerical("missing shearing cycle predecessor".into()))?;
            }
            let start = at;
            let mut cycle = vec![start];
            loop {
                at = predecessor[at]
                    .ok_or_else(|| Error::Numerical("missing shearing cycle predecessor".into()))?;
                cycle.push(at);
                if at == start {
                    break;
                }
                if cycle.len() > n + 1 {
                    return Err(Error::Numerical("invalid shearing cycle witness".into()));
                }
            }
            cycle.reverse();
            return Err(Error::Unsupported(format!(
                "diagonal power shearing has a negative constraint cycle {cycle:?}; a nondiagonal basis change may be required"
            )));
        }
    }
    let maximum = *shifts.iter().max().unwrap();
    let minimum = *shifts.iter().min().unwrap();
    if maximum
        .checked_sub(minimum)
        .is_none_or(|span| span > max_span)
    {
        return Err(Error::Limit(format!(
            "integer shearing exceeds {max_span} powers"
        )));
    }
    for shift in &mut shifts {
        *shift = shift
            .checked_sub(maximum)
            .ok_or_else(|| Error::Limit("integer shearing normalization overflow".into()))?;
    }
    // Check every original constraint after normalization, including loops.
    for (i, row) in orders.iter().enumerate() {
        for (j, order) in row.iter().enumerate() {
            if let Some(order) = order {
                let bound = shifts[j]
                    .checked_add(*order)
                    .and_then(|v| v.checked_add(allowed_pole))
                    .ok_or_else(|| Error::Limit("integer shearing validation overflow".into()))?;
                if shifts[i] > bound {
                    return Err(Error::Numerical(
                        "integer shearing failed its exact postcondition".into(),
                    ));
                }
            }
        }
    }
    Ok(shifts)
}
