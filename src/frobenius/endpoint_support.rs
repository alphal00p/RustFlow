//! Sufficient exact support test for finite projected coefficient limits.
//!
//! This is intentionally stronger than testing one numerically matched vector:
//! every direction in the rational lifted solution space must have a finite
//! physical projection. Root-constrained cancellations may therefore be refused.
//! No numerical threshold is used to decide whether a divergent sector exists.
use super::{PreparedFrobenius, valuation};
use crate::{Error, Result, RunContext};
use linnet::half_edge::{HedgeGraph, builder::HedgeGraphBuilder, involution::HedgePair};
use symbolica::prelude::*;

impl PreparedFrobenius {
    /// Establish finite physical projections from exact valuations and leading
    /// generalized eigenspaces. A possible nonpositive later recurrence term
    /// is refused rather than assumed to cancel in floating point.
    pub(crate) fn admit_endpoint_support(
        &self,
        physical_rows: usize,
        context: &RunContext,
    ) -> Result<()> {
        let n = self.normal.matrix.len();
        if physical_rows == 0 || physical_rows > n {
            return Err(Error::InvalidInput(
                "endpoint support projection dimensions".into(),
            ));
        }
        let z = Atom::var(self.source.variable);
        let mut builder = HedgeGraphBuilder::<i64, usize>::new();
        let nodes = (0..n).map(|i| builder.add_node(i)).collect::<Vec<_>>();
        for (i, row) in self.normal.matrix.iter().enumerate() {
            context.cancellation.check()?;
            for (j, entry) in row.iter().enumerate() {
                if entry.is_zero() {
                    continue;
                }
                let weight = valuation(&(entry * &z).together().cancel(), self.source.variable)?;
                if weight < 0 {
                    return Err(Error::Numerical(
                        "endpoint support requires a Fuchsian normalized connection".into(),
                    ));
                }
                // Native graph incidence is shared with integer shearing. A
                // zero-cost edge includes residue/inverse reachability; a
                // positive edge gives a lower bound on the recurrence degree.
                builder.add_edge(nodes[j], nodes[i], weight, true);
            }
        }
        let graph: HedgeGraph<i64, usize> = builder.build();
        let mut transform = vec![vec![None; n]; physical_rows];
        for (i, row) in self.transformation.iter().take(physical_rows).enumerate() {
            context.cancellation.check()?;
            for (j, entry) in row.iter().enumerate() {
                if entry.is_zero() {
                    continue;
                }
                let degree = valuation(entry, self.source.variable)?;
                // Dividing out its exact valuation leaves an ordinary Taylor
                // series; Symbolica owns extraction of the leading coefficient.
                let regular = (entry / z.clone().pow(degree)).together().cancel();
                let series = regular
                    .series(self.source.variable, 0, 1)
                    .map_err(|e| Error::Unsupported(e.to_string()))?;
                let leading = series
                    .terms()
                    .filter(|(power, _)| power.to_string() == "0")
                    .fold(Atom::new(), |sum, (_, value)| sum + value.to_owned());
                if leading.is_zero() {
                    return Err(Error::Numerical(
                        "endpoint transformation valuation has no exact leading coefficient".into(),
                    ));
                }
                transform[i][j] = Some((degree, leading));
            }
        }
        for space in &self.eigenspaces {
            let AtomView::Num(exponent) = space.exponent.as_view() else {
                return Err(Error::Unsupported("finite coefficient endpoint support requires exact rational exponents; dimensional-sector projection is separate".into()));
            };
            let symbolica::coefficient::Coefficient::Complex(exponent) =
                exponent.get_coeff_view().to_owned()
            else {
                return Err(Error::Unsupported(
                    "nonrational endpoint support exponent".into(),
                ));
            };
            if !exponent.im.is_zero() {
                return Err(Error::Unsupported(
                    "complex endpoint support exponent".into(),
                ));
            }
            for leading in &space.leading {
                context.cancellation.check()?;
                let mut degree = (0..n)
                    .map(|j| leading.iter().any(|row| !row[j].is_zero()).then_some(0_i64))
                    .collect::<Vec<_>>();
                // Domain-specific valuation propagation on native incidence.
                // All weights are nonnegative; n-1 rounds include every simple
                // path. Zero-cost cycles cannot lower a degree further.
                for _ in 1..n {
                    context.cancellation.check()?;
                    let mut changed = false;
                    for (pair, _, edge) in graph.iter_edges() {
                        let HedgePair::Paired { source, sink } = pair else {
                            return Err(Error::Numerical(
                                "internal endpoint support graph has an open edge".into(),
                            ));
                        };
                        let from = graph.node_id(source).0;
                        let to = graph.node_id(sink).0;
                        let Some(from_degree) = degree[from] else {
                            continue;
                        };
                        let candidate = from_degree.checked_add(*edge.data).ok_or_else(|| {
                            Error::Limit("endpoint support degree overflow".into())
                        })?;
                        if degree[to].is_none_or(|old| candidate < old) {
                            degree[to] = Some(candidate);
                            changed = true;
                        }
                    }
                    if !changed {
                        break;
                    }
                }
                for row in &transform {
                    context.cancellation.check()?;
                    let mut constant = vec![Atom::new(); leading.len()];
                    for (j, entry) in row.iter().enumerate() {
                        let (Some(recurrence_degree), Some((shift, coefficient))) =
                            (degree[j], entry)
                        else {
                            continue;
                        };
                        let shift = recurrence_degree.checked_add(*shift).ok_or_else(|| {
                            Error::Limit("endpoint projected degree overflow".into())
                        })?;
                        let power = &exponent.re + &Rational::from(shift);
                        if power < Rational::zero() {
                            return Err(Error::Accuracy("exact support permits a divergent physical endpoint power; constrained-sector evidence is required".into()));
                        }
                        if power.is_zero() {
                            if recurrence_degree != 0 {
                                return Err(Error::Accuracy("exact support permits a later nonpositive resonant endpoint sector; an exact recurrence or constrained-sector proof is required".into()));
                            }
                            for (log, vector) in leading.iter().enumerate() {
                                constant[log] += coefficient * &vector[j];
                            }
                        }
                    }
                    if constant
                        .iter()
                        .skip(1)
                        .any(|value| !value.together().cancel().is_zero())
                    {
                        return Err(Error::Accuracy("exact leading support contains a divergent endpoint logarithm; constrained-sector evidence is required".into()));
                    }
                }
            }
        }
        context.cancellation.check()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DifferentialSystem;

    fn prepared(matrix: &[&[&str]]) -> PreparedFrobenius {
        DifferentialSystem {
            variable: symbol!("endpoint_support_tests::z"),
            matrix: matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|value| {
                            Atom::parse(value, "endpoint_support_tests", Default::default())
                                .unwrap()
                        })
                        .collect()
                })
                .collect(),
        }
        .prepare_frobenius(&RunContext::default())
        .unwrap()
    }

    #[test]
    fn positive_power_logs_and_divergent_auxiliary_rows_are_distinct() {
        prepared(&[&["1/(2*z)", "1/z"], &["0", "1/(2*z)"]])
            .admit_endpoint_support(2, &RunContext::default())
            .unwrap();
        prepared(&[&["0", "0"], &["0", "-1/z"]])
            .admit_endpoint_support(1, &RunContext::default())
            .unwrap();
    }

    #[test]
    fn arbitrarily_tiny_leading_and_later_logs_remain_nonzero() {
        for matrix in [
            vec![vec!["0", "1/(10^300*z)"], vec!["0", "0"]],
            vec![vec!["0", "1/10^300"], vec!["0", "-1/z"]],
        ] {
            let rows = matrix.iter().map(Vec::as_slice).collect::<Vec<_>>();
            assert!(matches!(
                prepared(&rows).admit_endpoint_support(1, &RunContext::default()),
                Err(Error::Accuracy(_))
            ));
        }
    }
}
