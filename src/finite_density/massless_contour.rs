//! Exact positive-eta channel evidence, deliberately separate from numerical admission.
// This certificate does not replace positive_shells or enable an eta=0 projector.
use super::{DensityEdge, PreparedDensityInput};
use crate::{Error, Result};
use linnet::{
    half_edge::{HedgeGraph, NoData, builder::HedgeGraphBuilder},
    tree::{Forest, child_vec::ChildVecStore},
};
use serde::Serialize;
use symbolica::prelude::*;

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MasslessChannelKind {
    Zero,
    SingleLightlike,
    DifferenceOfFutureLightlike,
    Unproved,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct MasslessChannel {
    pub vertex_partition: [Vec<usize>; 2],
    /// In cut_slots order, multiplying future occupied momenta.
    pub future_cut_coefficients: Vec<i8>,
    pub kind: MasslessChannelKind,
}

#[derive(Clone, Debug, Serialize)]
pub struct MasslessChannelEvidence {
    pub input_identity: String,
    pub cut_slots: Vec<usize>,
    pub future_orientations: Vec<i8>,
    pub chemical_magnitudes: Vec<String>,
    pub shifted_physical_slots: Vec<usize>,
    pub physical_mass_squares: Vec<String>,
    pub inverse_routing: Vec<Vec<String>>,
    pub routing_determinant: String,
    pub candidate_partitions: usize,
    pub partition_budget: usize,
    pub channels: Vec<MasslessChannel>,
    pub certified: bool,
    pub scope: &'static str,
    pub endpoint_admitted: bool,
    pub numerical_evaluation_admitted: bool,
}

fn rational(atom: &Atom, what: &str) -> Result<Rational> {
    Rational::try_from(atom.as_view()).map_err(|_| {
        Error::Unsupported(format!(
            "massless channel certificate needs real rational {what}"
        ))
    })
}

fn induced_connected(
    vertices: usize,
    edges: &[DensityEdge],
    cuts: &[usize],
    side: &[bool],
    inside: bool,
) -> bool {
    if !side.iter().any(|&member| member == inside) {
        return false;
    }
    let mut builder = HedgeGraphBuilder::<(), ()>::new();
    let nodes = (0..vertices)
        .map(|v| (side[v] == inside).then(|| builder.add_node(())))
        .collect::<Vec<_>>();
    for (slot, edge) in edges.iter().enumerate() {
        if cuts.contains(&slot) {
            continue;
        }
        if let (Some(tail), Some(head)) = (nodes[edge.vertices[0]], nodes[edge.vertices[1]]) {
            builder.add_edge(tail, head, (), true);
        }
    }
    let graph: HedgeGraph<(), (), NoData, Forest<(), ChildVecStore<()>>> = builder.build();
    let isolated = graph
        .iter_nodes()
        .map(|(_, mut crown, _)| usize::from(crown.next().is_none()))
        .sum::<usize>();
    // Undirected native connectivity; directed SCCs would be incorrect here.
    graph.count_connected_components(&graph.full_filter()) + isolated == 1
}

/// Check every two-forest external channel of the admitted cut complement.
///
/// Every uncut physical slot must carry the same unit positive auxiliary mass;
/// `shifted_physical_slots` must come from the actual deformation descriptor.
/// This returns all unproved channels instead of calling them physical pinches.
/// It does not admit shell products, endpoint projection, or numerical flow.
pub fn massless_channel_evidence(
    input: &PreparedDensityInput,
    cuts: &[usize],
    shifted_physical_slots: &[usize],
    cut_budget: usize,
    partition_budget: usize,
) -> Result<MasslessChannelEvidence> {
    // Owns exact incidence/routing conservation, cut independence, complement
    // connectivity, future routing, nonzero chemical shifts and its determinant.
    let family = input.occupied_cut(cuts, cut_budget)?;
    let definition = input.input();
    let uncut = (0..definition.edges.len())
        .filter(|slot| !cuts.contains(slot))
        .collect::<Vec<_>>();
    if uncut.is_empty() {
        return Err(Error::Unsupported(
            "no uncut quadratic factors: use the compact polynomial terminal proof".into(),
        ));
    }
    if shifted_physical_slots != uncut {
        return Err(Error::Unsupported(
            "massless channel certificate requires precisely all uncut physical quadratics shifted by the same unit eta".into(),
        ));
    }
    for (slot, mass) in input.physical_masses().iter().enumerate() {
        let mass = rational(mass, "mass squares")?;
        if mass < Rational::zero() || (cuts.contains(&slot) && !mass.is_zero()) {
            return Err(Error::Unsupported(
                "massless channel certificate requires zero occupied masses and nonnegative uncut mass squares".into(),
            ));
        }
    }
    let potentials = definition
        .chemical_potentials
        .iter()
        .map(|text| {
            let atom = Atom::parse(text, "rustflow_density", Default::default())
                .map_err(|e| Error::InvalidInput(e.to_string()))?;
            rational(&atom, "chemical potentials")
        })
        .collect::<Result<Vec<_>>>()?;
    let future_orientations = cuts
        .iter()
        .map(|&slot| {
            let mu = definition.edges[slot]
                .charges
                .iter()
                .zip(&potentials)
                .fold(Rational::zero(), |sum, (&charge, mu)| {
                    sum + mu * &Rational::from(charge)
                });
            if mu < Rational::zero() { -1_i8 } else { 1_i8 }
        })
        .collect::<Vec<_>>();
    let vertices = definition.vertices;
    let exponent = u32::try_from(vertices - 1)
        .map_err(|_| Error::Limit("massless partition enumeration exponent".into()))?;
    let masks = 1_usize
        .checked_shl(exponent)
        .ok_or_else(|| Error::Limit("massless partition enumeration overflow".into()))?;
    let candidate_partitions = masks - 1;
    if candidate_partitions > partition_budget {
        return Err(Error::Limit(format!(
            "massless channel certificate needs {candidate_partitions} vertex partitions; budget {partition_budget}"
        )));
    }
    let mut channels = Vec::new();
    // Fix vertex zero inside. The last mask has no complement and is omitted.
    for mask in 0..candidate_partitions {
        let side = (0..vertices)
            .map(|v| v == 0 || (mask & (1_usize << (v - 1))) != 0)
            .collect::<Vec<_>>();
        if !induced_connected(vertices, &definition.edges, cuts, &side, true)
            || !induced_connected(vertices, &definition.edges, cuts, &side, false)
        {
            continue;
        }
        let coefficients = cuts
            .iter()
            .zip(&future_orientations)
            .map(|(&slot, &sign)| {
                let [tail, head] = definition.edges[slot].vertices;
                sign * (i8::from(side[head]) - i8::from(side[tail]))
            })
            .collect::<Vec<_>>();
        let nonzero = coefficients
            .iter()
            .copied()
            .filter(|&c| c != 0)
            .collect::<Vec<_>>();
        let kind = match nonzero.as_slice() {
            [] => MasslessChannelKind::Zero,
            [-1] | [1] => MasslessChannelKind::SingleLightlike,
            [-1, 1] | [1, -1] => MasslessChannelKind::DifferenceOfFutureLightlike,
            _ => MasslessChannelKind::Unproved,
        };
        channels.push(MasslessChannel {
            vertex_partition: [
                (0..vertices).filter(|&v| side[v]).collect(),
                (0..vertices).filter(|&v| !side[v]).collect(),
            ],
            future_cut_coefficients: coefficients,
            kind,
        });
    }
    let certified = channels
        .iter()
        .all(|c| c.kind != MasslessChannelKind::Unproved);
    Ok(MasslessChannelEvidence {
        input_identity: input.identity().to_owned(),
        cut_slots: cuts.to_vec(),
        future_orientations,
        chemical_magnitudes: family
            .shells()
            .iter()
            .map(|s| s.chemical_potential.to_string())
            .collect(),
        shifted_physical_slots: uncut,
        physical_mass_squares: input
            .physical_masses()
            .iter()
            .map(Atom::to_canonical_string)
            .collect(),
        inverse_routing: family
            .inverse_routing()
            .iter()
            .map(|row| row.iter().map(Atom::to_canonical_string).collect())
            .collect(),
        routing_determinant: family.routing_determinant().to_string(),
        candidate_partitions,
        partition_budget,
        channels,
        certified,
        scope: "only if certified: all-uncut common eta with Re(eta)>0; regulated complete-amplitude common external-energy continuation; F>=eta*U for real eta and projective positive parameters; no eta=0 endpoint admission",
        endpoint_admitted: false,
        numerical_evaluation_admitted: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    fn input(text: &str) -> PreparedDensityInput {
        serde_json::from_str::<DensityInput>(text)
            .unwrap()
            .prepare()
            .unwrap()
    }

    fn evidence(input: &PreparedDensityInput, cuts: &[usize]) -> MasslessChannelEvidence {
        let shifted = (0..input.input().edges.len())
            .filter(|i| !cuts.contains(i))
            .collect::<Vec<_>>();
        massless_channel_evidence(input, cuts, &shifted, 8192, 8192).unwrap()
    }

    #[test]
    fn generic_incidence_certifies_all_e7_and_prism_channels() {
        for (text, expected) in [
            (
                include_str!("../../examples/finite_density/chain_of_three_parallel_pairs.json"),
                3,
            ),
            (
                include_str!("../../examples/finite_density/triangular_prism.json"),
                7,
            ),
        ] {
            let input = input(text);
            let mut count = 0;
            for cut in input.cut_decomposition(8192).unwrap() {
                let slots = cut["cut_slots"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_u64().unwrap() as usize)
                    .collect::<Vec<_>>();
                if slots.is_empty() {
                    continue;
                }
                let report = evidence(&input, &slots);
                assert!(report.certified);
                assert!(!report.endpoint_admitted && !report.numerical_evaluation_admitted);
                count += 1;
            }
            assert_eq!(count, expected);
        }
    }

    #[test]
    fn actual_sum_channel_is_reported_unproved_without_omitting_the_sector() {
        let input = input(include_str!(
            "../../examples/finite_density/five_vertex_eight_edge.json"
        ));
        let report = evidence(&input, &[0, 6]);
        assert!(!report.certified);
        assert!(report.channels.iter().any(|channel| {
            channel.vertex_partition == [vec![0, 2], vec![1, 3, 4]]
                && channel.future_cut_coefficients == vec![-1, -1]
                && channel.kind == MasslessChannelKind::Unproved
        }));
    }

    #[test]
    fn strict_deformation_and_partition_limits_are_not_silently_relaxed() {
        let input = input(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ));
        let shifted = vec![1, 2, 3, 5, 6];
        assert!(matches!(
            massless_channel_evidence(&input, &[0, 4], &shifted, 8192, 1),
            Err(Error::Limit(_))
        ));
        assert!(matches!(
            massless_channel_evidence(&input, &[0, 4], &[1, 2], 8192, 8192),
            Err(Error::Unsupported(_))
        ));
    }

    #[test]
    fn orientation_reversal_preserves_future_channels_and_positive_cut_mass_is_rejected() {
        let mut definition: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap();
        let before = evidence(&definition.prepare().unwrap(), &[1, 5, 7]);
        assert_eq!(before.future_orientations, vec![1, 1, -1]);
        definition.name = "arbitrary_name_is_not_a_contour_certificate".into();
        for slot in [2, 7] {
            let edge = &mut definition.edges[slot];
            edge.vertices.swap(0, 1);
            edge.routing = edge.routing.iter().map(|r| format!("-({r})")).collect();
            for charge in &mut edge.charges {
                *charge = -*charge;
            }
        }
        let after = evidence(&definition.prepare().unwrap(), &[1, 5, 7]);
        assert_eq!(before.channels, after.channels);
        assert_eq!(after.future_orientations, vec![1, 1, 1]);
        definition.edges[1].mass_squared = "1/4".into();
        let input = definition.prepare().unwrap();
        assert!(matches!(
            massless_channel_evidence(&input, &[1, 5, 7], &[0, 2, 3, 4, 6, 8], 8192, 8192),
            Err(Error::Unsupported(_))
        ));
    }
}
