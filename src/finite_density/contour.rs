//! Sufficient common-contour certificates for fixed occupied shells.
use super::{PreparedDensityInput, geometry::OccupiedCutFamily};
use crate::{Error, Result};
use symbolica::prelude::*;

/// A common simultaneous Wick rotation exists if every uncut physical line's
/// rest mass exceeds the largest energy offset from compact coordinates.
/// This proves a contour domain, not UV convergence or reduction closure.
/// Analytic regularization of UV behavior precedes meromorphic continuation.
pub(crate) fn heavy_edge_domain(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
) -> Result<String> {
    let loops = input.input().loops;
    let compact = family.shells().len();
    if family.loops() != loops
        || family.physical_slots() != input.input().edges.len()
        || family.inverse_routing().len() != loops
        || family
            .inverse_routing()
            .iter()
            .any(|row| row.len() != loops)
    {
        return Err(Error::InvalidInput(
            "input and occupied contour geometry disagree".into(),
        ));
    }
    let mut evidence = Vec::new();
    for (slot, edge) in input.input().edges.iter().enumerate() {
        if family.shells().iter().any(|s| s.physical_slot == slot) {
            continue;
        }
        let routing = edge
            .routing
            .iter()
            .map(|text| {
                let atom = Atom::parse(text, "rustflow_density", Default::default())
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                Rational::try_from(atom.as_view()).map_err(|_| {
                    Error::Unsupported(
                        "contour certification requires rational real routing".into(),
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut bound = Rational::zero();
        for j in 0..compact {
            let coefficient = (0..loops)
                .map(|i| {
                    let inverse = Rational::try_from(family.inverse_routing()[i][j].as_view())
                        .map_err(|_| {
                            Error::Unsupported("nonrational inverse cut routing".into())
                        })?;
                    Ok(&routing[i] * &inverse)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Rational::zero(), |sum, c| sum + c);
            let absolute = if coefficient < Rational::zero() {
                -coefficient
            } else {
                coefficient
            };
            bound += absolute * &family.shells()[j].chemical_potential;
        }
        let mass = Rational::try_from(input.physical_masses()[slot].as_view()).map_err(|_| {
            Error::Unsupported(
                "contour certification requires assigned real rational mass squares".into(),
            )
        })?;
        let margin = &mass - &(&bound * &bound);
        if margin <= Rational::zero() {
            return Err(Error::Unsupported(format!(
                "uncut edge {slot} has no strict heavy-edge contour margin: m²={mass}, occupied energy bound={bound}; another common-contour certificate is required"
            )));
        }
        evidence.push(format!(
            "edge={slot},energy_bound={bound},mass_squared={mass},strict_margin={margin}"
        ));
    }
    Ok(format!(
        "regulated thermal contour; simultaneous virtual Wick rotation and common cut-energy continuation lambda in [0,1]; all uncut physical lines satisfy strict heavy-edge bound; fixed shells; Re(eta)>=0; UV analytic regularization before meromorphic dimension; {}",
        evidence.join(";")
    ))
}

/// Cut distributions are currently admitted at strictly positive shell mass.
/// The heavy-edge proof and the three-edge invariant proof concern only the
/// uncut subamplitude; neither one defines a massless shell-endpoint product.
pub(crate) fn positive_shells(family: &OccupiedCutFamily) -> Result<()> {
    for shell in family.shells() {
        let mass = Rational::try_from(shell.mass_squared.as_view()).map_err(|_| {
            Error::Unsupported(
                "occupied shell mass squares must be assigned real rational numbers".into(),
            )
        })?;
        if mass <= Rational::zero() {
            return Err(Error::Unsupported("occupied numerical flow currently requires strictly positive shell masses; a massless or tachyonic endpoint is not admitted by an uncut-line contour certificate".into()));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    #[test]
    fn four_loop_heavy_edge_certificate_checks_every_cut_and_strict_margin() {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap();
        for edge in &mut input.edges {
            edge.mass_squared = if edge.charges.iter().any(|&c| c != 0) {
                "1/4"
            } else {
                "16"
            }
            .into();
        }
        // Use the two charged edge momenta as independent loop coordinates.
        // The sufficient bound is routing-dependent: the original neutral
        // completion can leave an avoidable occupied energy in a charged line.
        let routings = [
            [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
            [1, -1, 0, 0], [0, 1, 0, 0], [-1, 1, 1, 0], [-1, 1, 0, 1],
        ];
        for (edge, routing) in input.edges.iter_mut().zip(routings) {
            edge.routing = routing.into_iter().map(|value| value.to_string()).collect();
        }
        input.loop_charges = vec![vec![1], vec![1], vec![0], vec![0]];
        let prepared = input.prepare().unwrap();
        for cuts in [vec![0], vec![4], vec![0, 4]] {
            let family = prepared
                .occupied_cut(&cuts, 16)
                .unwrap()
                .at_physical_masses();
            positive_shells(&family).unwrap();
            assert!(
                heavy_edge_domain(&prepared, &family)
                    .unwrap()
                    .contains("strict_margin=")
            );
        }
        // At equality this sufficient contour proof fails; the test does not
        // mislabel that failure as a demonstrated physical pinch.
        input.edges[3].mass_squared = "4".into();
        let prepared = input.prepare().unwrap();
        let family = prepared
            .occupied_cut(&[0, 4], 16)
            .unwrap()
            .at_physical_masses();
        assert!(matches!(
            heavy_edge_domain(&prepared, &family),
            Err(Error::Unsupported(_))
        ));
        input.edges[0].mass_squared = "0".into();
        let prepared = input.prepare().unwrap();
        let family = prepared
            .occupied_cut(&[0], 16)
            .unwrap()
            .at_physical_masses();
        assert!(positive_shells(&family).is_err());
    }
}
