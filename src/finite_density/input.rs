use super::InversePropagatorBasis;
use crate::algebra::{determinant, rref};
use crate::family::{LinearCombination, substitute};
use crate::{Error, Result};
use linnet::{
    half_edge::{HedgeGraph, NoData, builder::HedgeGraphBuilder},
    tree::{Forest, child_vec::ChildVecStore},
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NumeratorConvention {
    #[default]
    ShiftedEuclidean,
    UnshiftedEuclidean,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityEdge {
    /// Directed incidence, independent of the squared denominator.
    pub vertices: [usize; 2],
    pub routing: Vec<String>,
    pub mass_squared: String,
    /// Signed conserved species charges in this edge's orientation.
    pub charges: Vec<i32>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityTarget {
    /// Uniform inverse quadratic powers. Negative powers mean that quadratic,
    /// never the sign-dependent legacy scalar-product convention of the oracle.
    pub powers: Vec<i16>,
    pub numerator: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityInput {
    pub name: String,
    pub loops: usize,
    pub vertices: usize,
    pub edges: Vec<DensityEdge>,
    pub loop_charges: Vec<Vec<i32>>,
    pub chemical_potentials: Vec<String>,
    pub targets: Vec<DensityTarget>,
    #[serde(default)]
    pub numerator_convention: NumeratorConvention,
    pub laurent_orders: [i32; 2],
    pub digits: u32,
}

/// Algebraically admitted input, retaining original graph, chemical assignments,
/// independent masses, and exact target conversion. This is not a solved flow.
#[derive(Clone, Debug)]
pub struct PreparedDensityInput {
    input: DensityInput,
    basis: InversePropagatorBasis,
    targets: Vec<LinearCombination>,
    independent_masses: Vec<Symbol>,
    physical_masses: Vec<Atom>,
    loop_chemical_potentials: Vec<Atom>,
    routings: Vec<Vec<Atom>>,
    identity: String,
}

fn parse(expression: &str) -> Result<Atom> {
    Atom::parse(expression, "rustflow_density", Default::default())
        .map_err(|e| Error::InvalidInput(format!("Symbolica expression: {e}")))
}

fn connected(vertices: usize, edges: &[DensityEdge], removed: &[usize]) -> bool {
    let mut builder = HedgeGraphBuilder::<(), ()>::new();
    let nodes = (0..vertices)
        .map(|_| builder.add_node(()))
        .collect::<Vec<_>>();
    for (i, edge) in edges.iter().enumerate() {
        if !removed.contains(&i) {
            builder.add_edge(nodes[edge.vertices[0]], nodes[edge.vertices[1]], (), true);
        }
    }
    let graph: HedgeGraph<(), (), NoData, Forest<(), ChildVecStore<()>>> = builder.build();
    // Physical connectivity ignores charge/momentum orientation. The native
    // edge-filter traversal excludes isolated vertices, so count those too.
    let isolated = graph
        .iter_nodes()
        .map(|(_, mut crown, _)| usize::from(crown.next().is_none()))
        .sum::<usize>();
    graph.count_connected_components(&graph.full_filter()) + isolated == 1
}

impl DensityInput {
    pub fn prepare(&self) -> Result<PreparedDensityInput> {
        let l = self.loops;
        let species = self.chemical_potentials.len();
        if l == 0
            || self.vertices == 0
            || self.edges.is_empty()
            || self.targets.is_empty()
            || self.digits == 0
            || self.laurent_orders[0] > self.laurent_orders[1]
        {
            return Err(Error::InvalidInput(
                "empty family, target, precision or Laurent range".into(),
            ));
        }
        if self
            .edges
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_sub(self.vertices))
            != Some(l)
        {
            return Err(Error::InvalidInput(
                "graph cycle dimension differs from loop count".into(),
            ));
        }
        if self.loop_charges.len() != l || self.loop_charges.iter().any(|c| c.len() != species) {
            return Err(Error::InvalidInput("loop charge dimensions".into()));
        }
        if self.edges.iter().any(|e| {
            e.vertices.iter().any(|&v| v >= self.vertices)
                || e.routing.len() != l
                || e.charges.len() != species
        }) {
            return Err(Error::InvalidInput(
                "edge incidence, routing or charge dimensions".into(),
            ));
        }
        if !connected(self.vertices, &self.edges, &[]) {
            return Err(Error::InvalidInput(
                "finite-density vacuum graph is disconnected".into(),
            ));
        }
        if (0..self.edges.len()).any(|i| !connected(self.vertices, &self.edges, &[i])) {
            return Err(Error::Unsupported(
                "cutting admission requires a one-particle-irreducible graph".into(),
            ));
        }
        let routings = self
            .edges
            .iter()
            .map(|e| {
                e.routing
                    .iter()
                    .map(|x| {
                        let a = parse(x)?;
                        Rational::try_from(a.as_view()).map_err(|_| {
                            Error::InvalidInput("routing must be exact real rational".into())
                        })?;
                        Ok(a)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        if rref(routings.clone()).1.len() != l {
            return Err(Error::InvalidInput(
                "routing does not span the graph cycle space".into(),
            ));
        }
        for vertex in 0..self.vertices {
            for loop_index in 0..l {
                let conservation =
                    self.edges
                        .iter()
                        .zip(&routings)
                        .fold(Atom::new(), |a, (edge, routing)| {
                            a + Atom::num(
                                i32::from(edge.vertices[1] == vertex)
                                    - i32::from(edge.vertices[0] == vertex),
                            ) * &routing[loop_index]
                        });
                if !conservation.together().cancel().is_zero() {
                    return Err(Error::InvalidInput(format!(
                        "momentum routing violates incidence at vertex {vertex}"
                    )));
                }
            }
            for charge in 0..species {
                let mut incoming = 0_i64;
                let mut outgoing = 0_i64;
                for edge in &self.edges {
                    let c = i64::from(edge.charges[charge]);
                    if edge.vertices[1] == vertex {
                        incoming += c.max(0);
                        outgoing += (-c).max(0);
                    }
                    if edge.vertices[0] == vertex {
                        outgoing += c.max(0);
                        incoming += (-c).max(0);
                    }
                }
                if incoming != outgoing || incoming > 1 {
                    return Err(Error::Unsupported(format!(
                        "species {charge} does not form directed nonbranching fermion cycles at vertex {vertex}"
                    )));
                }
            }
        }
        for (edge, routing) in self.edges.iter().zip(&routings) {
            if edge.charges.iter().filter(|&&q| q != 0).count() > 1 {
                return Err(Error::Unsupported(
                    "multiply charged lines need an extended cutting admission".into(),
                ));
            }
            for species_index in 0..species {
                let routed = routing
                    .iter()
                    .zip(&self.loop_charges)
                    .fold(Atom::new(), |a, (r, c)| a + r * Atom::num(c[species_index]));
                if !(routed - Atom::num(edge.charges[species_index]))
                    .together()
                    .cancel()
                    .is_zero()
                {
                    return Err(Error::InvalidInput(
                        "edge chemical charge disagrees with routed loop charge".into(),
                    ));
                }
            }
        }
        let potentials = self
            .chemical_potentials
            .iter()
            .map(|s| parse(s))
            .collect::<Result<Vec<_>>>()?;
        let loop_chemical_potentials = self
            .loop_charges
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&potentials)
                    .fold(Atom::new(), |sum, (&c, mu)| sum + Atom::num(c) * mu)
            })
            .collect::<Vec<_>>();
        let mut coordinates = Vec::new();
        for i in 0..l {
            for j in i..l {
                coordinates.push(parse(&format!("g{}_{}", i + 1, j + 1))?);
            }
        }
        let scalar_count = coordinates.len();
        for i in 0..l {
            coordinates.push(parse(&format!("u{}", i + 1))?);
        }
        let independent_masses = (0..self.edges.len())
            .map(|i| {
                parse(&format!("rustflow_density::line_mass_squared_{i}")).map(|a| {
                    if let AtomView::Var(v) = a.as_view() {
                        v.get_symbol()
                    } else {
                        unreachable!()
                    }
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let physical_masses = self
            .edges
            .iter()
            .map(|e| parse(&e.mass_squared))
            .collect::<Result<Vec<_>>>()?;
        for parameter in physical_masses.iter().chain(&potentials) {
            if coordinates
                .iter()
                .chain(
                    independent_masses
                        .iter()
                        .map(|s| Atom::var(*s))
                        .collect::<Vec<_>>()
                        .iter(),
                )
                .any(|x| {
                    if let AtomView::Var(v) = x.as_view() {
                        !parameter.derivative(v.get_symbol()).is_zero()
                    } else {
                        false
                    }
                })
            {
                return Err(Error::InvalidInput(
                    "physical parameter depends on a momentum or reserved mass coordinate".into(),
                ));
            }
        }
        let propagators = routings
            .iter()
            .zip(&independent_masses)
            .map(|(routing, &mass)| {
                let mut result = Atom::var(mass);
                let mut index = 0;
                for i in 0..l {
                    for j in i..l {
                        result += &routing[i]
                            * &routing[j]
                            * Atom::num(if i == j { 1 } else { 2 })
                            * &coordinates[index];
                        index += 1;
                    }
                }
                result
            })
            .collect::<Vec<_>>();
        let labels = (0..propagators.len() + coordinates.len())
            .map(|i| parse(&format!("rustflow_density::rho_{i}")))
            .collect::<Result<Vec<_>>>()?;
        let basis = InversePropagatorBasis::new(coordinates.clone(), propagators, labels)?;
        for parameter in physical_masses.iter().chain(&potentials) {
            if basis.labels().iter().any(|x| {
                let AtomView::Var(v) = x.as_view() else {
                    unreachable!()
                };
                !parameter.derivative(v.get_symbol()).is_zero()
            }) {
                return Err(Error::InvalidInput(
                    "physical parameter uses a reserved inverse-propagator label".into(),
                ));
            }
        }
        let mut shifts = BTreeMap::new();
        if matches!(
            self.numerator_convention,
            NumeratorConvention::UnshiftedEuclidean
        ) {
            let imaginary = Atom::num(symbolica::domains::float::Complex::new(
                Rational::from(0),
                Rational::from(1),
            ));
            let mut index = 0;
            for i in 0..l {
                for j in i..l {
                    shifts.insert(
                        coordinates[index].clone(),
                        &coordinates[index]
                            - &imaginary
                                * (&loop_chemical_potentials[i] * &coordinates[scalar_count + j]
                                    + &loop_chemical_potentials[j]
                                        * &coordinates[scalar_count + i])
                            - &loop_chemical_potentials[i] * &loop_chemical_potentials[j],
                    );
                    index += 1;
                }
            }
            for i in 0..l {
                shifts.insert(
                    coordinates[scalar_count + i].clone(),
                    &coordinates[scalar_count + i] - &imaginary * &loop_chemical_potentials[i],
                );
            }
        }
        let targets = self
            .targets
            .iter()
            .map(|t| {
                let numerator = parse(&t.numerator)?;
                if independent_masses
                    .iter()
                    .any(|&m| !numerator.derivative(m).is_zero())
                {
                    return Err(Error::InvalidInput(
                        "input numerator uses a reserved independent mass symbol".into(),
                    ));
                }
                basis.convert(&substitute(&numerator, &shifts), &t.powers)
            })
            .collect::<Result<Vec<_>>>()?;
        let identity = blake3::hash(
            serde_json::to_string(self)
                .map_err(|e| Error::InvalidInput(e.to_string()))?
                .as_bytes(),
        )
        .to_hex()
        .to_string();
        Ok(PreparedDensityInput {
            input: self.clone(),
            basis,
            targets,
            independent_masses,
            physical_masses,
            loop_chemical_potentials,
            routings,
            identity,
        })
    }
}

impl PreparedDensityInput {
    pub fn input(&self) -> &DensityInput {
        &self.input
    }
    pub fn basis(&self) -> &InversePropagatorBasis {
        &self.basis
    }
    pub fn targets(&self) -> &[LinearCombination] {
        &self.targets
    }
    pub fn independent_masses(&self) -> &[Symbol] {
        &self.independent_masses
    }
    pub fn physical_masses(&self) -> &[Atom] {
        &self.physical_masses
    }
    pub fn loop_chemical_potentials(&self) -> &[Atom] {
        &self.loop_chemical_potentials
    }
    pub fn identity(&self) -> &str {
        &self.identity
    }

    /// All connected complements of chemically charged cut sets, including the
    /// vacuum. A cut set's completion to loop coordinates retains the exact
    /// determinant. Charge signs specify its occupied energy orientation.
    /// This constructs combinatorics, not a contour or a numerical boundary.
    pub fn cut_decomposition(&self, budget: usize) -> Result<Vec<Value>> {
        let charged = self
            .input
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, e)| e.charges.iter().any(|&c| c != 0).then_some(i))
            .collect::<Vec<_>>();
        let count = 1_usize
            .checked_shl(
                u32::try_from(charged.len())
                    .map_err(|_| Error::Limit("cut subset count".into()))?,
            )
            .ok_or_else(|| Error::Limit("cut subset count".into()))?;
        if count > budget {
            return Err(Error::Limit(format!(
                "cut enumeration requires {count} subsets; budget {budget}"
            )));
        }
        let l = self.input.loops;
        let mut result = Vec::new();
        for mask in 0..count {
            let cuts = charged
                .iter()
                .enumerate()
                .filter_map(|(j, &slot)| ((mask >> j) & 1 != 0).then_some(slot))
                .collect::<Vec<_>>();
            if cuts.len() > l || !connected(self.input.vertices, &self.input.edges, &cuts) {
                continue;
            }
            let mut matrix = cuts
                .iter()
                .map(|&i| self.routings[i].clone())
                .collect::<Vec<_>>();
            if rref(matrix.clone()).1.len() != cuts.len() {
                return Err(Error::InvalidInput(
                    "connected cut complement has dependent cut routings".into(),
                ));
            }
            for i in 0..l {
                let mut candidate = matrix.clone();
                candidate.push((0..l).map(|j| Atom::num(i32::from(i == j))).collect());
                if rref(candidate.clone()).1.len() > matrix.len() {
                    matrix = candidate;
                }
            }
            let determinant = determinant(matrix.clone());
            result.push(json!({"cut_slots": cuts, "cut_sign": if cuts.len() % 2 == 0 { 1 } else { -1 },
                "loop_routing": matrix.iter().map(|r| r.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "determinant": determinant.to_string(), "jacobian_convention": "abs(determinant)^(-D) before energy-shell integration"}));
        }
        Ok(result)
    }

    pub fn preparation_report(&self, cut_budget: usize) -> Result<Value> {
        Ok(
            json!({"schema_version": 2, "operation": "finite-density-prepare", "status": "algebraic_preparation_only",
            "input_identity": self.identity, "input": self.input,
            "basis": {"coordinates": self.basis.coordinates().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "slots": self.basis.slots().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "physical_slots": self.basis.physical_slots(),
                "independent_masses": self.independent_masses.iter().map(|&s| Atom::var(s).to_string()).collect::<Vec<_>>()},
            "targets": self.targets.iter().map(|target| target.iter().map(|(i,c)| json!({"powers":i.0,"coefficient":c.to_string()})).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "cut_decomposition": self.cut_decomposition(cut_budget)?,
            "numerical_evaluation": {
                "performed": false,
                "rust_interface": "finite_density::assembly::PreparedDensityFlow",
                "cli_interfaces": ["finite-density", "finite-density-sample"],
                "python_interface": "evaluate_finite_density",
                "supported_scope": "docs/finite-density.md",
                "requirement": "evaluation separately requires admission and closure of every contributing cut sector"
            }}),
        )
    }
}
