//! Assembly at a common regulator sample, before Laurent reconstruction.
use super::flow::{OccupiedFlowEvaluation, PreparedOccupiedFlow};
use super::reduction::{WeightedClosureDiagnostics, WeightedClosureOptions};
use super::vacuum::PreparedDensityVacuum;
use super::{DensityInput, PreparedDensityInput};
use crate::{
    ComplexFloat, Error, FlowOptions, LaurentExpansion, Precision, Progress, Result, RunContext,
};
use symbolica::prelude::*;

trait OccupiedEvaluation {
    fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
    ) -> Result<OccupiedFlowEvaluation>;
    fn diagnostics(&self) -> &WeightedClosureDiagnostics;
}
impl<const N: usize> OccupiedEvaluation for PreparedOccupiedFlow<N> {
    fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
    ) -> Result<OccupiedFlowEvaluation> {
        self.evaluate_report(epsilon, options, context, start_scale)
    }
    fn diagnostics(&self) -> &WeightedClosureDiagnostics {
        self.closure_diagnostics()
    }
}

/// Unscaled Euclidean values, with each connected cut contribution retained.
pub struct DensityEvaluation {
    pub contributions: Vec<(Vec<usize>, Vec<ComplexFloat>)>,
    pub values: Vec<ComplexFloat>,
    pub occupied_reports: Vec<(Vec<usize>, OccupiedFlowEvaluation)>,
    pub empty_support: Vec<(Vec<usize>, String)>,
}

enum OccupiedSector {
    Flow(Box<dyn OccupiedEvaluation>),
    Compact(super::terminal::PreparedOccupiedTerminal),
    EmptySupport(String),
}

/// Complete vacuum plus occupied amplitude. Every contributing cut is admitted
/// and source-closed during preparation; a failed sector cannot be omitted.
pub struct PreparedDensityFlow {
    input: PreparedDensityInput,
    vacuum: PreparedDensityVacuum,
    occupied: Vec<(Vec<usize>, OccupiedSector)>,
}

fn empty_occupied_support(input: &PreparedDensityInput, cuts: &[usize]) -> Result<Option<String>> {
    let parse = |text: &str| -> Result<Rational> {
        let atom = Atom::parse(text, "rustflow_density", Default::default())
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        Rational::try_from(atom.as_view()).map_err(|_| {
            Error::Unsupported(
                "occupied support requires assigned real rational masses and potentials".into(),
            )
        })
    };
    let potentials = input
        .input()
        .chemical_potentials
        .iter()
        .map(|text| parse(text))
        .collect::<Result<Vec<_>>>()?;
    for &slot in cuts {
        let edge = &input.input().edges[slot];
        let mass = parse(&edge.mass_squared)?;
        let mu = edge
            .charges
            .iter()
            .zip(&potentials)
            .fold(Rational::zero(), |sum, (&charge, potential)| {
                sum + Rational::from(charge) * potential
            });
        if &mu * &mu < mass {
            return Ok(Some(format!(
                "edge {slot}: mu²={} < mass²={mass}; the shell and occupied interval have disjoint support with a strict gap, including all finite shell/endpoint derivatives",
                &mu * &mu
            )));
        }
    }
    Ok(None)
}

impl PreparedDensityFlow {
    pub fn prepare(
        input: &DensityInput,
        options: &FlowOptions,
        closure: WeightedClosureOptions,
        context: &RunContext,
    ) -> Result<Self> {
        super::flow::validate_options(options)?;
        let input = input.prepare()?;
        let vacuum = PreparedDensityVacuum::prepare(&input, options, context)?;
        let mut occupied = Vec::new();
        for certificate in input.cut_decomposition(65536)? {
            let cuts = certificate["cut_slots"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_u64().unwrap() as usize)
                .collect::<Vec<_>>();
            if cuts.is_empty() {
                continue;
            }
            if let Some(proof) = empty_occupied_support(&input, &cuts)? {
                occupied.push((cuts, OccupiedSector::EmptySupport(proof)));
                continue;
            }
            if let Some(terminal) =
                super::terminal::PreparedOccupiedTerminal::prepare(&input, &cuts, options)?
            {
                occupied.push((cuts, OccupiedSector::Compact(terminal)));
                continue;
            }
            let arity = input.basis().slots().len() + 2 * cuts.len();
            let mut limits = closure.clone();
            if let Some(root) = &closure.checkpoints {
                limits.checkpoints = Some(root.join(format!(
                        "cut-{}",
                        cuts.iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join("-")
                    )));
            }
            // Native guarded programs have const-generic storage. Preserve the
            // established exact capacities, then select the smallest compiled
            // bucket for larger or intermediate physical arities. Padding is
            // confined to zero-tail index storage, never physical factors.
            let flow: Box<dyn OccupiedEvaluation> = match arity {
                7 => Box::new(PreparedOccupiedFlow::<7>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                9 => Box::new(PreparedOccupiedFlow::<9>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                1..=12 => Box::new(PreparedOccupiedFlow::<12>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                13..=16 => Box::new(PreparedOccupiedFlow::<16>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                17..=20 => Box::new(PreparedOccupiedFlow::<20>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                21..=24 => Box::new(PreparedOccupiedFlow::<24>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                25..=32 => Box::new(PreparedOccupiedFlow::<32>::prepare(
                    &input, &cuts, options, limits, context,
                )?),
                _ => {
                    return Err(Error::Unsupported(format!(
                        "occupied physical arity {arity} exceeds the compiled native storage capacities [7, 9, 12, 16, 20, 24, 32]; use the const-generic Rust flow API to compile a larger capacity"
                    )));
                }
            };
            occupied.push((cuts, OccupiedSector::Flow(flow)));
        }
        Ok(Self {
            input,
            vacuum,
            occupied,
        })
    }

    pub fn input(&self) -> &PreparedDensityInput {
        &self.input
    }

    pub fn closure_diagnostics(&self) -> Vec<(&[usize], &WeightedClosureDiagnostics)> {
        self.occupied
            .iter()
            .filter_map(|(cuts, sector)| match sector {
                OccupiedSector::Flow(flow) => Some((cuts.as_slice(), flow.diagnostics())),
                OccupiedSector::Compact(_) | OccupiedSector::EmptySupport(_) => None,
            })
            .collect()
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
    ) -> Result<DensityEvaluation> {
        super::flow::validate_options(options)?;
        if epsilon.is_zero() || start_scale < 4 {
            return Err(Error::InvalidInput("density evaluation requires nonzero epsilon and occupied start scale at least four".into()));
        }
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("density precision overflow".into()))?,
        )?;
        let mut contributions = vec![(vec![], self.vacuum.evaluate(epsilon, options, context)?)];
        let mut occupied_reports = Vec::new();
        let mut empty_support = Vec::new();
        for (cuts, sector) in &self.occupied {
            context.emit(Progress::Stage {
                name: format!("evaluating occupied cut {cuts:?}"),
            })?;
            let values = match sector {
                OccupiedSector::Flow(flow) => {
                    let report = flow.evaluate(epsilon, options, context, start_scale)?;
                    let values = report.values.clone();
                    occupied_reports.push((cuts.clone(), report));
                    values
                }
                OccupiedSector::Compact(terminal) => {
                    let report = terminal.evaluate(epsilon, options, context)?;
                    let values = report.values.clone();
                    occupied_reports.push((cuts.clone(), report));
                    values
                }
                OccupiedSector::EmptySupport(proof) => {
                    empty_support.push((cuts.clone(), proof.clone()));
                    vec![p.zero(); self.input.input().targets.len()]
                }
            };
            contributions.push((cuts.clone(), values));
        }
        let mut values = vec![p.zero(); self.input.input().targets.len()];
        for (_, part) in &contributions {
            if part.len() != values.len() {
                return Err(Error::InvalidInput(
                    "density target count changed after preparation".into(),
                ));
            }
            for (sum, value) in values.iter_mut().zip(part) {
                *sum = p.add(sum, value);
            }
        }
        Ok(DensityEvaluation {
            contributions,
            values,
            occupied_reports,
            empty_support,
        })
    }

    /// The whole amplitude is assembled at each exact epsilon before fitting.
    /// The existing owner independently refines its grid, precision and order.
    /// Numerical sample values contain no oracle coefficients or constraints.
    pub fn solve(
        &self,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        self.solve_with_start_scale(options, context, 8)
    }

    pub fn solve_with_start_scale(
        &self,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
    ) -> Result<Vec<LaurentExpansion>> {
        self.solve_with_sampling_grid(options, context, start_scale, 1000)
    }

    /// Independently change the initial exact epsilon grid while keeping the
    /// auxiliary start scale and numerical options fixed. The shared Laurent
    /// owner still requires a further independent grid/precision refinement.
    pub fn solve_with_sampling_grid(
        &self,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
        epsilon_grid_denominator: i64,
    ) -> Result<Vec<LaurentExpansion>> {
        let input = self.input.input();
        let mut pole_bound = i32::try_from(input.loops)
            .ok()
            .and_then(|l| l.checked_mul(-2))
            .ok_or_else(|| Error::Limit("density Laurent pole bound overflow".into()))?;
        // The inherited 2L bound is conservative in this massive admitted
        // domain. Explicit regulator poles in original numerator coefficients
        // extend it; they must never become implicit fit constraints.
        let epsilon = symbol!("rustflow_occupied::epsilon");
        let geometry = self.input.continued_cut(&[], 1)?.at_physical_masses();
        let mut coefficient_floor = 0;
        let allowed = std::collections::BTreeSet::from([Atom::var(epsilon)]);
        for coefficient in geometry.targets().iter().flat_map(|target| target.values()) {
            coefficient_floor = coefficient_floor.min(crate::engine::projection_weight_valuation(
                coefficient,
                epsilon,
                &allowed,
            )?);
        }
        pole_bound = pole_bound
            .checked_add(coefficient_floor)
            .ok_or_else(|| Error::Limit("density numerator pole bound overflow".into()))?;
        let mut expansions = crate::engine::fit_samples_refined_leading_with_grid(
            input.targets.len(),
            input.laurent_orders[0].min(pole_bound),
            input.laurent_orders[1],
            options,
            epsilon_grid_denominator,
            |samples, refined| {
                samples
                    .iter()
                    .enumerate()
                    .map(|(index, epsilon)| {
                        context.emit(Progress::Sample {
                            index,
                            total: samples.len(),
                        })?;
                        self.evaluate(epsilon, refined, context, start_scale)
                            .map(|v| v.values)
                    })
                    .collect()
            },
        )?;
        for expansion in &mut expansions {
            expansion
                .coefficients
                .retain(|power, _| *power >= input.laurent_orders[0]);
            expansion
                .comparison_errors
                .retain(|power, _| *power >= input.laurent_orders[0]);
        }
        Ok(expansions)
    }
}
