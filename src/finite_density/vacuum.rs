//! The zero-cut contribution uses the ordinary native AMF owner.
use super::PreparedDensityInput;
use super::normalization::native_measure_to_euclidean;
use crate::family::LinearCombination;
use crate::recursive::{RecursiveBoundary, RecursiveTerminalPolicy};
use crate::{
    ComplexFloat, Error, FlowOptions, KinematicPoint, Precision, PreparedFlow, Result, RunContext,
    RustRedBackend,
};
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub struct PreparedDensityVacuum {
    flow: Option<PreparedFlow>,
    loops: usize,
    dimension: i64,
    targets: usize,
}

impl PreparedDensityVacuum {
    pub fn prepare(
        input: &PreparedDensityInput,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        let geometry = input.continued_cut(&[], 1)?.at_physical_masses();
        let family =
            geometry.region_family(symbol!("rustflow_occupied::epsilon"), options.dimension)?;
        let requested = geometry
            .targets()
            .iter()
            .flat_map(|t| t.keys().cloned())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let flow = if requested.is_empty() {
            None
        } else {
            let backend = RustRedBackend {
                bubble_subloops: false,
                ..Default::default()
            };
            let mut flow = PreparedFlow::new(
                &family,
                &requested,
                &KinematicPoint::default(),
                &backend,
                options,
                context,
            )?;
            let indexed = requested
                .iter()
                .cloned()
                .zip(flow.reduced.targets.clone())
                .collect::<BTreeMap<_, _>>();
            // Compose the original fixed polynomial coefficients with the
            // native reductions BEFORE rational endpoint projection. This
            // retains cancellation between singular reduction weights.
            flow.reduced.targets = geometry
                .targets()
                .iter()
                .map(|target| {
                    let mut out = LinearCombination::new();
                    for (integral, coefficient) in target {
                        for (basis, value) in &indexed[integral] {
                            *out.entry(basis.clone()).or_default() += coefficient * value;
                        }
                    }
                    for value in out.values_mut() {
                        *value = value.together().cancel();
                    }
                    out.retain(|_, c| !c.is_zero());
                    out
                })
                .collect();
            Some(flow)
        };
        Ok(Self {
            flow,
            loops: geometry.loops(),
            dimension: options.dimension,
            targets: geometry.targets().len(),
        })
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        options.validate()?;
        if options.dimension != self.dimension {
            return Err(Error::InvalidInput(
                "zero-cut evaluation dimension changed after preparation".into(),
            ));
        }
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("vacuum precision overflow".into()))?,
        )?;
        let Some(flow) = &self.flow else {
            return Ok(vec![p.zero(); self.targets]);
        };
        let backend = RustRedBackend {
            bubble_subloops: false,
            ..Default::default()
        };
        let boundary = RecursiveBoundary::new(&backend, options, context)
            .with_terminal_policy(RecursiveTerminalPolicy::TadpolesOnly);
        let native = flow.evaluate(epsilon, options, &boundary, context)?;
        let dimension = Atom::num(self.dimension) - Atom::num(2) * Atom::num(epsilon.clone());
        let measure = p.eval(
            &native_measure_to_euclidean(self.loops, 0, &dimension)?,
            &ahash::HashMap::default(),
        )?;
        Ok(native.iter().map(|v| p.mul(v, &measure)).collect())
    }
}
