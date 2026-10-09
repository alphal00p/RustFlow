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
    target_conditions: Vec<Atom>,
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
        let target_conditions = crate::physical_conditions::rational_denominator_conditions(
            &geometry
                .targets()
                .iter()
                .flat_map(|target| target.values().cloned())
                .collect::<Vec<_>>(),
            &std::collections::BTreeSet::from([family.epsilon]),
        )?;
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
            flow.reduced
                .nonzero_conditions
                .extend(target_conditions.clone());
            Some(flow)
        };
        Ok(Self {
            flow,
            loops: geometry.loops(),
            dimension: options.dimension,
            targets: geometry.targets().len(),
            target_conditions,
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
        crate::physical_conditions::validate_conditions_at(
            &self.target_conditions,
            symbol!("rustflow_occupied::epsilon"),
            &BTreeMap::from([(
                symbol!("rustflow_occupied::epsilon"),
                Atom::num(epsilon.clone()),
            )]),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    #[test]
    fn zero_vacuum_projection_retains_original_coefficient_poles() {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_one_loop_tadpole.json"
        ))
        .unwrap();
        input.targets.truncate(1);
        input.targets[0].numerator = "u1/(rustflow_occupied::epsilon-1/2)".into();
        let input = input.prepare().unwrap();
        let options = FlowOptions::default();
        let context = RunContext::default();
        let prepared = PreparedDensityVacuum::prepare(&input, &options, &context).unwrap();
        assert!(
            prepared
                .evaluate(&Rational::from((1, 2)), &options, &context)
                .is_err()
        );
        assert!(
            prepared
                .evaluate(&Rational::from((1, 3)), &options, &context)
                .unwrap()[0]
                .is_zero()
        );
    }
}
