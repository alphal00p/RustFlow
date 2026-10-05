//! Analytic, log-free initialization at a pure-epsilon regular singular origin.
use super::*;
use crate::fixed_series::{coefficients, fixed_series};

fn maximum(a: Float, b: Float) -> Float {
    if a > b { a } else { b }
}

/// Local series settings; satisfying these checks does not certify unknown input errors.
#[derive(Clone, Debug)]
pub struct AnalyticOriginOptions {
    /// Highest retained nonnegative power of the local coordinate.
    pub order: usize,
    /// Scaled decimal tolerance for tails, root charts and differential defects.
    pub check_digits: u32,
    /// Shared cancellation token, checked between exact operations and series work.
    pub cancellation: crate::CancellationToken,
}

impl Default for AnalyticOriginOptions {
    fn default() -> Self {
        Self {
            order: 80,
            check_digits: 20,
            cancellation: Default::default(),
        }
    }
}

/// Derived regular boundary and diagnostics for the selected analytic sector.
pub struct AnalyticOriginSeed {
    pub boundary: EpsilonBoundary,
    pub seeds: BTreeMap<Symbol, RootSeed>,
    pub last_terms: Vec<Float>,
    pub maximum_defect: Float,
}

impl AlgebraicSystem {
    /// Derive the analytic, log-free sector from exact encoded boundary constants.
    /// Boundary rows begin at epsilon^0; the returned boundary has leading=0.
    /// This deliberately excludes a root branch point, a non-pure-epsilon matrix,
    /// higher-order poles and nonzero residue action. No uncertainty certificate
    /// for an unconstrained boundary-error box near the singular origin is implied.
    /// Compatibility is checked formally with registered root symbols independent;
    /// a cancellation valid only on one selected origin sheet is conservatively
    /// rejected. Exact encoded finite decimals do not acquire extra physical accuracy.
    pub fn analytic_origin(
        &self,
        exact_boundary: &[Vec<Atom>],
        constants: &ahash::HashMap<Atom, C>,
        seeds: &BTreeMap<Symbol, RootSeed>,
        p: Precision,
        offset: &C,
        options: &AnalyticOriginOptions,
    ) -> Result<AnalyticOriginSeed> {
        let system = self;
        let order = options.order;
        let check_digits = options.check_digits;
        let cancellation = &options.cancellation;
        cancellation.check()?;
        system.validate()?;
        let n = system.system.matrices[0].len();
        let count = exact_boundary.len();
        if system.system.matrices.len() < 2
            || system
                .system
                .matrices
                .iter()
                .enumerate()
                .any(|(k, m)| k != 1 && m.iter().flatten().any(|a| !a.is_zero()))
        {
            return Err(Error::Unsupported(
                "analytic origin requires a pure-epsilon connection".into(),
            ));
        }
        if count != system.system.matrices.len()
            || count > 128
            || exact_boundary.iter().any(|r| r.len() != n)
            || !(2..=512).contains(&order)
            || check_digits == 0
            || !p.finite(offset)
            || *offset == p.zero()
        {
            return Err(Error::InvalidInput(
                "invalid analytic-origin boundary, offset or order".into(),
            ));
        }
        if constants.values().any(|value| !p.finite(value)) {
            return Err(Error::InvalidInput(
                "nonfinite analytic-origin constants".into(),
            ));
        }
        for value in exact_boundary.iter().flatten() {
            cancellation.check()?;
            let mut symbols = BTreeSet::new();
            crate::family::scalar_symbols(value.as_view(), &mut symbols)?;
            if symbols.contains(&Atom::var(system.system.variable))
                || system
                    .roots
                    .iter()
                    .any(|root| symbols.contains(&Atom::var(root.symbol)))
            {
                return Err(Error::Unsupported("analytic boundary must contain constant values, independent of the coordinate and registered roots".into()));
            }
        }
        let x = Atom::var(system.system.variable);
        let at_origin = BTreeMap::from([(x.clone(), Atom::new())]);
        for root in &system.roots {
            cancellation.check()?;
            let value = substitute(&root.radicand, &at_origin).together().cancel();
            if value.is_zero() {
                return Err(Error::Unsupported(
                    "analytic origin cannot be a root branch point".into(),
                ));
            }
        }
        // G=x*A is a different regular coefficient function. Its canceled origin
        // poles are removed only for constructing this local solution, while the
        // original system and its complete domains are used at the nonzero offset.
        let g = system.system.matrices[1]
            .iter()
            .map(|row| {
                cancellation.check()?;
                Ok(row
                    .iter()
                    .map(|a| (&x * a).together().cancel())
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?;
        let regular = AlgebraicSystem {
            system: EpsilonSystem {
                variable: system.system.variable,
                matrices: vec![g.clone()],
            },
            roots: system.roots.clone(),
            nonzero_conditions: Vec::new(),
        };
        cancellation.check()?;
        let compiled = regular.compile(p)?;
        compiled.check_domain(&p.zero())?;
        cancellation.check()?;
        let original = system.compile(p)?;
        cancellation.check()?;
        original.check_domain(offset)?;
        let residue = g
            .iter()
            .map(|row| {
                cancellation.check()?;
                Ok(row
                    .iter()
                    .map(|a| substitute(a, &at_origin).together().cancel())
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?;
        for (k, boundary) in exact_boundary.iter().enumerate() {
            for (i, row) in residue.iter().enumerate() {
                cancellation.check()?;
                let sum = row
                    .iter()
                    .zip(boundary)
                    .fold(Atom::new(), |a, (r, y)| a + r * y)
                    .together()
                    .cancel();
                if !sum.is_zero() {
                    return Err(Error::InvalidInput(format!(
                        "exact encoded residue compatibility failed at epsilon{k}, row{i}: {sum}"
                    )));
                }
            }
        }
        let radius = compiled
            .poles
            .iter()
            .chain(original.poles.iter().filter(|v| **v != p.zero()))
            .map(|v| p.norm(v))
            .min_by(|a, b| a.partial_cmp(b).unwrap());
        if radius.as_ref().is_some_and(|r| p.norm(offset) * 2 >= *r) {
            return Err(Error::InvalidInput(
                "analytic offset lies outside half the nearest regularized pole radius".into(),
            ));
        }
        let run = AlgebraicRun {
            compiled: &compiled,
            seeds,
        };
        let initial = BoundaryData {
            point: p.zero(),
            values: vec![p.zero(); n],
        };
        cancellation.check()?;
        let state = run.initial_state(&initial)?;
        // Reuse the existing rational square-root ODE Taylor facility.
        cancellation.check()?;
        let root_coefficients = if let Some(root_system) = &compiled.root_system {
            root_system.taylor(
                &p.zero(),
                &state.roots.iter().map(|r| r.1.clone()).collect::<Vec<_>>(),
                order,
            )?
        } else {
            vec![Vec::new(); order + 1]
        };
        let root_series = (0..compiled.roots.len())
            .map(|i| {
                cancellation.check()?;
                fixed_series(
                    p,
                    system.system.variable,
                    &root_coefficients
                        .iter()
                        .map(|r| r[i].clone())
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let kernels = compiled
            .terms
            .iter()
            .map(|term| {
                cancellation.check()?;
                let mut kernel = fixed_series(
                    p,
                    system.system.variable,
                    &term.coefficient.series(p, &p.zero(), order)?,
                )?;
                for &root in &term.roots {
                    kernel = &kernel * &root_series[root];
                }
                coefficients(&kernel, order + 1)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut series = vec![vec![p.zero(); count * n]; order + 1];
        series[0] = exact_boundary
            .iter()
            .flatten()
            .map(|a| {
                cancellation.check()?;
                p.eval(a, constants)
            })
            .collect::<Result<Vec<_>>>()?;
        if series[0].iter().any(|v| !p.finite(v)) {
            return Err(Error::InvalidInput("nonfinite analytic boundary".into()));
        }
        let zero = p.zero();
        // The residue term involves epsilon k-1 at the same x-order, so k must be
        // solved upward inside each n. This selects the analytic sector explicitly.
        for power in 1..=order {
            cancellation.check()?;
            for k in 1..count {
                cancellation.check()?;
                for (term, kernel) in compiled.terms.iter().zip(&kernels) {
                    cancellation.check()?;
                    for (l, coefficient) in kernel.iter().enumerate().take(power + 1) {
                        let previous = &series[power - l][(k - 1) * n + term.column];
                        if *previous == zero || *coefficient == zero {
                            continue;
                        }
                        let contribution = p.mul(coefficient, previous);
                        let output = &mut series[power][k * n + term.row];
                        *output = p.add(output, &contribution);
                    }
                }
                for value in &mut series[power][k * n..(k + 1) * n] {
                    *value = p.scale(value, 1, power as i64);
                }
            }
        }
        if series.iter().flatten().any(|v| !p.finite(v)) {
            return Err(Error::Numerical("nonfinite analytic series".into()));
        }
        let chart = RootChart {
            coordinate: TaylorCoordinate::Identity,
            local_system: None,
            residual: Some(AlgebraicResidualChart::roots_only(
                &compiled,
                &p.zero(),
                &root_coefficients,
            )?),
            center: p.zero(),
            coefficients: root_coefficients,
        };
        let tolerance = p.tolerance(check_digits);
        let branches = run
            .accepted_state(&chart, offset, &tolerance)?
            .ok_or_else(|| {
                Error::Accuracy(
                    "analytic-origin root chart failed its existing acceptance checks".into(),
                )
            })?;
        cancellation.check()?;
        let (values, tails) = evaluate_taylor(p, &series, offset);
        for (value, tail) in values.iter().zip(&tails) {
            let scale = maximum(p.norm(value), p.real(1));
            if *tail > tolerance.clone() * scale {
                return Err(Error::Accuracy(
                    "analytic-origin Taylor tail did not meet tolerance".into(),
                ));
            }
        }
        let mut maximum_defect = p.real(0);
        for point in [offset.clone(), p.scale(offset, 1, 2)] {
            cancellation.check()?;
            let (local, _) = evaluate_taylor(p, &series, &point);
            let derivative = (0..count * n)
                .map(|i| {
                    let mut v = p.zero();
                    for (power, coefficient) in series.iter().enumerate().skip(1).rev() {
                        v = p.add(
                            &p.mul(&v, &point),
                            &p.scale(&coefficient[i], power as i64, 1),
                        );
                    }
                    v
                })
                .collect::<Vec<_>>();
            for k in 0..count {
                cancellation.check()?;
                let rhs = if k == 0 {
                    vec![p.zero(); n]
                } else {
                    run.rhs(&point, &local[(k - 1) * n..k * n], &chart)?
                };
                for (i, value) in rhs.iter().enumerate() {
                    let lhs = p.mul(&point, &derivative[k * n + i]);
                    let defect = p.norm(&p.sub(&lhs, value));
                    if !defect.is_finite() {
                        return Err(Error::Numerical("nonfinite analytic defect".into()));
                    }
                    maximum_defect = maximum(maximum_defect, defect.clone());
                    let scale = maximum(maximum(p.norm(&lhs), p.norm(value)), p.real(1));
                    if defect > tolerance.clone() * scale {
                        return Err(Error::Accuracy(
                            "analytic-origin differential defect failed".into(),
                        ));
                    }
                }
            }
        }
        // Original source domains are mandatory away from the singular origin.
        cancellation.check()?;
        Ok(AnalyticOriginSeed {
            boundary: EpsilonBoundary {
                point: offset.clone(),
                leading: 0,
                coefficients: values.chunks(n).map(<[C]>::to_vec).collect(),
            },
            seeds: branches
                .roots
                .iter()
                .map(|(s, v)| (*s, RootSeed::Value(v.clone())))
                .collect(),
            last_terms: tails,
            maximum_defect,
        })
    }
}
