//! Optional numerical propagation of incoming errors in an exact diagonal gauge.
//! MPFR disk norms and refinement are estimates, not directed certificates.
use crate::algebraic::{AlgebraicSystem, SquareRoot, normalized_root_entry};
use crate::diffexp::{EpsilonSolution, EpsilonSystem, amplification_from_integral};
use crate::kinematics::{KinematicPath, KinematicSystem};
use crate::{ComplexFloat, DifferentialSystem, Error, Precision, Result, RunContext};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use symbolica::prelude::*;

/// Mathematical hint only: g_i = factors[i] * sqrt(radicands[i]).
/// Both arrays contain exact rational functions of the physical coordinates.
/// The owner checks the original connection; no error bound or success flag
/// can be supplied here. The default physical transport does not use a hint.
#[derive(Clone, Debug)]
pub struct DiagonalGaugeEstimate {
    pub factors: Vec<Atom>,
    pub radicands: Vec<Atom>,
}

pub(super) struct PreparedGauge {
    system: AlgebraicSystem,
    factors: Vec<Atom>,
    radicands: Vec<Atom>,
    longest: usize,
}

fn exact(a: &Atom) -> Atom {
    a.replace(Atom::var(crate::family::imaginary_parameter()))
        .with(Atom::i())
        .together()
        .cancel()
}

fn expansion(raw: &DifferentialSystem, epsilon: Symbol, count: usize) -> Result<EpsilonSystem> {
    let system = EpsilonSystem::from_differential_system(raw, epsilon, count - 1)?;
    for i in 0..raw.matrix.len() {
        for j in 0..raw.matrix.len() {
            let reconstructed = system
                .matrices
                .iter()
                .enumerate()
                .fold(Atom::zero(), |sum, (k, matrix)| {
                    sum + &matrix[i][j] * Atom::var(epsilon).pow(k as i64)
                });
            if !exact(&(&raw.matrix[i][j] - reconstructed)).is_zero() {
                return Err(Error::Unsupported(
                    "diagonal gauge estimate requires a complete polynomial epsilon connection"
                        .into(),
                ));
            }
        }
    }
    Ok(system)
}

fn longest_path(nodes: usize, edges: &BTreeSet<(usize, usize)>) -> Result<usize> {
    let mut next = vec![Vec::new(); nodes];
    let mut indegree = vec![0usize; nodes];
    for &(to, from) in edges {
        next[from].push(to);
        indegree[to] += 1;
    }
    let mut queue = (0..nodes)
        .filter(|&i| indegree[i] == 0)
        .collect::<VecDeque<_>>();
    let mut distance = vec![0usize; nodes];
    let mut visited = 0;
    while let Some(from) = queue.pop_front() {
        visited += 1;
        for &to in &next[from] {
            distance[to] = distance[to].max(distance[from] + 1);
            indegree[to] -= 1;
            if indegree[to] == 0 {
                queue.push_back(to);
            }
        }
    }
    if visited != nodes {
        return Err(Error::Unsupported(
            "diagonal gauge residual has a cyclic common epsilon support".into(),
        ));
    }
    Ok(distance.into_iter().max().unwrap_or(0))
}

pub(super) fn prepare(
    hint: &DiagonalGaugeEstimate,
    original: &KinematicSystem,
    path: &KinematicPath,
    count: usize,
    context: &RunContext,
) -> Result<PreparedGauge> {
    context.cancellation.check()?;
    original.validate()?;
    let size = original.derivatives.values().next().unwrap().len();
    // A small optional estimator, not a general algebraic gauge search.
    if count == 0
        || size.checked_mul(count).is_none_or(|n| n > 32)
        || hint.factors.len() != size
        || hint.radicands.len() != size
        || hint.radicands.iter().filter(|r| !exact(r).is_one()).count() > 2
    {
        return Err(Error::Unsupported(
            "diagonal gauge estimate needs <=32 channels and <=2 nontrivial square roots".into(),
        ));
    }
    let allowed = original
        .derivatives
        .keys()
        .map(|&s| Atom::var(s))
        .chain([Atom::var(crate::family::imaginary_parameter())])
        .collect::<BTreeSet<_>>();
    for value in hint.factors.iter().chain(&hint.radicands) {
        if exact(value).is_zero() {
            return Err(Error::InvalidInput(
                "zero diagonal gauge factor/radicand".into(),
            ));
        }
        normalized_root_entry(value, &[], &allowed)?;
    }
    let mut edges = BTreeSet::new();
    // Check all partial connections and their union, not only one sampled path.
    for (&variable, matrix) in &original.derivatives {
        context.cancellation.check()?;
        let blocks = expansion(
            &DifferentialSystem {
                variable,
                matrix: matrix.clone(),
            },
            original.epsilon,
            count,
        )?;
        for i in 0..size {
            let dlog = hint.factors[i].derivative(variable) / &hint.factors[i]
                + hint.radicands[i].derivative(variable) / (&hint.radicands[i] * 2);
            if !exact(&(&blocks.matrices[0][i][i] - dlog)).is_zero() {
                return Err(Error::InvalidInput(
                    "diagonal gauge logarithmic derivative does not match the original connection"
                        .into(),
                ));
            }
        }
        for order in 0..count {
            for shift in 0..=order {
                for i in 0..size {
                    for j in 0..size {
                        if (shift != 0 || i != j) && !exact(&blocks.matrices[shift][i][j]).is_zero()
                        {
                            edges.insert((order * size + i, (order - shift) * size + j));
                        }
                    }
                }
            }
        }
    }
    let longest = longest_path(size * count, &edges)?;
    let substitutions = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect::<BTreeMap<_, _>>();
    let pull = |a: &Atom| exact(&crate::family::substitute(a, &substitutions));
    let factors = hint.factors.iter().map(pull).collect::<Vec<_>>();
    let radicands = hint.radicands.iter().map(pull).collect::<Vec<_>>();
    if factors.iter().chain(&radicands).any(Atom::is_zero) {
        return Err(Error::InvalidInput(
            "path lies in a diagonal gauge zero/pole".into(),
        ));
    }
    let mut roots = Vec::new();
    let mut g = Vec::new();
    let mut inverse = Vec::new();
    for (i, (factor, radicand)) in factors.iter().zip(&radicands).enumerate() {
        let root = if radicand.is_one() {
            Atom::one()
        } else {
            let symbol = symbol!(&format!("symbolica_amflow::diagonal_error_root_{i}"));
            roots.push(SquareRoot {
                symbol,
                radicand: radicand.clone(),
            });
            Atom::var(symbol)
        };
        g.push(factor * &root);
        inverse.push(root / (factor * radicand));
    }
    let allowed = roots
        .iter()
        .map(|r| Atom::var(r.symbol))
        .chain([
            Atom::var(path.parameter),
            Atom::var(crate::family::imaginary_parameter()),
        ])
        .collect();
    let reduce = |value: &Atom| -> Result<Atom> {
        Ok(normalized_root_entry(value, &roots, &allowed)?.expression(&roots))
    };
    for i in 0..size {
        if !reduce(&(&g[i] * &inverse[i] - Atom::one()))?.is_zero() {
            return Err(Error::InvalidInput(
                "diagonal gauge root inverse relation failed".into(),
            ));
        }
    }
    let raw = original.pullback(path)?;
    let mut transformed = expansion(&raw, original.epsilon, count)?;
    for k in 0..count {
        for i in 0..size {
            for j in 0..size {
                let dlog = if k == 0 && i == j {
                    factors[i].derivative(path.parameter) / &factors[i]
                        + radicands[i].derivative(path.parameter) / (&radicands[i] * 2)
                } else {
                    Atom::zero()
                };
                transformed.matrices[k][i][j] =
                    reduce(&(&inverse[i] * &transformed.matrices[k][i][j] * &g[j] - dlog))?;
            }
        }
    }
    // Also verify the complete pulled connection identity in the root quotient.
    for i in 0..size {
        for j in 0..size {
            let h = transformed
                .matrices
                .iter()
                .enumerate()
                .fold(Atom::zero(), |sum, (k, m)| {
                    sum + &m[i][j] * Atom::var(original.epsilon).pow(k as i64)
                });
            let dlog = if i == j {
                factors[i].derivative(path.parameter) / &factors[i]
                    + radicands[i].derivative(path.parameter) / (&radicands[i] * 2)
            } else {
                Atom::zero()
            };
            // Epsilon is an exact coefficient variable for this identity only.
            let mut full_allowed = allowed.clone();
            full_allowed.insert(Atom::var(original.epsilon));
            let residual = &g[i] * h - &raw.matrix[i][j] * &g[j] + dlog * &g[i];
            if !normalized_root_entry(&residual, &roots, &full_allowed)?
                .expression(&roots)
                .is_zero()
            {
                return Err(Error::InvalidInput(
                    "complete pulled diagonal gauge identity failed".into(),
                ));
            }
        }
    }
    let conditions = factors.iter().chain(&radicands).cloned().collect();
    let system = AlgebraicSystem {
        system: transformed,
        roots,
        nonzero_conditions: conditions,
    };
    system.validate()?;
    Ok(PreparedGauge {
        system,
        factors,
        radicands,
        longest,
    })
}

impl PreparedGauge {
    fn moduli(&self, p: Precision, point: &ComplexFloat) -> Result<Vec<Float>> {
        let values =
            ahash::HashMap::from_iter([(Atom::var(self.system.system.variable), point.clone())]);
        self.factors
            .iter()
            .zip(&self.radicands)
            .map(|(f, r)| {
                let value = p.norm(&p.eval(f, &values)?) * p.norm(&p.eval(r, &values)?).sqrt();
                if !value.is_finite() || value <= p.real(0) {
                    return Err(Error::Accuracy(
                        "nonpositive/nonfinite diagonal gauge modulus".into(),
                    ));
                }
                Ok(value)
            })
            .collect()
    }

    pub(super) fn increments(
        &self,
        solution: &EpsilonSolution,
        weights: &[Float],
        relative: &Float,
        credit: u32,
        subdivisions: i64,
        context: &RunContext,
    ) -> Result<(Vec<Vec<Float>>, String)> {
        let p = Precision {
            bits: solution.diagnostics.working_bits,
        };
        let first = solution.segments.first().ok_or_else(|| {
            Error::Unsupported("missing trajectory for diagonal gauge estimate".into())
        })?;
        if solution.segments.last().unwrap().end != solution.point {
            return Err(Error::InvalidInput(
                "diagonal gauge endpoint/trajectory mismatch".into(),
            ));
        }
        let source_moduli = self.moduli(p, &first.center)?;
        let size = source_moduli.len();
        let what = weights
            .iter()
            .enumerate()
            .map(|(i, w)| w.clone() / &source_moduli[i % size])
            .collect::<Vec<_>>();
        let compiled = self
            .system
            .prepare_with_context(context)?
            .compile_with_context(p, context)?;
        let mut integral = p.real(0);
        let mut endpoints = Vec::new();
        let mut calls = 0usize;
        for segment in &solution.segments {
            context.cancellation.check()?;
            let delta = p.sub(&segment.end, &segment.center);
            let mut start = segment.center.clone();
            for part in 1..=subdivisions {
                context.cancellation.check()?;
                let end = if part == subdivisions {
                    segment.end.clone()
                } else {
                    p.add(&segment.center, &p.scale(&delta, part, subdivisions))
                };
                integral += compiled.error_norm_integral_weighted(&start, &end, &what)?;
                calls += 1;
                start = end;
            }
            let amplification = amplification_from_integral(p, &integral, Some(self.longest + 1))?;
            let moduli = self.moduli(p, &segment.end)?;
            let increments = what
                .iter()
                .enumerate()
                .map(|(i, w)| moduli[i % size].clone() * w * &amplification * relative)
                .collect::<Vec<_>>();
            if increments.iter().any(|v| !v.is_finite() || *v <= p.real(0)) {
                return Err(Error::Accuracy(
                    "diagonal gauge error estimate exceeds numerical range".into(),
                ));
            }
            endpoints.push(increments);
        }
        let report = serde_json::json!({"method":"diagonal gauge COMMON-DAG NUMERICAL ESTIMATE; not directed/certified", "longest_path":self.longest,
            "incoming_credit":credit,"working_bits":p.bits,"subdivisions":subdivisions,"norm_calls":calls,
            "fixed_weights":what.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "source_g_moduli":source_moduli.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "relative_with_actual_credit_floor":relative.to_string(),"whole_integral":integral.to_string()});
        Ok((endpoints, report.to_string()))
    }
}
