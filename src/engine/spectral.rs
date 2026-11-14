//! Spectral export of an auxiliary-mass flow: besides the physical endpoint
//! values, the `eta -> infinity` expansion of every target and the target
//! values on an extra path ("arm") of the auxiliary-mass plane. These are the
//! inputs of a dispersive gluing (see [`crate::spectral::gluing`]): the target
//! is integrated against a spectral weight along a contour in `eta`.
use super::*;
use crate::ode::{BoundaryData, TaylorSegment};

/// `eta -> infinity` expansion of one target in `x = 1/eta` for one exponent
/// class: `sum_n coefficients[n] x^(base + slope * eps + n)`.
#[derive(Clone, Debug)]
pub struct SpectralExpansion {
    pub base: f64,
    pub slope: f64,
    pub coefficients: Vec<ComplexFloat>,
}

#[derive(Clone, Debug)]
pub struct SpectralSample {
    /// Physical endpoint values, identical to [`PreparedFlow::evaluate`].
    pub endpoint: Vec<ComplexFloat>,
    pub nodes: Vec<ComplexFloat>,
    /// `[node][target]`
    pub node_values: Vec<Vec<ComplexFloat>>,
    /// `[target]` -> exponent classes
    pub expansions: Vec<Vec<SpectralExpansion>>,
    /// Largest `|coefficient * constant|` of a logarithmic (log power >= 1)
    /// term of a single infinity column. The expansions keep only the power
    /// series, so these terms must cancel between the columns of a class; this
    /// per-column size does not show whether they do (see `log_terms_residual`).
    pub log_terms_max: f64,
    /// Largest dropped logarithmic coefficient of a target expansion (class sum
    /// over the columns, projected with the reduction weights), relative to the
    /// largest kept coefficient of the same target and class. The expansions
    /// are complete only when this is negligible at the working precision.
    pub log_terms_residual: f64,
    pub arm_segments: usize,
    pub arm_steps: usize,
    pub arm_rejected_steps: usize,
    pub arm_end_values: Vec<ComplexFloat>,
}

/// What `evaluate_spectral` adds to the sample pipeline.
pub(super) struct SpectralRequest<'a> {
    pub(super) arm_path: &'a [ComplexFloat],
    pub(super) nodes: &'a [ComplexFloat],
    pub(super) expansion_order: usize,
}

/// Spectral data formed before the infinity basis is released.
pub(super) struct SpectralPrepared<'a> {
    request: SpectralRequest<'a>,
    weights: Vec<TargetWeights>,
    expansions: Vec<Vec<SpectralExpansion>>,
    log_terms_max: f64,
    log_terms_residual: f64,
    arm_start: BoundaryData,
}

pub(super) struct SpectralArmValues {
    node_values: Vec<Vec<ComplexFloat>>,
    expansions: Vec<Vec<SpectralExpansion>>,
    log_terms_max: f64,
    log_terms_residual: f64,
    arm_segments: usize,
    arm_steps: usize,
    arm_rejected_steps: usize,
    arm_end_values: Vec<ComplexFloat>,
}

/// Reduction weight of one target: per master, numerator and denominator
/// coefficients in `eta` (ascending powers), `None` for a zero weight.
pub(super) struct TargetWeights {
    per_master: Vec<Option<(Vec<ComplexFloat>, Vec<ComplexFloat>)>>,
}

impl TargetWeights {
    fn evaluate_dot(
        &self,
        p: Precision,
        eta: &ComplexFloat,
        masters: &[ComplexFloat],
    ) -> ComplexFloat {
        let horner = |c: &[ComplexFloat]| {
            c.iter()
                .rev()
                .fold(p.zero(), |s, a| p.add(&p.mul(&s, eta), a))
        };
        let mut out = p.zero();
        for (w, m) in self.per_master.iter().zip(masters) {
            if let Some((num, den)) = w {
                out = p.add(&out, &p.mul(&p.div(&horner(num), &horner(den)), m));
            }
        }
        out
    }
}

fn spectral_weights(
    reduced: &ReducedSystem,
    variable: Symbol,
    p: Precision,
    parameters: &ahash::HashMap<Atom, ComplexFloat>,
) -> Result<Vec<TargetWeights>> {
    use symbolica::domains::integer::IntegerRing;
    use symbolica::domains::{integer::Z, rational::Q, rational_polynomial::RationalPolynomial};
    let mut out = Vec::new();
    for target in &reduced.targets {
        let mut per_master = Vec::new();
        for basis in &reduced.basis {
            let Some(w) = target.get(basis).filter(|w| !w.is_zero()) else {
                per_master.push(None);
                continue;
            };
            let fraction: RationalPolynomial<IntegerRing, u16> = w
                .try_to_rational_polynomial(&Q, &Z, None)
                .map_err(|e| Error::Unsupported(format!("spectral weight is not rational: {e}")))?;
            let num = crate::ode::polynomial_coefficients(
                &fraction.numerator.to_expression(),
                variable,
                p,
                parameters,
            )?;
            let den = crate::ode::polynomial_coefficients(
                &fraction.denominator.to_expression(),
                variable,
                p,
                parameters,
            )?;
            per_master.push(Some((num, den)));
        }
        out.push(TargetWeights { per_master });
    }
    Ok(out)
}

/// Series of `num(1/x) / den(1/x) = x^shift * sum_m q_m x^m`.
fn inverse_series(
    p: Precision,
    num: &[ComplexFloat],
    den: &[ComplexFloat],
    order: usize,
) -> Result<(i64, Vec<ComplexFloat>)> {
    let a = num.len() - 1;
    let b = den.len() - 1;
    if den[b] == p.zero() {
        return Err(Error::Numerical(
            "spectral weight: vanishing leading denominator coefficient".into(),
        ));
    }
    let nt: Vec<ComplexFloat> = (0..=order)
        .map(|k| if k <= a { num[a - k].clone() } else { p.zero() })
        .collect();
    let dt: Vec<ComplexFloat> = (0..=order)
        .map(|k| if k <= b { den[b - k].clone() } else { p.zero() })
        .collect();
    let mut q: Vec<ComplexFloat> = Vec::with_capacity(order + 1);
    for m in 0..=order {
        let mut s = nt[m].clone();
        for k in 1..=m.min(b) {
            s = p.sub(&s, &p.mul(&dt[k], &q[m - k]));
        }
        q.push(p.div(&s, &dt[0]));
    }
    Ok((b as i64 - a as i64, q))
}

/// Target expansions at `eta -> infinity`, grouped by exponent class
/// `(slope, base mod 1)`. Only the power series (log power 0) is kept; also
/// returns `log_terms_max` and `log_terms_residual` (see [`SpectralSample`]).
fn spectral_expansions(
    p: Precision,
    infinity: &crate::frobenius::FrobeniusBasis,
    constants: &[ComplexFloat],
    weights: &[TargetWeights],
    parameters: &ahash::HashMap<Atom, ComplexFloat>,
    eps_atom: Atom,
    order: usize,
) -> Result<(Vec<Vec<SpectralExpansion>>, f64, f64)> {
    use symbolica::domains::float::RealLike;
    // exponent classes: (slope, base mod 1) -> columns with integer base offsets
    let mut classes: Vec<(f64, f64, Vec<(usize, i64)>)> = Vec::new();
    let mut log_max = 0f64;
    let mut log_residual = 0f64;
    for (j, column) in infinity.columns.iter().enumerate() {
        if constants[j] == p.zero() {
            continue;
        }
        let mut at = parameters.clone();
        at.insert(eps_atom.clone(), p.i(0));
        let e0 = p.eval(&column.exponent, &at)?;
        at.insert(eps_atom.clone(), p.i(1));
        let e1 = p.eval(&column.exponent, &at)?;
        let base = e0.re.to_f64();
        let slope = (e1.re.to_f64()) - base;
        for coeff in &column.coefficients {
            for row in coeff.iter().skip(1) {
                for v in row {
                    let m = p.norm(&p.mul(v, &constants[j])).to_f64();
                    if m > log_max {
                        log_max = m;
                    }
                }
            }
        }
        let frac = base - base.floor();
        let ib = base.floor() as i64;
        match classes
            .iter_mut()
            .find(|(s, f, _)| (s - slope).abs() < 1e-9 && (f - frac).abs() < 1e-9)
        {
            Some((_, _, cols)) => cols.push((j, ib)),
            None => classes.push((slope, frac, vec![(j, ib)])),
        }
    }
    let n = infinity
        .columns
        .first()
        .map_or(0, |c| c.coefficients[0][0].len());
    let mut out: Vec<Vec<SpectralExpansion>> = vec![Vec::new(); weights.len()];
    // weight series per target and master
    let wseries: Vec<Vec<Option<(i64, Vec<ComplexFloat>)>>> = weights
        .iter()
        .map(|w| {
            w.per_master
                .iter()
                .map(|x| {
                    x.as_ref()
                        .map(|(num, den)| inverse_series(p, num, den, order))
                        .transpose()
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    for (slope, frac, cols) in &classes {
        let bmin = cols.iter().map(|c| c.1).min().unwrap();
        // class series per master: S_i(x) = sum_j c_j x^(b_j - bmin) sum_n a_{n,i} x^n
        let mut s = vec![vec![p.zero(); order + 1]; n];
        // The same class sums of the logarithmic rows (log power >= 1). The
        // expansion keeps only the power series, so these must cancel.
        let log_rows = cols
            .iter()
            .flat_map(|&(j, _)| infinity.columns[j].coefficients.iter().map(Vec::len))
            .max()
            .unwrap_or(1)
            .saturating_sub(1);
        let mut logs = vec![vec![vec![p.zero(); order + 1]; n]; log_rows];
        for &(j, b) in cols {
            let shift = (b - bmin) as usize;
            for (k, coeff) in infinity.columns[j].coefficients.iter().enumerate() {
                if k + shift > order {
                    break;
                }
                for i in 0..n {
                    s[i][k + shift] = p.add(&s[i][k + shift], &p.mul(&coeff[0][i], &constants[j]));
                }
                for (l, row) in coeff.iter().enumerate().skip(1) {
                    for i in 0..n {
                        logs[l - 1][i][k + shift] =
                            p.add(&logs[l - 1][i][k + shift], &p.mul(&row[i], &constants[j]));
                    }
                }
            }
        }
        for (t, ws) in wseries.iter().enumerate() {
            let emin = ws.iter().flatten().map(|(e, _)| *e).min().unwrap_or(0);
            let project = |series: &[Vec<ComplexFloat>]| {
                let mut acc = vec![p.zero(); order + 1];
                for (i, w) in ws.iter().enumerate() {
                    let Some((e, q)) = w else { continue };
                    let off = (e - emin) as usize;
                    for m in 0..=order {
                        if m + off > order {
                            break;
                        }
                        if q[m] == p.zero() {
                            continue;
                        }
                        for k in 0..=(order - m - off) {
                            acc[m + k + off] =
                                p.add(&acc[m + k + off], &p.mul(&q[m], &series[i][k]));
                        }
                    }
                }
                acc
            };
            let largest = |series: &[ComplexFloat]| {
                series
                    .iter()
                    .map(|c| p.norm(c))
                    .fold(p.real(0), |m, v| if v > m { v } else { m })
            };
            let acc = project(&s);
            let kept = largest(&acc);
            for log_series in &logs {
                let dropped = largest(&project(log_series));
                if dropped > p.real(0) {
                    let relative = if kept > p.real(0) {
                        (dropped / &kept).to_f64()
                    } else {
                        f64::INFINITY
                    };
                    log_residual = log_residual.max(relative);
                }
            }
            out[t].push(SpectralExpansion {
                base: bmin as f64 + frac + emin as f64,
                slope: *slope,
                coefficients: acc,
            });
        }
    }
    Ok((out, log_max, log_residual))
}

/// Master values at a node lying on a saved segment `[center, end]`.
fn spectral_node_masters(
    p: Precision,
    segments: &[TaylorSegment],
    node: &ComplexFloat,
) -> Result<Vec<ComplexFloat>> {
    // Same interval test as diffexp's saved-segment lookup:
    // t = (node - center)/(end - center) in [0, 1].
    let tolerance = p.tolerance((p.bits / 4).min(30));
    let mut best: Option<(usize, Float)> = None;
    for (i, s) in segments.iter().enumerate() {
        let span = p.sub(&s.end, &s.center);
        if span == p.zero() {
            continue;
        }
        let t = p.div(&p.sub(node, &s.center), &span);
        if t.im > tolerance
            || t.im < -tolerance.clone()
            || t.re < -tolerance.clone()
            || t.re > p.real(1) + tolerance.clone()
        {
            continue;
        }
        if best.as_ref().is_none_or(|(_, b)| t.re < *b) {
            best = Some((i, t.re.clone()));
        }
    }
    let (i, _) =
        best.ok_or_else(|| Error::Numerical("spectral node not on a saved arm segment".into()))?;
    let s = &segments[i];
    let local = s.coordinate.local_point(p, &s.center, node)?;
    Ok(s.evaluate_local(&local)?.0)
}

impl<'a> SpectralPrepared<'a> {
    /// Transport the infinity boundary vector along the arm path with saved
    /// segments, read the masters off at the nodes and project onto the targets.
    pub(super) fn transport(
        self,
        p: Precision,
        compiled: &crate::ode::CompiledSystem,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<SpectralArmValues> {
        let mut segments = Vec::new();
        let arm = crate::ode::transport_series(
            compiled,
            &self.arm_start,
            self.request.arm_path,
            options,
            context,
            Some(&mut segments),
        )?;
        let mut node_values = Vec::with_capacity(self.request.nodes.len());
        for node in self.request.nodes {
            let masters = spectral_node_masters(p, &segments, node)?;
            node_values.push(
                self.weights
                    .iter()
                    .map(|w| w.evaluate_dot(p, node, &masters))
                    .collect::<Vec<_>>(),
            );
        }
        Ok(SpectralArmValues {
            node_values,
            expansions: self.expansions,
            log_terms_max: self.log_terms_max,
            log_terms_residual: self.log_terms_residual,
            arm_segments: segments.len(),
            arm_steps: arm.diagnostics.steps,
            arm_rejected_steps: arm.diagnostics.rejected_steps,
            arm_end_values: arm.values,
        })
    }
}

impl PreparedFlow {
    /// Spectral data that needs the infinity basis (see `evaluate_with`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn spectral_prepare<'a>(
        &self,
        request: SpectralRequest<'a>,
        p: Precision,
        infinity: &crate::frobenius::FrobeniusBasis,
        constants: &[ComplexFloat],
        parameters: &ahash::HashMap<Atom, ComplexFloat>,
        start: &ComplexFloat,
        initial: &[ComplexFloat],
    ) -> Result<SpectralPrepared<'a>> {
        let weights = spectral_weights(&self.reduced, self.system.variable, p, parameters)?;
        let (expansions, log_terms_max, log_terms_residual) = spectral_expansions(
            p,
            infinity,
            constants,
            &weights,
            parameters,
            Atom::var(self.family.epsilon),
            request.expansion_order,
        )?;
        Ok(SpectralPrepared {
            request,
            weights,
            expansions,
            log_terms_max,
            log_terms_residual,
            arm_start: BoundaryData {
                point: start.clone(),
                values: initial.to_vec(),
            },
        })
    }

    /// [`PreparedFlow::evaluate`] plus the spectral data of every target:
    /// (i) its `eta -> infinity` expansion (`x = 1/eta`, grouped by exponent
    /// class, to `expansion_order`) and (ii) its values at `nodes`, obtained by
    /// transporting the infinity boundary vector along `arm_path` (waypoints
    /// after the transport start point) and reading the saved Taylor segments.
    /// Every node must lie on the arm path. The endpoint values are bitwise
    /// those of `evaluate` with the same arguments.
    pub fn evaluate_spectral(
        &self,
        eps: &Rational,
        options: &FlowOptions,
        boundary: &dyn BoundaryProvider,
        context: &RunContext,
        arm_path: &[ComplexFloat],
        nodes: &[ComplexFloat],
        expansion_order: usize,
    ) -> Result<SpectralSample> {
        let request = SpectralRequest {
            arm_path,
            nodes,
            expansion_order,
        };
        let (endpoint, arm) = self.evaluate_with(eps, options, boundary, context, Some(request))?;
        let arm = arm.ok_or_else(|| Error::Numerical("spectral arm missing".into()))?;
        Ok(SpectralSample {
            endpoint,
            nodes: nodes.to_vec(),
            node_values: arm.node_values,
            expansions: arm.expansions,
            log_terms_max: arm.log_terms_max,
            log_terms_residual: arm.log_terms_residual,
            arm_segments: arm.arm_segments,
            arm_steps: arm.arm_steps,
            arm_rejected_steps: arm.arm_rejected_steps,
            arm_end_values: arm.arm_end_values,
        })
    }
}
