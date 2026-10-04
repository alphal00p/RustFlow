//! Direct epsilon-coefficient transport of supplied differential equations.
//!
//! This independently implemented series recurrence shares the AMF continuation
//! engine. For A = sum epsilon^j A_j and Y = epsilon^leading sum epsilon^k Y_k,
//! it solves Y_k' = sum_{j=0}^k A_j Y_{k-j}, without sampled-epsilon fitting or
//! construction of a dense augmented matrix. Matrix coefficients must be
//! regular in epsilon at zero. A common Laurent prefactor of the solution is
//! allowed and is carried explicitly.
use crate::ode::{
    CompiledSystem, PolynomialRow, SeriesSystem, TaylorSegment, compile_rows, evaluate_taylor,
    shift_polynomial, transport_series,
};
use crate::{
    BoundaryData, ComplexFloat as C, DifferentialSystem, Error, FlowDiagnostics, FlowOptions,
    Precision, Result, RunContext,
};
use symbolica::prelude::*;

/// Exact coefficient matrices A_0, ..., A_order. Omitted higher orders do not
/// affect solution coefficients through `order` when the boundary starts at
/// the declared common leading power.
#[derive(Clone, Debug)]
pub struct EpsilonSystem {
    pub variable: Symbol,
    pub matrices: Vec<Vec<Vec<Atom>>>,
}

/// Boundary coefficients indexed first by epsilon offset, then by component.
/// Components may have different leading powers by supplying exact zeros.
#[derive(Clone, Debug)]
pub struct EpsilonBoundary {
    pub point: C,
    pub leading: i32,
    pub coefficients: Vec<Vec<C>>,
}

#[derive(Clone, Debug)]
pub struct EpsilonSolution {
    pub point: C,
    pub leading: i32,
    pub coefficients: Vec<Vec<C>>,
    pub diagnostics: FlowDiagnostics,
    /// Present only when requested; epsilon offsets precede component indices
    /// in each segment's flattened coefficient vectors.
    pub segments: Vec<TaylorSegment>,
    /// Set by independent precision/order recomputation, not working precision.
    pub verified_digits: Option<u32>,
    pub comparison_errors: Vec<Vec<Float>>,
    /// Independently compared accepted endpoints, suitable for a physical
    /// boundary bank after mapping coordinates and accounting for input errors.
    pub checkpoints: Vec<VerifiedCheckpoint>,
}

#[derive(Clone, Debug)]
pub struct VerifiedCheckpoint {
    pub segment: usize,
    pub point: C,
    pub coefficients: Vec<Vec<C>>,
    pub comparison_errors: Vec<Vec<Float>>,
}

#[derive(Clone, Debug)]
pub struct CompiledEpsilonSystem {
    rows: CompiledSystem,
    count: usize,
}

impl EpsilonSystem {
    pub fn validate(&self) -> Result<()> {
        let n = self.matrices.first().map_or(0, Vec::len);
        if n == 0
            || self
                .matrices
                .iter()
                .any(|m| m.len() != n || m.iter().any(|r| r.len() != n))
        {
            return Err(Error::InvalidInput(
                "epsilon matrices must be nonempty and square with equal dimensions".into(),
            ));
        }
        if self.matrices.len() > 1024 {
            return Err(Error::Limit("at most 1024 epsilon orders supported".into()));
        }
        Ok(())
    }

    /// Expand a rational epsilon dependence exactly. Negative epsilon powers
    /// in the matrix require a basis change and are rejected explicitly.
    pub fn from_differential_system(
        system: &DifferentialSystem,
        epsilon: Symbol,
        order: usize,
    ) -> Result<Self> {
        system.validate()?;
        if epsilon == system.variable {
            return Err(Error::InvalidInput(
                "epsilon and path parameter must differ".into(),
            ));
        }
        if order >= 1024 {
            return Err(Error::Limit("at most 1024 epsilon orders supported".into()));
        }
        let n = system.matrix.len();
        let mut matrices = vec![vec![vec![Atom::new(); n]; n]; order + 1];
        for (i, row) in system.matrix.iter().enumerate() {
            for (j, a) in row.iter().enumerate() {
                // In particular, do not accept floating coefficients or
                // fractional epsilon powers as an ordinary Taylor series.
                let a = crate::family::encode_complex(a);
                let _: RationalPolynomial<IntegerRing, u16> =
                    a.try_to_rational_polynomial(&Q, &Z, None).map_err(|e| {
                        Error::Unsupported(format!("rational epsilon system required: {e}"))
                    })?;
                let expansion = a
                    .series(epsilon, 0, order as i64 + 1)
                    .map_err(|e| Error::Unsupported(e.to_string()))?;
                if !expansion.is_zero() && expansion.get_trailing_exponent() < 0 {
                    return Err(Error::Unsupported(
                        "matrix has an epsilon pole; supply an epsilon-regular basis".into(),
                    ));
                }
                for (k, matrix) in matrices.iter_mut().enumerate() {
                    matrix[i][j] = expansion
                        .coefficient(Rational::from(k as i64))
                        .unwrap_or_default()
                        .together()
                        .cancel();
                }
            }
        }
        Ok(Self {
            variable: system.variable,
            matrices,
        })
    }

    pub fn compile(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
    ) -> Result<CompiledEpsilonSystem> {
        self.validate()?;
        // One denominator LCM per physical row across all epsilon orders.
        // Reuse these finite polynomials at every coefficient order; there is
        // no repeated (n*orders)^2 matrix allocation or pole finding.
        let n = self.matrices[0].len();
        let rows = (0..n)
            .map(|i| {
                self.matrices
                    .iter()
                    .flat_map(|matrix| matrix[i].iter().cloned())
                    .collect()
            })
            .collect::<Vec<_>>();
        Ok(CompiledEpsilonSystem {
            rows: compile_rows(self.variable, &rows, p, values)?,
            count: self.matrices.len(),
        })
    }
}

impl CompiledEpsilonSystem {
    pub fn poles(&self) -> &[C] {
        &self.rows.poles
    }

    pub fn plan_path(&self, start: &C, end: &C, side: i64) -> Result<Vec<C>> {
        crate::ode::plan_path(self.rows.p, self.poles(), start, end, side)
    }

    pub fn transport(
        &self,
        boundary: &EpsilonBoundary,
        waypoints: &[C],
        options: &FlowOptions,
        context: &RunContext,
        save_segments: bool,
    ) -> Result<EpsilonSolution> {
        let n = self.rows.dimension();
        if boundary.coefficients.len() != self.count
            || boundary.coefficients.iter().any(|r| r.len() != n)
        {
            return Err(Error::InvalidInput("epsilon boundary dimensions".into()));
        }
        boundary
            .leading
            .checked_add(self.count as i32 - 1)
            .ok_or_else(|| Error::InvalidInput("epsilon power range overflows".into()))?;
        let flat = BoundaryData {
            point: boundary.point.clone(),
            values: boundary.coefficients.iter().flatten().cloned().collect(),
        };
        let mut segments = Vec::new();
        let result = transport_series(
            self,
            &flat,
            waypoints,
            options,
            context,
            save_segments.then_some(&mut segments),
        )?;
        Ok(EpsilonSolution {
            point: result.point,
            leading: boundary.leading,
            coefficients: result.values.chunks(n).map(<[C]>::to_vec).collect(),
            diagnostics: result.diagnostics,
            segments,
            verified_digits: None,
            comparison_errors: Vec::new(),
            checkpoints: Vec::new(),
        })
    }

    /// A norm-based estimate of amplification of boundary uncertainty along a
    /// regular segment. Numerators are bounded above on each disk and common
    /// denominators below by the triangle inequality; subdivide if that lower
    /// bound is inconclusive. MPFR rounding and initial error estimates prevent
    /// interpreting this as an interval-arithmetic proof.
    pub fn error_amplification(&self, start: &C, end: &C) -> Result<Float> {
        let p = self.rows.p;
        let mut intervals = vec![(start.clone(), end.clone(), 0usize)];
        let mut exponent = p.real(0);
        while let Some((start, end, depth)) = intervals.pop() {
            let radius = p.norm(&p.sub(&end, &start));
            if radius == p.real(0) {
                continue;
            }
            let mut norm = p.real(0);
            let mut valid = true;
            for row in &self.rows.polynomial_rows {
                let denominator =
                    shift_polynomial(p, &row.denominator, &start, row.denominator.len() - 1);
                let mut lower = p.norm(&denominator[0]);
                let mut power = radius.clone();
                for a in denominator.iter().skip(1) {
                    lower -= p.norm(a) * &power;
                    power *= &radius;
                }
                if lower <= p.real(0) {
                    valid = false;
                    break;
                }
                let mut upper = p.real(0);
                for (_, coefficients) in &row.entries {
                    let shifted = shift_polynomial(p, coefficients, &start, coefficients.len() - 1);
                    let mut power = p.real(1);
                    for a in shifted {
                        upper += p.norm(&a) * &power;
                        power *= &radius;
                    }
                }
                let row_norm = upper / lower;
                if row_norm > norm {
                    norm = row_norm;
                }
            }
            if valid {
                exponent += radius * norm;
            } else {
                if depth >= 32 {
                    return Err(Error::Accuracy(
                        "cannot bound cached boundary error along this path".into(),
                    ));
                }
                let middle = p.scale(&p.add(&start, &end), 1, 2);
                if middle == start || middle == end {
                    return Err(Error::Accuracy(
                        "error-bound subdivision lost to rounding".into(),
                    ));
                }
                intervals.push((middle.clone(), end, depth + 1));
                intervals.push((start, middle, depth + 1));
            }
        }
        let amplification = p.exp(&C::new(exponent, p.real(0))).re;
        if !amplification.is_finite() {
            return Err(Error::Accuracy(
                "boundary-error amplification exceeds numerical range".into(),
            ));
        }
        Ok(amplification)
    }
}

impl SeriesSystem for CompiledEpsilonSystem {
    fn precision(&self) -> Precision {
        self.rows.p
    }
    fn dimension(&self) -> usize {
        self.rows.dimension() * self.count
    }
    fn poles(&self) -> &[C] {
        self.poles()
    }

    fn taylor(&self, center: &C, values: &[C], order: usize) -> Result<Vec<Vec<C>>> {
        let p = self.rows.p;
        let n = self.rows.dimension();
        if values.len() != self.dimension() {
            return Err(Error::InvalidInput("epsilon boundary dimensions".into()));
        }
        let rows = self
            .rows
            .polynomial_rows
            .iter()
            .map(|row| {
                let mut denominator = shift_polynomial(
                    p,
                    &row.denominator,
                    center,
                    (row.denominator.len() - 1).min(order.saturating_sub(1)),
                );
                if denominator[0] == p.zero() {
                    return Err(Error::Numerical("expansion center is a pole".into()));
                }
                let divisor = denominator[0].clone();
                for a in &mut denominator {
                    *a = p.div(a, &divisor);
                }
                let entries = row
                    .entries
                    .iter()
                    .map(|(j, coefficients)| {
                        let mut coefficients = shift_polynomial(
                            p,
                            coefficients,
                            center,
                            (coefficients.len() - 1).min(order.saturating_sub(1)),
                        );
                        for a in &mut coefficients {
                            *a = p.div(a, &divisor);
                        }
                        (*j, coefficients)
                    })
                    .collect();
                Ok(PolynomialRow {
                    denominator,
                    entries,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let zero = p.zero();
        let mut y = vec![vec![zero.clone(); self.dimension()]; order + 1];
        let mut derivative = vec![vec![zero.clone(); self.dimension()]; order];
        y[0] = values.to_vec();
        for k in 0..order {
            for e in 0..self.count {
                for (i, row) in rows.iter().enumerate() {
                    let index = e * n + i;
                    let mut sum = zero.clone();
                    for (column, coefficients) in &row.entries {
                        let epsilon_shift = column / n;
                        if epsilon_shift > e {
                            continue;
                        }
                        let source = (e - epsilon_shift) * n + column % n;
                        for (l, a) in coefficients.iter().enumerate().take(k + 1) {
                            if *a != zero {
                                sum = p.add(&sum, &p.mul(a, &y[k - l][source]));
                            }
                        }
                    }
                    for (l, a) in row.denominator.iter().enumerate().skip(1).take(k) {
                        if *a != zero {
                            sum = p.sub(&sum, &p.mul(a, &derivative[k - l][index]));
                        }
                    }
                    derivative[k][index] = sum;
                    y[k + 1][index] = p.scale(&derivative[k][index], 1, (k + 1) as i64);
                }
            }
        }
        Ok(y)
    }

    fn rhs(&self, point: &C, values: &[C]) -> Result<Vec<C>> {
        let p = self.rows.p;
        let n = self.rows.dimension();
        let mut rhs = vec![p.zero(); self.dimension()];
        // Evaluate each nonzero rational coefficient only once, sharing it
        // across epsilon orders when checking the differential-equation defect.
        for (i, row) in self.rows.polynomial_rows.iter().enumerate() {
            let denominator = shift_polynomial(p, &row.denominator, point, 0).remove(0);
            if denominator == p.zero() {
                return Err(Error::Numerical("path reached a pole".into()));
            }
            for (column, polynomial) in &row.entries {
                let a = p.div(&shift_polynomial(p, polynomial, point, 0)[0], &denominator);
                let shift = column / n;
                for e in shift..self.count {
                    rhs[e * n + i] = p.add(
                        &rhs[e * n + i],
                        &p.mul(&a, &values[(e - shift) * n + column % n]),
                    );
                }
            }
        }
        Ok(rhs)
    }
}

impl EpsilonSolution {
    /// Evaluate only when the retained path identifies one value at the point.
    /// Self-intersections carrying different monodromies are rejected.
    pub fn evaluate_path(&self, point: &C) -> Result<Vec<Vec<C>>> {
        let p = Precision {
            bits: self.diagnostics.working_bits,
        };
        let mut value: Option<Vec<Vec<C>>> = None;
        for i in 0..self.segments.len() {
            if let Ok(candidate) = self.evaluate_segment(i, point) {
                if let Some(previous) = &value {
                    if !candidate
                        .iter()
                        .flatten()
                        .zip(previous.iter().flatten())
                        .all(|(a, b)| p.close(a, b, 30))
                    {
                        return Err(Error::InvalidInput(
                            "saved path has multiple branch values at this point; select a segment"
                                .into(),
                        ));
                    }
                } else {
                    value = Some(candidate);
                }
            }
        }
        value.ok_or_else(|| Error::InvalidInput("point is outside the saved path".into()))
    }
    /// Evaluate a saved segment only on its traversed straight interval. This
    /// avoids silently switching analytic branches at overlapping disks. At a
    /// self-intersection the caller chooses the desired segment explicitly.
    pub fn evaluate_segment(&self, index: usize, point: &C) -> Result<Vec<Vec<C>>> {
        let segment = self
            .segments
            .get(index)
            .ok_or_else(|| Error::InvalidInput("unknown saved segment".into()))?;
        let p = Precision {
            bits: segment.working_bits,
        };
        let n = self.coefficients.first().map_or(0, Vec::len);
        if n == 0
            || segment.center == segment.end
            || segment.coefficients.is_empty()
            || segment
                .coefficients
                .iter()
                .any(|r| r.len() != n * self.coefficients.len())
        {
            return Err(Error::InvalidInput("malformed saved segment".into()));
        }
        if !p.finite(point) {
            return Err(Error::InvalidInput("nonfinite evaluation point".into()));
        }
        let t = p.div(
            &p.sub(point, &segment.center),
            &p.sub(&segment.end, &segment.center),
        );
        let tolerance = p.tolerance((p.bits / 4).min(30));
        if t.im > tolerance
            || t.im < -tolerance.clone()
            || t.re < -tolerance.clone()
            || t.re > p.real(1) + tolerance
        {
            return Err(Error::InvalidInput(
                "point is outside the saved segment interval".into(),
            ));
        }
        let (flat, _) = evaluate_taylor(p, &segment.coefficients, &p.sub(point, &segment.center));
        Ok(flat.chunks(n).map(<[C]>::to_vec).collect())
    }
}

/// Recompute the complete transport with increasing precision and series order.
/// The provider must supply fresh boundary coefficients at the requested
/// precision; any uncertainty in those inputs also limits the result. Exact
/// waypoints avoid rounding a path through binary64.
pub fn transport_epsilon(
    system: &EpsilonSystem,
    boundary_provider: impl Fn(Precision) -> Result<EpsilonBoundary>,
    waypoints: &[Atom],
    options: &FlowOptions,
    context: &RunContext,
    save_segments: bool,
) -> Result<EpsilonSolution> {
    options.validate()?;
    let mut previous: Option<EpsilonSolution> = None;
    for attempt in 0..=options.max_precision_attempts {
        context.cancellation.check()?;
        let mut refined = options.clone();
        refined.guard_digits = refined
            .guard_digits
            .checked_add(20 * attempt as u32)
            .ok_or_else(|| Error::Limit("working precision overflow".into()))?;
        refined.series_order = refined
            .series_order
            .checked_add(32 * attempt)
            .ok_or_else(|| Error::Limit("series order overflow".into()))?;
        refined.validate()?;
        let p = Precision::decimal(refined.digits + refined.guard_digits)?;
        let compiled = system.compile(p, &Default::default())?;
        let points = waypoints
            .iter()
            .map(|a| p.eval(a, &Default::default()))
            .collect::<Result<Vec<_>>>()?;
        let mut result = compiled.transport(
            &boundary_provider(p)?,
            &points,
            &refined,
            context,
            save_segments,
        )?;
        if let Some(old) = &previous {
            if old.leading != result.leading {
                return Err(Error::InvalidInput(
                    "boundary provider changed the leading epsilon power".into(),
                ));
            }
            result.comparison_errors = result
                .coefficients
                .iter()
                .zip(&old.coefficients)
                .map(|(new, old)| {
                    new.iter()
                        .zip(old)
                        .map(|(a, b)| p.norm(&p.sub(a, b)))
                        .collect()
                })
                .collect();
            if result
                .coefficients
                .iter()
                .flatten()
                .zip(old.coefficients.iter().flatten())
                .all(|(a, b)| p.close(a, b, options.digits))
            {
                result.verified_digits = Some(options.digits);
                if save_segments {
                    for (index, segment) in result.segments.iter().enumerate() {
                        if let Ok(reference) = old.evaluate_path(&segment.end) {
                            let values = result.evaluate_segment(index, &segment.end)?;
                            if values
                                .iter()
                                .flatten()
                                .zip(reference.iter().flatten())
                                .all(|(a, b)| p.close(a, b, options.digits))
                            {
                                let errors = values
                                    .iter()
                                    .zip(reference)
                                    .map(|(a, b)| {
                                        a.iter()
                                            .zip(b)
                                            .map(|(a, b)| p.norm(&p.sub(a, &b)))
                                            .collect()
                                    })
                                    .collect();
                                result.checkpoints.push(VerifiedCheckpoint {
                                    segment: index,
                                    point: segment.end.clone(),
                                    coefficients: values,
                                    comparison_errors: errors,
                                });
                            }
                        }
                    }
                }
                return Ok(result);
            }
        }
        previous = Some(result);
    }
    Err(Error::Accuracy(
        "epsilon transport did not stabilize under precision/order refinement".into(),
    ))
}
