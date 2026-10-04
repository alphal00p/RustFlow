//! Rational linear differential systems and arbitrary-precision Taylor transport.
use crate::{ComplexFloat as C, Error, FlowOptions, Precision, Progress, Result, RunContext};
use std::sync::Arc;
use symbolica::domains::float::FloatField;
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct DifferentialSystem {
    pub variable: Symbol,
    pub matrix: Vec<Vec<Atom>>,
}

#[derive(Clone, Debug)]
pub struct BoundaryData {
    pub point: C,
    pub values: Vec<C>,
}

#[derive(Clone, Debug, Default)]
pub struct FlowDiagnostics {
    pub steps: usize,
    pub rejected_steps: usize,
    pub working_bits: u32,
    pub expansion_order: usize,
}

#[derive(Clone, Debug)]
pub struct FlowResult {
    pub point: C,
    pub values: Vec<C>,
    pub diagnostics: FlowDiagnostics,
}

#[derive(Clone, Debug)]
pub(crate) struct NumericRational {
    pub numerator: Vec<C>,
    pub denominator: Vec<C>,
}

pub(crate) fn polynomial_coefficients(
    a: &Atom,
    variable: Symbol,
    p: Precision,
    values: &ahash::HashMap<Atom, C>,
) -> Result<Vec<C>> {
    let x = Atom::var(variable);
    let mut result = vec![p.zero()];
    for (monomial, coefficient) in a.coefficient_list::<i32>(std::slice::from_ref(&x)) {
        let degree = if monomial.is_one() {
            0
        } else if monomial == x {
            1
        } else if let AtomView::Pow(power) = monomial.as_view() {
            let (base, exponent) = power.get_base_exp();
            if base != x.as_view() {
                return Err(Error::InvalidInput("non-polynomial coefficient".into()));
            }
            exponent
                .to_string()
                .parse::<usize>()
                .map_err(|_| Error::InvalidInput("noninteger polynomial power".into()))?
        } else {
            return Err(Error::InvalidInput(format!(
                "non-polynomial monomial {monomial}"
            )));
        };
        if degree > 100_000 {
            return Err(Error::Limit("polynomial degree exceeds 100000".into()));
        }
        result.resize(result.len().max(degree + 1), p.zero());
        result[degree] = p.add(&result[degree], &p.eval(&coefficient, values)?);
    }
    while result.len() > 1 && result.last() == Some(&p.zero()) {
        result.pop();
    }
    Ok(result)
}

pub(crate) fn shift_polynomial(
    p: Precision,
    coefficients: &[C],
    center: &C,
    order: usize,
) -> Vec<C> {
    let mut out = vec![p.zero(); order + 1];
    for a in coefficients.iter().rev() {
        for k in (1..=order).rev() {
            out[k] = p.add(&p.mul(&out[k], center), &out[k - 1]);
        }
        out[0] = p.add(&p.mul(&out[0], center), a);
    }
    out
}

pub(crate) fn quotient_series(
    p: Precision,
    numerator: &[C],
    denominator: &[C],
    order: usize,
) -> Result<Vec<C>> {
    if denominator.is_empty() || denominator[0] == p.zero() {
        return Err(Error::Numerical("expansion center is a pole".into()));
    }
    let mut out = vec![p.zero(); order + 1];
    for n in 0..=order {
        let mut v = numerator.get(n).cloned().unwrap_or_else(|| p.zero());
        for j in 1..=n.min(denominator.len() - 1) {
            v = p.sub(&v, &p.mul(&denominator[j], &out[n - j]));
        }
        out[n] = p.div(&v, &denominator[0]);
    }
    Ok(out)
}

fn polynomial_roots(p: Precision, coefficients: &[C], variable: Symbol) -> Result<Vec<C>> {
    let degree = coefficients.len() - 1;
    if degree == 1 {
        return Ok(vec![p.neg(&p.div(&coefficients[0], &coefficients[1]))]);
    }
    // Fujiwara's bound, rounded upward to a power of two. With x = scale*z,
    // every root lies in the unit disk and the monic coefficients stay bounded.
    // Symbolica's Aberth solver uses an absolute polynomial residual; passing
    // raw integer coefficients can make its stopping tolerance unattainable.
    let leading = &coefficients[degree];
    let mut exponent = None;
    for (i, coefficient) in coefficients[..degree].iter().enumerate() {
        let magnitude = p.norm(&p.div(coefficient, leading));
        if !magnitude.is_finite() {
            return Err(Error::Numerical("nonfinite pole polynomial".into()));
        }
        if let Some(e) = magnitude.as_raw().get_exp() {
            let n = (degree - i) as i64;
            let bound = (i64::from(e) + n - 1).div_euclid(n) + 1;
            exponent = Some(exponent.map_or(bound, |old: i64| old.max(bound)));
        }
    }
    let scale = p.powi(&p.i(2), exponent.unwrap_or(0));
    if !p.finite(&scale) || scale == p.zero() {
        return Err(Error::Numerical(
            "pole scaling exceeds numerical range".into(),
        ));
    }
    let mut scaled = coefficients
        .iter()
        .enumerate()
        .map(|(i, c)| p.div(&p.div(c, leading), &p.powi(&scale, (degree - i) as i64)))
        .collect::<Vec<_>>();
    scaled[degree] = p.i(1);
    let poly = UnivariatePolynomial::from_coefficients(
        &FloatField::from_rep(p.zero()),
        scaled.clone(),
        Arc::new(PolyVariable::Symbol(variable)),
    );
    let roots = poly.roots(1000, &p.tolerance(p.bits / 4)).map_err(|_| {
        Error::Numerical(format!(
            "degree-{degree} scaled pole root finding did not converge"
        ))
    })?;
    if roots.len() != degree || roots.iter().any(|root| !p.finite(root)) {
        return Err(Error::Numerical("incomplete pole root set".into()));
    }
    // Verify Newton corrections as well as the polynomial reconstructed from
    // all roots. A small residual alone is weak near clustered roots.
    let tolerance = p.tolerance(p.bits / 5);
    let mut reconstructed = vec![p.i(1)];
    for (i, root) in roots.iter().enumerate() {
        let mut value = p.i(1);
        let mut derivative = p.zero();
        for coefficient in scaled[..degree].iter().rev() {
            derivative = p.add(&p.mul(&derivative, root), &value);
            value = p.add(&p.mul(&value, root), coefficient);
        }
        if derivative == p.zero()
            || p.norm(&p.div(&value, &derivative)) > tolerance
            || roots[..i]
                .iter()
                .any(|other| p.norm(&p.sub(root, other)) <= tolerance)
        {
            return Err(Error::Numerical(
                "unresolved or ill-conditioned pole roots".into(),
            ));
        }
        let mut next = vec![p.zero(); reconstructed.len() + 1];
        for (j, coefficient) in reconstructed.iter().enumerate() {
            next[j] = p.sub(&next[j], &p.mul(root, coefficient));
            next[j + 1] = p.add(&next[j + 1], coefficient);
        }
        reconstructed = next;
    }
    if reconstructed
        .iter()
        .zip(&scaled)
        .any(|(a, b)| p.norm(&p.sub(a, b)) > tolerance)
    {
        return Err(Error::Numerical(
            "pole roots fail polynomial reconstruction".into(),
        ));
    }
    Ok(roots.iter().map(|root| p.mul(root, &scale)).collect())
}

impl NumericRational {
    pub(crate) fn series(&self, p: Precision, center: &C, order: usize) -> Result<Vec<C>> {
        quotient_series(
            p,
            &shift_polynomial(p, &self.numerator, center, order),
            &shift_polynomial(p, &self.denominator, center, order),
            order,
        )
    }
}

type ExactPolynomial = MultivariatePolynomial<IntegerRing, u16>;

fn multiply_polynomials(a: &ExactPolynomial, b: &ExactPolynomial) -> Result<ExactPolynomial> {
    for (i, variable) in a.variables().iter().enumerate() {
        if let Some(j) = b.variables().iter().position(|v| v == variable)
            && a.degree(i).checked_add(b.degree(j)).is_none()
        {
            return Err(Error::Limit(
                "cleared differential polynomial degree overflow".into(),
            ));
        }
    }
    Ok(a * b)
}

/// A row of D(x) y'(x) = A(x) y(x), with an exact common denominator
/// cleared before numerical specialization. Structural zero entries are absent.
#[derive(Clone, Debug)]
struct PolynomialRow {
    denominator: Vec<C>,
    entries: Vec<(usize, Vec<C>)>,
}

#[derive(Clone, Debug)]
pub struct CompiledSystem {
    pub(crate) p: Precision,
    pub(crate) matrix: Vec<Vec<NumericRational>>,
    pub poles: Vec<C>,
    polynomial_rows: Vec<PolynomialRow>,
}

impl DifferentialSystem {
    /// Change basis using old_y = transformation * new_y, preserving the
    /// derivative term: B = T^-1 A T - T^-1 dT/dx.
    pub fn change_basis(&self, transformation: &[Vec<Atom>]) -> Result<Self> {
        self.validate()?;
        if transformation.len() != self.matrix.len()
            || transformation.iter().any(|r| r.len() != self.matrix.len())
        {
            return Err(Error::InvalidInput(
                "basis transformation dimensions".into(),
            ));
        }
        let inverse = crate::algebra::inverse(transformation)?;
        let mut product = crate::algebra::matmul(&self.matrix, transformation);
        for (i, row) in product.iter_mut().enumerate() {
            for (j, a) in row.iter_mut().enumerate() {
                *a = (&*a - transformation[i][j].derivative(self.variable))
                    .together()
                    .cancel();
            }
        }
        Ok(Self {
            variable: self.variable,
            matrix: crate::algebra::matmul(&inverse, &product),
        })
    }
    pub fn validate(&self) -> Result<()> {
        let n = self.matrix.len();
        if n == 0 || self.matrix.iter().any(|r| r.len() != n) {
            return Err(Error::InvalidInput(
                "nonempty square differential matrix required".into(),
            ));
        }
        Ok(())
    }
    pub fn compile(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
    ) -> Result<CompiledSystem> {
        self.validate()?;
        let mut matrix = Vec::new();
        let mut poles = Vec::new();
        let mut seen_factors = ahash::HashSet::default();
        let mut polynomial_rows = Vec::new();
        for row in &self.matrix {
            let mut out = Vec::new();
            let mut exact = Vec::new();
            let mut common_denominator: Option<ExactPolynomial> = None;
            for a in row {
                let rational: RationalPolynomial<IntegerRing, u16> = a
                    .try_to_rational_polynomial(&Q, &Z, None)
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                let numerator = polynomial_coefficients(
                    &rational.numerator.to_expression(),
                    self.variable,
                    p,
                    values,
                )?;
                let denominator = polynomial_coefficients(
                    &rational.denominator.to_expression(),
                    self.variable,
                    p,
                    values,
                )?;
                if denominator.iter().all(|v| *v == p.zero()) {
                    return Err(Error::Numerical(
                        "identically zero specialized denominator".into(),
                    ));
                }
                let factored = rational.denominator.to_expression().factor();
                let factors = if let AtomView::Mul(m) = factored.as_view() {
                    m.iter().map(|v| v.to_owned()).collect::<Vec<_>>()
                } else {
                    vec![factored]
                };
                for factor in factors {
                    let base = if let AtomView::Pow(v) = factor.as_view() {
                        v.get_base_exp().0.to_owned()
                    } else {
                        factor
                    };
                    if !seen_factors.insert(base.clone()) {
                        continue;
                    }
                    let coefficients = polynomial_coefficients(&base, self.variable, p, values)?;
                    if coefficients.len() <= 1 {
                        continue;
                    }
                    let roots = polynomial_roots(p, &coefficients, self.variable)?;
                    for root in roots {
                        if !poles.iter().any(|v| {
                            let a = p.norm(v);
                            let b = p.norm(&root);
                            let scale = if a < b { a } else { b };
                            p.norm(&p.sub(v, &root)) <= p.tolerance(p.bits / 5) * scale
                        }) {
                            poles.push(root);
                        }
                    }
                }
                out.push(NumericRational {
                    numerator,
                    denominator,
                });
                common_denominator = Some(if let Some(previous) = common_denominator {
                    let quotient = previous
                        .try_div(&previous.gcd(&rational.denominator))
                        .ok_or_else(|| {
                            Error::Numerical("exact denominator LCM division failed".into())
                        })?;
                    multiply_polynomials(&quotient, &rational.denominator)?
                } else {
                    rational.denominator.clone()
                });
                exact.push(rational);
            }
            let common_denominator = common_denominator.unwrap();
            let denominator = polynomial_coefficients(
                &common_denominator.to_expression(),
                self.variable,
                p,
                values,
            )?;
            let mut entries = Vec::new();
            for (j, rational) in exact.iter().enumerate() {
                if rational.numerator.is_zero() {
                    continue;
                }
                let multiplier = common_denominator
                    .try_div(&rational.denominator)
                    .ok_or_else(|| Error::Numerical("exact denominator clearing failed".into()))?;
                let cleared = multiply_polynomials(&rational.numerator, &multiplier)?;
                let coefficients =
                    polynomial_coefficients(&cleared.to_expression(), self.variable, p, values)?;
                if coefficients.iter().any(|c| *c != p.zero()) {
                    entries.push((j, coefficients));
                }
            }
            polynomial_rows.push(PolynomialRow {
                denominator,
                entries,
            });
            matrix.push(out);
        }
        Ok(CompiledSystem {
            p,
            matrix,
            poles,
            polynomial_rows,
        })
    }
    /// Strongly connected blocks in dependency order, independent of input ordering.
    pub fn blocks(&self) -> Result<Vec<Vec<usize>>> {
        self.validate()?;
        let n = self.matrix.len();
        let mut reach = vec![vec![false; n]; n];
        for (i, row) in self.matrix.iter().enumerate() {
            for (j, a) in row.iter().enumerate() {
                reach[i][j] = i == j || !a.is_zero();
            }
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    reach[i][j] |= reach[i][k] && reach[k][j];
                }
            }
        }
        let mut remaining = (0..n).collect::<std::collections::BTreeSet<_>>();
        let mut blocks = Vec::new();
        while !remaining.is_empty() {
            let i = *remaining
                .iter()
                .find(|&&i| remaining.iter().all(|&j| !reach[i][j] || reach[j][i]))
                .unwrap();
            let block = remaining
                .iter()
                .copied()
                .filter(|&j| reach[i][j] && reach[j][i])
                .collect::<Vec<_>>();
            for j in &block {
                remaining.remove(j);
            }
            blocks.push(block);
        }
        Ok(blocks)
    }
}

impl CompiledSystem {
    /// A real radius strictly inside the local convergence disk around an
    /// endpoint, excluding only a pole exactly at that endpoint.
    pub fn endpoint_radius(&self, endpoint: &C) -> Float {
        let p = self.p;
        self.poles
            .iter()
            .map(|pole| p.norm(&p.sub(pole, endpoint)))
            .filter(|r| *r > p.real(0))
            .fold(p.real(1), |a, b| if a < b { a } else { b })
            / 8
    }
    pub fn dimension(&self) -> usize {
        self.matrix.len()
    }
    pub(crate) fn taylor(&self, center: &C, values: &[C], order: usize) -> Result<Vec<Vec<C>>> {
        let p = self.p;
        let n = self.dimension();
        if values.len() != n {
            return Err(Error::InvalidInput("boundary dimension".into()));
        }
        let zero = p.zero();
        // Translate only the finite polynomial degrees needed at this order;
        // no rational matrix entry is expanded into an order-long power series.
        let rows = self
            .polynomial_rows
            .iter()
            .map(|row| {
                let mut denominator = shift_polynomial(
                    p,
                    &row.denominator,
                    center,
                    (row.denominator.len() - 1).min(order.saturating_sub(1)),
                );
                if denominator[0] == zero {
                    return Err(Error::Numerical("expansion center is a pole".into()));
                }
                let divisor = denominator[0].clone();
                for c in &mut denominator {
                    *c = p.div(c, &divisor);
                }
                let entries = row
                    .entries
                    .iter()
                    .map(|(j, coefficients)| {
                        let mut shifted = shift_polynomial(
                            p,
                            coefficients,
                            center,
                            (coefficients.len() - 1).min(order.saturating_sub(1)),
                        );
                        for c in &mut shifted {
                            *c = p.div(c, &divisor);
                        }
                        (*j, shifted)
                    })
                    .collect::<Vec<_>>();
                Ok(PolynomialRow {
                    denominator,
                    entries,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut y = vec![vec![zero.clone(); n]; order + 1];
        // Store y' coefficients as well as y. Then the denominator convolution
        // does not repeatedly multiply y_k by its integer derivative weight.
        let mut derivative = vec![vec![zero.clone(); n]; order];
        y[0] = values.to_vec();
        // With D_0 normalized to one:
        // y'_{k,i} = sum_{j,l} A_{ij,l} y_{k-l,j}
        //             - sum_{l>=1} D_{i,l} y'_{k-l,i}.
        for k in 0..order {
            for (i, row) in rows.iter().enumerate() {
                let mut sum = zero.clone();
                for (j, coefficients) in &row.entries {
                    for (l, coefficient) in coefficients.iter().enumerate().take(k + 1) {
                        if *coefficient != zero {
                            sum = p.add(&sum, &p.mul(coefficient, &y[k - l][*j]));
                        }
                    }
                }
                for (l, coefficient) in row.denominator.iter().enumerate().skip(1).take(k) {
                    if *coefficient != zero {
                        sum = p.sub(&sum, &p.mul(coefficient, &derivative[k - l][i]));
                    }
                }
                derivative[k][i] = sum;
                y[k + 1][i] = p.scale(&derivative[k][i], 1, (k + 1) as i64);
            }
        }
        Ok(y)
    }
    /// Plan a polygonal contour around disjoint pole disks. `side` is +1
    /// for a counterclockwise detour relative to the segment, or -1 for the
    /// opposite side. The choice fixes the continuation homotopy.
    pub fn plan_path(&self, start: &C, end: &C, side: i64) -> Result<Vec<C>> {
        if ![-1, 1].contains(&side) || !self.p.finite(start) || !self.p.finite(end) {
            return Err(Error::InvalidInput(
                "invalid contour endpoints or side".into(),
            ));
        }
        let p = self.p;
        let mut disks = Vec::new();
        for (i, pole) in self.poles.iter().enumerate() {
            let mut radius = p.norm(&p.sub(start, pole));
            let end_distance = p.norm(&p.sub(end, pole));
            if end_distance < radius {
                radius = end_distance;
            }
            if radius == p.real(0) {
                return Err(Error::InvalidInput("contour endpoint is a pole".into()));
            }
            for (j, other) in self.poles.iter().enumerate() {
                if i != j && pole != other {
                    let distance = p.norm(&p.sub(pole, other));
                    if distance < radius {
                        radius = distance;
                    }
                }
            }
            disks.push((pole, radius / 4));
        }
        let mut segments = vec![(start.clone(), end.clone())];
        let mut output = Vec::new();
        let mut iterations = 0;
        while let Some((a, b)) = segments.pop() {
            iterations += 1;
            if iterations > 4096 {
                return Err(Error::Limit(
                    "contour planning exceeds 4096 subdivisions".into(),
                ));
            }
            let delta = p.sub(&b, &a);
            let length = p.norm(&delta);
            if length == p.real(0) {
                output.push(b);
                continue;
            }
            let obstruction = disks
                .iter()
                .filter_map(|(pole, radius)| {
                    let t = p.div(&p.sub(pole, &a), &delta).re;
                    if t <= p.real(0) || t >= p.real(1) {
                        return None;
                    }
                    let projection = p.add(&a, &p.mul(&delta, &C::new(t.clone(), p.real(0))));
                    (p.norm(&p.sub(pole, &projection)) < *radius).then_some((t, *pole, radius))
                })
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            if let Some((_, pole, radius)) = obstruction {
                let normal = p.div(
                    &p.mul(&delta, &p.complex(0, side)),
                    &C::new(length, p.real(0)),
                );
                let waypoint = p.add(
                    pole,
                    &p.mul(&normal, &C::new(radius.clone() * 2, p.real(0))),
                );
                segments.push((waypoint.clone(), b));
                segments.push((a, waypoint));
            } else {
                output.push(b);
            }
        }
        Ok(output)
    }

    fn residual_is_small(
        &self,
        center: &C,
        step: &C,
        coefficients: &[Vec<C>],
        values: &[C],
        tolerance: &Float,
    ) -> Result<bool> {
        let p = self.p;
        let point = p.add(center, step);
        for (i, row) in self.matrix.iter().enumerate() {
            let mut derivative = p.zero();
            for k in (1..coefficients.len()).rev() {
                derivative = p.add(
                    &p.mul(&derivative, step),
                    &p.scale(&coefficients[k][i], k as i64, 1),
                );
            }
            let mut rhs = p.zero();
            for (entry, value) in row.iter().zip(values) {
                rhs = p.add(&rhs, &p.mul(&entry.series(p, &point, 0)?[0], value));
            }
            let defect = p.mul(step, &p.sub(&derivative, &rhs));
            let magnitude = p.norm(&values[i]);
            let scale = if magnitude > p.real(1) {
                magnitude
            } else {
                p.real(1)
            };
            if !p.finite(&defect) || p.norm(&defect) > tolerance.clone() * scale {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub fn transport(
        &self,
        boundary: &BoundaryData,
        waypoints: &[C],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<FlowResult> {
        options.validate()?;
        let p = self.p;
        if boundary.values.len() != self.dimension()
            || !p.finite(&boundary.point)
            || boundary.values.iter().any(|v| !p.finite(v))
            || waypoints.iter().any(|v| !p.finite(v))
        {
            return Err(Error::InvalidInput("invalid boundary or path".into()));
        }
        let mut center = boundary.point.clone();
        let mut values = boundary.values.clone();
        let mut diagnostics = FlowDiagnostics {
            working_bits: p.bits,
            expansion_order: options.series_order,
            ..Default::default()
        };
        for target in waypoints {
            while center != *target {
                context.emit(Progress::Step {
                    index: diagnostics.steps,
                })?;
                if diagnostics.steps + diagnostics.rejected_steps >= options.max_steps {
                    return Err(Error::Limit(
                        "Taylor continuation step budget exhausted".into(),
                    ));
                }
                let delta = p.sub(target, &center);
                let distance = p.norm(&delta);
                let radius = self
                    .poles
                    .iter()
                    .map(|s| p.norm(&p.sub(&center, s)))
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap_or_else(|| distance.clone() * 2);
                if radius == p.real(0) {
                    return Err(Error::Numerical(
                        "path reached a differential-equation pole".into(),
                    ));
                }
                let mut step = delta.clone();
                let safe = radius / 3;
                if distance > safe {
                    step = p.mul(
                        &delta,
                        &p.div(&C::new(safe, p.real(0)), &C::new(distance, p.real(0))),
                    );
                }
                let coefficients = self.taylor(&center, &values, options.series_order)?;
                let mut accepted = None;
                for _ in 0..32 {
                    let (v, tail) = evaluate_taylor(p, &coefficients, &step);
                    // Regulator fitting and endpoint matching consume guard
                    // digits too. Increasing arithmetic precision must also
                    // tighten truncation, even at a fixed expansion order.
                    let truncation_digits = options
                        .digits
                        .saturating_add(8)
                        .max((options.digits + options.guard_digits).saturating_sub(10));
                    let available_digits = (u64::from(p.bits) * 1000 / 3322) as u32;
                    let tolerance =
                        p.tolerance(truncation_digits.min(available_digits.saturating_sub(3)));
                    let good = v.iter().zip(&tail).all(|(v, t)| {
                        let scale = p.norm(v);
                        let scale = if scale > p.real(1) { scale } else { p.real(1) };
                        p.finite(v) && *t <= tolerance.clone() * scale
                    });
                    if good
                        && self.residual_is_small(&center, &step, &coefficients, &v, &tolerance)?
                    {
                        // Vanishing final Taylor terms do not bound omitted
                        // terms for sparse systems such as y' = x^20 y.
                        // Check the differential equation at the step endpoint
                        // and midpoint as an independent defect test.
                        let half_step = p.scale(&step, 1, 2);
                        let (middle, _) = evaluate_taylor(p, &coefficients, &half_step);
                        if self.residual_is_small(
                            &center,
                            &half_step,
                            &coefficients,
                            &middle,
                            &tolerance,
                        )? {
                            accepted = Some(v);
                            break;
                        }
                    }
                    step = p.scale(&step, 1, 2);
                    diagnostics.rejected_steps += 1;
                    if diagnostics.steps + diagnostics.rejected_steps >= options.max_steps {
                        break;
                    }
                }
                values = accepted
                    .ok_or_else(|| Error::Accuracy("Taylor tail did not meet tolerance".into()))?;
                let next = if step == delta {
                    target.clone()
                } else {
                    p.add(&center, &step)
                };
                if next == center {
                    return Err(Error::Accuracy("continuation step lost to rounding".into()));
                }
                center = next;
                diagnostics.steps += 1;
            }
        }
        Ok(FlowResult {
            point: center,
            values,
            diagnostics,
        })
    }
}

fn evaluate_taylor(p: Precision, coefficients: &[Vec<C>], h: &C) -> (Vec<C>, Vec<Float>) {
    let n = coefficients[0].len();
    let mut out = vec![p.zero(); n];
    let mut tail = vec![p.real(0); n];
    let start = coefficients.len().saturating_sub(6);
    for k in (0..coefficients.len()).rev() {
        for i in 0..n {
            out[i] = p.add(&p.mul(&out[i], h), &coefficients[k][i]);
        }
    }
    let mut power = p.powi(h, start as i64);
    for row in coefficients.iter().skip(start) {
        for (i, c) in row.iter().enumerate() {
            tail[i] += p.norm(&p.mul(c, &power));
        }
        power = p.mul(&power, h);
    }
    (out, tail)
}

#[cfg(test)]
mod recurrence_tests {
    use super::*;

    // Independent reference recurrence: first expand each rational matrix entry,
    // then match coefficients of y' = M y without clearing denominators.
    fn rational_series_taylor(
        system: &CompiledSystem,
        center: &C,
        values: &[C],
        order: usize,
    ) -> Result<Vec<Vec<C>>> {
        let p = system.p;
        let matrix = system
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|entry| entry.series(p, center, order.saturating_sub(1)))
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut y = vec![vec![p.zero(); values.len()]; order + 1];
        y[0] = values.to_vec();
        for k in 0..order {
            for (i, row) in matrix.iter().enumerate() {
                let mut sum = p.zero();
                for (j, entry) in row.iter().enumerate() {
                    for (l, a) in entry.iter().enumerate().take(k + 1) {
                        sum = p.add(&sum, &p.mul(a, &y[k - l][j]));
                    }
                }
                y[k + 1][i] = p.scale(&sum, 1, (k + 1) as i64);
            }
        }
        Ok(y)
    }

    #[test]
    fn cleared_recurrence_matches_rational_series_at_complex_centers() {
        let p = Precision::decimal(80).unwrap();
        let x = Atom::var(symbol!("recurrence_x"));
        let epsilon = Atom::var(symbol!("recurrence_eps"));
        let first = &x - Atom::num(1);
        let second = &x + Atom::num(2);
        let system = DifferentialSystem {
            variable: symbol!("recurrence_x"),
            matrix: vec![
                vec![
                    (x.clone().pow(3) + Atom::num(2)) / first.clone().pow(2),
                    &epsilon / &second,
                    Atom::new(),
                ],
                vec![Atom::num(1) / (&first * &second), Atom::new(), x.pow(20)],
                vec![Atom::new(); 3],
            ],
        }
        .compile(
            p,
            &ahash::HashMap::from_iter([(epsilon, p.rational(&Rational::from((1, 37))))]),
        )
        .unwrap();
        assert_eq!(system.polynomial_rows[0].denominator.len(), 4);
        assert_eq!(system.polynomial_rows[0].entries.len(), 2);
        assert!(system.polynomial_rows[2].entries.is_empty());
        let values = vec![p.i(2), p.complex(3, 1), p.i(7)];
        for center in [p.zero(), p.parse("0.375", "0.25").unwrap()] {
            for order in [0, 1, 12, 32] {
                let finite = system.taylor(&center, &values, order).unwrap();
                let reference = rational_series_taylor(&system, &center, &values, order).unwrap();
                for (a, b) in finite.iter().flatten().zip(reference.iter().flatten()) {
                    assert!(p.close(a, b, 65), "order {order}: {a} != {b}");
                }
                assert!(finite.iter().skip(1).all(|row| row[2] == p.zero()));
            }
        }
        assert!(matches!(
            system.taylor(&p.i(1), &values, 12),
            Err(Error::Numerical(_))
        ));
    }

    #[test]
    fn polynomial_clearing_reports_degree_overflow() {
        let a: RationalPolynomial<IntegerRing, u16> = parse!("overflow_x^40000")
            .try_to_rational_polynomial(&Q, &Z, None)
            .unwrap();
        assert!(matches!(
            multiply_polynomials(&a.numerator, &a.numerator),
            Err(Error::Limit(_))
        ));
    }
}
