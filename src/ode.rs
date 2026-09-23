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

#[derive(Clone, Debug)]
pub struct CompiledSystem {
    pub(crate) p: Precision,
    pub(crate) matrix: Vec<Vec<NumericRational>>,
    pub poles: Vec<C>,
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
        for row in &self.matrix {
            let mut out = Vec::new();
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
                    let coefficients = polynomial_coefficients(&base, self.variable, p, values)?;
                    if coefficients.len() <= 1 {
                        continue;
                    }
                    let field = FloatField::from_rep(p.zero());
                    let poly = UnivariatePolynomial::from_coefficients(
                        &field,
                        coefficients,
                        Arc::new(PolyVariable::Symbol(self.variable)),
                    );
                    let roots = poly.roots(1000, &p.tolerance(p.bits / 4)).map_err(|_| {
                        Error::Numerical("pole root finding did not converge".into())
                    })?;
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
            }
            matrix.push(out);
        }
        Ok(CompiledSystem { p, matrix, poles })
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
        let a = self
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| v.series(p, center, order - 1))
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut y = vec![vec![p.zero(); n]; order + 1];
        y[0] = values.to_vec();
        for k in 0..order {
            for (i, row) in a.iter().enumerate() {
                let mut sum = p.zero();
                for j in 0..n {
                    for l in 0..=k {
                        sum = p.add(&sum, &p.mul(&row[j][l], &y[k - l][j]));
                    }
                }
                y[k + 1][i] = p.scale(&sum, 1, (k + 1) as i64);
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
                    let tolerance = p.tolerance(options.digits + 8);
                    let good = v.iter().zip(&tail).all(|(v, t)| {
                        let scale = p.norm(v);
                        let scale = if scale > p.real(1) { scale } else { p.real(1) };
                        p.finite(v) && *t <= tolerance.clone() * scale
                    });
                    if good {
                        accepted = Some(v);
                        break;
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
