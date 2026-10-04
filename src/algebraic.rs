//! Regular-contour transport with explicitly registered square-root sheets.
//!
//! Roots are independent algebraic generators satisfying r² = R(x). The first
//! implementation accepts root monomials over rational functions and rejects
//! sums of roots in denominators. It reuses Symbolica's exact series and fixed-
//! precision numeric series multiplication. Branch state belongs to each call;
//! only accepted continuation steps can advance it. Singular endpoints require
//! a separate Frobenius treatment and are outside this regular-contour API.
use crate::diffexp::{EpsilonBoundary, EpsilonSolution, EpsilonSystem};
use crate::family::{encode_complex, imaginary_parameter, scalar_symbols, substitute};
use crate::fixed_series::{coefficients, evaluate_series, fixed_series};
use crate::ode::{
    NumericRational, SeriesSystem, compile_rows, evaluate_taylor, transport_series_with_state,
};
use crate::{
    BoundaryData, ComplexFloat as C, DifferentialSystem, Error, FlowOptions, Precision,
    Prescription, Result, RunContext,
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct SquareRoot {
    pub symbol: Symbol,
    pub radicand: Atom,
}

#[derive(Clone, Debug)]
pub struct AlgebraicSystem {
    pub system: EpsilonSystem,
    pub roots: Vec<SquareRoot>,
}

/// Explicit initial sheet choice. `Value` selects only a discrete sign: its
/// magnitude is recomputed from the exact radicand at the working precision.
/// A zero, nonfinite, or ambiguous sign hint is rejected. `I0` means the side
/// of the *radicand* at a real starting point, not the side of the path variable.
#[derive(Clone, Debug)]
pub enum RootSeed {
    Principal,
    Opposite,
    I0(Prescription),
    Value(C),
}

/// Roots retain declaration order, identity and endpoint for branch comparisons.
/// Values have working precision; this type does not claim verified accuracy.
#[derive(Clone, Debug)]
pub struct BranchState {
    pub point: C,
    pub roots: Vec<(Symbol, C)>,
}

#[derive(Clone, Debug)]
pub struct AlgebraicSolution {
    pub solution: EpsilonSolution,
    pub branches: BranchState,
}

#[derive(Clone, Debug)]
struct RootKernel {
    definition: SquareRoot,
    value: NumericRational,
    logarithmic_derivative: NumericRational,
}
#[derive(Clone, Debug)]
struct KernelTerm {
    shift: usize,
    row: usize,
    column: usize,
    coefficient: NumericRational,
    roots: Vec<usize>,
}
#[derive(Clone, Debug)]
pub struct CompiledAlgebraicSystem {
    variable: Symbol,
    p: Precision,
    size: usize,
    count: usize,
    roots: Vec<RootKernel>,
    terms: Vec<KernelTerm>,
    poles: Vec<C>,
}

fn exact_point(p: Precision, value: &C) -> Result<Atom> {
    if !p.finite(value) {
        return Err(Error::InvalidInput(
            "nonfinite algebraic chart coordinate".into(),
        ));
    }
    Ok(Atom::num(Complex::new(
        value.re.to_rational(),
        value.im.to_rational(),
    )))
}
fn decoded(a: &Atom) -> Atom {
    substitute(
        a,
        &BTreeMap::from([(
            Atom::var(imaginary_parameter()),
            Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
        )]),
    )
    .together()
    .cancel()
}
fn rational(a: &Atom, allowed: &BTreeSet<Atom>) -> Result<RationalPolynomial<IntegerRing, u16>> {
    let mut symbols = BTreeSet::new();
    scalar_symbols(a.as_view(), &mut symbols)?;
    if !symbols.is_subset(allowed) {
        return Err(Error::Unsupported(
            "algebraic coefficients contain undeclared symbols".into(),
        ));
    }
    encode_complex(a)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| {
            Error::Unsupported(format!(
                "rational functions of the declared root generators are required: {e}"
            ))
        })
}
fn powers(monomial: &Atom, roots: &[Atom]) -> Result<Vec<i64>> {
    roots
        .iter()
        .map(|root| {
            let symbol = if let AtomView::Var(v) = root.as_view() {
                v.get_symbol()
            } else {
                unreachable!()
            };
            // Native exact differentiation extracts the degree of a monomial,
            // including its constant case, without parsing product structure.
            (monomial.derivative(symbol) * root / monomial)
                .cancel()
                .to_string()
                .parse::<i64>()
                .map_err(|_| Error::Unsupported("noninteger root monomial".into()))
        })
        .collect()
}

fn normalized_polynomial(
    expression: &Atom,
    root_atoms: &[Atom],
    roots: &[SquareRoot],
) -> Result<BTreeMap<Vec<usize>, Atom>> {
    let mut terms = BTreeMap::<Vec<usize>, Atom>::new();
    for (monomial, mut coefficient) in expression.coefficient_list::<i32>(root_atoms) {
        let mut odd = Vec::new();
        for (index, exponent) in powers(&monomial, root_atoms)?.into_iter().enumerate() {
            coefficient *= encode_complex(&roots[index].radicand).pow(exponent.div_euclid(2));
            if exponent.rem_euclid(2) == 1 {
                odd.push(index);
            }
        }
        let previous = terms.entry(odd).or_default();
        *previous = encode_complex(&decoded(&(&*previous + coefficient)))
            .together()
            .cancel();
    }
    terms.retain(|_, coefficient| !coefficient.is_zero());
    Ok(terms)
}

impl AlgebraicSystem {
    pub fn ordinary(system: DifferentialSystem, roots: Vec<SquareRoot>) -> Self {
        Self {
            system: EpsilonSystem {
                variable: system.variable,
                matrices: vec![system.matrix],
            },
            roots,
        }
    }

    pub fn validate(&self) -> Result<()> {
        self.system.validate()?;
        let mut declared = BTreeSet::from([
            Atom::var(self.system.variable),
            Atom::var(imaginary_parameter()),
        ]);
        for root in &self.roots {
            if !declared.insert(Atom::var(root.symbol)) {
                return Err(Error::InvalidInput("root symbols must be unique and differ from the path parameter and reserved imaginary unit".into()));
            }
        }
        let allowed = BTreeSet::from([
            Atom::var(self.system.variable),
            Atom::var(imaginary_parameter()),
        ]);
        for root in &self.roots {
            rational(&root.radicand, &allowed)?;
            if decoded(&root.radicand).is_zero() {
                return Err(Error::InvalidInput(
                    "a square-root radicand is identically zero".into(),
                ));
            }
        }
        for matrix in &self.system.matrices {
            for row in matrix {
                for a in row {
                    rational(a, &declared)?;
                }
            }
        }
        Ok(())
    }

    pub fn compile(&self, p: Precision) -> Result<CompiledAlgebraicSystem> {
        if p.bits < 2 {
            return Err(Error::InvalidInput(
                "algebraic precision must be at least two bits".into(),
            ));
        }
        self.validate()?;
        let root_atoms = self
            .roots
            .iter()
            .map(|r| Atom::var(r.symbol))
            .collect::<Vec<_>>();
        let mut allowed = BTreeSet::from([
            Atom::var(self.system.variable),
            Atom::var(imaginary_parameter()),
        ]);
        allowed.extend(root_atoms.iter().cloned());
        let mut entries = Vec::new();
        let mut normalized = Vec::new();
        for (shift, matrix) in self.system.matrices.iter().enumerate() {
            for (row, values) in matrix.iter().enumerate() {
                for (column, value) in values.iter().enumerate() {
                    let fraction = rational(value, &allowed)?;
                    // Normalize even powers before checking the denominator:
                    // 1/(r²+1) is rational, while 1/(r+1) needs rationalization.
                    let denominator = normalized_polynomial(
                        &fraction.denominator.to_expression(),
                        &root_atoms,
                        &self.roots,
                    )?;
                    if denominator.is_empty() {
                        return Err(Error::InvalidInput(
                            "denominator vanishes under the declared square-root relations".into(),
                        ));
                    }
                    if denominator.len() != 1 {
                        return Err(Error::Unsupported("sums involving square roots in denominators require algebraic rationalization, which is not supported by this regular-contour adapter".into()));
                    }
                    let (denominator_roots, denominator_coefficient) =
                        denominator.first_key_value().unwrap();
                    let mut terms = BTreeMap::<Vec<usize>, Atom>::new();
                    for (numerator_roots, numerator_coefficient) in normalized_polynomial(
                        &fraction.numerator.to_expression(),
                        &root_atoms,
                        &self.roots,
                    )? {
                        let mut coefficient = numerator_coefficient / denominator_coefficient;
                        let mut odd = Vec::new();
                        for (index, root) in self.roots.iter().enumerate() {
                            let numerator = numerator_roots.contains(&index);
                            let denominator = denominator_roots.contains(&index);
                            if numerator != denominator {
                                odd.push(index);
                            }
                            if denominator && !numerator {
                                coefficient /= encode_complex(&root.radicand);
                            }
                        }
                        let previous = terms.entry(odd).or_default();
                        *previous = (&*previous + coefficient).together().cancel();
                    }
                    for (roots, coefficient) in terms {
                        if !coefficient.is_zero() {
                            entries.push(vec![coefficient]);
                            normalized.push((shift, row, column, roots));
                        }
                    }
                }
            }
        }
        let term_count = entries.len();
        for root in &self.roots {
            let r = encode_complex(&root.radicand).together().cancel();
            entries.push(vec![r.clone()]);
            entries.push(vec![Atom::num(1) / &r]); // zeros are branch singularities
            entries.push(vec![
                (r.derivative(self.system.variable) / (&r * 2)).cancel(),
            ]);
        }
        // A zero differential system with no roots still needs a valid native
        // compilation row, while it retains no artificial coupling terms.
        if entries.is_empty() {
            entries.push(vec![Atom::new()]);
        }
        let compiled = compile_rows(self.system.variable, &entries, p, &Default::default())?;
        let terms = normalized
            .into_iter()
            .zip(&compiled.matrix)
            .map(|((shift, row, column, roots), values)| KernelTerm {
                shift,
                row,
                column,
                roots,
                coefficient: values[0].clone(),
            })
            .collect();
        let roots = self
            .roots
            .iter()
            .enumerate()
            .map(|(i, root)| RootKernel {
                definition: SquareRoot {
                    symbol: root.symbol,
                    radicand: decoded(&root.radicand),
                },
                value: compiled.matrix[term_count + 3 * i][0].clone(),
                logarithmic_derivative: compiled.matrix[term_count + 3 * i + 2][0].clone(),
            })
            .collect();
        Ok(CompiledAlgebraicSystem {
            variable: self.system.variable,
            p,
            size: self.system.matrices[0].len(),
            count: self.system.matrices.len(),
            roots,
            terms,
            poles: compiled.poles,
        })
    }
}

impl CompiledAlgebraicSystem {
    /// Includes rational poles and zeros/poles of every registered radicand.
    pub fn singularities(&self) -> &[C] {
        &self.poles
    }
    pub fn transport(
        &self,
        boundary: &EpsilonBoundary,
        waypoints: &[C],
        seeds: &BTreeMap<Symbol, RootSeed>,
        options: &FlowOptions,
        context: &RunContext,
        save_segments: bool,
    ) -> Result<AlgebraicSolution> {
        if boundary.coefficients.len() != self.count
            || boundary.coefficients.iter().any(|r| r.len() != self.size)
        {
            return Err(Error::InvalidInput(
                "algebraic epsilon boundary dimensions".into(),
            ));
        }
        boundary
            .leading
            .checked_add(self.count as i32 - 1)
            .ok_or_else(|| Error::InvalidInput("epsilon power range overflows".into()))?;
        let run = AlgebraicRun {
            compiled: self,
            seeds,
        };
        let flat = BoundaryData {
            point: boundary.point.clone(),
            values: boundary.coefficients.iter().flatten().cloned().collect(),
        };
        let mut segments = Vec::new();
        let (result, branches) = transport_series_with_state(
            &run,
            &flat,
            waypoints,
            options,
            context,
            save_segments.then_some(&mut segments),
        )?;
        Ok(AlgebraicSolution {
            solution: EpsilonSolution {
                point: result.point,
                leading: boundary.leading,
                coefficients: result.values.chunks(self.size).map(<[C]>::to_vec).collect(),
                diagnostics: result.diagnostics,
                segments,
                verified_digits: None,
                comparison_errors: Vec::new(),
                checkpoints: Vec::new(),
            },
            branches,
        })
    }
    pub fn transport_ordinary(
        &self,
        boundary: &BoundaryData,
        waypoints: &[C],
        seeds: &BTreeMap<Symbol, RootSeed>,
        options: &FlowOptions,
        context: &RunContext,
        save_segments: bool,
    ) -> Result<AlgebraicSolution> {
        if self.count != 1 {
            return Err(Error::InvalidInput(
                "ordinary transport requires one epsilon matrix".into(),
            ));
        }
        self.transport(
            &EpsilonBoundary {
                point: boundary.point.clone(),
                leading: 0,
                coefficients: vec![boundary.values.clone()],
            },
            waypoints,
            seeds,
            options,
            context,
            save_segments,
        )
    }
}

struct AlgebraicRun<'a> {
    compiled: &'a CompiledAlgebraicSystem,
    seeds: &'a BTreeMap<Symbol, RootSeed>,
}
struct RootChart {
    center: C,
    coefficients: Vec<Vec<C>>,
}
fn principal_sqrt(p: Precision, value: &C) -> Result<C> {
    if !p.finite(value) || *value == p.zero() {
        return Err(Error::Numerical(
            "square-root branch point or nonfinite radicand".into(),
        ));
    }
    // Exact real-axis choices avoid numerical imaginary remnants from exp(log).
    let value = if value.im == p.real(0) {
        if value.re > p.real(0) {
            C::new(value.re.clone().sqrt(), p.real(0))
        } else {
            C::new(p.real(0), (-value.re.clone()).sqrt())
        }
    } else {
        p.pow(value, &p.rational(&Rational::from((1, 2))))
    };
    if !p.finite(&value) {
        return Err(Error::Numerical("nonfinite square root".into()));
    }
    Ok(p.round(&value))
}
fn select_sign(p: Precision, principal: &C, hint: &C) -> Result<C> {
    if !p.finite(hint) || *hint == p.zero() {
        return Err(Error::InvalidInput(
            "nonfinite or zero square-root sign hint".into(),
        ));
    }
    let minus = p.neg(principal);
    let a = p.norm(&p.sub(principal, hint));
    let b = p.norm(&p.sub(&minus, hint));
    let (best, other, value) = if a < b {
        (a, b, principal.clone())
    } else {
        (b, a, minus)
    };
    if best * 4 >= other {
        return Err(Error::Accuracy("ambiguous square-root sign hint".into()));
    }
    Ok(value)
}
impl AlgebraicRun<'_> {
    fn project(&self, chart: &RootChart, point: &C) -> Result<Vec<C>> {
        let p = self.compiled.p;
        let step = p.sub(point, &chart.center);
        let (predicted, _) = evaluate_taylor(p, &chart.coefficients, &step);
        self.compiled
            .roots
            .iter()
            .zip(&predicted)
            .map(|(root, hint)| {
                if !p.finite(hint) || *hint == p.zero() {
                    return Err(Error::Accuracy(
                        "local root chart cannot select a sheet at this trial point".into(),
                    ));
                }
                let value = &root.value.series(p, point, 0)?[0];
                select_sign(p, &principal_sqrt(p, value)?, hint)
            })
            .collect()
    }
    fn roots_accurate(
        &self,
        chart: &RootChart,
        point: &C,
        tolerance: &Float,
    ) -> Result<Option<Vec<C>>> {
        let p = self.compiled.p;
        let step = p.sub(point, &chart.center);
        let (predicted, tails) = evaluate_taylor(p, &chart.coefficients, &step);
        let roots = match self.project(chart, point) {
            Ok(r) => r,
            Err(Error::Accuracy(_)) => return Ok(None),
            Err(e) => return Err(e),
        };
        for (i, root) in roots.iter().enumerate() {
            let scale = p.norm(root);
            if tails[i] > tolerance.clone() * &scale
                || p.norm(&p.sub(root, &predicted[i])) > tolerance.clone() * &scale
            {
                return Ok(None);
            }
            let mut derivative = p.zero();
            for k in (1..chart.coefficients.len()).rev() {
                derivative = p.add(
                    &p.mul(&derivative, &step),
                    &p.scale(&chart.coefficients[k][i], k as i64, 1),
                );
            }
            let rhs = p.mul(
                &self.compiled.roots[i]
                    .logarithmic_derivative
                    .series(p, point, 0)?[0],
                root,
            );
            let defect = p.mul(&step, &p.sub(&derivative, &rhs));
            let radicand = &self.compiled.roots[i].value.series(p, point, 0)?[0];
            let square_residual = p.sub(&p.mul(root, root), radicand);
            if !p.finite(&defect)
                || !p.finite(&square_residual)
                || p.norm(&defect) > tolerance.clone() * scale
                || p.norm(&square_residual) > tolerance.clone() * p.norm(radicand)
            {
                return Ok(None);
            }
        }
        Ok(Some(roots))
    }
}
impl SeriesSystem for AlgebraicRun<'_> {
    type State = BranchState;
    type Chart = RootChart;
    fn precision(&self) -> Precision {
        self.compiled.p
    }
    fn dimension(&self) -> usize {
        self.compiled.size * self.compiled.count
    }
    fn poles(&self) -> &[C] {
        &self.compiled.poles
    }
    fn initial_state(&self, boundary: &BoundaryData) -> Result<BranchState> {
        let p = self.compiled.p;
        if self.seeds.len() != self.compiled.roots.len()
            || self
                .compiled
                .roots
                .iter()
                .any(|r| !self.seeds.contains_key(&r.definition.symbol))
        {
            return Err(Error::InvalidInput(
                "provide exactly one seed per registered root".into(),
            ));
        }
        let roots = self
            .compiled
            .roots
            .iter()
            .map(|root| {
                let value = &root.value.series(p, &boundary.point, 0)?[0];
                let principal = principal_sqrt(p, value)?;
                let chosen = match &self.seeds[&root.definition.symbol] {
                    RootSeed::Principal => principal,
                    RootSeed::Opposite => p.neg(&principal),
                    RootSeed::Value(hint) => select_sign(p, &principal, hint)?,
                    RootSeed::I0(side) => {
                        if value.im != p.real(0) {
                            return Err(Error::InvalidInput(
                                "a radicand i0 seed requires a real starting radicand".into(),
                            ));
                        }
                        if value.re < p.real(0) && *side == Prescription::MinusI0 {
                            p.neg(&principal)
                        } else {
                            principal
                        }
                    }
                };
                Ok((root.definition.symbol, chosen))
            })
            .collect::<Result<_>>()?;
        Ok(BranchState {
            point: boundary.point.clone(),
            roots,
        })
    }
    fn local_chart(
        &self,
        center: &C,
        values: &[C],
        order: usize,
        state: &BranchState,
    ) -> Result<(Vec<Vec<C>>, RootChart)> {
        let c = self.compiled;
        let p = c.p;
        if state.point != *center || state.roots.len() != c.roots.len() {
            return Err(Error::Numerical(
                "algebraic branch state and expansion center differ".into(),
            ));
        }
        let center_exact = exact_point(p, center)?;
        let mut roots = Vec::new();
        for (index, root) in c.roots.iter().enumerate() {
            if state.roots[index].0 != root.definition.symbol {
                return Err(Error::Numerical("algebraic branch identity changed".into()));
            }
            let at_center = substitute(
                &root.definition.radicand,
                &BTreeMap::from([(Atom::var(c.variable), center_exact.clone())]),
            )
            .together()
            .cancel();
            if at_center.is_zero() {
                return Err(Error::Numerical(
                    "root chart starts at a branch singularity".into(),
                ));
            }
            let exact = (&root.definition.radicand / &at_center)
                .series(c.variable, &center_exact, order as i64 + 1)
                .map_err(|e| Error::Numerical(e.to_string()))?
                .rpow(Rational::from((1, 2)))
                .map_err(|e| Error::Numerical(e.to_string()))?;
            if exact.coefficient(Rational::from(0)) != Some(Atom::num(1)) {
                return Err(Error::Numerical(
                    "normalized root series does not start at exact one".into(),
                ));
            }
            let numeric = evaluate_series(&exact, p, &Default::default())?;
            let coefficients = coefficients(&numeric, order + 1)?
                .into_iter()
                .map(|v| p.mul(&v, &state.roots[index].1))
                .collect::<Vec<_>>();
            roots.push(fixed_series(p, c.variable, &coefficients)?);
        }
        let root_coefficients = (0..=order)
            .map(|k| {
                roots
                    .iter()
                    .map(|root| root.coefficient(Rational::from(k as i64)).unwrap())
                    .collect()
            })
            .collect();
        let mut kernels = Vec::new();
        for term in &c.terms {
            let mut kernel = fixed_series(
                p,
                c.variable,
                &term.coefficient.series(p, center, order - 1)?,
            )?;
            for &index in &term.roots {
                kernel = &kernel * &roots[index];
            }
            kernels.push(coefficients(&kernel, order)?);
        }
        let mut result = vec![vec![p.zero(); self.dimension()]; order + 1];
        result[0] = values.to_vec();
        for k in 0..order {
            for (term, kernel) in c.terms.iter().zip(&kernels) {
                for epsilon in term.shift..c.count {
                    let output = epsilon * c.size + term.row;
                    let input = (epsilon - term.shift) * c.size + term.column;
                    for (l, coefficient) in kernel.iter().take(k + 1).enumerate() {
                        result[k + 1][output] = p.add(
                            &result[k + 1][output],
                            &p.mul(coefficient, &result[k - l][input]),
                        );
                    }
                }
            }
            for value in &mut result[k + 1] {
                *value = p.scale(value, 1, k as i64 + 1);
            }
        }
        Ok((
            result,
            RootChart {
                center: center.clone(),
                coefficients: root_coefficients,
            },
        ))
    }
    fn rhs(&self, point: &C, values: &[C], chart: &RootChart) -> Result<Vec<C>> {
        let c = self.compiled;
        let p = c.p;
        let roots = self.project(chart, point)?;
        let mut result = vec![p.zero(); self.dimension()];
        for term in &c.terms {
            let mut coefficient = term.coefficient.series(p, point, 0)?[0].clone();
            for &root in &term.roots {
                coefficient = p.mul(&coefficient, &roots[root]);
            }
            for epsilon in term.shift..c.count {
                let output = epsilon * c.size + term.row;
                let input = (epsilon - term.shift) * c.size + term.column;
                result[output] = p.add(&result[output], &p.mul(&coefficient, &values[input]));
            }
        }
        Ok(result)
    }
    fn accepted_state(
        &self,
        chart: &RootChart,
        point: &C,
        tolerance: &Float,
    ) -> Result<Option<BranchState>> {
        let p = self.compiled.p;
        let midpoint = p.scale(&p.add(&chart.center, point), 1, 2);
        if self.roots_accurate(chart, &midpoint, tolerance)?.is_none() {
            return Ok(None);
        }
        Ok(self
            .roots_accurate(chart, point, tolerance)?
            .map(|roots| BranchState {
                point: point.clone(),
                roots: self
                    .compiled
                    .roots
                    .iter()
                    .zip(roots)
                    .map(|(root, value)| (root.definition.symbol, value))
                    .collect(),
            }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejected_chart_does_not_advance_branch_state() {
        let x = symbol!("branch_rejection::x");
        let r = symbol!("branch_rejection::r");
        let p = Precision::decimal(40).unwrap();
        let compiled = AlgebraicSystem::ordinary(
            DifferentialSystem {
                variable: x,
                matrix: vec![vec![Atom::new()]],
            },
            vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(x),
            }],
        )
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
        let run = AlgebraicRun {
            compiled: &compiled,
            seeds: &seeds,
        };
        let boundary = BoundaryData {
            point: p.i(1),
            values: vec![p.i(1)],
        };
        let state = run.initial_state(&boundary).unwrap();
        let (_, chart) = run
            .local_chart(&boundary.point, &boundary.values, 8, &state)
            .unwrap();
        assert!(
            run.accepted_state(&chart, &p.i(2), &p.tolerance(25))
                .unwrap()
                .is_none()
        );
        assert_eq!(state.point, p.i(1));
        assert_eq!(state.roots[0].1, p.i(1));
        let nearby = p.add(&p.i(1), &p.parse("1e-8", "0").unwrap());
        let next = run
            .accepted_state(&chart, &nearby, &p.tolerance(25))
            .unwrap()
            .unwrap();
        assert!(
            p.norm(&p.sub(&p.mul(&next.roots[0].1, &next.roots[0].1), &nearby)) < p.tolerance(35)
        );
        assert_eq!(state.point, p.i(1));
    }
    #[test]
    fn sparse_root_tail_requires_independent_defect_checks() {
        let x = symbol!("branch_sparse::x");
        let r = symbol!("branch_sparse::r");
        let p = Precision::decimal(40).unwrap();
        let compiled = AlgebraicSystem::ordinary(
            DifferentialSystem {
                variable: x,
                matrix: vec![vec![Atom::new()]],
            },
            vec![SquareRoot {
                symbol: r,
                radicand: Atom::num(1) + Atom::var(x).pow(20),
            }],
        )
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
        let run = AlgebraicRun {
            compiled: &compiled,
            seeds: &seeds,
        };
        let boundary = BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        };
        let state = run.initial_state(&boundary).unwrap();
        let (_, chart) = run
            .local_chart(&boundary.point, &boundary.values, 8, &state)
            .unwrap();
        assert!(chart.coefficients[1..].iter().all(|row| row[0] == p.zero()));
        assert!(
            run.accepted_state(&chart, &p.parse("0.2", "0").unwrap(), &p.tolerance(25))
                .unwrap()
                .is_none()
        );
    }
}
