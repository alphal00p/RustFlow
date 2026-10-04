//! Regular-contour transport with explicitly registered square-root sheets.
//!
//! Roots retain their named generators and independent signs. Formal denominator
//! inversion uses native Symbolica quotient arithmetic. Source domains survive
//! path and kernel cancellation; their norm conditions can conservatively exclude
//! points regular on a particular sheet. Only accepted continuation steps advance
//! branch state. Singular endpoints require a separate Frobenius treatment.
mod canonical;
mod quotient;
use crate::diffexp::{EpsilonBoundary, EpsilonSolution, EpsilonSystem};
use crate::family::{encode_complex, imaginary_parameter, scalar_symbols, substitute};
use crate::fixed_series::{coefficients, fixed_series};
use crate::ode::{
    CompiledSystem, NumericRational, SeriesSystem, compile_rows, evaluate_taylor,
    transport_series_with_state,
};
use crate::{
    BoundaryData, ComplexFloat as C, DifferentialSystem, Error, FlowOptions, Precision,
    Prescription, Result, RunContext,
};
pub use canonical::CanonicalAlgebraicSystem;
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
    /// Exact regular-domain conditions retained before path or kernel cancellation.
    pub nonzero_conditions: Vec<Atom>,
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

/// Exact differential matrices in physical invariants, with explicitly declared
/// square-root generators. Root radicands depend on physical invariants only;
/// epsilon remains a separate symbol until the exact path is pulled back.
#[derive(Clone, Debug)]
pub struct AlgebraicKinematicSystem {
    pub epsilon: Symbol,
    pub derivatives: BTreeMap<Symbol, Vec<Vec<Atom>>>,
    pub roots: Vec<SquareRoot>,
}

impl AlgebraicKinematicSystem {
    /// Exact physical-domain conditions before any path/Jacobian cancellation.
    /// Formal root norms can exclude a regular point on one particular sheet;
    /// this interface conservatively requires invertibility on every sheet.
    pub fn nonzero_conditions(&self) -> Result<Vec<Atom>> {
        self.validate()?;
        registered_domain_conditions(
            &self
                .derivatives
                .values()
                .flatten()
                .flatten()
                .cloned()
                .collect::<Vec<_>>(),
            &self.roots,
            &self
                .derivatives
                .keys()
                .copied()
                .chain([self.epsilon])
                .collect(),
        )
    }

    pub fn validate(&self) -> Result<()> {
        let Some((&first_variable, first)) = self.derivatives.first_key_value() else {
            return Err(Error::InvalidInput(
                "an algebraic kinematic system needs physical derivative matrices".into(),
            ));
        };
        DifferentialSystem {
            variable: first_variable,
            matrix: first.clone(),
        }
        .validate()?;
        let n = first.len();
        if self.epsilon == imaginary_parameter()
            || self.derivatives.contains_key(&self.epsilon)
            || self.derivatives.contains_key(&imaginary_parameter())
        {
            return Err(Error::InvalidInput(
                "physical variables must differ from epsilon and the reserved imaginary unit"
                    .into(),
            ));
        }
        let mut physical = self
            .derivatives
            .keys()
            .map(|&v| Atom::var(v))
            .collect::<BTreeSet<_>>();
        physical.insert(Atom::var(imaginary_parameter()));
        let mut allowed = physical.clone();
        allowed.insert(Atom::var(self.epsilon));
        for root in &self.roots {
            if !allowed.insert(Atom::var(root.symbol)) {
                return Err(Error::InvalidInput(
                    "root symbols must be distinct from epsilon and physical variables".into(),
                ));
            }
            rational(&root.radicand, &physical)?;
            if decoded(&root.radicand).is_zero() {
                return Err(Error::InvalidInput(
                    "a square-root radicand is identically zero".into(),
                ));
            }
        }
        for matrix in self.derivatives.values() {
            if matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
                return Err(Error::InvalidInput(
                    "algebraic kinematic matrices must be square with equal dimensions".into(),
                ));
            }
            for a in matrix.iter().flatten() {
                rational(a, &allowed)?;
            }
        }
        Ok(())
    }

    /// Construct dY = epsilon * sum M_letter dlog(letter) Y. Matrices are exact
    /// constants. Native differentiation includes dr/dv = (dR/dv)/(2r), so
    /// registered roots are not treated as constants under the chain rule.
    pub fn canonical_dlog(
        epsilon: Symbol,
        variables: &[Symbol],
        letters: &[Atom],
        matrices: &[Vec<Vec<Atom>>],
        roots: Vec<SquareRoot>,
    ) -> Result<Self> {
        if variables.is_empty()
            || variables.iter().copied().collect::<BTreeSet<_>>().len() != variables.len()
            || letters.is_empty()
            || letters.len() != matrices.len()
        {
            return Err(Error::InvalidInput("canonical algebraic data needs distinct variables and equal nonempty letters/matrices".into()));
        }
        let n = matrices[0].len();
        let mut result = Self {
            epsilon,
            derivatives: variables
                .iter()
                .map(|&v| (v, vec![vec![Atom::new(); n]; n]))
                .collect(),
            roots,
        };
        result.validate()?;
        let constants = BTreeSet::from([Atom::var(imaginary_parameter())]);
        let mut allowed = variables
            .iter()
            .map(|&v| Atom::var(v))
            .collect::<BTreeSet<_>>();
        allowed.extend(result.roots.iter().map(|r| Atom::var(r.symbol)));
        allowed.insert(Atom::var(imaginary_parameter()));
        let root_atoms = result
            .roots
            .iter()
            .map(|r| Atom::var(r.symbol))
            .collect::<Vec<_>>();
        for (letter, matrix) in letters.iter().zip(matrices) {
            if matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
                return Err(Error::InvalidInput(
                    "canonical algebraic matrices must have identical square dimensions".into(),
                ));
            }
            for a in matrix.iter().flatten() {
                rational(a, &constants)?;
            }
            let fraction = rational(letter, &allowed)?;
            if normalized_polynomial(
                &fraction.numerator.to_expression(),
                &root_atoms,
                &result.roots,
            )?
            .is_empty()
                || normalized_polynomial(
                    &fraction.denominator.to_expression(),
                    &root_atoms,
                    &result.roots,
                )?
                .is_empty()
            {
                return Err(Error::InvalidInput("an algebraic dlog letter cannot vanish or have zero denominator under its root relations".into()));
            }
            for (&variable, output) in &mut result.derivatives {
                let mut derivative = letter.derivative(variable);
                for root in &result.roots {
                    derivative += letter.derivative(root.symbol)
                        * root.radicand.derivative(variable)
                        / (Atom::var(root.symbol) * 2);
                }
                let dlog = (derivative / letter).together().cancel();
                for (row, source) in output.iter_mut().zip(matrix) {
                    for (entry, coefficient) in row.iter_mut().zip(source) {
                        *entry += Atom::var(epsilon) * coefficient * &dlog;
                    }
                }
            }
        }
        for entry in result.derivatives.values_mut().flatten().flatten() {
            *entry = entry.together().cancel();
        }
        result.validate()?;
        Ok(result)
    }

    /// Pull the physical system back to one exact rational path and retain
    /// epsilon coefficients through `order`. Coordinate substitution is
    /// simultaneous, including in every root definition. Endpoints and sheet
    /// seeds remain explicit inputs to the compiled regular-contour transport.
    pub fn pullback(
        &self,
        path: &crate::kinematics::KinematicPath,
        order: usize,
    ) -> Result<AlgebraicSystem> {
        self.validate()?;
        path.validate()?;
        if path.parameter == self.epsilon
            || path.parameter == imaginary_parameter()
            || self.derivatives.contains_key(&path.parameter)
            || self.roots.iter().any(|r| r.symbol == path.parameter)
            || !self.derivatives.keys().eq(path.coordinates.keys())
        {
            return Err(Error::InvalidInput("algebraic path must specify exactly the physical variables and use a distinct parameter".into()));
        }
        let rules = path
            .coordinates
            .iter()
            .map(|(&symbol, value)| (Atom::var(symbol), value.clone()))
            .collect::<BTreeMap<_, _>>();
        let original_symbols = self
            .derivatives
            .keys()
            .copied()
            .chain([self.epsilon, imaginary_parameter()])
            .chain(self.roots.iter().map(|r| r.symbol))
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        let restricted_symbols = [path.parameter, self.epsilon, imaginary_parameter()]
            .into_iter()
            .chain(self.roots.iter().map(|r| r.symbol))
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        let path_symbols =
            BTreeSet::from([Atom::var(path.parameter), Atom::var(imaginary_parameter())]);
        for coordinate in path.coordinates.values() {
            rational(coordinate, &path_symbols)?;
        }
        // Restrict numerator and denominator separately, before a zero
        // Jacobian or another term can erase an undefined source entry.
        let restrict = |expression: &Atom| -> Result<Atom> {
            let fraction = rational(expression, &original_symbols)?;
            let denominator = substitute(&fraction.denominator.to_expression(), &rules)
                .together()
                .cancel();
            rational(&denominator, &restricted_symbols)?;
            if decoded(&denominator).is_zero() {
                return Err(Error::InvalidInput(
                    "algebraic path lies identically on a source denominator pole".into(),
                ));
            }
            let numerator = substitute(&fraction.numerator.to_expression(), &rules)
                .together()
                .cancel();
            rational(&numerator, &restricted_symbols)?;
            let restricted = (numerator / denominator).together().cancel();
            rational(&restricted, &restricted_symbols)?;
            Ok(restricted)
        };
        // Preserve every source domain before the Jacobian or another row
        // cancels it. Regular epsilon transport must also retain the generic
        // epsilon valuation at this physical point.
        let mut nonzero_conditions = crate::physical_conditions::rational_denominator_conditions(
            &path.coordinates.values().cloned().collect::<Vec<_>>(),
            &BTreeSet::from([path.parameter, imaginary_parameter()]),
        )?;
        for guard in self.nonzero_conditions()? {
            // Select the generic epsilon valuation before physical restriction.
            // Otherwise an exceptional path can silently remove its leading term.
            let fraction = rational(&guard, &original_symbols)?;
            for part in [
                fraction.numerator.to_expression(),
                fraction.denominator.to_expression(),
            ] {
                let coefficient = crate::physical_conditions::epsilon_leading_coefficient(
                    &decoded(&part),
                    self.epsilon,
                )?;
                let coefficient = restrict(&coefficient)?;
                if decoded(&coefficient).is_zero() {
                    return Err(Error::InvalidInput(
                        "path changes the generic epsilon valuation of a source domain condition"
                            .into(),
                    ));
                }
                if !matches!(coefficient.as_view(), AtomView::Num(_)) {
                    nonzero_conditions.push(coefficient);
                }
            }
        }
        nonzero_conditions.sort();
        nonzero_conditions.dedup();
        let roots = self
            .roots
            .iter()
            .map(|root| {
                let radicand = restrict(&root.radicand)?;
                if decoded(&radicand).is_zero() {
                    return Err(Error::InvalidInput(
                        "algebraic path lies identically on a root branch locus".into(),
                    ));
                }
                Ok(SquareRoot {
                    symbol: root.symbol,
                    radicand,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let root_atoms = roots
            .iter()
            .map(|r| Atom::var(r.symbol))
            .collect::<Vec<_>>();
        let n = self.derivatives.first_key_value().unwrap().1.len();
        let mut matrix = vec![vec![Atom::new(); n]; n];
        for (variable, source) in &self.derivatives {
            let jacobian = path.coordinates[variable].derivative(path.parameter);
            rational(&jacobian, &path_symbols)?;
            for (row, source) in matrix.iter_mut().zip(source) {
                for (entry, coefficient) in row.iter_mut().zip(source) {
                    let restricted = restrict(coefficient)?;
                    let fraction = rational(&restricted, &restricted_symbols)?;
                    let denominator = normalized_polynomial(
                        &fraction.denominator.to_expression(),
                        &root_atoms,
                        &roots,
                    )?;
                    if denominator.is_empty() {
                        return Err(Error::InvalidInput("algebraic path lies identically on a matrix pole under its root relations".into()));
                    }
                    if denominator.len() != 1 {
                        quotient::invert_denominator(
                            &quotient::recompose(&denominator, &roots),
                            &roots,
                            &restricted_symbols,
                        )?;
                    }
                    *entry += &jacobian * restricted;
                }
            }
        }
        for entry in matrix.iter_mut().flatten() {
            *entry = entry.together().cancel();
        }
        let result = AlgebraicSystem {
            system: EpsilonSystem::from_differential_system(
                &DifferentialSystem {
                    variable: path.parameter,
                    matrix,
                },
                self.epsilon,
                order,
            )?,
            roots,
            nonzero_conditions,
        };
        result.validate()?;
        Ok(result)
    }
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
    root_system: Option<CompiledSystem>,
    terms: Vec<KernelTerm>,
    poles: Vec<C>,
    pole_polynomials: Vec<Atom>,
    domain_guards: Vec<Atom>,
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
    let fraction: RationalPolynomial<IntegerRing, u16> = encode_complex(a)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| {
            Error::Unsupported(format!(
                "rational functions of the declared root generators are required: {e}"
            ))
        })?;
    if fraction.numerator.variables().iter().any(|variable| {
        !matches!(variable, PolyVariable::Symbol(symbol) if allowed.contains(&Atom::var(*symbol)))
    }) {
        return Err(Error::Unsupported(
            "nonrational powers require an explicitly registered supported root".into(),
        ));
    }
    Ok(fraction)
}
fn normalized_polynomial(
    expression: &Atom,
    root_atoms: &[Atom],
    roots: &[SquareRoot],
) -> Result<BTreeMap<Vec<usize>, Atom>> {
    let mut terms = BTreeMap::<Vec<usize>, Atom>::new();
    // Retain native exact rational coefficients. The general AtomField view
    // performs statistical zero tests on large symbolic coefficients; a native
    // polynomial-of-polynomials groups the same root powers using exact Q.
    let polynomial: MultivariatePolynomial<_, i32> = encode_complex(expression)
        .try_to_polynomial(&Q, root_atoms)
        .map_err(|e| Error::Unsupported(format!("exact root polynomial required: {e}")))?;
    let root_indices = (0..root_atoms.len()).collect::<Vec<_>>();
    let grouped = polynomial.to_polynomial_in(&root_indices);
    for term in grouped.into_iter() {
        let mut coefficient = term.coefficient.to_expression();
        let mut odd = Vec::new();
        for (index, &exponent) in term.exponents.iter().enumerate() {
            coefficient *=
                encode_complex(&roots[index].radicand).pow(i64::from(exponent.div_euclid(2)));
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

/// Preserve raw denominator occurrences before native rational cancellation.
/// Named-root inversion then removes root generators from these exact domains.
fn registered_domain_conditions(
    expressions: &[Atom],
    roots: &[SquareRoot],
    physical: &BTreeSet<Symbol>,
) -> Result<Vec<Atom>> {
    let mut variables = physical.clone();
    variables.insert(imaginary_parameter());
    let mut conditions = crate::physical_conditions::canonical_conditions(
        &roots.iter().map(|r| r.radicand.clone()).collect::<Vec<_>>(),
        &variables,
    )?;
    variables.extend(roots.iter().map(|r| r.symbol));
    let allowed = variables.iter().map(|&s| Atom::var(s)).collect();
    let root_atoms = roots
        .iter()
        .map(|r| Atom::var(r.symbol))
        .collect::<Vec<_>>();
    for condition in
        crate::physical_conditions::rational_denominator_conditions(expressions, &variables)?
    {
        let fraction = rational(&condition, &allowed)?;
        let normalized =
            normalized_polynomial(&fraction.numerator.to_expression(), &root_atoms, roots)?;
        if normalized.is_empty() {
            return Err(Error::InvalidInput(
                "source denominator vanishes under root relations".into(),
            ));
        }
        for coefficient in
            quotient::invert_denominator(&quotient::recompose(&normalized, roots), roots, &allowed)?
                .values()
        {
            conditions.push(rational(coefficient, &allowed)?.denominator.to_expression());
        }
    }
    conditions = conditions.into_iter().map(|a| decoded(&a)).collect();
    conditions.sort();
    conditions.dedup();
    conditions.retain(|a| !matches!(a.as_view(), AtomView::Num(_)));
    Ok(conditions)
}

impl AlgebraicSystem {
    pub fn ordinary(system: DifferentialSystem, roots: Vec<SquareRoot>) -> Self {
        Self {
            system: EpsilonSystem {
                variable: system.variable,
                matrices: vec![system.matrix],
            },
            roots,
            nonzero_conditions: Vec::new(),
        }
    }

    pub fn validate(&self) -> Result<()> {
        self.system.validate()?;
        if self.system.variable == imaginary_parameter() {
            return Err(Error::InvalidInput(
                "path parameter collides with the reserved imaginary unit".into(),
            ));
        }
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
        for condition in &self.nonzero_conditions {
            rational(condition, &allowed)?;
            if decoded(condition).is_zero() {
                return Err(Error::InvalidInput(
                    "algebraic domain contains an identically zero condition".into(),
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
        let mut domain_entries = Vec::new();
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
                        let inverse = quotient::invert_denominator(
                            &quotient::recompose(&denominator, &self.roots),
                            &self.roots,
                            &allowed,
                        )?;
                        domain_entries.extend(
                            inverse
                                .values()
                                .map(|coefficient| vec![coefficient.clone()]),
                        );
                        let numerator = normalized_polynomial(
                            &fraction.numerator.to_expression(),
                            &root_atoms,
                            &self.roots,
                        )?;
                        let product = quotient::recompose(&numerator, &self.roots)
                            * quotient::recompose(&inverse, &self.roots);
                        for (roots, coefficient) in
                            normalized_polynomial(&product, &root_atoms, &self.roots)?
                        {
                            if !coefficient.is_zero() {
                                entries.push(vec![coefficient]);
                                normalized.push((shift, row, column, roots));
                            }
                        }
                        continue;
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
        let mut root_matrix = vec![vec![Atom::new(); self.roots.len()]; self.roots.len()];
        for (index, root) in self.roots.iter().enumerate() {
            let r = encode_complex(&root.radicand).together().cancel();
            entries.push(vec![r.clone()]);
            entries.push(vec![Atom::num(1) / &r]); // zeros are branch singularities
            let logarithmic_derivative = (r.derivative(self.system.variable) / (&r * 2)).cancel();
            entries.push(vec![logarithmic_derivative.clone()]);
            root_matrix[index][index] = logarithmic_derivative;
        }
        // Retain norm poles even when multiplication by a numerator cancels them.
        let mut domain_guards = domain_entries
            .iter()
            .map(|row| rational(&row[0], &allowed).map(|r| r.denominator.to_expression()))
            .collect::<Result<Vec<_>>>()?;
        let mut source_conditions = registered_domain_conditions(
            &self
                .system
                .matrices
                .iter()
                .flatten()
                .flatten()
                .cloned()
                .collect::<Vec<_>>(),
            &self.roots,
            &BTreeSet::from([self.system.variable]),
        )?;
        source_conditions.extend(self.nonzero_conditions.iter().cloned());
        for condition in &source_conditions {
            let condition = rational(condition, &allowed)?;
            for part in [
                condition.numerator.to_expression(),
                condition.denominator.to_expression(),
            ] {
                domain_entries.push(vec![Atom::one() / &part]);
                domain_guards.push(part);
            }
        }
        domain_guards.sort();
        domain_guards.dedup();
        domain_guards.retain(|a| !matches!(a.as_view(), AtomView::Num(_)));
        entries.extend(domain_entries);
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
        let root_system = if root_matrix.is_empty() {
            None
        } else {
            Some(compile_rows(
                self.system.variable,
                &root_matrix,
                p,
                &Default::default(),
            )?)
        };
        Ok(CompiledAlgebraicSystem {
            variable: self.system.variable,
            p,
            size: self.system.matrices[0].len(),
            count: self.system.matrices.len(),
            roots,
            root_system,
            terms,
            poles: compiled.poles,
            pole_polynomials: compiled.pole_polynomials,
            domain_guards,
        })
    }
}

impl CompiledAlgebraicSystem {
    fn check_domain(&self, point: &C) -> Result<()> {
        if self.domain_guards.is_empty() {
            return Ok(());
        }
        let coordinate = exact_point(self.p, point)?;
        let rules = BTreeMap::from([(Atom::var(self.variable), coordinate)]);
        if self
            .domain_guards
            .iter()
            .any(|a| decoded(&substitute(a, &rules)).is_zero())
        {
            return Err(Error::InvalidInput(
                "point lies on a retained root-denominator norm pole".into(),
            ));
        }
        Ok(())
    }
    /// Includes rational poles and zeros/poles of every registered radicand.
    pub fn singularities(&self) -> &[C] {
        &self.poles
    }
    /// Exact denominator factors, including zeros and poles of radicands.
    /// These retain exact root identities for prescribed-contour planning;
    /// nearby numerical root estimates are insufficient to identify aliases.
    pub fn singularity_polynomials(&self) -> &[Atom] {
        &self.pole_polynomials
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
        self.compiled.check_domain(&boundary.point)?;
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
                    RootSeed::Value(hint) => {
                        if !p.finite(hint) || *hint == p.zero() {
                            return Err(Error::InvalidInput(
                                "nonfinite or zero square-root sign hint".into(),
                            ));
                        }
                        let magnitude = p.norm(hint);
                        if !magnitude.is_finite() || magnitude == p.real(0) {
                            return Err(Error::InvalidInput(
                                "invalid square-root sign hint magnitude".into(),
                            ));
                        }
                        let normalized = p.mul(
                            hint,
                            &p.div(
                                &C::new(p.norm(&principal), p.real(0)),
                                &C::new(magnitude, p.real(0)),
                            ),
                        );
                        select_sign(p, &principal, &normalized)?
                    }
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
        self.compiled.check_domain(center)?;
        let c = self.compiled;
        let p = c.p;
        if state.point != *center || state.roots.len() != c.roots.len() {
            return Err(Error::Numerical(
                "algebraic branch state and expansion center differ".into(),
            ));
        }
        let center_exact = exact_point(p, center)?;
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
            // Preserve exact-center domain validation independently of the
            // numeric logarithmic-derivative recurrence.
            p.eval(&at_center, &Default::default())?;
        }
        // A square root satisfies a rational scalar differential equation:
        // r'=(R'/2R)r. Reuse the already-tested finite-polynomial recurrence;
        // accepted endpoint sheets still supply its initial values.
        let root_coefficients = if let Some(system) = &c.root_system {
            system.taylor(
                center,
                &state.roots.iter().map(|r| r.1.clone()).collect::<Vec<_>>(),
                order,
            )?
        } else {
            vec![Vec::new(); order + 1]
        };
        let roots = (0..c.roots.len())
            .map(|index| {
                fixed_series(
                    p,
                    c.variable,
                    &root_coefficients
                        .iter()
                        .map(|row| row[index].clone())
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Result<Vec<_>>>()?;
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
        self.compiled.check_domain(point)?;
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

#[cfg(test)]
mod root_taylor_comparisons {
    use super::*;

    #[test]
    fn rational_root_ode_matches_native_exact_rpow_on_difficult_charts() {
        let p = Precision::decimal(90).unwrap();
        let x = symbol!("root_chart_comparison::x");
        let r = symbol!("root_chart_comparison::r");
        let variable = Atom::var(x);
        let near = &variable - 1;
        // The longer optional cases reproduce the retained chart comparison
        // report; ordinary regression runs keep the costly exact oracle short.
        let extended = std::env::var_os("RUSTFLOW_EXTENDED_ROOT_CHART_CHECK").is_some();
        let cases = [
            (
                variable.clone().pow(4) + variable.clone().pow(3) * 7 - &variable * 11 + 5,
                p.zero(),
                80,
            ),
            (
                (near.clone().pow(3) / (variable.clone() + 2).pow(2)).cancel(),
                p.parse("1.000001", "0.0000001").unwrap(),
                if extended { 24 } else { 16 },
            ),
            (
                near.pow(2) + (Atom::num(1) / Atom::num(Integer::from(10).pow(30))),
                p.parse("1.000000000000001", "0.0000000000000001").unwrap(),
                if extended { 20 } else { 16 },
            ),
            (
                Atom::num(Complex::new(Rational::from(1), Rational::from(1)))
                    * (variable.clone().pow(2) + 2)
                    / (variable.pow(2) + 3),
                p.parse(
                    "0.123456789123456789123456789",
                    "-0.078901234567890123456789",
                )
                .unwrap(),
                if extended { 24 } else { 12 },
            ),
        ];
        for (index, (radicand, center, order)) in cases.into_iter().enumerate() {
            let compiled = AlgebraicSystem::ordinary(
                DifferentialSystem {
                    variable: x,
                    matrix: vec![vec![Atom::num(1) / Atom::var(r)]],
                },
                vec![SquareRoot {
                    symbol: r,
                    radicand,
                }],
            )
            .compile(p)
            .unwrap();
            let seeds = BTreeMap::from([(r, RootSeed::Opposite)]);
            let run = AlgebraicRun {
                compiled: &compiled,
                seeds: &seeds,
            };
            let boundary = BoundaryData {
                point: center.clone(),
                values: vec![p.i(1)],
            };
            let state = run.initial_state(&boundary).unwrap();
            let started = std::time::Instant::now();
            let actual = compiled
                .root_system
                .as_ref()
                .unwrap()
                .taylor(&center, &[state.roots[0].1.clone()], order)
                .unwrap();
            let numeric_time = started.elapsed();
            let center_exact = exact_point(p, &center).unwrap();
            let definition = &compiled.roots[0].definition.radicand;
            let at_center = substitute(
                definition,
                &BTreeMap::from([(Atom::var(x), center_exact.clone())]),
            )
            .together()
            .cancel();
            let started = std::time::Instant::now();
            let exact = (definition / &at_center)
                .series(x, &center_exact, order as i64 + 1)
                .unwrap()
                .rpow(Rational::from((1, 2)))
                .unwrap();
            assert_eq!(exact.coefficient(Rational::from(0)), Some(Atom::num(1)));
            let series =
                crate::fixed_series::evaluate_series(&exact, p, &Default::default()).unwrap();
            let expected = coefficients(&series, order + 1).unwrap();
            let exact_time = started.elapsed();
            for (k, expected) in expected.iter().enumerate() {
                let expected = p.mul(expected, &state.roots[0].1);
                let magnitude = p.norm(&expected);
                let scale = if magnitude > p.real(1) {
                    magnitude
                } else {
                    p.real(1)
                };
                let error = p.norm(&p.sub(&actual[k][0], &expected));
                assert!(
                    error < p.tolerance(50) * scale,
                    "case{index} coefficient{k}: {error}"
                );
            }
            println!("case={index} order={order} numeric={numeric_time:?} exact={exact_time:?}");
        }
    }
}
impl CompiledAlgebraicSystem {
    /// Recompute root magnitudes from the exact system at this working precision;
    /// seed values specify only discrete sheets and supply no accuracy evidence.
    pub fn branch_state_at(
        &self,
        point: &C,
        seeds: &BTreeMap<Symbol, RootSeed>,
    ) -> Result<BranchState> {
        AlgebraicRun {
            compiled: self,
            seeds,
        }
        .initial_state(&BoundaryData {
            point: point.clone(),
            values: Vec::new(),
        })
    }

    /// Weighted Gronwall estimate with rational disk bounds and |sqrt(R)| =
    /// sqrt(|R|). It bounds inherited boundary uncertainty, not arithmetic error;
    /// independently refined transport checks supply the latter evidence.
    pub fn error_amplification_weighted(
        &self,
        start: &C,
        end: &C,
        weights: &[Float],
    ) -> Result<Float> {
        let p = self.p;
        if weights.len() != self.size * self.count
            || weights.iter().any(|w| !w.is_finite() || *w <= p.real(0))
            || !p.finite(start)
            || !p.finite(end)
        {
            return Err(Error::InvalidInput("algebraic error amplification requires finite endpoints and positive coefficient weights".into()));
        }
        crate::diffexp::integrate_error_amplification(p, start, end, |center, radius| {
            let mut roots = Vec::new();
            for root in &self.roots {
                let Some(upper) =
                    crate::diffexp::rational_disk_bound(p, &root.value, center, radius)?
                else {
                    return Ok(None);
                };
                // Separate every root from zero on each disk as well as from
                // its poles; no branch crossing can be hidden by a finite upper bound.
                let inverse = NumericRational {
                    numerator: root.value.denominator.clone(),
                    denominator: root.value.numerator.clone(),
                };
                if crate::diffexp::rational_disk_bound(p, &inverse, center, radius)?.is_none() {
                    return Ok(None);
                }
                roots.push(upper.sqrt());
            }
            let mut rows = vec![p.real(0); self.size * self.count];
            for term in &self.terms {
                let Some(mut upper) =
                    crate::diffexp::rational_disk_bound(p, &term.coefficient, center, radius)?
                else {
                    return Ok(None);
                };
                for &root in &term.roots {
                    upper *= &roots[root];
                }
                for order in term.shift..self.count {
                    rows[order * self.size + term.row] += upper.clone()
                        * &weights[(order - term.shift) * self.size + term.column]
                        / &weights[order * self.size + term.row];
                }
            }
            Ok(Some(
                rows.into_iter()
                    .fold(p.real(0), |a, b| if a > b { a } else { b }),
            ))
        })
    }
}
