//! Exact multivariate differential equations restricted to kinematic paths.
//!
//! The pullback is the chain rule `A_x = sum_s (ds/dx) A_s(s(x), epsilon)`.
//! Algebraic powers are retained symbolically; numerical support for a pulled-back
//! matrix is determined separately by the differential-equation solver.
use crate::{DifferentialSystem, Error, Result, family::substitute};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;

/// Partial-derivative matrices `dY/ds = A_s Y` in a common ordered basis.
///
/// Every physical variable appearing in an entry must have a matrix, including
/// variables whose derivative matrix is zero. The regulator is independent of
/// these variables. Validation checks the representation, not flatness of the
/// multivariate connection.
#[derive(Clone, Debug)]
pub struct KinematicSystem {
    pub epsilon: Symbol,
    pub derivatives: BTreeMap<Symbol, Vec<Vec<Atom>>>,
}

/// An exact coordinate map from one parameter into kinematic space.
///
/// Coordinates may be rational or algebraic expressions in `parameter` and
/// exact rational-complex constants. They cannot contain undeclared parameters
/// or refer to another coordinate by name.
#[derive(Clone, Debug)]
pub struct KinematicPath {
    pub parameter: Symbol,
    pub coordinates: BTreeMap<Symbol, Atom>,
}

fn exact_expression(a: AtomView<'_>, symbols: &BTreeSet<Symbol>) -> Result<()> {
    match a {
        AtomView::Num(n) => match n.get_coeff_view().to_owned() {
            Coefficient::Complex(_) => Ok(()),
            _ => Err(Error::InvalidInput(
                "kinematic expressions require finite exact rational-complex numbers".into(),
            )),
        },
        AtomView::Var(v) if symbols.contains(&v.get_symbol()) => Ok(()),
        AtomView::Var(_) => Err(Error::InvalidInput(format!(
            "undeclared kinematic variable {a}"
        ))),
        AtomView::Add(v) => v.iter().try_for_each(|a| exact_expression(a, symbols)),
        AtomView::Mul(v) => v.iter().try_for_each(|a| exact_expression(a, symbols)),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            exact_expression(base, symbols)?;
            if let AtomView::Num(n) = exponent
                && let Coefficient::Complex(value) = n.get_coeff_view().to_owned()
                && value.im.is_zero()
            {
                return Ok(());
            }
            Err(Error::Unsupported(
                "kinematic powers require exact rational exponents".into(),
            ))
        }
        _ => Err(Error::Unsupported(
            "kinematic expressions support arithmetic and algebraic powers, not function calls"
                .into(),
        )),
    }
}

fn square_dimension(matrix: &[Vec<Atom>]) -> Result<usize> {
    let n = matrix.len();
    if n == 0 || matrix.iter().any(|row| row.len() != n) {
        return Err(Error::InvalidInput(
            "kinematic derivative matrices must be nonempty and square".into(),
        ));
    }
    Ok(n)
}

impl KinematicSystem {
    pub fn validate(&self) -> Result<()> {
        let Some(first) = self.derivatives.values().next() else {
            return Err(Error::InvalidInput(
                "no kinematic derivative matrices".into(),
            ));
        };
        if self.derivatives.contains_key(&self.epsilon) {
            return Err(Error::InvalidInput(
                "epsilon cannot also be a physical variable".into(),
            ));
        }
        let n = square_dimension(first)?;
        let mut symbols: BTreeSet<_> = self.derivatives.keys().copied().collect();
        symbols.insert(self.epsilon);
        for matrix in self.derivatives.values() {
            if square_dimension(matrix)? != n {
                return Err(Error::InvalidInput(
                    "kinematic derivative matrices have different dimensions".into(),
                ));
            }
            for entry in matrix.iter().flatten() {
                exact_expression(entry.as_view(), &symbols)?;
            }
        }
        Ok(())
    }

    /// Restrict the system using simultaneous substitution and the exact Jacobian.
    /// A path contained in a pole of a matrix entry is rejected, including when
    /// that coordinate is stationary; a limiting prescription is then required.
    pub fn pullback(&self, path: &KinematicPath) -> Result<DifferentialSystem> {
        self.validate()?;
        path.validate()?;
        if path.parameter == self.epsilon || path.coordinates.contains_key(&self.epsilon) {
            return Err(Error::InvalidInput(
                "the path parameter and physical coordinates must differ from epsilon".into(),
            ));
        }
        if !self.derivatives.keys().eq(path.coordinates.keys()) {
            return Err(Error::InvalidInput(
                "path coordinates must match all physical derivative variables exactly".into(),
            ));
        }
        let substitutions = path
            .coordinates
            .iter()
            .map(|(&symbol, value)| (Atom::var(symbol), value.clone()))
            .collect();
        let n = self.derivatives.values().next().unwrap().len();
        let mut matrix = vec![vec![Atom::new(); n]; n];
        let symbols = BTreeSet::from([path.parameter, self.epsilon]);
        for (variable, derivative) in &self.derivatives {
            let jacobian = path.coordinates[variable].derivative(path.parameter);
            exact_expression(jacobian.as_view(), &symbols)?;
            for (row, partial_row) in matrix.iter_mut().zip(derivative) {
                for (entry, partial) in row.iter_mut().zip(partial_row) {
                    let restricted = substitute(partial, &substitutions);
                    exact_expression(restricted.as_view(), &symbols)?;
                    *entry += &jacobian * restricted;
                }
            }
        }
        for entry in matrix.iter_mut().flatten() {
            *entry = entry.together().cancel();
            exact_expression(entry.as_view(), &symbols)?;
        }
        Ok(DifferentialSystem {
            variable: path.parameter,
            matrix,
        })
    }

    /// Construct the canonical connection `dY = epsilon sum_k C_k dlog(l_k) Y`.
    ///
    /// `matrices[k]` is a constant exact matrix multiplying `letters[k]`.
    /// No logarithms are introduced: each derivative is formed as `(dl_k/ds)/l_k`.
    /// Branch prescriptions for algebraic letters belong to numerical transport.
    pub fn canonical_dlog(
        epsilon: Symbol,
        variables: &[Symbol],
        letters: &[Atom],
        matrices: &[Vec<Vec<Atom>>],
    ) -> Result<Self> {
        let physical: BTreeSet<_> = variables.iter().copied().collect();
        if physical.is_empty() || physical.len() != variables.len() || physical.contains(&epsilon) {
            return Err(Error::InvalidInput(
                "canonical physical variables must be nonempty, distinct, and different from epsilon"
                    .into(),
            ));
        }
        if letters.is_empty() || letters.len() != matrices.len() {
            return Err(Error::InvalidInput(
                "canonical letters and matrices must have equal nonzero lengths".into(),
            ));
        }
        let n = square_dimension(&matrices[0])?;
        let constant_symbols = BTreeSet::new();
        for (letter, matrix) in letters.iter().zip(matrices) {
            exact_expression(letter.as_view(), &physical)?;
            if letter.together().cancel().is_zero() {
                return Err(Error::InvalidInput("a dlog letter cannot be zero".into()));
            }
            if square_dimension(matrix)? != n {
                return Err(Error::InvalidInput(
                    "canonical matrices have different dimensions".into(),
                ));
            }
            for entry in matrix.iter().flatten() {
                exact_expression(entry.as_view(), &constant_symbols)?;
            }
        }
        let regulator = Atom::var(epsilon);
        let mut derivatives = BTreeMap::new();
        for variable in variables {
            let mut matrix = vec![vec![Atom::new(); n]; n];
            for (letter, coefficients) in letters.iter().zip(matrices) {
                let logarithmic_derivative = letter.derivative(*variable) / letter;
                for (row, constants) in matrix.iter_mut().zip(coefficients) {
                    for (entry, coefficient) in row.iter_mut().zip(constants) {
                        *entry += &regulator * coefficient * &logarithmic_derivative;
                    }
                }
            }
            for entry in matrix.iter_mut().flatten() {
                *entry = entry.together().cancel();
            }
            derivatives.insert(*variable, matrix);
        }
        let result = Self {
            epsilon,
            derivatives,
        };
        result.validate()?;
        Ok(result)
    }
}

impl KinematicPath {
    pub fn validate(&self) -> Result<()> {
        if self.coordinates.is_empty() || self.coordinates.contains_key(&self.parameter) {
            return Err(Error::InvalidInput(
                "path coordinates must be nonempty and different from the parameter".into(),
            ));
        }
        let symbols = BTreeSet::from([self.parameter]);
        for coordinate in self.coordinates.values() {
            exact_expression(coordinate.as_view(), &symbols)?;
        }
        Ok(())
    }

    /// Interpolate exact points: `coordinate(x) = start + x * (end - start)`.
    /// Endpoints must have identical nonempty coordinate sets and no free symbols.
    pub fn straight_line(
        parameter: Symbol,
        start: &BTreeMap<Symbol, Atom>,
        end: &BTreeMap<Symbol, Atom>,
    ) -> Result<Self> {
        if start.is_empty() || !start.keys().eq(end.keys()) || start.contains_key(&parameter) {
            return Err(Error::InvalidInput(
                "straight-line endpoints need identical coordinates distinct from the parameter"
                    .into(),
            ));
        }
        let constant_symbols = BTreeSet::new();
        for value in start.values().chain(end.values()) {
            exact_expression(value.as_view(), &constant_symbols)?;
        }
        let x = Atom::var(parameter);
        let result = Self {
            parameter,
            coordinates: start
                .iter()
                .map(|(&variable, initial)| {
                    (
                        variable,
                        (initial + &x * (&end[&variable] - initial)).expand(),
                    )
                })
                .collect(),
        };
        result.validate()?;
        Ok(result)
    }
}
