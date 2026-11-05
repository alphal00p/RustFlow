use super::*;
use crate::diffexp::EpsilonSystem;
use crate::physical_conditions::{PhysicalPhaseTiming, rational_denominator_conditions};
use std::sync::Arc;
use symbolica::domains::{float::FloatField, rational_polynomial::FromNumeratorAndDenominator};

/// Immutable epsilon expansion of a rational physical connection, independent
/// of any requested point or path. Each partial connection must be regular at
/// epsilon=0. Path-dependent cancellation of epsilon poles remains available
/// through the legacy `KinematicSystem::pullback` route.
#[derive(Clone, Debug)]
pub struct PreparedKinematicEpsilon {
    coefficients: Vec<KinematicSystem>,
    // Parsed once, before any path or numerical precision is selected.
    polynomial_coefficients: Vec<BTreeMap<Symbol, Vec<Vec<ExactRational>>>>,
    original_poles: Vec<PreparedPole>,
    expansion_poles: Vec<PreparedPole>,
}

fn decode_complex(value: &Atom) -> Atom {
    value
        .replace(Atom::var(crate::family::imaginary_parameter()))
        .with(Atom::i())
}

// FloatField is Symbolica's field adapter, including for exact Gaussian
// rationals. No floating conversion is performed by these types.
type GaussianField = FloatField<Complex<Rational>>;
type ExactPolynomial = MultivariatePolynomial<GaussianField, u16>;
type ExactRational = RationalPolynomial<GaussianField, u16>;

#[derive(Clone, Debug)]
struct PreparedPole {
    // Retain the checked symbolic guard for non-affine paths or a polynomial
    // conversion the native owner cannot represent.
    expression: Atom,
    // Physical coordinates in the same order as polynomial_coefficients,
    // followed by formal epsilon. Never specialize epsilon in an original pole.
    polynomial: Option<ExactPolynomial>,
}

fn gaussian_field() -> GaussianField {
    FloatField::from_rep(Complex::new(Rational::from(0), Rational::from(0)))
}

// FloatField's generic PolynomialGCD is a stub returning one, and its
// from_num_den ignores do_gcd. Use the owner's explicit univariate Euclidean
// GCD only on the completed row entries, never on expanded Atom expressions.
fn normalize_line_entry(mut entry: ExactRational) -> ExactRational {
    let field = gaussian_field();
    if entry.is_zero() {
        return ExactRational::new(&field, entry.get_variables().clone());
    }
    if !entry.numerator.is_constant() && !entry.denominator.is_constant() {
        let numerator = entry.numerator.to_univariate_from_univariate(0);
        let denominator = entry.denominator.to_univariate_from_univariate(0);
        let gcd = numerator.gcd_euclidean(&denominator);
        if !gcd.is_constant() {
            let gcd = gcd.to_multivariate::<u16>();
            entry.numerator = &entry.numerator / &gcd;
            entry.denominator = &entry.denominator / &gcd;
        }
    }
    ExactRational::from_num_den(entry.numerator, entry.denominator, &field, false)
}

impl KinematicSystem {
    /// Expand all rational partial matrices once through `max_order` inclusive.
    /// For a solution requested on [leading,last], prepare at least last-leading.
    /// No numerical precision, boundary credit or cache identity is changed.
    pub fn prepare_epsilon(&self, max_order: usize) -> Result<PreparedKinematicEpsilon> {
        let _timing = PhysicalPhaseTiming::start("transport.epsilon_prepare", 0, 0, None);
        self.validate()?;
        if max_order >= 1024 {
            return Err(Error::Limit("at most 1024 epsilon orders supported".into()));
        }
        let mut variables: BTreeSet<_> = self.derivatives.keys().copied().collect();
        variables.insert(self.epsilon);
        let expressions = self
            .derivatives
            .values()
            .flatten()
            .flatten()
            .cloned()
            .collect::<Vec<_>>();
        // Keep original denominators even above the retained order or in a
        // stationary coordinate; expansion/cancellation must not erase them.
        let original_poles = rational_denominator_conditions(&expressions, &variables)?;
        let mut coefficients = vec![
            KinematicSystem {
                epsilon: self.epsilon,
                derivatives: BTreeMap::new(),
            };
            max_order + 1
        ];
        let mut expansion_poles = Vec::new();
        for (&variable, matrix) in &self.derivatives {
            let expanded = EpsilonSystem::from_differential_system(
                &DifferentialSystem {
                    variable,
                    matrix: matrix.clone(),
                },
                self.epsilon,
                max_order,
            )?;
            for (coefficient, mut matrix) in coefficients.iter_mut().zip(expanded.matrices) {
                for entry in matrix.iter_mut().flatten() {
                    *entry = decode_complex(entry);
                    exact_expression(entry.as_view(), &variables)?;
                    if entry.contains_symbol(self.epsilon) {
                        return Err(Error::Unsupported(
                            "epsilon remained in a prepared coefficient".into(),
                        ));
                    }
                }
                coefficient.derivatives.insert(variable, matrix);
            }
            for entry in matrix.iter().flatten() {
                let rational: RationalPolynomial<IntegerRing, u16> =
                    crate::family::encode_complex(entry)
                        .try_to_rational_polynomial(&Q, &Z, None)
                        .map_err(|e| {
                            Error::Unsupported(format!("rational epsilon system required: {e}"))
                        })?;
                // A zero constant denominator after specialization can change
                // epsilon valuation, even when the specialized legacy function
                // is regular (e.g. epsilon/(s+epsilon) at s=0). Reject it rather
                // than silently using the generic Taylor coefficients.
                let constant = decode_complex(&rational.denominator.to_expression())
                    .replace(Atom::var(self.epsilon))
                    .with(Atom::zero())
                    .expand();
                if constant.is_zero() {
                    return Err(Error::Unsupported(
                        "partial matrix has an epsilon pole; supply an epsilon-regular basis"
                            .into(),
                    ));
                }
                if !matches!(constant.as_view(), AtomView::Num(_)) {
                    expansion_poles.push(constant);
                }
            }
        }
        expansion_poles.sort_by_cached_key(AtomCore::to_canonical_string);
        expansion_poles.dedup();
        let polynomial_variables = Arc::new(
            self.derivatives
                .keys()
                .copied()
                .map(PolyVariable::Symbol)
                .collect::<Vec<_>>(),
        );
        let polynomial_coefficients = {
            let _timing =
                PhysicalPhaseTiming::start("transport.epsilon_prepare.polynomials", 0, 0, None);
            let field = gaussian_field();
            let variables = polynomial_variables.clone();
            let mut polynomials = Vec::with_capacity(coefficients.len());
            for coefficient in &coefficients {
                let mut partials = BTreeMap::new();
                for (&variable, matrix) in &coefficient.derivatives {
                    let mut polynomial_matrix = Vec::with_capacity(matrix.len());
                    for row in matrix {
                        let mut polynomial_row = Vec::with_capacity(row.len());
                        for entry in row {
                            let rational: ExactRational = entry
                                .try_to_rational_polynomial(&field, &field, variables.clone())
                                .map_err(|e| {
                                    Error::Unsupported(format!(
                                        "rational prepared coefficient required: {e}"
                                    ))
                                })?;
                            if rational.get_variables() != &variables {
                                return Err(Error::Unsupported(
                                    "unexpected variable in prepared polynomial coefficient".into(),
                                ));
                            }
                            polynomial_row.push(rational);
                        }
                        polynomial_matrix.push(polynomial_row);
                    }
                    partials.insert(variable, polynomial_matrix);
                }
                polynomials.push(partials);
            }
            polynomials
        };
        let (original_poles, expansion_poles) = {
            let _timing =
                PhysicalPhaseTiming::start("transport.epsilon_prepare.guards", 0, 0, None);
            let field = gaussian_field();
            let mut variables = polynomial_variables.as_ref().clone();
            variables.push(PolyVariable::Symbol(self.epsilon));
            let variables = Arc::new(variables);
            let prepare_poles = |poles: Vec<Atom>| {
                poles
                    .into_iter()
                    .map(|expression| {
                        let polynomial: Option<ExactPolynomial> = expression
                            .try_to_polynomial(&field, variables.clone())
                            .ok()
                            .filter(|p: &ExactPolynomial| p.variables() == &variables);
                        PreparedPole {
                            expression,
                            polynomial,
                        }
                    })
                    .collect()
            };
            (
                prepare_poles(original_poles),
                prepare_poles(expansion_poles),
            )
        };
        Ok(PreparedKinematicEpsilon {
            coefficients,
            polynomial_coefficients,
            original_poles,
            expansion_poles,
        })
    }
}

impl PreparedKinematicEpsilon {
    /// Maximum retained connection order (inclusive), independent of the
    /// solution's common Laurent leading power.
    pub fn max_order(&self) -> usize {
        self.coefficients.len() - 1
    }

    pub(crate) fn check_order(&self, order: usize) -> Result<()> {
        if order > self.max_order() {
            return Err(Error::InvalidInput(format!(
                "prepared epsilon depth {} does not cover requested order {order}",
                self.max_order()
            )));
        }
        Ok(())
    }

    // Return None explicitly for non-affine or non-rational-complex polynomial
    // coordinates. Those paths keep the existing symbolic pullback below.
    fn affine_coordinates(&self, path: &KinematicPath) -> Option<Vec<ExactPolynomial>> {
        let field = gaussian_field();
        let variables = Arc::new(vec![PolyVariable::Symbol(path.parameter)]);
        path.coordinates
            .values()
            .map(|coordinate| {
                let polynomial: ExactPolynomial = coordinate
                    .try_to_polynomial(&field, variables.clone())
                    .ok()?;
                (polynomial.variables() == &variables && polynomial.degree(0) <= 1)
                    .then_some(polynomial)
            })
            .collect()
    }

    fn check_poles(
        &self,
        path: &KinematicPath,
        coordinates: Option<&[ExactPolynomial]>,
    ) -> Result<()> {
        let _timing =
            PhysicalPhaseTiming::start("transport.pullback.guard_restriction", 0, 0, None);
        let field = gaussian_field();
        let ring = PolynomialRing::new(field.clone());
        // The coefficient pullback stays univariate. Only the guards need a
        // second variable so epsilon remains independent of the path parameter.
        let guard_coordinates = if let Some(coordinates) = coordinates {
            let variables = Arc::new(vec![
                PolyVariable::Symbol(path.parameter),
                PolyVariable::Symbol(self.coefficients[0].epsilon),
            ]);
            let parameter: ExactPolynomial = Atom::var(path.parameter)
                .try_to_polynomial(&field, variables.clone())
                .map_err(|e| {
                    Error::Unsupported(format!("polynomial path parameter required: {e}"))
                })?;
            let epsilon: ExactPolynomial = Atom::var(self.coefficients[0].epsilon)
                .try_to_polynomial(&field, variables)
                .map_err(|e| Error::Unsupported(format!("polynomial epsilon required: {e}")))?;
            let mut lifted = coordinates
                .iter()
                .map(|p| {
                    p.evaluate_with_coeff_map(
                        |c| parameter.constant(c.clone()),
                        std::slice::from_ref(&parameter),
                        &ring,
                    )
                })
                .collect::<Vec<_>>();
            lifted.push(epsilon);
            Some(lifted)
        } else {
            None
        };
        // Allocate symbolic substitutions only for a path or guard that needs
        // the fallback. Ordinary affine requests never substitute/expand Atoms.
        let substitutions = (coordinates.is_none()
            || self
                .original_poles
                .iter()
                .chain(&self.expansion_poles)
                .any(|pole| pole.polynomial.is_none()))
        .then(|| {
            path.coordinates
                .iter()
                .map(|(&s, a)| (Atom::var(s), a.clone()))
                .collect()
        });
        for (poles, expansion, phase) in [
            (
                &self.original_poles,
                false,
                "transport.pullback.original_guard_restriction",
            ),
            (
                &self.expansion_poles,
                true,
                "transport.pullback.expansion_guard_restriction",
            ),
        ] {
            for (index, pole) in poles.iter().enumerate() {
                let _timing = PhysicalPhaseTiming::start(phase, index, 0, None);
                let is_zero = if let (Some(polynomial), Some(coordinates)) =
                    (&pole.polynomial, &guard_coordinates)
                {
                    polynomial
                        .evaluate_with_coeff_map(
                            |c| coordinates[0].constant(c.clone()),
                            coordinates,
                            &ring,
                        )
                        .is_zero()
                } else {
                    let _timing = PhysicalPhaseTiming::start(
                        "transport.pullback.symbolic_guard_fallback",
                        index,
                        usize::from(expansion),
                        None,
                    );
                    substitute(&pole.expression, substitutions.as_ref().unwrap())
                        .expand()
                        .is_zero()
                };
                if is_zero {
                    return Err(if expansion {
                        Error::Unsupported("path has an exceptional epsilon specialization; use the legacy pullback".into())
                    } else {
                        Error::InvalidInput("path lies in an original kinematic matrix pole".into())
                    });
                }
            }
        }
        Ok(())
    }

    fn pullback_affine(
        &self,
        path: &KinematicPath,
        order: usize,
        coordinates: &[ExactPolynomial],
    ) -> Result<EpsilonSystem> {
        let field = gaussian_field();
        let ring = PolynomialRing::new(field.clone());
        let template = &coordinates[0];
        let jacobians = coordinates
            .iter()
            .map(|p| p.derivative(0).get_constant())
            .collect::<Vec<_>>();
        let n = self.coefficients[0]
            .derivatives
            .values()
            .next()
            .unwrap()
            .len();
        let mut matrices = Vec::with_capacity(order + 1);
        for (coefficient, partials) in self.polynomial_coefficients[..=order].iter().enumerate() {
            let mut rows =
                vec![vec![ExactRational::new(&field, template.variables().clone()); n]; n];
            {
                let _timing = PhysicalPhaseTiming::start(
                    "transport.pullback.precomputed_restriction",
                    coefficient,
                    0,
                    None,
                );
                for (partials, jacobian) in partials.values().zip(&jacobians) {
                    for (row, partial_row) in rows.iter_mut().zip(partials) {
                        for (entry, partial) in row.iter_mut().zip(partial_row) {
                            let compose = |p: &ExactPolynomial| {
                                p.evaluate_with_coeff_map(
                                    |c| template.constant(c.clone()),
                                    coordinates,
                                    &ring,
                                )
                            };
                            // Always restrict the denominator, even for zero Jacobians.
                            let denominator = compose(&partial.denominator);
                            if denominator.is_zero() {
                                return Err(Error::InvalidInput(
                                    "path lies in a prepared coefficient pole".into(),
                                ));
                            }
                            if field.is_zero(jacobian) || partial.is_zero() {
                                continue;
                            }
                            let numerator = compose(&partial.numerator).mul_coeff(jacobian.clone());
                            let restricted =
                                ExactRational::from_num_den(numerator, denominator, &field, false);
                            // Native rational addition combines equal denominators directly.
                            // Defer polynomial cancellation until the weighted sum is complete.
                            *entry = if entry.is_zero() {
                                restricted
                            } else {
                                &*entry + &restricted
                            };
                        }
                    }
                }
            }
            let mut matrix = Vec::with_capacity(n);
            for (row_index, row) in rows.into_iter().enumerate() {
                let _timing = PhysicalPhaseTiming::start(
                    "transport.pullback.coefficient_row_normalization",
                    coefficient,
                    row_index,
                    None,
                );
                matrix.push(
                    row.into_iter()
                        .map(|entry| {
                            normalize_line_entry(entry).to_expression_with_coeff_map(|_, c, out| {
                                // Emit the existing formal-i convention once, at the boundary.
                                *out = Atom::num(c.re.clone())
                                    + Atom::var(crate::family::imaginary_parameter())
                                        * Atom::num(c.im.clone());
                            })
                        })
                        .collect(),
                );
            }
            let pulled = DifferentialSystem {
                variable: path.parameter,
                matrix,
            };
            pulled.validate()?;
            matrices.push(pulled.matrix);
        }
        Ok(EpsilonSystem {
            variable: path.parameter,
            matrices,
        })
    }

    /// Restrict rational-complex affine paths with native polynomial composition;
    /// other paths explicitly fall back to the symbolic coefficient pullback.
    /// Original poles
    /// and exceptional epsilon specializations are checked before Jacobians or
    /// truncation can hide them. Exceptional lines are rejected; callers may
    /// explicitly try the legacy route to obtain their specialized expansion.
    pub fn pullback(&self, path: &KinematicPath, order: usize) -> Result<EpsilonSystem> {
        self.check_order(order)?;
        {
            let _timing =
                PhysicalPhaseTiming::start("transport.pullback.path_validation", 0, 0, None);
            self.coefficients[0].validate_path(path)?;
        }
        let coordinates = {
            let _timing =
                PhysicalPhaseTiming::start("transport.pullback.affine_classification", 0, 0, None);
            self.affine_coordinates(path)
        };
        self.check_poles(path, coordinates.as_deref())?;
        if let Some(coordinates) = coordinates {
            return self.pullback_affine(path, order, &coordinates);
        }
        // Explicit non-affine fallback: preserve the legacy symbolic behavior.
        let _timing =
            PhysicalPhaseTiming::start("transport.pullback.non_affine_fallback", 0, 0, None);
        let mut matrices = Vec::with_capacity(order + 1);
        for coefficient in &self.coefficients[..=order] {
            let pulled = coefficient.pullback_validated(path)?;
            pulled.validate()?;
            let mut matrix = pulled.matrix;
            // Match EpsilonSystem's formal complex-coefficient convention.
            for entry in matrix.iter_mut().flatten() {
                *entry = crate::family::encode_complex(entry).together().cancel();
            }
            matrices.push(matrix);
        }
        Ok(EpsilonSystem {
            variable: path.parameter,
            matrices,
        })
    }
}
