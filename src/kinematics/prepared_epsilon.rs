use super::*;
use crate::diffexp::EpsilonSystem;
use crate::physical_conditions::{PhysicalPhaseTiming, rational_denominator_conditions};
use std::collections::HashMap;
use std::sync::Arc;
use symbolica::domains::{float::FloatField, rational_polynomial::FromNumeratorAndDenominator};

/// Immutable epsilon expansion of a rational physical connection, independent
/// of any requested point or path. Each partial connection must be regular at
/// epsilon=0. Path-dependent cancellation of epsilon poles remains available
/// through the legacy `KinematicSystem::pullback` route.
#[derive(Clone, Debug)]
pub struct PreparedKinematicEpsilon {
    coefficients: Vec<KinematicSystem>,
    // One denominator per coefficient/entry, shared by all derivative variables.
    // Only kernel indices are repeated across rows and epsilon orders.
    polynomial_coefficients: Vec<Vec<Vec<PreparedEntry>>>,
    polynomial_kernels: Vec<Arc<ExactPolynomial>>,
    affine_monomials: PreparedAffineMonomials,
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

type IntegerPolynomial = MultivariatePolynomial<IntegerRing, u16>;
type IntegerRational = RationalPolynomial<IntegerRing, u16>;

#[derive(Clone, Debug)]
struct PreparedEntry {
    denominator: usize,
    // Same BTreeMap order as the physical coordinates and their Jacobians.
    numerators: Vec<usize>,
}

// Symbolica's evaluate_with_coeff_map has no shared substitution workspace.
// The generic ring evaluator has no native polynomial coefficient-mapping
// bridge (PolynomialRing lacks ConvertToRing; polynomials lack EvaluationDomain).
// Retain only term indices and coefficient-free monomials; all arithmetic below
// is still owned by ExactPolynomial. Preparation's interning table is dropped.
#[derive(Clone, Debug)]
struct PreparedAffineMonomials {
    factors: Vec<Vec<(usize, u16)>>,
    kernels: Vec<Vec<usize>>,
}

impl PreparedAffineMonomials {
    fn new(kernels: &[Arc<ExactPolynomial>]) -> Self {
        let mut indices: HashMap<Vec<u16>, usize> = HashMap::new();
        let mut factors = Vec::new();
        let kernels = kernels
            .iter()
            .map(|kernel| {
                kernel
                    .as_ref()
                    .into_iter()
                    .map(|term| {
                        if let Some(&index) = indices.get(term.exponents) {
                            return index;
                        }
                        let index = factors.len();
                        factors.push(
                            term.exponents
                                .iter()
                                .enumerate()
                                .filter_map(|(coordinate, &power)| {
                                    (power != 0).then_some((coordinate, power))
                                })
                                .collect(),
                        );
                        indices.insert(term.exponents.to_vec(), index);
                        index
                    })
                    .collect()
            })
            .collect();
        Self { factors, kernels }
    }
}

// One line only. Coordinate powers are lazy and shared even by distinct
// monomials. Cache a monomial product only if another active kernel needs it,
// then release it at its last use. Unrequested orders and stationary numerators
// do not keep products alive. Zero is a valid cached result, never a sentinel.
struct AffineMonomialRestriction<'a> {
    prepared: &'a PreparedAffineMonomials,
    coordinates: &'a [ExactPolynomial],
    powers: HashMap<(usize, u16), ExactPolynomial>,
    products: Vec<Option<ExactPolynomial>>,
    remaining: Vec<usize>,
    product_count: usize,
}

impl<'a> AffineMonomialRestriction<'a> {
    fn new(
        prepared: &'a PreparedAffineMonomials,
        coordinates: &'a [ExactPolynomial],
        active: &[bool],
    ) -> Self {
        let mut remaining = vec![0; prepared.factors.len()];
        for (terms, &active) in prepared.kernels.iter().zip(active) {
            if active {
                for &monomial in terms {
                    remaining[monomial] += 1;
                }
            }
        }
        Self {
            prepared,
            coordinates,
            powers: HashMap::new(),
            products: vec![None; prepared.factors.len()],
            remaining,
            product_count: 0,
        }
    }

    fn restrict(&mut self, kernel: usize, polynomial: &ExactPolynomial) -> ExactPolynomial {
        let mut result = self.coordinates[0].zero();
        for (term, &monomial) in polynomial.into_iter().zip(&self.prepared.kernels[kernel]) {
            let remaining = &mut self.remaining[monomial];
            debug_assert!(*remaining > 0);
            *remaining -= 1;
            let product = if *remaining == 0 && self.products[monomial].is_some() {
                self.products[monomial].take().unwrap()
            } else if let Some(product) = &self.products[monomial] {
                product.clone()
            } else {
                let mut product: Option<ExactPolynomial> = None;
                let factors = &self.prepared.factors[monomial];
                for &(coordinate, power) in factors {
                    let value = self
                        .powers
                        .entry((coordinate, power))
                        .or_insert_with(|| self.coordinates[coordinate].pow(usize::from(power)));
                    product = Some(match product {
                        None => value.clone(),
                        Some(product) => product * &*value,
                    });
                }
                let product = product.unwrap_or_else(|| self.coordinates[0].one());
                if factors.len() > 1 {
                    self.product_count += 1;
                    if *remaining != 0 {
                        self.products[monomial] = Some(product.clone());
                    }
                }
                product
            };
            result = result + product.mul_coeff(term.coefficient.clone());
        }
        result
    }
}

fn intern_polynomial<P: Eq + std::hash::Hash>(
    polynomial: P,
    polynomials: &mut Vec<Arc<P>>,
    indices: &mut HashMap<Arc<P>, usize>,
) -> usize {
    if let Some(&index) = indices.get(&polynomial) {
        return index;
    }
    let index = polynomials.len();
    let polynomial = Arc::new(polynomial);
    indices.insert(polynomial.clone(), index);
    polynomials.push(polynomial);
    index
}

fn exact_quotient(
    numerator: &IntegerPolynomial,
    denominator: &IntegerPolynomial,
) -> Result<IntegerPolynomial> {
    let (quotient, remainder) = numerator.quot_rem(denominator, true);
    if !remainder.is_zero() {
        return Err(Error::Unsupported(
            "native prepared denominator quotient was not exact".into(),
        ));
    }
    Ok(quotient)
}

// Preparation-only native forms. Integer contents (including cleared rational
// constants) participate in LCMs, so every cofactor is an exact Z quotient.
// Never form a*b before dividing by gcd(a,b). All these tables are dropped
// before returning the prepared connection; only the Gaussian kernels survive.
#[derive(Default)]
struct DenominatorPreparation {
    polynomials: Vec<Arc<IntegerPolynomial>>,
    indices: HashMap<Arc<IntegerPolynomial>, usize>,
    lcms: HashMap<(usize, usize), usize>,
    cofactors: HashMap<(usize, usize), usize>,
}

impl DenominatorPreparation {
    fn intern(&mut self, polynomial: IntegerPolynomial) -> usize {
        intern_polynomial(polynomial, &mut self.polynomials, &mut self.indices)
    }

    fn lcm(&mut self, a: usize, b: usize) -> Result<usize> {
        if a == b {
            return Ok(a);
        }
        let key = (a.min(b), a.max(b));
        if let Some(&index) = self.lcms.get(&key) {
            return Ok(index);
        }
        let a = &*self.polynomials[key.0];
        let b = &*self.polynomials[key.1];
        let gcd = a.gcd(b);
        let common = exact_quotient(a, &gcd)? * b;
        let index = self.intern(common);
        self.lcms.insert(key, index);
        Ok(index)
    }

    fn cofactor(&mut self, common: usize, partial: usize) -> Result<usize> {
        let key = (common, partial);
        if let Some(&index) = self.cofactors.get(&key) {
            return Ok(index);
        }
        let quotient = exact_quotient(&self.polynomials[common], &self.polynomials[partial])?;
        let index = self.intern(quotient);
        self.cofactors.insert(key, index);
        Ok(index)
    }
}

// Specialize only formal_i at setup, using the owner's native polynomial map.
// Physical variables stay independent until a path is requested.
fn gaussian_polynomial(
    polynomial: &IntegerPolynomial,
    coordinates: &[ExactPolynomial],
) -> ExactPolynomial {
    let field = gaussian_field();
    polynomial.evaluate_with_coeff_map(
        |c| coordinates[0].constant(Complex::new(Rational::from(c.clone()), Rational::from(0))),
        coordinates,
        &PolynomialRing::new(field),
    )
}

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
// from_num_den ignores do_gcd. Real physical lines can use the native rational
// polynomial GCD, avoiding exact Euclidean coefficient swell. Retain the
// Gaussian Euclidean route for genuinely complex coefficients.
fn normalize_line_entry(mut entry: ExactRational) -> ExactRational {
    let field = gaussian_field();
    if entry.is_zero() {
        return ExactRational::new(&field, entry.get_variables().clone());
    }
    if !entry.numerator.is_constant() && !entry.denominator.is_constant() {
        let real = entry
            .numerator
            .coefficients
            .iter()
            .chain(&entry.denominator.coefficients)
            .all(|c| c.im.is_zero());
        let gcd = if real {
            let numerator = entry.numerator.map_coeff(|c| c.re.clone(), Q);
            let denominator = entry.denominator.map_coeff(|c| c.re.clone(), Q);
            numerator.gcd(&denominator).map_coeff(
                |c| Complex::new(c.clone(), Rational::from(0)),
                field.clone(),
            )
        } else {
            let numerator = entry.numerator.to_univariate_from_univariate(0);
            let denominator = entry.denominator.to_univariate_from_univariate(0);
            numerator
                .gcd_euclidean(&denominator)
                .to_multivariate::<u16>()
        };
        if !gcd.is_constant() {
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
        let (polynomial_coefficients, polynomial_kernels) = {
            let _timing = PhysicalPhaseTiming::start(
                "transport.epsilon_prepare.kernel_preparation",
                0,
                0,
                None,
            );
            let field = gaussian_field();
            let template = ExactPolynomial::new(&field, None, polynomial_variables.clone());
            let mut native_variables = polynomial_variables.as_ref().clone();
            native_variables.push(PolyVariable::Symbol(crate::family::imaginary_parameter()));
            let native_variables = Arc::new(native_variables);
            let mut native_coordinates = polynomial_variables
                .iter()
                .map(|variable| template.variable(variable).map_err(Error::Unsupported))
                .collect::<Result<Vec<_>>>()?;
            native_coordinates
                .push(template.constant(Complex::new(Rational::from(0), Rational::from(1))));
            let mut denominators = DenominatorPreparation::default();
            let mut kernels = Vec::new();
            let mut kernel_indices = HashMap::new();
            let mut denominator_kernels = HashMap::new();
            let mut polynomials = Vec::with_capacity(coefficients.len());
            for (order, coefficient) in coefficients.iter().enumerate() {
                let n = coefficient.derivatives.values().next().unwrap().len();
                let mut matrix = Vec::with_capacity(n);
                for row in 0..n {
                    let _timing = PhysicalPhaseTiming::start(
                        "transport.epsilon_prepare.coefficient_row_kernels",
                        order,
                        row,
                        None,
                    );
                    let mut entries = Vec::with_capacity(n);
                    for column in 0..n {
                        // Parse one matrix entry across partials at a time, not
                        // another full collection of all coefficient matrices.
                        let mut partials = Vec::with_capacity(coefficient.derivatives.len());
                        for partial in coefficient.derivatives.values() {
                            // The Q -> Z conversion clears nonzero constant
                            // rational denominators into the native quotient.
                            let rational: IntegerRational =
                                crate::family::encode_complex(&partial[row][column])
                                    .try_to_rational_polynomial(&Q, &Z, native_variables.clone())
                                    .map_err(|e| {
                                        Error::Unsupported(format!(
                                            "rational prepared coefficient required: {e}"
                                        ))
                                    })?;
                            if rational.get_variables() != &native_variables {
                                return Err(Error::Unsupported(
                                    "unexpected variable in prepared polynomial coefficient".into(),
                                ));
                            }
                            let denominator = denominators.intern(rational.denominator);
                            partials.push((rational.numerator, denominator));
                        }
                        let common = {
                            let _timing = PhysicalPhaseTiming::start(
                                "transport.epsilon_prepare.common_denominator_lcm",
                                order,
                                row * n + column,
                                None,
                            );
                            let mut common = partials[0].1;
                            for &(_, denominator) in &partials[1..] {
                                common = denominators.lcm(common, denominator)?;
                            }
                            common
                        };
                        let denominator = if let Some(&kernel) = denominator_kernels.get(&common) {
                            kernel
                        } else {
                            let polynomial = gaussian_polynomial(
                                &denominators.polynomials[common],
                                &native_coordinates,
                            );
                            if polynomial.is_zero() {
                                return Err(Error::Unsupported(
                                    "prepared common denominator vanishes at formal i".into(),
                                ));
                            }
                            let kernel =
                                intern_polynomial(polynomial, &mut kernels, &mut kernel_indices);
                            denominator_kernels.insert(common, kernel);
                            kernel
                        };
                        let mut numerators = Vec::with_capacity(partials.len());
                        for (numerator, partial_denominator) in partials {
                            let cofactor = denominators.cofactor(common, partial_denominator)?;
                            let cleared = numerator * &*denominators.polynomials[cofactor];
                            let polynomial = gaussian_polynomial(&cleared, &native_coordinates);
                            numerators.push(intern_polynomial(
                                polynomial,
                                &mut kernels,
                                &mut kernel_indices,
                            ));
                        }
                        entries.push(PreparedEntry {
                            denominator,
                            numerators,
                        });
                    }
                    matrix.push(entries);
                }
                polynomials.push(matrix);
            }
            (polynomials, kernels)
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
        let affine_monomials = {
            let _timing = PhysicalPhaseTiming::start(
                "transport.epsilon_prepare.affine_monomials",
                0,
                0,
                None,
            );
            PreparedAffineMonomials::new(&polynomial_kernels)
        };
        Ok(PreparedKinematicEpsilon {
            coefficients,
            polynomial_coefficients,
            polynomial_kernels,
            affine_monomials,
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
        // An affine coordinate's constant term is its exact value at path
        // parameter zero. Keep epsilon last, as in each retained guard.
        let witness_points = coordinates.map(|coordinates| {
            let mut point = coordinates
                .iter()
                .map(|p| p.get_constant())
                .collect::<Vec<_>>();
            point.push(field.zero());
            let mut alternate = point.clone();
            *alternate.last_mut().unwrap() = field.one();
            [point, alternate]
        });
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
                    let witnessed_nonzero = {
                        let _timing = PhysicalPhaseTiming::start(
                            "transport.pullback.native_guard_witness",
                            index,
                            usize::from(expansion),
                            None,
                        );
                        witness_points.as_ref().unwrap().iter().any(|point| {
                            !field.is_zero(&polynomial.evaluate_with_coeff_map(
                                |c| c.clone(),
                                point,
                                &field,
                            ))
                        })
                    };
                    // One exact nonzero value proves the formal (parameter,
                    // epsilon) restriction is not identically zero. Zero at
                    // both points proves nothing: a source may be on a guard
                    // while its line leaves it. Preserve the full exact check.
                    if witnessed_nonzero {
                        false
                    } else {
                        let _timing = PhysicalPhaseTiming::start(
                            "transport.pullback.native_guard_fallback",
                            index,
                            usize::from(expansion),
                            None,
                        );
                        polynomial
                            .evaluate_with_coeff_map(
                                |c| coordinates[0].constant(c.clone()),
                                coordinates,
                                &ring,
                            )
                            .is_zero()
                    }
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

    fn restrict_kernel<'a>(
        &self,
        index: usize,
        monomials: &mut AffineMonomialRestriction<'_>,
        restrictions: &'a mut [Option<ExactPolynomial>],
        unique_count: &mut usize,
    ) -> &'a ExactPolynomial {
        restrictions[index].get_or_insert_with(|| {
            *unique_count += 1;
            monomials.restrict(index, &self.polynomial_kernels[index])
        })
    }

    fn pullback_affine(
        &self,
        path: &KinematicPath,
        order: usize,
        coordinates: &[ExactPolynomial],
    ) -> Result<EpsilonSystem> {
        let field = gaussian_field();
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
        let mut restrictions = vec![None; self.polynomial_kernels.len()];
        let mut monomials = {
            let _timing = PhysicalPhaseTiming::start(
                "transport.pullback.affine_monomial_setup",
                order,
                0,
                None,
            );
            let mut active = vec![false; self.polynomial_kernels.len()];
            for entry in self.polynomial_coefficients[..=order]
                .iter()
                .flatten()
                .flatten()
            {
                active[entry.denominator] = true;
                for (&kernel, jacobian) in entry.numerators.iter().zip(&jacobians) {
                    if !field.is_zero(jacobian) && !self.polynomial_kernels[kernel].is_zero() {
                        active[kernel] = true;
                    }
                }
            }
            AffineMonomialRestriction::new(&self.affine_monomials, coordinates, &active)
        };
        let mut unique_count = 0;
        let mut matrices = Vec::with_capacity(order + 1);
        for (coefficient, entries) in self.polynomial_coefficients[..=order].iter().enumerate() {
            let mut rows = Vec::with_capacity(n);
            {
                let _timing = PhysicalPhaseTiming::start(
                    "transport.pullback.precomputed_restriction",
                    coefficient,
                    0,
                    None,
                );
                for prepared_row in entries {
                    let mut row = Vec::with_capacity(n);
                    for entry in prepared_row {
                        // The common denominator contains every partial's
                        // denominator, including stationary coordinates. Check
                        // it before skipping any zero Jacobian or numerator.
                        if self
                            .restrict_kernel(
                                entry.denominator,
                                &mut monomials,
                                &mut restrictions,
                                &mut unique_count,
                            )
                            .is_zero()
                        {
                            return Err(Error::InvalidInput(
                                "path lies in a prepared coefficient pole".into(),
                            ));
                        }
                        let mut numerator = template.zero();
                        for (&kernel, jacobian) in entry.numerators.iter().zip(&jacobians) {
                            if field.is_zero(jacobian) || self.polynomial_kernels[kernel].is_zero()
                            {
                                continue;
                            }
                            let partial = self.restrict_kernel(
                                kernel,
                                &mut monomials,
                                &mut restrictions,
                                &mut unique_count,
                            );
                            numerator = numerator + partial.clone().mul_coeff(jacobian.clone());
                        }
                        // Exactly one fraction; no per-request fraction algebra
                        // can multiply repeated Gram factors across partials.
                        row.push(ExactRational::from_num_den(
                            numerator,
                            restrictions[entry.denominator].as_ref().unwrap().clone(),
                            &field,
                            false,
                        ));
                    }
                    rows.push(row);
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
        // Timing metadata only: index is the number actually restricted on
        // this line, subindex is the number retained across all epsilon orders.
        let _counts = PhysicalPhaseTiming::start(
            "transport.pullback.unique_restriction_count",
            unique_count,
            self.polynomial_kernels.len(),
            None,
        );
        // Counts only; no coordinate or coefficient values enter timings.
        let _monomial_counts = PhysicalPhaseTiming::start(
            "transport.pullback.affine_power_product_count",
            monomials.powers.len(),
            monomials.product_count,
            None,
        );
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

#[cfg(test)]
mod affine_monomial_tests {
    use super::*;

    #[test]
    fn distinct_kernels_share_exact_powers_and_products_until_last_use() {
        let (s, t, x) = symbol!(
            "affine_monomial_test::s",
            "affine_monomial_test::t",
            "affine_monomial_test::x"
        );
        let (sv, tv, xv) = (Atom::var(s), Atom::var(t), Atom::var(x));
        let field = gaussian_field();
        let variables = Arc::new(vec![PolyVariable::Symbol(s), PolyVariable::Symbol(t)]);
        let m = sv.pow(3) * tv.pow(2);
        let n = sv.pow(2) * tv.pow(3);
        let expressions: [Atom; 2] = [m.clone() + 2 * &n + 3, 5 * &m - &n + sv.pow(3)];
        let kernels = expressions
            .iter()
            .map(|a| Arc::new(a.try_to_polynomial(&field, variables.clone()).unwrap()))
            .collect::<Vec<Arc<ExactPolynomial>>>();
        let prepared = PreparedAffineMonomials::new(&kernels);
        for source in [3 + &xv / 7, Atom::zero()] {
            let variables = Arc::new(vec![PolyVariable::Symbol(x)]);
            let coordinates = [source, 2 - &xv / 11]
                .iter()
                .map(|a| a.try_to_polynomial(&field, variables.clone()).unwrap())
                .collect::<Vec<ExactPolynomial>>();
            let mut line = AffineMonomialRestriction::new(&prepared, &coordinates, &[true, true]);
            for (index, kernel) in kernels.iter().enumerate() {
                assert_eq!(
                    line.restrict(index, kernel),
                    kernel.evaluate_with_coeff_map(
                        |c| coordinates[0].constant(c.clone()),
                        &coordinates,
                        &PolynomialRing::new(field.clone()),
                    )
                );
                // Both full products and all four coordinate powers are
                // computed once, including when their exact value is zero.
                assert_eq!(line.product_count, 2);
                assert_eq!(line.powers.len(), 4);
            }
            assert!(line.products.iter().all(Option::is_none));
            assert!(line.remaining.iter().all(|&uses| uses == 0));
        }
    }
}
