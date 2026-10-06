//! Optional signed fundamental matrices for supplied boundary uncertainty.
//!
//! Polynomial candidates remain numerical. Native ball products certify a
//! whole-chart inverse and the retained exact-source defect; a small cumulative
//! defect/jump budget gives an arithmetic-only Neumann bound. No global accuracy
//! claim replaces the supplied evidence or the independent central refinement.
use super::source::{
    BallPolynomial, ExactSourceResidual, ExactSpecialization, complex_upper, dyadic_ball,
    exact_ball, polynomial_lower,
};
use super::{CompiledSystem, TaylorSegment};
use crate::{ComplexFloat as C, Error, Precision, Result, RunContext};
use std::sync::Arc;
use symbolica::domains::float::{ComplexBall, FloatField, RoundingDirection};
use symbolica::poly::univariate::{UnivariatePolynomial, UnivariatePolynomialRing};
use symbolica::prelude::*;

type BallMatrix = Matrix<FloatField<ComplexBall>>;
type PolynomialMatrix = Matrix<UnivariatePolynomialRing<FloatField<ComplexBall>>>;

// Admission estimates, not a hard allocator/time quota. A full epsilon
// hierarchy must not be silently flattened into this dense owner.
const MAX_DIMENSION: usize = 16;
const MAX_ORDER: usize = 256;
const MAX_BITS: u32 = 4096;
const MAX_COEFFICIENTS: u64 = 1_000_000;
const MAX_WORK: u64 = 100_000_000;
const MAX_STORAGE_ESTIMATE: u64 = 256 * 1024 * 1024;
const MAX_EXACT_BITS: u64 = 65_536;

#[derive(Debug)]
pub(crate) struct FundamentalErrors {
    /// One absolute source-error vector per supplied physical segment.
    pub(crate) endpoints: Vec<Vec<Float>>,
    #[cfg(test)]
    pub(crate) jump_budget: Float,
}

struct Arithmetic {
    p: Precision,
    weights: Vec<Float>,
    field: FloatField<ComplexBall>,
    polynomial_field: UnivariatePolynomialRing<FloatField<ComplexBall>>,
    variable: Arc<PolyVariable>,
}
impl Arithmetic {
    fn new(p: Precision, weights: &[Float]) -> Self {
        let field = FloatField::from_rep(dyadic_ball(p, &p.zero()));
        let variable = Arc::new(PolyVariable::Temporary(0));
        let polynomial_field = UnivariatePolynomialRing::new(field.clone(), variable.clone());
        Self {
            p,
            weights: weights.to_vec(),
            field,
            polynomial_field,
            variable,
        }
    }
    fn add(&self, a: &Float, b: &Float) -> Float {
        a.add_round(b, self.p.bits, RoundingDirection::Up)
    }
    fn mul(&self, a: &Float, b: &Float) -> Float {
        a.mul_round(b, self.p.bits, RoundingDirection::Up)
    }
    fn div(&self, a: &Float, b: &Float) -> Float {
        a.div_round(b, self.p.bits, RoundingDirection::Up)
    }
    fn complement(&self, a: &Float) -> Float {
        self.p
            .real(1)
            .sub_round(a, self.p.bits, RoundingDirection::Down)
    }
    fn polynomial(&self, coefficients: Vec<ComplexBall>) -> BallPolynomial {
        UnivariatePolynomial::from_coefficients(&self.field, coefficients, self.variable.clone())
    }
    /// Channels are columns for P and rows for Q (the latter solves -A^T).
    fn matrix(&self, coefficients: &[Vec<C>], transpose: bool) -> Result<PolynomialMatrix> {
        let n = self.weights.len();
        let entries = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| {
                let index = if transpose { i * n + j } else { j * n + i };
                self.polynomial(
                    coefficients
                        .iter()
                        .map(|row| dyadic_ball(self.p, &row[index]))
                        .collect(),
                )
            })
            .collect();
        Matrix::from_linear(entries, n as u32, n as u32, self.polynomial_field.clone())
            .map_err(Error::InvalidInput)
    }
    fn evaluate(&self, matrix: &PolynomialMatrix, point: &ComplexBall) -> BallMatrix {
        matrix.map(|a| a.evaluate(point), self.field.clone())
    }
    fn polynomial_upper(&self, a: &BallPolynomial, radius: &Float) -> Float {
        let mut power = self.p.real(1);
        let mut sum = self.p.real(0);
        for coefficient in a.coefficients() {
            sum = self.add(&sum, &self.mul(&complex_upper(coefficient, self.p), &power));
            power = self.mul(&power, radius);
        }
        sum
    }
    fn norm(&self, mut entry: impl FnMut(usize, usize) -> Float) -> Float {
        let n = self.weights.len();
        (0..n)
            .map(|i| {
                (0..n).fold(self.p.real(0), |sum, j| {
                    self.add(
                        &sum,
                        &self.div(&self.mul(&entry(i, j), &self.weights[j]), &self.weights[i]),
                    )
                })
            })
            .reduce(|a, b| if a > b { a } else { b })
            .unwrap_or_else(|| self.p.real(0))
    }
    fn matrix_norm(&self, matrix: &BallMatrix) -> Float {
        self.norm(|i, j| complex_upper(&matrix[(i as u32, j as u32)], self.p))
    }
    fn disk_norm(&self, matrix: &PolynomialMatrix, radius: &Float) -> Float {
        self.norm(|i, j| self.polynomial_upper(&matrix[(i as u32, j as u32)], radius))
    }
    fn center_channels(&self, matrix: &BallMatrix, transpose: bool) -> Vec<C> {
        let n = self.weights.len();
        (0..n)
            .flat_map(|channel| (0..n).map(move |component| (channel, component)))
            .map(|(channel, component)| {
                let (i, j) = if transpose {
                    (channel, component)
                } else {
                    (component, channel)
                };
                let value = &matrix[(i as u32, j as u32)];
                C::new(value.re.center.clone(), value.im.center.clone())
            })
            .collect()
    }
    fn propagated_errors(
        &self,
        matrix: &BallMatrix,
        input: &[Float],
        theta: &Float,
    ) -> Result<Vec<Float>> {
        let rho = input
            .iter()
            .zip(&self.weights)
            .map(|(error, weight)| self.div(error, weight))
            .reduce(|a, b| if a > b { a } else { b })
            .unwrap_or_else(|| self.p.real(0));
        let correction = self.mul(&rho, &self.div(theta, &self.complement(theta)));
        let n = input.len();
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let mut direct = self.p.real(0);
            let mut scaled_row = self.p.real(0);
            for (j, error) in input.iter().enumerate() {
                let magnitude = complex_upper(&matrix[(i as u32, j as u32)], self.p);
                direct = self.add(&direct, &self.mul(&magnitude, error));
                scaled_row = self.add(&scaled_row, &self.mul(&magnitude, &self.weights[j]));
            }
            let bound = self.add(&direct, &self.mul(&scaled_row, &correction));
            if !bound.is_finite() {
                return Err(Error::Accuracy(
                    "nonfinite fundamental source-error bound".into(),
                ));
            }
            result.push(bound);
        }
        Ok(result)
    }
}

fn preflight(
    system: &CompiledSystem,
    order: usize,
    charts: usize,
    context: &RunContext,
) -> Result<()> {
    let n = system.dimension();
    if n == 0 || n > MAX_DIMENSION || order > MAX_ORDER || system.p.bits > MAX_BITS {
        return Err(Error::Limit(format!(
            "fundamental proof admits dimension <= {MAX_DIMENSION}, order <= {MAX_ORDER}, precision <= {MAX_BITS} bits"
        )));
    }
    let n = n as u64;
    let order = order as u64 + 1;
    let source_degree = system
        .exact_source_rows
        .iter()
        .flat_map(|row| {
            std::iter::once(row.denominator.len()).chain(row.entries.iter().map(|(_, a)| a.len()))
        })
        .max()
        .unwrap_or(1) as u64;
    // Transposition can combine n distinct original row denominators.
    let source_degree = source_degree.saturating_mul(n);
    let guard_cells = (system.pole_polynomials.len() as u64).checked_mul(MAX_ORDER as u64 + 1);
    let coefficients = n
        .checked_mul(n)
        .and_then(|a| a.checked_mul(8 * order + 4 * source_degree))
        .and_then(|a| a.checked_add(guard_cells?.checked_mul(4)?));
    let work = n
        .checked_pow(3)
        .and_then(|a| a.checked_mul(order))
        .and_then(|a| a.checked_mul(order + source_degree))
        .and_then(|a| a.checked_mul(8))
        .and_then(|a| a.checked_add(guard_cells?.checked_mul(4 * (MAX_ORDER as u64 + 1))?))
        .and_then(|a| a.checked_mul(charts as u64));
    // Four MPFR values per rectangular complex ball, plus headers/limb slack.
    let bytes_per_coefficient =
        std::mem::size_of::<ComplexBall>() as u64 + 4 * (u64::from(system.p.bits).div_ceil(8) + 32);
    let endpoint_cells = n.checked_mul(charts as u64);
    let float_bytes =
        std::mem::size_of::<Float>() as u64 + u64::from(system.p.bits).div_ceil(8) + 32;
    let endpoint_storage = endpoint_cells
        .and_then(|a| a.checked_mul(float_bytes))
        .and_then(|a| {
            a.checked_add((charts as u64).checked_mul(std::mem::size_of::<Vec<Float>>() as u64)?)
        });
    // Include the newly compiled transpose's exact coefficient arrays, whose
    // rational limb storage can exceed that of the working-precision balls.
    let exact_coefficient_bytes = std::mem::size_of::<super::source::Gaussian>() as u64
        + 4 * (MAX_EXACT_BITS.div_ceil(8) + 32);
    let exact_storage = n
        .checked_mul(n)
        .and_then(|a| a.checked_mul(2 * source_degree))
        .and_then(|a| a.checked_mul(exact_coefficient_bytes));
    let storage = coefficients
        .and_then(|a| a.checked_mul(bytes_per_coefficient))
        .and_then(|a| a.checked_add(endpoint_storage?))
        .and_then(|a| a.checked_add(exact_storage?));
    if coefficients.is_none_or(|a| a > MAX_COEFFICIENTS)
        || endpoint_cells.is_none_or(|a| a > MAX_COEFFICIENTS)
        || work.is_none_or(|a| a > MAX_WORK)
        || storage.is_none_or(|a| a > MAX_STORAGE_ESTIMATE)
    {
        return Err(Error::Limit(
            "fundamental proof coefficient/work/storage estimate exhausted".into(),
        ));
    }
    // Reuse the shared metadata-only rational-operation admission before the
    // native transpose compilation and guard conversion. It bounds exponent,
    // source height/depth and polynomial growth without performing arithmetic.
    let limits = crate::frobenius::exact::ExactFrobeniusLimits {
        max_dimension: MAX_DIMENSION,
        max_order: MAX_ORDER,
        max_coefficient_bits: MAX_EXACT_BITS,
        max_scalar_cells: MAX_COEFFICIENTS as usize,
    };
    let mut source_bytes = 0_usize;
    if system.pole_polynomials.len() > 256 || !system.source.values.is_empty() {
        return Err(Error::Unsupported(
            "fundamental proof requires a specialized rational source and at most 256 guards"
                .into(),
        ));
    }
    for expression in system
        .source
        .rows
        .iter()
        .flatten()
        .chain(&system.pole_polynomials)
    {
        context.cancellation.check()?;
        source_bytes = source_bytes.saturating_add(expression.as_view().get_byte_size());
        if source_bytes > MAX_COEFFICIENTS as usize {
            return Err(Error::Limit(
                "fundamental source expression budget exhausted".into(),
            ));
        }
        crate::frobenius::exact::preflight_operation(
            &[expression],
            system.source.variable,
            1,
            1,
            &limits,
            context,
        )?;
    }
    for column in 0..system.dimension() {
        context.cancellation.check()?;
        let expressions = system
            .source
            .rows
            .iter()
            .map(|row| &row[column])
            .collect::<Vec<_>>();
        crate::frobenius::exact::preflight_operation(
            &expressions,
            system.source.variable,
            1,
            system.dimension(),
            &limits,
            context,
        )?;
    }
    Ok(())
}

/// Compute all endpoint contributions transactionally. Optional callers may
/// retain the existing scalar estimate on Limit/Unsupported/Accuracy, while
/// cancellation and malformed inputs remain typed failures.
pub(crate) fn propagate(
    system: &CompiledSystem,
    segments: &[TaylorSegment],
    input_errors: &[Float],
    weights: &[Float],
    order: usize,
    context: &RunContext,
) -> Result<FundamentalErrors> {
    context.cancellation.check()?;
    let p = system.p;
    let n = system.dimension();
    if weights.len() != n
        || input_errors.len() != n
        || order == 0
        || weights.iter().any(|a| !a.is_finite() || *a <= p.real(0))
        || input_errors
            .iter()
            .any(|a| !a.is_finite() || *a < p.real(0))
        || segments
            .iter()
            .any(|segment| !p.finite(&segment.center) || !p.finite(&segment.end))
        || segments
            .windows(2)
            .any(|pair| pair[0].end != pair[1].center)
        || system.source.rows.iter().any(|row| row.len() != n)
    {
        return Err(Error::InvalidInput(
            "fundamental source-error dimensions, path or bounds".into(),
        ));
    }
    preflight(system, order, segments.len(), context)?;
    context.emit(crate::Progress::Stage {
        name: "propagate supplied errors with signed fundamental matrices".into(),
    })?;
    let arithmetic = Arithmetic::new(p, weights);
    let inverse_rows = (0..n)
        .map(|i| (0..n).map(|j| -&system.source.rows[j][i]).collect())
        .collect::<Vec<_>>();
    context.cancellation.check()?;
    let inverse =
        super::prepared::PreparedRows::new(system.source.variable, &inverse_rows, context)?
            .compile(p, &system.source.values, context)?;
    context.cancellation.check()?;
    for row in &inverse.exact_source_rows {
        for polynomial in std::iter::once(&row.denominator)
            .chain(row.entries.iter().map(|(_, a)| a))
            .chain(row.denominator_factors.iter().map(|(a, _)| a))
        {
            context.cancellation.check()?;
            if polynomial.len() > MAX_ORDER + 1
                || polynomial
                    .iter()
                    .any(|a| crate::frobenius::exact::height(a) > MAX_EXACT_BITS)
            {
                return Err(Error::Limit(
                    "fundamental observed transpose degree/height exceeds admission".into(),
                ));
            }
        }
    }
    let specialization = ExactSpecialization::new(p, &system.source.values)?;
    let guards = system
        .pole_polynomials
        .iter()
        .map(|a| {
            context.cancellation.check()?;
            let coefficients = specialization.polynomial(a, system.source.variable)?;
            if coefficients.len() > MAX_ORDER + 1 {
                return Err(Error::Limit(
                    "fundamental source guard degree exceeded".into(),
                ));
            }
            Ok(arithmetic.polynomial(coefficients.iter().map(|a| exact_ball(a, p)).collect()))
        })
        .collect::<Result<Vec<_>>>()?;
    let identity = Matrix::identity(n as u32, arithmetic.field.clone());
    let mut left_endpoint = identity.clone();
    let mut p_initial = arithmetic.center_channels(&identity, false);
    let mut q_initial = arithmetic.center_channels(&identity, true);
    let mut theta = p.real(0);
    #[cfg(test)]
    let mut jump_budget = p.real(0);
    let mut endpoints = Vec::with_capacity(segments.len());
    for segment in segments {
        context.cancellation.check()?;
        // Native subtraction encloses the EXACT difference of stored endpoints.
        let delta = dyadic_ball(p, &segment.end) - dyadic_ball(p, &segment.center);
        let radius = complex_upper(&delta, p);
        if !radius.is_finite() {
            return Err(Error::Accuracy("nonfinite fundamental chart radius".into()));
        }
        let center_ball = dyadic_ball(p, &segment.center);
        for guard in &guards {
            context.cancellation.check()?;
            if polynomial_lower(&guard.shift_var(&center_ball), &radius, p) <= p.real(0) {
                return Err(Error::Accuracy(
                    "fundamental proof disk meets an original source exclusion".into(),
                ));
            }
        }
        let p_coefficients = system.taylor_channels(&segment.center, &p_initial, order, n)?;
        context.cancellation.check()?;
        let q_coefficients = inverse.taylor_channels(&segment.center, &q_initial, order, n)?;
        context.cancellation.check()?;
        let candidate = arithmetic.matrix(&p_coefficients, false)?;
        let inverse_candidate = arithmetic.matrix(&q_coefficients, true)?;
        let identity_polynomial = Matrix::identity(n as u32, arithmetic.polynomial_field.clone());
        let inverse_defect = &identity_polynomial - &(&inverse_candidate * &candidate);
        context.cancellation.check()?;
        let eta = arithmetic.disk_norm(&inverse_defect, &radius);
        if !eta.is_finite() || eta >= p.real(1) / p.real(2) {
            return Err(Error::Accuracy(
                "fundamental whole-chart inverse is inconclusive".into(),
            ));
        }
        let denominator = arithmetic.complement(&eta);
        let q_upper = arithmetic.disk_norm(&inverse_candidate, &radius);
        let candidate_start = arithmetic.evaluate(&candidate, &dyadic_ball(p, &p.zero()));
        let difference = &left_endpoint - &candidate_start;
        let jump = arithmetic.div(
            &arithmetic.mul(&q_upper, &arithmetic.matrix_norm(&difference)),
            &denominator,
        );
        #[cfg(test)]
        {
            jump_budget = arithmetic.add(&jump_budget, &jump);
        }
        let residual = ExactSourceResidual::new(
            p,
            &system.exact_source_rows,
            &segment.center,
            &p_coefficients,
            n,
        )?;
        context.cancellation.check()?;
        let defects = residual.defect_bounds(p, &C::new(radius.clone(), p.real(0)), None)?;
        let defect_norm = arithmetic.norm(|i, j| defects[j * n + i].clone());
        let local = arithmetic.div(&arithmetic.mul(&q_upper, &defect_norm), &denominator);
        theta = arithmetic.add(&theta, &arithmetic.add(&local, &jump));
        if !theta.is_finite() || theta > p.real(1) / p.real(4) {
            return Err(Error::Accuracy(
                "fundamental cumulative defect/jump budget is inconclusive".into(),
            ));
        }
        left_endpoint = arithmetic.evaluate(&candidate, &delta);
        endpoints.push(arithmetic.propagated_errors(&left_endpoint, input_errors, &theta)?);
        p_initial = arithmetic.center_channels(&left_endpoint, false);
        q_initial =
            arithmetic.center_channels(&arithmetic.evaluate(&inverse_candidate, &delta), true);
        context.cancellation.check()?;
    }
    Ok(FundamentalErrors {
        endpoints,
        #[cfg(test)]
        jump_budget,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::local_coordinates::TaylorCoordinate;
    use crate::ode::compile_rows;

    fn system(input: &[&[&str]], p: Precision) -> CompiledSystem {
        let rows = input
            .iter()
            .map(|row| {
                row.iter()
                    .map(|a| Atom::parse(a, "fundamental_tests", Default::default()).unwrap())
                    .collect()
            })
            .collect::<Vec<_>>();
        compile_rows(
            symbol!("fundamental_tests::x"),
            &rows,
            p,
            &Default::default(),
        )
        .unwrap()
    }
    fn segment(p: Precision, start: C, end: C) -> TaylorSegment {
        TaylorSegment {
            center: start,
            end,
            coefficients: Vec::new(),
            pade: None,
            working_bits: p.bits,
            coordinate: TaylorCoordinate::Identity,
            conditioning_digits: 20,
        }
    }
    #[test]
    fn nilpotent_signed_mixing_retains_polynomial_error_growth() {
        for bits in [160, 282, 400] {
            let p = Precision { bits };
            let system = system(&[&["5800", "5800"], &["-5800", "-5800"]], p);
            let error = p.tolerance(30);
            let result = propagate(
                &system,
                &[segment(p, p.zero(), p.i(1))],
                &[error.clone(), error.clone()],
                &[p.real(1), p.real(1)],
                112,
                &RunContext::default(),
            )
            .unwrap();
            for (bound, factor) in result.endpoints[0].iter().zip([11601, 11599]) {
                let exact = error.clone() * p.real(factor);
                assert!(*bound >= exact);
                assert!(*bound < exact.clone() * p.real(2));
            }
        }
    }
    #[test]
    fn noncommuting_complex_polynomial_fundamental_retains_chart_jumps() {
        let p = Precision { bits: 201 };
        // Exact F=[[1+x^2,x],[x,1]], det(F)=1. A(x) does not commute at different x.
        let system = system(&[&["x", "1-x^2"], &["1", "-x"]], p);
        let middle = p.rational(&Rational::from((1, 3)));
        let endpoint = p.complex(1, 1);
        let segments = [
            segment(p, p.zero(), middle.clone()),
            segment(p, middle, endpoint.clone()),
        ];
        let error = p.tolerance(30);
        let result = propagate(
            &system,
            &segments,
            &[error.clone(), error.clone()],
            &[p.real(1), p.real(1)],
            24,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(result.endpoints.len(), 2);
        assert!(
            result.jump_budget > p.real(0),
            "rounded chart reinitialization must be enclosed"
        );
        // The arithmetic-only complex one-norm is 3+2=5 in row0 and 2+1=3 in row1.
        for (bound, factor) in result.endpoints[1].iter().zip([5, 3]) {
            assert!(*bound >= error.clone() * p.real(factor));
            assert!(*bound < error.clone() * p.real(factor + 1));
        }
    }
    #[test]
    fn genuinely_exponential_errors_are_not_suppressed() {
        let p = Precision { bits: 282 };
        let system = system(&[&["20"]], p);
        let error = p.tolerance(30);
        let result = propagate(
            &system,
            &[segment(p, p.zero(), p.i(1))],
            std::slice::from_ref(&error),
            &[p.real(1)],
            112,
            &RunContext::default(),
        )
        .unwrap();
        let expected = p.exp(&p.i(20)).re * error;
        assert!(result.endpoints[0][0] >= expected);
        assert!(result.endpoints[0][0] < expected * p.real(2));
    }
    #[test]
    fn source_holes_limits_invalid_data_and_cancellation_are_explicit() {
        let p = Precision { bits: 160 };
        let mut system = system(&[&["0"]], p);
        system
            .exclude_polynomials(&[
                Atom::parse("x-1/2", "fundamental_tests", Default::default()).unwrap(),
            ])
            .unwrap();
        let segments = [segment(p, p.zero(), p.i(1))];
        assert!(matches!(
            propagate(
                &system,
                &segments,
                &[p.tolerance(30)],
                &[p.real(1)],
                24,
                &RunContext::default()
            ),
            Err(Error::Accuracy(_))
        ));
        assert!(matches!(
            propagate(
                &system,
                &segments,
                &[],
                &[p.real(1)],
                24,
                &RunContext::default()
            ),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            propagate(
                &system,
                &segments,
                &[p.tolerance(30)],
                &[p.real(1)],
                257,
                &RunContext::default()
            ),
            Err(Error::Limit(_))
        ));
        let context = RunContext::default();
        context.cancellation.cancel();
        assert!(matches!(
            propagate(
                &system,
                &segments,
                &[p.tolerance(30)],
                &[p.real(1)],
                24,
                &context
            ),
            Err(Error::Cancelled)
        ));
    }

    #[test]
    fn insufficient_inverse_order_is_rejected_and_rational_source_refines() {
        let p = Precision { bits: 160 };
        let exponential = system(&[&["20"]], p);
        let failed = propagate(
            &exponential,
            &[segment(p, p.zero(), p.i(1))],
            &[p.tolerance(30)],
            &[p.real(1)],
            8,
            &RunContext::default(),
        );
        assert!(matches!(failed,Err(Error::Accuracy(ref message)) if message.contains("inverse")));
        let mut excesses = Vec::new();
        for bits in [160, 282, 400] {
            let p = Precision { bits };
            let rational = system(&[&["1/3", "1/3"], &["-1/3", "-1/3"]], p);
            let error = p.tolerance(30);
            let result = propagate(
                &rational,
                &[segment(p, p.zero(), p.i(1))],
                &[error.clone(), error.clone()],
                &[p.real(1), p.real(1)],
                24,
                &RunContext::default(),
            )
            .unwrap();
            let excess =
                result.endpoints[0][0].to_rational() / error.to_rational() - Rational::from((5, 3));
            assert!(excess >= Rational::zero());
            excesses.push(excess);
        }
        assert!(excesses[1] < excesses[0] && excesses[2] < excesses[1]);
    }

    #[test]
    fn transpose_growth_is_admitted_before_native_preparation() {
        let p = Precision { bits: 160 };
        let source = system(&[&["x^130", "0"], &["x^130", "0"]], p);
        let result = propagate(
            &source,
            &[segment(p, p.zero(), p.i(1))],
            &[p.tolerance(30), p.tolerance(30)],
            &[p.real(1), p.real(1)],
            24,
            &RunContext::default(),
        );
        assert!(
            matches!(result,Err(Error::Limit(ref message)) if message.contains("rational-operation"))
        );
    }

    #[test]
    fn retained_endpoints_are_charged_at_their_storage_precision() {
        let ordinary = system(&[&["0"]], Precision { bits: 160 });
        let large = system(&[&["0"]], Precision { bits: 4096 });
        assert!(preflight(&ordinary, 1, 600_000, &RunContext::default()).is_ok());
        assert!(matches!(
            preflight(&large, 1, 600_000, &RunContext::default()),
            Err(Error::Limit(_))
        ));
        assert!(matches!(
            preflight(&ordinary, 1, 2_000_000, &RunContext::default()),
            Err(Error::Limit(_))
        ));
    }
}
