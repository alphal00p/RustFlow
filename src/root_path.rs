//! Exact registered-square-root sheet transport along a regular real parameter.
//! Polynomial arithmetic, factorization, root isolation and interval refinement
//! remain native Symbolica operations. This module only counts principal-cut
//! intersections and handles the principal convention at negative-axis endpoints.
use crate::{Error, Precision, Result, RunContext};
use std::{cmp::Ordering, collections::BTreeSet, sync::Arc};
use symbolica::domains::algebraic::{AlgebraicExtension, AlgebraicNumber};
use symbolica::domains::float::{FloatField, RealBall};
use symbolica::domains::rational::RationalField;
use symbolica::poly::univariate::{ExactComplexPolynomial, UnivariatePolynomial};
use symbolica::prelude::*;

type RealPolynomial = UnivariatePolynomial<RationalField>;

pub(crate) fn exact_complex_rational(a: &Atom) -> Result<()> {
    let value = a.together().cancel();
    if let AtomView::Num(n) = value.as_view()
        && matches!(
            n.get_coeff_view().to_owned(),
            symbolica::coefficient::Coefficient::Complex(_)
        )
    {
        return Ok(());
    }
    Err(Error::Unsupported("algebraic cache coordinates and specialized radicands must be exact rational-complex constants".into()))
}

fn complex_coefficient(c: &AlgebraicNumber<RationalField>) -> Complex<Rational> {
    let mut value = Complex::new(Rational::zero(), Rational::zero());
    // Gaussian extension elements are natively reduced to degree < 2.
    for term in c.poly() {
        if term.exponents.iter().all(|&e| e == 0) {
            value.re = term.coefficient.clone();
        } else {
            debug_assert_eq!(term.exponents, &[1]);
            value.im = term.coefficient.clone();
        }
    }
    value
}

pub(crate) fn rational_parts(a: &Atom, parameter: Symbol) -> Result<[ExactComplexPolynomial; 2]> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(a.as_view(), &mut symbols)?;
    // scalar_symbols marks native complex coefficients with the reducer's
    // formal imaginary symbol; the native Gaussian field absorbs that marker.
    symbols.remove(&Atom::var(crate::family::imaginary_parameter()));
    if !symbols.is_subset(&BTreeSet::from([Atom::var(parameter)])) {
        return Err(Error::Unsupported(
            "root paths require rational-complex coefficients in one parameter".into(),
        ));
    }
    let gaussian = AlgebraicExtension::complex(Q);
    let variables = Arc::new(vec![PolyVariable::Symbol(parameter)]);
    let fraction: RationalPolynomial<AlgebraicExtension<RationalField>, u16> = a
        .try_to_rational_polynomial(&gaussian, &gaussian, Some(variables))
        .map_err(|e| Error::Unsupported(format!("rational root path required: {e}")))?;
    for poly in [&fraction.numerator, &fraction.denominator] {
        if poly
            .variables()
            .iter()
            .any(|variable| *variable != PolyVariable::Symbol(parameter))
        {
            return Err(Error::Unsupported(
                "literal foreign variable in the Gaussian root-path polynomial".into(),
            ));
        }
    }
    Ok([fraction.numerator, fraction.denominator].map(|poly| {
        poly.to_univariate_from_univariate(0).map_coeff(
            complex_coefficient,
            FloatField::from_rep(Complex::new(Rational::one(), Rational::zero())),
        )
    }))
}

fn real_parts(a: &ExactComplexPolynomial) -> [RealPolynomial; 2] {
    [
        a.map_coeff(|c| c.re.clone(), Q),
        a.map_coeff(|c| c.im.clone(), Q),
    ]
}

fn common_zero_polynomial(a: &ExactComplexPolynomial) -> RealPolynomial {
    let [re, im] = real_parts(a);
    if re.is_zero() {
        return im;
    }
    if im.is_zero() {
        return re;
    }
    re.to_multivariate::<u16>()
        .gcd(&im.to_multivariate::<u16>())
        .to_univariate_from_univariate(0)
}

// Returns a strictly interior isolating interval or None for an exterior root.
// Endpoint roots are removed by callers before isolation. No midpoint sampling
// is used to assert a root's location.
fn interior_interval(
    polynomial: &RealPolynomial,
    mut interval: (Rational, Rational),
    context: &RunContext,
) -> Result<Option<(Rational, Rational)>> {
    let mut tolerance = Rational::from((1, 16));
    for _ in 0..64 {
        context.cancellation.check()?;
        if interval.1 <= 0 || interval.0 >= 1 {
            return Ok(None);
        }
        if interval.0 > 0 && interval.1 < 1 {
            return Ok(Some(interval));
        }
        if (&interval.0 + &interval.1).is_zero() {
            return Err(Error::Accuracy(
                "root interval endpoint separation has zero midpoint".into(),
            ));
        }
        interval = polynomial.refine_root_interval(interval, &tolerance);
        tolerance *= Rational::from((1, 16));
    }
    Err(Error::Accuracy(
        "root interval could not be separated from path endpoints".into(),
    ))
}

fn has_unit_interval_zero(polynomial: &RealPolynomial, context: &RunContext) -> Result<bool> {
    if polynomial.is_zero() {
        return Ok(true);
    }
    if polynomial.is_constant() {
        return Ok(false);
    }
    if polynomial.evaluate(&Rational::zero()).is_zero()
        || polynomial.evaluate(&Rational::one()).is_zero()
    {
        return Ok(true);
    }
    context.cancellation.check()?;
    // The native isolator itself is synchronous; cancellation is checked before
    // and after each native call and during bounded certificate refinement.
    let intervals = polynomial.isolate_real_root_intervals();
    context.cancellation.check()?;
    for (lower, upper, _) in intervals {
        if interior_interval(polynomial, (lower, upper), context)?.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

fn negative_at_root(
    real: &RealPolynomial,
    imaginary: &RealPolynomial,
    mut interval: (Rational, Rational),
    p: Precision,
    context: &RunContext,
) -> Result<bool> {
    // Exact native remainder reduces cancellation without altering values at
    // roots of Im(H). The centered native RealBall evaluates the certificate.
    let value_poly = real.rem(imaginary);
    let mut tolerance = Rational::from((1, 16));
    let mut bits = p.bits;
    for _ in 0..64 {
        context.cancellation.check()?;
        let center = (&interval.0 + &interval.1) / Rational::from(2);
        let radius = (&interval.1 - &interval.0) / Rational::from(2);
        let offset = RealBall::from_rational_bounds(&(-radius.clone()), &radius, bits);
        let value = value_poly
            .shift_var(&center)
            .map_coeff(
                |c| RealBall::from_rational_bounds(c, c, bits),
                FloatField::from_rep(RealBall::exact(Float::with_val(bits, 0))),
            )
            .evaluate(&offset);
        if value.is_strictly_negative() {
            return Ok(true);
        }
        if value.is_strictly_positive() {
            return Ok(false);
        }
        if center.is_zero() {
            return Err(Error::Accuracy(
                "principal-cut sign certificate has zero midpoint".into(),
            ));
        }
        interval = imaginary.refine_root_interval(interval, &tolerance);
        tolerance *= Rational::from((1, 16));
        bits = bits
            .checked_add(32)
            .ok_or_else(|| Error::Limit("root-sheet sign precision overflow".into()))?;
    }
    Err(Error::Accuracy(
        "could not certify the real sign at a principal-cut intersection".into(),
    ))
}

fn lower_side_at_endpoint(imaginary: &RealPolynomial, endpoint: i64) -> bool {
    // At an endpoint on the negative axis, principal sqrt uses +i sqrt(|R|).
    // Its interior lower-half-plane limit differs by a sign. Determine that
    // one-sided sign from the first nonzero exact Taylor coefficient.
    imaginary
        .shift_var(&Rational::from(endpoint))
        .coefficients()
        .iter()
        .enumerate()
        .find(|(_, c)| !c.is_zero())
        .is_some_and(|(order, coefficient)| {
            let negative = coefficient < &Rational::zero();
            negative ^ (endpoint == 1 && order % 2 == 1)
        })
}

/// None means a genuine branch point/pole on the path; Some(flip) is an exact
/// transition relative to endpoint principal roots. Unresolved certificates are
/// typed errors, never permission to continue.
pub(crate) fn principal_flip(
    radicand: &Atom,
    parameter: Symbol,
    p: Precision,
    context: &RunContext,
) -> Result<Option<bool>> {
    context.cancellation.check()?;
    if p.bits < 2 {
        return Err(Error::InvalidInput(
            "root-path precision must be at least two bits".into(),
        ));
    }
    let [numerator, denominator] = rational_parts(radicand, parameter)?;
    for part in [&numerator, &denominator] {
        if has_unit_interval_zero(&common_zero_polynomial(part), context)? {
            return Ok(None);
        }
    }
    let conjugate = denominator.map_coeff(
        Complex::conj,
        FloatField::from_rep(Complex::new(Rational::one(), Rational::zero())),
    );
    let product = numerator * &conjugate;
    let [real, imaginary] = real_parts(&product);
    if imaginary.is_zero() {
        return Ok(Some(false));
    }
    let mut flip = false;
    for endpoint in [0, 1] {
        let value = Rational::from(endpoint);
        if imaginary.evaluate(&value).is_zero()
            && real.evaluate(&value).cmp(&Rational::zero()) == Ordering::Less
            && lower_side_at_endpoint(&imaginary, endpoint)
        {
            flip = !flip;
        }
    }
    // Remove endpoint factors before native isolation so no arbitrary broad
    // interval at an exact endpoint can obstruct classification.
    let mut interior = imaginary.clone();
    let x = UnivariatePolynomial::from_coefficients(
        &Q,
        vec![Rational::zero(), Rational::one()],
        imaginary.get_vars(),
    );
    for endpoint in [0, 1] {
        let factor = x.clone() - x.constant(Rational::from(endpoint));
        while interior.evaluate(&Rational::from(endpoint)).is_zero() {
            interior = interior.quot_rem(&factor).0;
        }
    }
    context.cancellation.check()?;
    let intervals = interior.isolate_real_root_intervals();
    context.cancellation.check()?;
    for (lower, upper, multiplicity) in intervals {
        let Some(interval) = interior_interval(&interior, (lower, upper), context)? else {
            continue;
        };
        if multiplicity % 2 == 1 && negative_at_root(&real, &interior, interval, p, context)? {
            flip = !flip;
        }
    }
    Ok(Some(flip))
}

#[cfg(test)]
mod tests {
    use crate::{Error, Precision, Result, RunContext};
    use symbolica::prelude::*;

    fn x() -> Symbol {
        symbol!("complex_path_probe::x")
    }
    fn parse(text: &str) -> Atom {
        Atom::parse(text, "complex_path_probe", Default::default()).unwrap()
    }
    fn flip(expression: &str) -> Result<Option<bool>> {
        super::principal_flip(
            &parse(expression),
            x(),
            Precision::decimal(50)?,
            &RunContext::default(),
        )
    }
    #[test]
    fn exact_cut_crossings_tangencies_and_endpoint_sides() -> Result<()> {
        for (expression, expected) in [
            ("-1+𝑖*(1-2*x)", true),
            ("-1+𝑖*(2*x-1)", true),
            ("-1+𝑖*(x-1/2)^2", false),
            ("-1-𝑖*(x-1/2)^2", false),
            ("-1+𝑖*(x-1/2)^3", true),
            ("-1+𝑖*x", false),
            ("-1-𝑖*x", true),
            ("-1+𝑖*(1-x)", false),
            ("-1-𝑖*(1-x)", true),
            ("-1-𝑖*x*(1-x)", false),
            ("-1", false),
            ("1+𝑖*(1-2*x)", false),
            ("1/(-1+𝑖*(1-2*x))", true),
            ("(-1+𝑖*(1-2*x))/(2+𝑖*x)", true),
        ] {
            assert_eq!(flip(expression)?, Some(expected), "{expression}");
        }
        Ok(())
    }
    #[test]
    fn exact_real_parameter_branch_points_and_poles_rejected() -> Result<()> {
        for expression in ["x", "1-x", "x-1/2", "1/(x-1/2)", "(1+𝑖)*(x-1/2)", "0"] {
            assert_eq!(flip(expression)?, None, "{expression}");
        }
        assert_eq!(flip("x-1/2+𝑖/1000")?, Some(false));
        assert!(matches!(flip("sqrt(x)"), Err(Error::Unsupported(_))));
        assert!(matches!(
            flip("symbolica_amflow::imaginary_unit*x"),
            Err(Error::Unsupported(_))
        ));
        Ok(())
    }
}
