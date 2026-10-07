//! Exact sufficient certificate for an auxiliary-mass path on the positive axis.
use super::*;

fn nonnegative_rational(value: &Atom) -> bool {
    let AtomView::Num(number) = value.as_view() else {
        return false;
    };
    let symbolica::coefficient::Coefficient::Complex(value) = number.get_coeff_view().to_owned()
    else {
        return false;
    };
    value.im.is_zero() && value.re >= Rational::zero()
}

/// Nonnegative coefficients of U and the Euclidean F certify that positive
/// auxiliary mass stays on the Euclidean sheet. Test the complete physical
/// sector, including every pinched subsector; ISPs are never denominators.
/// Unknown, complex, or non-Euclidean coefficients use the causal complex path.
pub(super) fn positive_mass_contour(family: &IntegralFamily, eta: Symbol) -> Result<bool> {
    // Avoid adding a large symbolic preparation to unrelated high-loop inputs.
    if family.loops.len() > 2 || family.physical_propagators > 8 {
        return Ok(false);
    }
    let converted = family.convert()?;
    let Ok(polynomials) =
        rustred::family::symanzik::SymanzikPolynomials::try_from_family_with_limits(
            &converted.family,
            Default::default(),
        )
    else {
        return Ok(false);
    };
    let at_zero = BTreeMap::from([(Atom::var(eta), Atom::num(0))]);
    let at_one = BTreeMap::from([(Atom::var(eta), Atom::num(1))]);
    for (polynomial, sign) in [(polynomials.u(), 1), (polynomials.f(), -1)] {
        for (coefficient, powers) in polynomial.terms() {
            if powers[family.physical_propagators..]
                .iter()
                .any(|&power| power != 0)
            {
                continue;
            }
            let expression = Atom::num(sign)
                * crate::family::substitute(&coefficient.to_expression(), &converted.reverse);
            let constant = crate::family::substitute(&expression, &at_zero)
                .together()
                .cancel();
            let slope = (crate::family::substitute(&expression, &at_one) - &constant)
                .together()
                .cancel();
            if !nonnegative_rational(&constant)
                || !nonnegative_rational(&slope)
                || !(expression - &constant - &slope * Atom::var(eta))
                    .expand()
                    .is_zero()
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bubble(momentum: Atom, mass: Atom) -> IntegralFamily {
        let gram = vec![vec![momentum]];
        IntegralFamily {
            name: "euclidean_contour_test".into(),
            loops: vec!["k".into()],
            external: vec!["p".into()],
            propagators: vec![
                Propagator::quadratic(&[1], &[0], mass.clone(), &gram).unwrap(),
                Propagator::quadratic(&[1], &[-1], mass, &gram).unwrap(),
            ],
            external_gram: gram,
            physical_propagators: 2,
            epsilon: symbol!("contour_test::epsilon"),
            dimension: 4,
        }
    }

    #[test]
    fn exact_certificate_keeps_complex_and_timelike_paths() {
        let eta = symbol!("contour_test::eta");
        for (momentum, mass, expected) in [
            (Atom::num(-2), Atom::num(1), true),
            (Atom::num(-2), Atom::num(0), true),
            (Atom::num(2), Atom::num(1), true),
            (Atom::num(6), Atom::num(1), false),
            (Atom::num(-2), Atom::num(-1), false),
            (
                Atom::num(-2),
                Atom::var(symbol!("contour_test::unknown")),
                false,
            ),
            (
                Atom::num(-2),
                Atom::num(Complex::new(Rational::from(1), Rational::from(1))),
                false,
            ),
        ] {
            let (family, _) = bubble(momentum, mass).deform(eta, &MassMode::All).unwrap();
            assert_eq!(
                positive_mass_contour(&family, eta).unwrap(),
                expected,
                "{family:?}"
            );
        }
    }
}
