//! Symbolica-only reproducer: nonzero coefficients disappear above f64 range.
use std::sync::Arc;
use symbolica::prelude::*;

fn exact_collection(expression: &Atom, variable: Symbol) -> Vec<(Atom, Atom)> {
    let field = AtomField {
        statistical_zero_test: false,
        cancel_check_on_division: true,
        custom_normalization: None,
    };
    let marker = PolyVariable::Symbol(variable);
    let polynomial: MultivariatePolynomial<_, i32> = expression
        .try_to_polynomial(&field, Arc::new(vec![marker.clone()]))
        .unwrap();
    let index = polynomial
        .variables()
        .iter()
        .position(|v| *v == marker)
        .unwrap();
    // Both conversion and collection use existing Symbolica polynomial APIs.
    polynomial
        .to_univariate_polynomial_list(index)
        .into_iter()
        .map(|(coefficient, degree)| (Atom::var(variable).pow(degree), coefficient.flatten(false)))
        .collect()
}

fn reconstruct(terms: &[(Atom, Atom)]) -> Atom {
    terms
        .iter()
        .fold(Atom::zero(), |sum, (monomial, coefficient)| {
            sum + monomial * coefficient
        })
}

fn main() {
    // "collect" asserts the failing public API, "exact" the passing control.
    // With no argument, report every case without stopping at the first failure.
    let mode = std::env::args().nth(1).unwrap_or_else(|| "report".into());
    assert!(matches!(mode.as_str(), "report" | "collect" | "exact"));
    let variable = symbol!("coefficient_list_mre::x");
    let x = Atom::var(variable);
    let epsilon = Atom::var(symbol!("coefficient_list_mre::epsilon"));
    for power in [20, 300, 308, 309, 400] {
        let big = Atom::num(Integer::from(10).pow(power));
        let coefficient = &big + &big * &epsilon;
        let expression = &coefficient * &x;
        assert!(!expression.is_zero());
        let public = expression.coefficient_list::<i32>(std::slice::from_ref(&x));
        let exact = exact_collection(&expression, variable);
        let public_correct = (reconstruct(&public) - &expression).expand().is_zero();
        let exact_correct = (reconstruct(&exact) - &expression).expand().is_zero();
        println!(
            "power={power} zero_test={:?} public_terms={} public_exact={public_correct} exact_field_terms={} exact_field_exact={exact_correct}",
            coefficient.zero_test(10, f64::EPSILON),
            public.len(),
            exact.len(),
        );
        if mode == "collect" {
            assert!(
                public_correct,
                "coefficient_list discarded a nonzero coefficient"
            );
        }
        if mode == "exact" {
            assert!(
                exact_correct,
                "exact native polynomial collection changed the expression"
            );
        }
    }
}
