use symbolica::prelude::*;
use symbolica_amflow::{integrand, *};

fn template(loops: usize, external: usize) -> IntegralFamily {
    IntegralFamily {
        name: "native_boundary_owner".into(),
        loops: (0..loops).map(|i| format!("k{i}")).collect(),
        external: (0..external).map(|i| format!("p{i}")).collect(),
        external_gram: vec![vec![Atom::num(7); external]; external],
        propagators: vec![],
        physical_propagators: 0,
        epsilon: symbol!("boundary_owner::epsilon"),
        dimension: 6,
    }
}

fn reconstruct(terms: &[integrand::IntegralTerm], variables: &[Atom]) -> Atom {
    terms.iter().fold(Atom::new(), |sum, term| {
        let value = term.family.propagators.iter().zip(&term.integral.0).fold(
            term.coefficient.clone(),
            |value, (d, &power)| {
                let denominator = d
                    .scalar_products
                    .iter()
                    .zip(variables)
                    .fold(d.constant.clone(), |sum, (c, x)| sum + c * x);
                value * denominator.pow(-i64::from(power))
            },
        );
        sum + value
    })
}

#[test]
fn native_completion_preserves_interleaved_products_and_exact_family_metadata() {
    let variables = [
        parse!("a"),
        parse!("b"),
        parse!("c"),
        parse!("d"),
        parse!("e"),
    ];
    let family = template(2, 1);
    let expression = parse!("(b+2*d-3*e)^2/((a-1)^2*(c-2)*(a+c-3))");
    let terms = integrand::to_integrals(&expression, &variables, &family, 100).unwrap();
    assert!(!terms.is_empty());
    for term in &terms {
        term.family.validate_integral(&term.integral).unwrap();
        assert_eq!(term.family.loops, family.loops);
        assert_eq!(term.family.external, family.external);
        assert_eq!(term.family.external_gram, family.external_gram);
        assert_eq!(term.family.epsilon, family.epsilon);
        assert_eq!(term.family.dimension, family.dimension);
        assert_eq!(term.family.propagators.len(), variables.len());
    }
    assert!(
        (reconstruct(&terms, &variables) - expression)
            .together()
            .cancel()
            .is_zero()
    );
}

#[test]
fn native_partial_fractions_preserve_symbolic_degeneracy_and_gaussian_numerators() {
    let variables = [parse!("x")];
    let family = template(1, 0);
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    for expression in [
        (parse!("x^3") + &i * parse!("x")) / parse!("(x-1)^2*(x-2)"),
        parse!("(x+1)^2/((x-a)*(x-b))"),
        parse!("(x+1)^2/(x-2)^2"),
        (Atom::one() + &i) * parse!("x^3"),
    ] {
        let terms = integrand::to_integrals(&expression, &variables, &family, 100).unwrap();
        assert!(
            (reconstruct(&terms, &variables) - expression)
                .together()
                .cancel()
                .is_zero()
        );
    }
    let error =
        integrand::to_integrals(&parse!("1/((x-1)*(x-2))"), &variables, &family, 1).unwrap_err();
    assert!(
        matches!(error, Error::Limit(ref message) if message.contains("native boundary partial fractions"))
    );
}

#[test]
fn gaussian_affine_dependencies_remain_exact_without_formal_imaginary_rank() {
    let variables = [parse!("x"), parse!("y")];
    let family = template(1, 1);
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let d1 = &variables[0] + &i * &variables[1] - 1;
    let d2 = -&i * &variables[0] + &variables[1] - 2;
    // The rows are proportional over Q(i), but independent if i is encoded as
    // a free parameter t (their determinant would be 1+t^2).
    for expression in [
        Atom::one() / (&d1 * &d2),
        (&variables[0] + &i * &variables[1]).pow(2) / (d1.pow(2) * &d2),
    ] {
        let terms = integrand::to_integrals(&expression, &variables, &family, 100).unwrap();
        assert!(!terms.is_empty());
        assert!(terms.iter().all(|t| t.family.physical_propagators == 1));
        assert!(
            (reconstruct(&terms, &variables) - expression)
                .together()
                .cancel()
                .is_zero()
        );
    }
    let expression = Atom::one() / ((&variables[0] - &i) * (&variables[0] - 2));
    let terms = integrand::to_integrals(&expression, &variables, &family, 100).unwrap();
    assert!(
        (reconstruct(&terms, &variables) - expression)
            .together()
            .cancel()
            .is_zero()
    );
}

#[test]
fn invalid_boundary_coordinates_return_typed_errors_before_native_conversion() {
    for variables in [vec![parse!("x+y")], vec![], vec![parse!("x"), parse!("x")]] {
        assert!(matches!(
            integrand::to_integrals(&Atom::one(), &variables, &template(1, 0), 100),
            Err(Error::InvalidInput(_))
        ));
    }
}

#[test]
fn internal_denominator_labels_cannot_capture_scalar_parameters() {
    let x = parse!("x");
    let parameter = Atom::var(symbol!("symbolica_amflow::boundary_d_0"));
    let nested = Atom::var(symbol!("symbolica_amflow::boundary_d_0_fresh_1"));
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    for mass in [Atom::num(2), i] {
        let expression = (&parameter + &nested * &x) / (&x - mass);
        let terms =
            integrand::to_integrals(&expression, std::slice::from_ref(&x), &template(1, 0), 100)
                .unwrap();
        assert!(
            (reconstruct(&terms, std::slice::from_ref(&x)) - expression)
                .together()
                .cancel()
                .is_zero()
        );
        assert!(
            terms
                .iter()
                .any(|t| t.coefficient.contains(parameter.as_view()))
        );
    }
}

#[test]
fn internal_momenta_cannot_capture_scalar_parameters_or_coordinates() {
    let variables = [
        Atom::var(symbol!("symbolica_amflow::boundary_native_loop_0")),
        parse!("y"),
    ];
    let scalar = Atom::var(symbol!("symbolica_amflow::boundary_native_loop_0_fresh_1"));
    let mass = Atom::var(symbol!("symbolica_amflow::boundary_native_external_0"));
    let dimension = Atom::var(symbol!("symbolica_amflow::boundary_native_dimension"));
    let expression =
        (&scalar * &variables[1] + dimension) / ((&variables[0] - mass) * (&variables[1] - 3));
    let terms = integrand::to_integrals(&expression, &variables, &template(1, 1), 100).unwrap();
    assert!(
        (reconstruct(&terms, &variables) - expression)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        terms
            .iter()
            .any(|term| term.coefficient.contains(scalar.as_view()))
    );
}

#[test]
fn native_extreme_symbolic_coefficients_stay_scalar_and_reconstruct_exactly() {
    let variables = [parse!("x"), parse!("y")];
    let x = &variables[0];
    let y = &variables[1];
    for q in [
        parse!("(a+b)/10^1000"),
        parse!("10^1000*(a+b)"),
        parse!("a+(a+b)/10^1000"),
    ] {
        let expression = Atom::one() / ((x + q * y - 1) * (x - 2) * (y - 3));
        let terms = integrand::to_integrals(&expression, &variables, &template(1, 1), 100).unwrap();
        assert!(terms.iter().all(|term| {
            variables
                .iter()
                .all(|v| !term.coefficient.contains(v.as_view()))
        }));
        assert!(
            (reconstruct(&terms, &variables) - expression)
                .together()
                .cancel()
                .is_zero()
        );
    }
}
