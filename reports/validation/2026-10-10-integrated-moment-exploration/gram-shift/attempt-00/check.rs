//! Generic vacuum dimension shift, tested using the independent native Gaussian terminal.
use symbolica::prelude::*;
use symbolica_amflow::{gaussian, integrand, *};

fn determinant(a: &[Vec<Atom>]) -> Atom {
    if a.len() == 1 { return a[0][0].clone(); }
    (0..a.len()).fold(Atom::zero(), |sum, j| {
        let minor = a[1..].iter().map(|r| r.iter().enumerate()
            .filter(|(k, _)| *k != j).map(|(_, x)| x.clone()).collect()).collect::<Vec<_>>();
        sum + Atom::num(if j % 2 == 0 { 1 } else { -1 }) * &a[0][j] * determinant(&minor)
    })
}

fn integrate(expr: &Atom, variables: &[Atom], h: usize, dimension: i64, eps: &Rational, p: Precision) -> Result<(ComplexFloat, usize)> {
    let family = IntegralFamily {
        name: "gram_shift_generic_check".into(), loops: (0..h).map(|i| format!("k{i}")).collect(),
        external: vec![], external_gram: vec![], propagators: vec![], physical_propagators: 0,
        epsilon: symbol!("gram_shift_check::epsilon"), dimension,
    };
    let terms = integrand::to_integrals(expr, variables, &family, 200000)?;
    let mut sum = p.zero();
    for term in &terms {
        let c = Rational::try_from(term.coefficient.as_view()).expect("rational coefficient");
        sum = p.add(&sum, &p.mul(&p.rational(&c), &gaussian::terminal(&term.family, &term.integral, eps, p)?));
    }
    Ok((sum, terms.len()))
}

fn main() -> Result<()> {
    let mut results = vec![];
    for digits in [50, 80] {
        let p = Precision::decimal(digits)?;
        for (h, n) in [(1usize, 1usize), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (3, 1), (4, 1)] {
            let mut variables = vec![];
            let mut gram = vec![vec![Atom::zero(); h]; h];
            let mut quadratic = Atom::num(-3);
            for i in 0..h { for j in i..h {
                let x = Atom::var(symbol!(format!("gram_shift_check::g{i}_{j}")));
                variables.push(x.clone()); gram[i][j] = x.clone(); gram[j][i] = x.clone();
                // A = 2 I + 11^T is positive definite and has off-diagonal entries.
                quadratic += Atom::num(if i == j { 3 } else { 2 }) * x;
            }}
            let eps = Rational::from((1, 7));
            let d = Rational::from(6) - &eps * 2;
            let power = 40;
            let denominator = quadratic.pow(-power);
            let mut factor = Rational::one();
            for r in 0..n { for j in 0..h {
                factor *= Rational::from(-2) / (&d + Rational::from((2*r) as i64) - Rational::from(j as i64));
            }}
            let (shifted, _) = integrate(&denominator, &variables, h, 6 + (2*n) as i64, &eps, p)?;
            let (inserted, terms) = integrate(&(determinant(&gram).pow(n as i64) * denominator), &variables, h, 6, &eps, p)?;
            let bridged = p.mul(&p.rational(&factor), &inserted);
            let relative = p.div(&p.sub(&bridged, &shifted), &shifted);
            let error = p.norm(&relative);
            let passed = error < p.parse("1e-35", "0")?.re;
            results.push(serde_json::json!({"loops":h,"shift_count":n,"digits":digits,"base_dimension":"40/7","quadratic":"A=2I+11^T; mass_squared=3; denominator_power=40","native_bridge_factor":factor.to_string(),"converted_terms":terms,"shifted_value":format!("{}",shifted.re),"relative_error":error.to_string(),"passed":passed}));
        }
    }
    let passed = results.iter().all(|r| r["passed"].as_bool() == Some(true));
    let output = serde_json::json!({"scope":"Generic vacuum dimension-shift normalization checks, no physical finite-density prediction", "passed":passed,"cases":results});
    std::fs::write(std::env::args().nth(1).expect("output path"), serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    println!("checks={} passed={passed}",output["cases"].as_array().unwrap().len());
    assert!(passed); Ok(())
}
