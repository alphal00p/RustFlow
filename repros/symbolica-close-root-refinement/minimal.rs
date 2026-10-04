use std::sync::Arc;
use std::time::Instant;
use symbolica::domains::integer::Integer;
use symbolica::domains::rational::{Q, Rational};
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::PolyVariable;
use symbolica::symbol;

fn main() {
    let n: u64 = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "100".into())
        .parse()
        .unwrap();
    let centered = std::env::args().nth(2).as_deref() == Some("centered");
    let start = Instant::now();
    let small = Rational::from((Integer::from(1), Integer::from(10).pow(n)));
    let coefficients = if centered {
        vec![small, Rational::from(0), Rational::from(1)]
    } else {
        vec![
            Rational::from((1, 4)) + small,
            Rational::from(-1),
            Rational::from(1),
        ]
    };
    let p = UnivariatePolynomial::from_coefficients(
        &Q,
        coefficients,
        Arc::new(PolyVariable::Symbol(symbol!("close_root_mre::x"))),
    );
    eprintln!(
        "begin isolate_roots: exponent={n}, centered={centered}, construction={:?}",
        start.elapsed()
    );
    let mut roots = p.isolate_roots();
    eprintln!(
        "finished isolate_roots: elapsed={:?}, roots={}",
        start.elapsed(),
        roots.len()
    );
    let bits: u64 = std::env::args()
        .nth(3)
        .unwrap_or_else(|| "232".into())
        .parse()
        .unwrap();
    let tolerance = Rational::from((Integer::from(1), Integer::from(2).pow(bits)));
    for (index, (root, _)) in roots.iter_mut().enumerate() {
        eprintln!(
            "refining root={index} bits={bits} elapsed={:?}",
            start.elapsed()
        );
        *root = root.clone().refined(&tolerance);
        eprintln!("refined root={index} elapsed={:?}", start.elapsed());
        let location = root.classify_location();
        eprintln!(
            "classified root={index} location={location:?} elapsed={:?}",
            start.elapsed()
        );
    }
}
