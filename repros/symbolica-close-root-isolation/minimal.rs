use std::{sync::Arc, time::Instant};
use symbolica::{
    domains::{
        integer::Integer,
        rational::{Q, Rational},
    },
    poly::univariate::UnivariatePolynomial,
    prelude::PolyVariable,
    symbol,
};

fn main() {
    let delta = Rational::from((Integer::from(16), Integer::from(10).pow(43)))
        + Rational::from((Integer::from(1), Integer::from(10).pow(68)));
    let polynomial = UnivariatePolynomial::from_coefficients(
        &Q,
        vec![
            Rational::from((1, 4)) + &delta * &delta,
            Rational::from(-1),
            Rational::from(1),
        ],
        Arc::new(PolyVariable::Symbol(symbol!("close_root_isolation_mre::x"))),
    );
    let started = Instant::now();
    eprintln!("begin isolate_roots: (x-1/2)^2 + (16/10^43 + 1/10^68)^2");
    let roots = polynomial.isolate_roots();
    eprintln!(
        "finished: elapsed={:?}, roots={}",
        started.elapsed(),
        roots.len()
    );
    assert_eq!(roots.len(), 2);
}
