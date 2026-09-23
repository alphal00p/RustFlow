//! Three-loop single-mass vacuum, using recursive FT boundary evaluation.
use symbolica::prelude::*;
use symbolica_amflow::*;
fn main() -> Result<()> {
    let mut propagators = [
        ([1, 0, 0], 1),
        ([0, 1, 0], 0),
        ([0, 0, 1], 0),
        ([1, 1, 1], 0),
    ]
    .into_iter()
    .map(|(loops, m)| Propagator::quadratic(&loops, &[], Atom::num(m), &[]))
    .collect::<Result<Vec<_>>>()?;
    for coordinate in [1, 2] {
        let mut scalar_products = vec![Atom::new(); 6];
        scalar_products[coordinate] = Atom::num(1);
        propagators.push(Propagator {
            constant: Atom::new(),
            scalar_products,
        });
    }
    let family = IntegralFamily {
        name: "single_mass_banana".into(),
        loops: vec!["l1".into(), "l2".into(), "l3".into()],
        external: vec![],
        external_gram: vec![],
        propagators,
        physical_propagators: 4,
        epsilon: symbol!("banana_eps"),
        dimension: 4,
    };
    let options = FlowOptions::default();
    let backend = cache::CachedBackend {
        backend: RustRedBackend::default(),
        directory: ".amflow-cache".into(),
    };
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|e| {
            if !matches!(e, Progress::Step { .. }) {
                eprintln!("{e:?}");
            }
        })),
        ..Default::default()
    };
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let p = Precision::decimal(60)?;
    let value = provider.evaluate(
        &family,
        &Integral(vec![1, 1, 1, 1, 0, 0]),
        &Rational::from((3, 4)),
        p,
    )?;
    let gamma = |n, d| p.gamma_real(&p.rational(&Rational::from((n, d))).re);
    let expected = p.div(
        &p.mul(&gamma(1, 2)?, &p.powi(&gamma(1, 4)?, 4)),
        &gamma(5, 4)?,
    );
    if !p.close(&value, &expected, 20) {
        return Err(Error::Accuracy(format!("{value} != {expected}")));
    }
    println!("{value}\nIndependent gamma-function comparison agrees to 20 digits.");
    Ok(())
}
