use symbolica::prelude::*;
use symbolica_amflow::*;

fn main() -> Result<()> {
    let family = IntegralFamily {
        name: "massless_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(1)]],
        propagators: vec![
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            Propagator {
                constant: Atom::num(1),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    let options = FlowOptions::default();
    let result = solve_integrals(
        &family,
        &[Integral(vec![1, 1])],
        &KinematicPoint::default(),
        0,
        &options,
        &RustRedBackend::default(),
        &RunContext::default(),
    )?;
    for (power, value) in &result[0].coefficients {
        println!("eps^{power}: {value}");
    }
    println!(
        "independently verified digits: {:?}",
        result[0].verified_digits
    );
    Ok(())
}
