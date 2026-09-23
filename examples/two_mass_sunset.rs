//! Recursive AMF (or FT with --ft) for a genuine two-loop vacuum family.
use symbolica::prelude::*;
use symbolica_amflow::*;
fn main() -> Result<()> {
    let family = IntegralFamily {
        name: "two_mass_sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(1), &[])?,
            Propagator::quadratic(&[0, 1], &[], Atom::num(1), &[])?,
            Propagator::quadratic(&[1, 1], &[], Atom::new(), &[])?,
        ],
        physical_propagators: 3,
        epsilon: symbol!("eps"),
        dimension: 4,
    };
    let options = FlowOptions {
        recursion: if std::env::args().any(|s| s == "--ft") {
            RecursionMode::Ft
        } else {
            RecursionMode::Amf
        },
        cache_directory: Some(".amflow-cache".into()),
        ..Default::default()
    };
    let values = solve_integrals(
        &family,
        &[Integral(vec![1, 1, 1])],
        &KinematicPoint::default(),
        0,
        &options,
        &RustRedBackend::default(),
        &RunContext::default(),
    )?;
    for (power, coefficient) in &values[0].coefficients {
        println!("epsilon^{power}: {coefficient}");
    }
    println!("verified digits: {:?}", values[0].verified_digits);
    Ok(())
}
