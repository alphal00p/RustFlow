use symbolica::prelude::*;
use symbolica_amflow::*;
#[test]
fn isotropic_projection_through_rank_six() {
    let d = parse!("D");
    let mut projector = tensor::TensorProjector::new(d.clone());
    for (rank, numerator, denominator) in [
        (2, 1, parse!("D")),
        (4, 3, parse!("D*(D+2)")),
        (6, 15, parse!("D*(D+2)*(D+4)")),
    ] {
        let hard = vec![vec![parse!("k2"); rank]; rank];
        let external = vec![vec![parse!("p2"); rank]; rank];
        let expected =
            Atom::num(numerator) * (parse!("k2*p2")).pow((rank / 2) as i64) / denominator;
        assert!(
            (projector.project(&hard, &external).unwrap() - expected)
                .together()
                .cancel()
                .is_zero()
        );
    }
    assert!(
        projector
            .project(&[vec![parse!("k2")]], &[vec![parse!("p2")]])
            .unwrap()
            .is_zero()
    );
}
#[test]
fn analytic_vacuum_scaling_and_mass_derivative() {
    let p = Precision::decimal(70).unwrap();
    let dimension = Rational::from((18, 5));
    let one = vacuum::single_mass_sunset([1, 1, 1], &p.i(1), &dimension, p).unwrap();
    let two = vacuum::single_mass_sunset([2, 1, 1], &p.i(1), &dimension, p).unwrap();
    assert!(p.close(
        &two,
        &p.mul(&one, &p.rational(&(&dimension - &Rational::from(3)))),
        55
    ));
    let scaled = vacuum::single_mass_sunset([1, 1, 1], &p.i(4), &dimension, p).unwrap();
    assert!(p.close(
        &scaled,
        &p.mul(
            &one,
            &p.pow(&p.i(4), &p.rational(&(&dimension - &Rational::from(3))))
        ),
        55
    ));
    let tadpole = vacuum::tadpole(1, &p.i(1), &dimension, p).unwrap();
    let dotted = vacuum::tadpole(2, &p.i(1), &dimension, p).unwrap();
    assert!(p.close(
        &dotted,
        &p.mul(
            &tadpole,
            &p.rational(&(&dimension / &Rational::from(2) - Rational::from(1)))
        ),
        55
    ));
}
