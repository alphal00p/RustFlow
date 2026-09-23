use symbolica::prelude::*;
use symbolica_amflow::*;
#[test]
fn isotropic_projection_through_rank_eight() {
    let d = parse!("D");
    let mut projector = tensor::TensorProjector::new(d.clone());
    for (rank, numerator, denominator) in [
        (2, 1, parse!("D")),
        (4, 3, parse!("D*(D+2)")),
        (6, 15, parse!("D*(D+2)*(D+4)")),
        (8, 105, parse!("D*(D+2)*(D+4)*(D+6)")),
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
fn repeated_vector_orbits_preserve_mixed_tensor_contractions() {
    let mut projector = tensor::TensorProjector::new(parse!("D"));
    for labels in [vec![0, 0, 0, 0, 1, 1, 1, 1], vec![0, 1, 0, 1, 0, 1, 0, 1]] {
        let hard = labels
            .iter()
            .map(|&i| {
                labels
                    .iter()
                    .map(|&j| {
                        if i != j {
                            parse!("uv")
                        } else if i == 0 {
                            parse!("u2")
                        } else {
                            parse!("v2")
                        }
                    })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let external = vec![vec![parse!("p2"); 8]; 8];
        let expected = parse!("p2^4*(9*u2^2*v2^2+72*u2*v2*uv^2+24*uv^4)/(D*(D+2)*(D+4)*(D+6))");
        assert!(
            (projector.project(&hard, &external).unwrap() - expected)
                .together()
                .cancel()
                .is_zero()
        );
    }
    let hard = vec![
        vec![parse!("u2"), parse!("u2"), parse!("uv"), parse!("uv")],
        vec![parse!("u2"), parse!("u2"), parse!("uv"), parse!("uv")],
        vec![parse!("uv"), parse!("uv"), parse!("v2"), parse!("v2")],
        vec![parse!("uv"), parse!("uv"), parse!("v2"), parse!("v2")],
    ];
    let external = [[1, 1, 2, 3], [1, 1, 4, 5], [2, 4, 1, 6], [3, 5, 6, 1]]
        .map(|row| row.map(Atom::num).to_vec())
        .to_vec();
    let expected = parse!("(6*((D+1)*u2*v2-2*uv^2)+22*(D*uv^2-u2*v2))/(D*(D-1)*(D+2))");
    assert!(
        (projector.project(&hard, &external).unwrap() - expected)
            .together()
            .cancel()
            .is_zero()
    );
}

#[test]
fn higher_rank_repeated_tensors_use_exact_angular_moments() {
    let dimension = parse!("D");
    let mut projector = tensor::TensorProjector::new(dimension.clone());
    for rank in [10, 16, 32] {
        let hard = vec![vec![parse!("k2"); rank]; rank];
        let external = vec![vec![parse!("p2"); rank]; rank];
        let numerator = (1..rank)
            .step_by(2)
            .fold(Atom::num(1), |a, k| a * Atom::num(k as i64));
        let denominator = (0..rank / 2).fold(Atom::num(1), |a, k| {
            a * (&dimension + Atom::num(2 * k as i64))
        });
        let expected = numerator * parse!("k2*p2").pow((rank / 2) as i64) / denominator;
        assert!(
            (projector.project(&hard, &external).unwrap() - expected)
                .together()
                .cancel()
                .is_zero()
        );
    }
    for (u, v, numerator) in [
        (4, 6, parse!("45*u2^2*v2^3+540*u2*v2^2*uv^2+360*v2*uv^4")),
        (
            6,
            6,
            parse!("225*u2^3*v2^3+4050*u2^2*v2^2*uv^2+5400*u2*v2*uv^4+720*uv^6"),
        ),
    ] {
        let rank = u + v;
        let hard = (0..rank)
            .map(|i| {
                (0..rank)
                    .map(|j| {
                        if (i < u) != (j < u) {
                            parse!("uv")
                        } else if i < u {
                            parse!("u2")
                        } else {
                            parse!("v2")
                        }
                    })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let external = vec![vec![parse!("p2"); rank]; rank];
        let denominator = (0..rank / 2).fold(Atom::num(1), |a, k| {
            a * (&dimension + Atom::num(2 * k as i64))
        });
        let expected = numerator * parse!("p2").pow((rank / 2) as i64) / denominator;
        assert!(
            (projector.project(&hard, &external).unwrap() - expected)
                .together()
                .cancel()
                .is_zero()
        );
    }
    let null = vec![vec![Atom::new(); 16]; 16];
    assert!(
        projector
            .project(&vec![vec![parse!("k2"); 16]; 16], &null)
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
