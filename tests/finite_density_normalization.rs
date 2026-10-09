use ahash::HashMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{compact::CompactShell, normalization::*};
use symbolica_amflow::finite_density::boundary::{
    IntegratedOccupiedBoundary, OccupiedBoundaryDistribution, OccupiedBoundaryLimits,
};
use symbolica_amflow::*;

fn equal(left: Atom, right: Atom) {
    assert!((left - right).together().cancel().is_zero());
}

#[test]
fn small_graph_cut_phases_and_medium_wick_factors_are_applied_once() {
    let d = Atom::num(4);
    let pi = Atom::var(Symbol::PI);
    let virtual_measure = native_measure_to_euclidean(2, 0, &d).unwrap();
    for cuts in 0..=2 {
        equal(
            native_measure_to_euclidean(2, cuts, &d).unwrap(),
            &virtual_measure * (Atom::num(2) * &pi).pow(cuts as i64),
        );
    }
    // Sunset simple powers have the SAME total quadratic phase for zero,
    // one, and two occupied lines. Adding a separate (-1)^cuts is wrong.
    assert_eq!(quadratic_index_phase(&[1, 1, 1]), Atom::num(-1));
    assert_eq!(quadratic_index_phase(&[2, 1, 1]), Atom::num(1));
    assert_eq!(quadratic_index_phase(&[2, 1, 1, -1]), Atom::num(-1));
    assert!(native_measure_to_euclidean(2, 3, &d).is_err());
    let g = parse!("normalization_g12");
    let u = parse!("normalization_u1");
    let v = parse!("normalization_u2");
    let chemical = parse!("normalization_mu");
    let n = &u * &g + v.clone().pow(3) + &chemical * &u;
    equal(
        wick_polynomial(&n, std::slice::from_ref(&g), &[u.clone(), v.clone()]).unwrap(),
        -Atom::i() * &u * &g - Atom::i() * v.clone().pow(3) + Atom::i() * chemical * &u,
    );
    assert!(wick_polynomial(&(Atom::num(1) / &u), &[], std::slice::from_ref(&u)).is_err());
}

#[test]
fn unexpanded_msbar_conversion_matches_every_cut_multiplicity() {
    let p = Precision::decimal(60).unwrap();
    let epsilon = Atom::num((1, 7));
    let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
    let lambda = Atom::num((5, 3));
    for loops in 1..=4 {
        for cuts in 0..=loops {
            let direct = native_measure_to_euclidean(loops, cuts, &dimension).unwrap()
                * msbar_measure(loops, &epsilon, &lambda).unwrap()
                / scale_factored_normalization(loops, &epsilon, &lambda).unwrap();
            let reduced = native_measure_to_scale_factored_msbar(loops, cuts, &epsilon).unwrap();
            let actual = p.eval(&direct, &HashMap::default()).unwrap();
            let expected = p.eval(&reduced, &HashMap::default()).unwrap();
            assert!(
                p.close(&actual, &expected, 50),
                "L={loops}, cuts={cuts}: {actual} != {expected}"
            );
        }
    }
}

#[test]
fn massive_one_loop_vacuum_and_occupied_seeds_assemble_in_euclidean_units() {
    let mut results = Vec::new();
    for digits in [40, 60] {
        let p = Precision::decimal(digits + 20).unwrap();
        let epsilon = Rational::from((1, 2));
        let family = IntegralFamily {
            name: "finite_density_normalization_massive_tadpole".into(),
            loops: vec!["k".into()],
            external: vec![],
            external_gram: vec![],
            propagators: vec![Propagator::quadratic(&[1], &[], Atom::num((1, 4)), &[]).unwrap()],
            physical_propagators: 1,
            epsilon: symbol!("normalization_epsilon"),
            dimension: 4,
        };
        let backend = RustRedBackend {
            bubble_subloops: false,
            ..Default::default()
        };
        let options = FlowOptions {
            digits,
            ..Default::default()
        };
        let context = RunContext::default();
        let flow = PreparedFlow::new(
            &family,
            &[Integral(vec![1]), Integral(vec![2])],
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )
        .unwrap();
        let boundary = recursive::RecursiveBoundary::new(&backend, &options, &context)
            .with_terminal_policy(recursive::RecursiveTerminalPolicy::TadpolesOnly);
        let native = flow
            .evaluate(&epsilon, &options, &boundary, &context)
            .unwrap();
        let raw = native_measure_to_euclidean(1, 0, &Atom::num(3)).unwrap();
        let measure = p.eval(&raw, &HashMap::default()).unwrap();
        let vacuum = p.neg(&p.mul(&measure, &native[0]));
        let raised_vacuum = p.mul(&measure, &native[1]);
        let pi = p.eval(&Atom::var(Symbol::PI), &HashMap::default()).unwrap();
        let expected_vacuum = p.neg(&p.div(&p.i(1), &p.scale(&pi, 8, 1)));
        assert!(p.close(&vacuum, &expected_vacuum, 30));
        assert!(p.norm(&vacuum) > p.tolerance(10));
        let mut values = Vec::new();
        for mu in [
            Rational::from((1, 4)),
            Rational::from((1, 2)),
            Rational::from(1),
        ] {
            let shell = CompactShell {
                mass_squared: Rational::from((1, 4)),
                chemical_potential: mu.clone(),
            };
            let density = shell
                .raised_moment(&Rational::from(2), 1, 0, 0, digits, 10000, p)
                .unwrap();
            let assembled = p.add(&vacuum, &density);
            // Independently integrate the regulated T->0 energy residue in
            // d=2 spatial dimensions: the finite part is -max(m,mu)/(4*pi).
            let lower = if mu < Rational::from((1, 2)) {
                Rational::from((1, 2))
            } else {
                mu.clone()
            };
            let expected = p.neg(&p.div(&p.rational(&lower), &p.scale(&pi, 4, 1)));
            assert!(p.close(&assembled, &expected, 30));
            if mu > Rational::from((1, 2)) {
                let raised_density = shell
                    .raised_moment(&Rational::from(2), 2, 0, 0, digits, 10000, p)
                    .unwrap();
                assert!(p.norm(&p.add(&raised_vacuum, &raised_density)) < p.tolerance(30));
                // Independent native delta-shell angular integral in D=3 is
                // (mu-m)/sqrt(pi). This checks the 2*pi cut measure factor.
                let native_shell = p.div(
                    &p.rational(&(mu - Rational::from((1, 2)))),
                    &p.pow(&pi, &p.rational(&Rational::from((1, 2)))),
                );
                let cut_factor = p
                    .eval(
                        &native_measure_to_euclidean(1, 1, &Atom::num(3)).unwrap(),
                        &HashMap::default(),
                    )
                    .unwrap();
                let from_native = p.neg(&p.mul(&cut_factor, &native_shell));
                assert!(p.close(&density, &from_native, 30));
            }
            values.push(assembled);
        }
        results.push(values);
    }
    let p = Precision::decimal(80).unwrap();
    for (baseline, refined) in results[0].iter().zip(&results[1]) {
        assert!(p.close(baseline, refined, 30));
    }
}

#[test]
fn massive_sunset_integrated_boundaries_convert_to_euclidean_seed_products() {
    let gram = vec![vec![Atom::num(1)]];
    let family = IntegralFamily {
        name: "normalization_occupied_sunset".into(),
        loops: vec!["q1".into(), "q2".into()],
        external: vec!["u".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[0], Atom::num((1, 4)), &gram).unwrap(),
            Propagator::quadratic(&[0, 1], &[0], Atom::num((1, 4)), &gram).unwrap(),
            Propagator::quadratic(&[1, -1], &[0], Atom::num(1), &gram).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("normalization_sunset_epsilon"),
        dimension: 4,
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let owner = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let p = Precision::decimal(60).unwrap();
    let shell = CompactShell {
        mass_squared: Rational::from((1, 4)),
        chemical_potential: Rational::from(1),
    };
    let occupied_euclidean = shell
        .raised_moment(&Rational::from(2), 1, 0, 0, 40, 10000, p)
        .unwrap();
    let pi = p.eval(&Atom::var(Symbol::PI), &HashMap::default()).unwrap();
    // The hard unit-mass Gaussian seed with squared denominator in D=3.
    // This independent scalar value is used only to check the native owner's
    // recursively integrated coefficient and whole-measure conversion.
    let hard_euclidean = p.div(&p.i(1), &p.scale(&pi, 8, 1));
    let epsilon = Rational::from((1, 2));
    for (cuts, indices, shifted, hard, expected) in [
        (
            1,
            vec![0, 1, 1],
            vec![false, true, true],
            vec![false, true],
            p.mul(&occupied_euclidean, &hard_euclidean),
        ),
        (
            2,
            vec![0, 0, 1],
            vec![false, false, true],
            vec![false, false],
            p.mul(&occupied_euclidean, &occupied_euclidean),
        ),
    ] {
        let region = regions::LoopRegion {
            transformation: vec![
                vec![Atom::num(1), Atom::num(0)],
                vec![Atom::num(0), Atom::num(1)],
            ],
            hard: hard.clone(),
            hard_branches: vec![],
            jacobian_determinant: Atom::num(1),
        };
        let expansion = regions::expand_region(
            &family,
            &Integral(indices),
            &shifted,
            &region,
            0,
        )
        .unwrap();
        let projected = integrand::projected_factor_region(
            &expansion.coefficients[0],
            &expansion.coordinates,
            &family,
            &hard,
            100,
        )
        .unwrap();
        let distributions = (0..cuts)
            .map(|source_loop_index| OccupiedBoundaryDistribution {
                source_loop_index,
                shell: shell.clone(),
                cut_index: 1,
                upper_index: 0,
                lower_index: 0,
            })
            .collect::<Vec<_>>();
        let native = owner
            .evaluate_projected(
                &projected,
                &distributions,
                &epsilon,
                &HashMap::default(),
                p,
            )
            .unwrap();
        // All original physical indices, including the required occupied
        // shells, enter the common phase. There is no extra (-1)^cuts.
        let conversion = p
            .eval(
                &(quadratic_index_phase(&[1, 1, 1])
                    * native_measure_to_euclidean(2, cuts, &Atom::num(3)).unwrap()),
                &HashMap::default(),
            )
            .unwrap();
        let physical = p.mul(&conversion, &native);
        assert!(
            p.close(&physical, &expected, 35),
            "sunset cuts={cuts}: {physical} != {expected}"
        );
    }
}
