//! Isolated exact-arithmetic prototype and comparison with the production path.
//! The production constructor and its precision/rejection policy are unchanged.
use super::*;
use std::time::Instant;

type ExactPolynomial = UnivariatePolynomial<FloatField<Gaussian>>;

fn gaussian(value: &C) -> Gaussian {
    Gaussian::new(value.re.to_rational(), value.im.to_rational())
}

fn exact_constructor(
    p: Precision,
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
) -> ExactSourceResidual {
    exact_constructor_observed(p, rows, center, coefficients, channels, |_, _| {})
}

fn exact_constructor_observed(
    p: Precision,
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
    mut observe: impl FnMut(&ExactPolynomial, &BallPolynomial),
) -> ExactSourceResidual {
    // Test inputs have already passed the production constructor's validation.
    let n = rows.len();
    let size = n * channels;
    let zero = Gaussian::new(Rational::zero(), Rational::zero());
    let ring = FloatField::from_rep(zero);
    let ball_ring = FloatField::from_rep(dyadic_ball(p, &p.zero()));
    let variable = Arc::new(PolyVariable::Temporary(0));
    let polynomial = |a| UnivariatePolynomial::from_coefficients(&ring, a, variable.clone());
    let enclose = |a: &ExactPolynomial| {
        UnivariatePolynomial::from_coefficients(
            &ball_ring,
            a.coefficients().iter().map(|a| exact_ball(a, p)).collect(),
            variable.clone(),
        )
    };
    let solutions = (0..size)
        .map(|i| polynomial(coefficients.iter().map(|r| gaussian(&r[i])).collect()))
        .collect::<Vec<_>>();
    let derivatives = solutions.iter().map(|a| a.derivative()).collect::<Vec<_>>();
    let center = gaussian(center);
    let mut residuals = Vec::with_capacity(size);
    residuals.resize_with(size, || enclose(&polynomial(vec![ring.zero()])));
    let mut denominators = Vec::with_capacity(n);
    for (i, row) in rows.iter().enumerate() {
        let denominator = polynomial(row.denominator.clone()).shift_var(&center);
        let entries = row
            .entries
            .iter()
            .map(|(j, a)| (*j, polynomial(a.clone()).shift_var(&center)))
            .collect::<Vec<_>>();
        for channel in 0..channels {
            let index = channel * n + i;
            let mut numerator = &denominator * &derivatives[index];
            for (column, a) in &entries {
                let shift = column / n;
                if shift <= channel {
                    let source = (channel - shift) * n + column % n;
                    numerator = &numerator - &(a * &solutions[source]);
                }
            }
            let enclosed = enclose(&numerator);
            observe(&numerator, &enclosed);
            residuals[index] = enclosed;
        }
        let expanded = enclose(&denominator);
        observe(&denominator, &expanded);
        denominators.push(DenominatorEnclosure {
            expanded,
            factors: row
                .denominator_factors
                .iter()
                .map(|(a, multiplicity)| {
                    let shifted = polynomial(a.clone()).shift_var(&center);
                    let enclosed = enclose(&shifted);
                    observe(&shifted, &enclosed);
                    (enclosed, *multiplicity)
                })
                .collect(),
        });
    }
    ExactSourceResidual {
        residuals,
        denominators,
    }
}

fn contains(ball: &RealBall, exact: &Rational) -> bool {
    ball.lower_bound().to_rational() <= *exact && *exact <= ball.upper_bound().to_rational()
}

#[test]
fn exact_gaussian_polynomial_defect_encloses_independent_coefficients() {
    let p = Precision { bits: 64 };
    let a = Gaussian::new(Rational::from((3, 7)), Rational::from((2, 9)));
    let one = Gaussian::new(Rational::one(), Rational::zero());
    let rows = vec![ExactPolynomialRow {
        denominator: vec![one],
        denominator_factors: Vec::new(),
        entries: vec![(0, vec![a.clone()])],
    }];
    let coefficients = vec![vec![p.complex(2, -1)], vec![p.parse("0.1", "0.2").unwrap()]];
    let expected = [
        gaussian(&coefficients[1][0]) - &a * gaussian(&coefficients[0][0]),
        -(&a * gaussian(&coefficients[1][0])),
    ];
    let chart = exact_constructor(p, &rows, &p.complex(1, 2), &coefficients, 1);
    for (ball, exact) in chart.residuals[0].coefficients().iter().zip(expected) {
        assert!(contains(&ball.re, &exact.re));
        assert!(contains(&ball.im, &exact.im));
    }
}

#[test]
fn exact_polynomial_constructor_preserves_inconclusive_denominator_rejection() {
    let p = Precision::decimal(60).unwrap();
    let x = symbol!("exact_residual_prototype::x");
    let system = crate::ode::DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::one() / (Atom::one() - Atom::var(x)).pow(10)]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let coefficients = system.taylor(&p.zero(), &[p.i(1)], 16).unwrap();
    let chart = exact_constructor(p, &system.exact_source_rows, &p.zero(), &coefficients, 1);
    assert!(
        chart
            .defect_bounds(p, &p.rational(&Rational::from((1, 3))), None)
            .is_ok()
    );
    assert!(matches!(
        chart.defect_bounds(p, &p.i(2), None),
        Err(Error::Accuracy(_))
    ));
}

#[test]
fn exact_constructor_coupled_channels_match_an_independent_augmented_system() {
    let p = Precision::decimal(60).unwrap();
    let x = symbol!("exact_residual_channels::x");
    let z = Atom::zero();
    let a = Atom::one() / (Atom::one() - Atom::var(x));
    let b = Atom::num(2) / (Atom::num(2) + Atom::var(x));
    // Both current-channel entries are present, so every ordinary row retains
    // the same denominator as its corresponding shared epsilon row.
    let coupled = crate::ode::compile_rows(
        x,
        &[
            vec![a.clone(), b.clone(), z.clone(), b.clone()],
            vec![a.clone(), b.clone(), a.clone(), z.clone()],
        ],
        p,
        &Default::default(),
    )
    .unwrap();
    let ordinary = crate::ode::compile_rows(
        x,
        &[
            vec![
                a.clone(),
                b.clone(),
                z.clone(),
                z.clone(),
                z.clone(),
                z.clone(),
            ],
            vec![
                a.clone(),
                b.clone(),
                z.clone(),
                z.clone(),
                z.clone(),
                z.clone(),
            ],
            vec![
                z.clone(),
                b.clone(),
                a.clone(),
                b.clone(),
                z.clone(),
                z.clone(),
            ],
            vec![
                a.clone(),
                z.clone(),
                a.clone(),
                b.clone(),
                z.clone(),
                z.clone(),
            ],
            vec![
                z.clone(),
                z.clone(),
                z.clone(),
                b.clone(),
                a.clone(),
                b.clone(),
            ],
            vec![z.clone(), z.clone(), a.clone(), z, a, b],
        ],
        p,
        &Default::default(),
    )
    .unwrap();
    let center = p.parse("0.3", "0.125").unwrap();
    let values = [2, 3, 5, 7, 11, 13].map(|a| p.i(a));
    let coefficients = ordinary.taylor(&center, &values, 24).unwrap();
    let combined = exact_constructor(p, &coupled.exact_source_rows, &center, &coefficients, 3);
    let independent = exact_constructor(p, &ordinary.exact_source_rows, &center, &coefficients, 1);
    for step in [p.parse("0.01", "0").unwrap(), p.parse("0", "0.1").unwrap()] {
        assert_eq!(
            combined.defect_bounds(p, &step, None).unwrap(),
            independent.defect_bounds(p, &step, None).unwrap()
        );
    }
}

fn coefficient_height(rows: &[ExactPolynomialRow]) -> u64 {
    rows.iter()
        .flat_map(|row| {
            row.denominator
                .iter()
                .chain(row.entries.iter().flat_map(|(_, a)| a))
        })
        .flat_map(|a| [&a.re, &a.im])
        .map(|a| {
            a.numerator_ref()
                .significant_bits()
                .max(a.denominator_ref().significant_bits())
        })
        .max()
        .unwrap_or(0)
}

fn benchmark_case(
    name: &str,
    system: &crate::ode::CompiledSystem,
    center: &C,
    values: &[C],
    order: usize,
    repetitions: usize,
) -> serde_json::Value {
    eprintln!("benchmarking {name}, bits={}, order={order}", system.p.bits);
    let p = system.p;
    let started = Instant::now();
    let coefficients = system.taylor(center, values, order).unwrap();
    let taylor_ns = started.elapsed().as_nanos();
    let mut old_ns = Vec::new();
    let mut exact_ns = Vec::new();
    let mut old = None;
    let mut exact = None;
    // Alternate ordering to reduce systematic cache/warm-up bias.
    for repeat in 0..repetitions {
        for method in if repeat % 2 == 0 { [0, 1] } else { [1, 0] } {
            let started = Instant::now();
            let chart = if method == 0 {
                ExactSourceResidual::new(p, &system.exact_source_rows, center, &coefficients, 1)
                    .unwrap()
            } else {
                exact_constructor(p, &system.exact_source_rows, center, &coefficients, 1)
            };
            let elapsed = started.elapsed().as_nanos();
            if method == 0 {
                old_ns.push(elapsed);
                old = Some(chart);
            } else {
                exact_ns.push(elapsed);
                exact = Some(chart);
            }
        }
    }
    let old = old.unwrap();
    let exact = exact.unwrap();
    let mut old_polynomials = Vec::new();
    for (i, denominator) in old.denominators.iter().enumerate() {
        old_polynomials.push(&old.residuals[i]);
        old_polynomials.push(&denominator.expanded);
        old_polynomials.extend(denominator.factors.iter().map(|(a, _)| a));
    }
    let mut old_polynomials = old_polynomials.into_iter();
    let mut components_checked = 0;
    let mut exact_interval_not_nested = 0;
    let mut exact_width_larger = 0;
    let mut maximum_exact_bits = 0;
    exact_constructor_observed(
        p,
        &system.exact_source_rows,
        center,
        &coefficients,
        1,
        |raw, new| {
            let old = old_polynomials.next().unwrap();
            let zero_ball = dyadic_ball(p, &p.zero());
            let zero_exact = Gaussian::new(Rational::zero(), Rational::zero());
            for index in 0..old
                .coefficients()
                .len()
                .max(new.coefficients().len())
                .max(raw.coefficients().len())
            {
                let a = old.coefficients().get(index).unwrap_or(&zero_ball);
                let b = new.coefficients().get(index).unwrap_or(&zero_ball);
                let exact = raw.coefficients().get(index).unwrap_or(&zero_exact);
                for (a, b, exact) in [(&a.re, &b.re, &exact.re), (&a.im, &b.im, &exact.im)] {
                    assert!(
                        contains(a, exact),
                        "old enclosure does not contain exact coefficient in {name}, bits={}",
                        p.bits
                    );
                    assert!(
                        contains(b, exact),
                        "new enclosure does not contain exact coefficient in {name}, bits={}",
                        p.bits
                    );
                    components_checked += 1;
                    exact_interval_not_nested += usize::from(
                        a.lower_bound().to_rational() > b.lower_bound().to_rational()
                            || b.upper_bound().to_rational() > a.upper_bound().to_rational(),
                    );
                    exact_width_larger += usize::from(b.radius > a.radius);
                    maximum_exact_bits = maximum_exact_bits
                        .max(exact.numerator_ref().significant_bits())
                        .max(exact.denominator_ref().significant_bits());
                }
            }
        },
    );
    assert!(old_polynomials.next().is_none());
    let trials = ["0.001", "0.01", "0.1", "1"].map(|step| {
        let step = p.parse(step, "0").unwrap();
        let evaluated = crate::ode::evaluate_taylor(p, &coefficients, &step).0;
        let format = |chart: &ExactSourceResidual, budget| {
            match chart.defect_bounds(p, &step, budget) {
                Ok(bounds) => serde_json::json!({"bounds": bounds.iter().map(|v| v.as_raw().to_string()).collect::<Vec<_>>()}),
                Err(e) => serde_json::json!({"error": format!("{e:?}")}),
            }
        };
        let tolerance = p.tolerance(p.bits / 4);
        serde_json::json!({
            "step": step.re.as_raw().to_string(),
            "ball": format(&old, None), "exact": format(&exact, None),
            "ball_with_budget": format(&old, Some((&evaluated, &tolerance))),
            "exact_with_budget": format(&exact, Some((&evaluated, &tolerance))),
        })
    });
    serde_json::json!({
        "name": name, "bits": p.bits, "order": order, "dimension": system.dimension(),
        "center": [center.re.as_raw().to_string(), center.im.as_raw().to_string()],
        "coefficient_height_bits": coefficient_height(&system.exact_source_rows),
        "maximum_denominator_degree": system.exact_source_rows.iter().map(|r| r.denominator.len()-1).max(),
        "maximum_entry_degree": system.exact_source_rows.iter().flat_map(|r| &r.entries).map(|(_, a)| a.len()-1).max(),
        "taylor_ns": taylor_ns, "ball_constructor_ns": old_ns, "exact_constructor_ns": exact_ns,
        "exact_components_contained_in_both": components_checked,
        "exact_interval_not_nested": exact_interval_not_nested,
        "exact_width_larger": exact_width_larger,
        "maximum_exact_coefficient_bits": maximum_exact_bits,
        "trials": trials,
    })
}

#[test]
#[ignore = "explicit performance experiment; writes JSON to EXACT_RESIDUAL_BENCH_DIR"]
fn exact_source_residual_performance_comparison() {
    let directory = std::path::PathBuf::from(std::env::var("EXACT_RESIDUAL_BENCH_DIR").unwrap());
    let mut results = Vec::new();
    for filename in ["upstream12.json", "paper27.json"] {
        let input: serde_json::Value =
            serde_json::from_slice(&std::fs::read(directory.join(filename)).unwrap()).unwrap();
        for bits in [201, 333] {
            let p = Precision { bits };
            let parse =
                |text: &str| Atom::parse(text, "exact_residual_bench", Default::default()).unwrap();
            let epsilon = parse("eps");
            let value = parse(input["epsilon"].as_str().unwrap());
            let matrix = input["matrix"]
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|a| {
                            parse(a.as_str().unwrap())
                                .replace(epsilon.clone())
                                .with(value.clone())
                                .together()
                                .cancel()
                        })
                        .collect()
                })
                .collect();
            let system = crate::ode::DifferentialSystem {
                variable: symbol!("exact_residual_bench::eta"),
                matrix,
            }
            .compile(p, &Default::default())
            .unwrap();
            let pair = |a: &serde_json::Value| {
                p.parse(a[0].as_str().unwrap(), a[1].as_str().unwrap())
                    .unwrap()
            };
            let center = pair(&input["start"]);
            let values = input["boundary"]
                .as_array()
                .unwrap()
                .iter()
                .map(pair)
                .collect::<Vec<_>>();
            results.push(benchmark_case(filename, &system, &center, &values, 80, 3));
            std::fs::write(
                directory.join("comparison.json"),
                serde_json::to_vec_pretty(&results).unwrap(),
            )
            .unwrap();
        }
    }
    for (name, text, center, order) in [
        (
            "nondyadic_complex",
            "(1/3+2*ⅈ/7+x/11)/(1-x/5)",
            ["0.3", "0.125"],
            80,
        ),
        ("sparse_high_degree", "x^257", ["0", "0"], 80),
        (
            "shifted_repeated_degree64",
            "1/(1-x)^64",
            ["0.3", "0.125"],
            80,
        ),
        (
            "rational_growth_degree32",
            "1/(1-(2^512+1)/(2^512+3)*x)^32",
            ["0.3", "0.125"],
            80,
        ),
    ] {
        let p = Precision { bits: 201 };
        let system = crate::ode::DifferentialSystem {
            variable: symbol!("exact_residual_counterexample::x"),
            matrix: vec![vec![
                Atom::parse(text, "exact_residual_counterexample", Default::default())
                    .unwrap()
                    .replace(
                        Atom::parse("ⅈ", "exact_residual_counterexample", Default::default())
                            .unwrap(),
                    )
                    .with(Atom::num(Gaussian::new(Rational::zero(), Rational::one()))),
            ]],
        }
        .compile(p, &Default::default())
        .unwrap();
        results.push(benchmark_case(
            name,
            &system,
            &p.parse(center[0], center[1]).unwrap(),
            &[p.i(1)],
            order,
            3,
        ));
        std::fs::write(
            directory.join("comparison.json"),
            serde_json::to_vec_pretty(&results).unwrap(),
        )
        .unwrap();
    }
    for bits in [216, 415, 548] {
        let p = Precision { bits };
        let large = Rational::from(Integer::from(2).pow(400));
        let x = symbol!("exact_residual_cancellation::x");
        for (label, offset) in [
            ("dyadic_cancellation", Rational::one()),
            ("nondyadic_cancellation", Rational::from((1, 3))),
        ] {
            let system = crate::ode::DifferentialSystem {
                variable: x,
                matrix: vec![
                    vec![Atom::zero(), Atom::num(offset - &large)],
                    vec![Atom::zero(), Atom::zero()],
                ],
            }
            .compile(p, &Default::default())
            .unwrap();
            results.push(benchmark_case(
                label,
                &system,
                &p.zero(),
                &[p.rational(&large), p.i(1)],
                80,
                3,
            ));
            std::fs::write(
                directory.join("comparison.json"),
                serde_json::to_vec_pretty(&results).unwrap(),
            )
            .unwrap();
        }
    }
}
