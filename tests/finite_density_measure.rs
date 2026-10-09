//! Independent distribution tests of generated sources, before guarded replay.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::{guarded::IndexRole, measure::WeightedMeasure};

fn choose(n: i64, k: i64) -> Rational {
    if k < 0 || n < k {
        return Rational::from(0);
    }
    (0..k).fold(Rational::from(1), |value, j| {
        value * Rational::from(n - j) / Rational::from(j + 1)
    })
}

// Integrate x^a y^b C_n(y-2) H_upper(1-x) H_lower(x) directly.
// The distribution has compact support. Distinct endpoint supports avoid
// undefined products of distributions, and the opposite theta is locally one.
fn distribution_integral(indices: [i64; 5]) -> Rational {
    let [nx, ny, cut, upper, lower] = indices;
    if cut <= 0 {
        return Rational::from(0);
    }
    assert!(nx <= 0 && ny <= 0, "test functions must remain polynomial");
    assert!(upper >= 0 && lower >= 0, "occupation index must be defined");
    let a = -nx;
    let b = -ny;
    let order = cut - 1;
    if order > b {
        return Rational::from(0);
    }
    // C_n = (-1)^(n-1) delta^(n-1)/(n-1)!: testing it against
    // a polynomial gives its derivative divided by (n-1)!.
    let shell = choose(b, order) * Rational::from(1_i64 << (b - order));
    let spatial = match (upper, lower) {
        (0, 0) => Rational::from((1, a + 1)),
        (u, 0) => {
            let order = u - 1;
            choose(a, order) * Rational::from(if order % 2 == 0 { 1 } else { -1 })
        }
        (0, l) => Rational::from(i64::from(a == l - 1)),
        _ => Rational::from(0),
    };
    shell * spatial
}

#[test]
fn generated_sources_annihilate_exact_polynomial_distribution_tests() {
    let x = parse!("fd_distribution_x");
    let y = parse!("fd_distribution_y");
    let measure = WeightedMeasure::new(
        vec![x.clone(), y.clone()],
        [
            x.clone(),
            y.clone(),
            &y - Atom::num(2),
            Atom::num(1) - &x,
            x.clone(),
        ],
        [
            IndexRole::Ordinary,
            IndexRole::Ordinary,
            IndexRole::RequiredCut,
            IndexRole::Occupation,
            IndexRole::Occupation,
        ],
    )
    .unwrap();
    let indices = [
        symbol!("fd_distribution_nx"),
        symbol!("fd_distribution_ny"),
        symbol!("fd_distribution_nc"),
        symbol!("fd_distribution_nu"),
        symbol!("fd_distribution_nl"),
    ];
    let mut sources = measure.multiplication_sources().unwrap();
    sources.extend(
        measure
            .ibp(
                "constant-field",
                &indices,
                &[Atom::num(1), Atom::num(1)],
                &Atom::new(),
                4,
            )
            .unwrap(),
    );
    sources.extend(
        measure
            .ibp(
                "polynomial-field",
                &indices,
                &[&x * (Atom::num(1) - &x), y],
                &(Atom::num(2) - Atom::num(2) * &x),
                4,
            )
            .unwrap(),
    );
    let mut tested = 0;
    let mut nonzero_surface_terms = 0;
    for a in 0..=3 {
        for b in 0..=3 {
            for cut in -2..=4 {
                for upper in 0..=3 {
                    for lower in 0..=3 {
                        let point = [-a, -b, cut, upper, lower];
                        let substitutions = indices
                            .iter()
                            .zip(point)
                            .map(|(&s, value)| (Atom::var(s), Atom::num(value)))
                            .collect::<BTreeMap<_, _>>();
                        for source in &sources {
                            if !source.domain.contains(&point) {
                                continue;
                            }
                            let mut value = Rational::from(0);
                            for term in &source.terms {
                                let coefficient = symbolica_amflow::family::substitute(
                                    &term.coefficient,
                                    &substitutions,
                                )
                                .together()
                                .cancel();
                                let coefficient =
                                    Rational::try_from(coefficient.as_view()).unwrap();
                                if coefficient.is_zero() {
                                    continue;
                                }
                                let shifted =
                                    std::array::from_fn(|i| point[i] + i64::from(term.shift[i]));
                                let contribution = coefficient * distribution_integral(shifted);
                                if !contribution.is_zero() && (shifted[3] > 0 || shifted[4] > 0) {
                                    nonzero_surface_terms += 1;
                                }
                                value += contribution;
                            }
                            assert!(
                                value.is_zero(),
                                "invalid source {} at {point:?}: {value}",
                                source.id
                            );
                            tested += 1;
                        }
                    }
                }
            }
        }
    }
    assert!(tested > 1000);
    assert!(
        nonzero_surface_terms > 0,
        "endpoint identities must be nonvacuous"
    );
}
