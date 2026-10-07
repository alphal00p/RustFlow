//! Shared residual polynomial construction in real or complex directed balls.
use super::*;

type Polynomial<T> = UnivariatePolynomial<FloatField<T>>;

#[allow(clippy::type_complexity)] // One recurrence, with either coefficient field.
pub(super) fn construct<T>(
    ring: FloatField<T>,
    rows: &[ExactPolynomialRow],
    center: T,
    coefficients: &[Vec<C>],
    channels: usize,
    dyadic: impl Fn(&C) -> T,
    exact: impl Fn(&Gaussian) -> T,
) -> (
    Vec<Polynomial<T>>,
    Vec<(Polynomial<T>, Vec<(Polynomial<T>, usize)>)>,
)
where
    T: SingleFloat + std::hash::Hash + Eq + symbolica::domains::InternalOrdering,
{
    let n = rows.len();
    let size = n * channels;
    let variable = Arc::new(PolyVariable::Temporary(0));
    let polynomial = |a| Polynomial::from_coefficients(&ring, a, variable.clone());
    let exact_polynomial = |a: &[Gaussian]| polynomial(a.iter().map(&exact).collect());
    let solutions = (0..size)
        .map(|i| polynomial(coefficients.iter().map(|r| dyadic(&r[i])).collect()))
        .collect::<Vec<_>>();
    let mut residuals = vec![polynomial(vec![ring.zero()]); size];
    let mut denominators = Vec::with_capacity(n);
    for (i, row) in rows.iter().enumerate() {
        let denominator = exact_polynomial(&row.denominator).shift_var(&center);
        let entries = row
            .entries
            .iter()
            .map(|(j, a)| (*j, exact_polynomial(a).shift_var(&center)))
            .collect::<Vec<_>>();
        for channel in 0..channels {
            let index = channel * n + i;
            let mut numerator = &denominator * &solutions[index].derivative();
            for (column, a) in &entries {
                let shift = column / n;
                if shift <= channel {
                    let source = (channel - shift) * n + column % n;
                    numerator = &numerator - &(a * &solutions[source]);
                }
            }
            residuals[index] = numerator;
        }
        let factors = row
            .denominator_factors
            .iter()
            .map(|(a, multiplicity)| (exact_polynomial(a).shift_var(&center), *multiplicity))
            .collect();
        denominators.push((denominator, factors));
    }
    (residuals, denominators)
}

pub(super) fn all_real(rows: &[ExactPolynomialRow], center: &C, coefficients: &[Vec<C>]) -> bool {
    center.im.is_zero()
        && coefficients.iter().flatten().all(|a| a.im.is_zero())
        && rows.iter().all(|r| {
            r.denominator
                .iter()
                .chain(r.entries.iter().flat_map(|(_, a)| a))
                .chain(r.denominator_factors.iter().flat_map(|(a, _)| a))
                .all(|a| a.im.is_zero())
        })
}

pub(super) fn promote(p: Precision, polynomial: Polynomial<RealBall>) -> BallPolynomial {
    let ring = FloatField::from_rep(dyadic_ball(p, &p.zero()));
    BallPolynomial::from_coefficients(
        &ring,
        polynomial
            .coefficients()
            .iter()
            .map(|a| ComplexBall::new(a.clone(), RealBall::exact(p.real(0))))
            .collect(),
        Arc::new(PolyVariable::Temporary(0)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_field_encloses_the_exact_rational_residual_in_each_channel() {
        for bits in [53, 201, 400] {
            let p = Precision { bits };
            let gaussian = |n, d| Gaussian::new(Rational::from((n, d)), Rational::zero());
            let rows = vec![ExactPolynomialRow {
                denominator: vec![gaussian(1, 7), gaussian(-2, 3)],
                denominator_factors: vec![(vec![gaussian(1, 7), gaussian(-2, 3)], 1)],
                entries: vec![
                    (0, vec![gaussian(-4, 5), gaussian(1, 11)]),
                    (1, vec![gaussian(2, 13)]),
                ],
            }];
            let center = p.rational(&Rational::from((1, 3)));
            let coefficients = (0..9)
                .map(|k| {
                    vec![
                        p.rational(&Rational::from((k - 3, 7))),
                        p.rational(&Rational::from((2 * k - 5, 11))),
                    ]
                })
                .collect::<Vec<_>>();
            assert!(all_real(&rows, &center, &coefficients));
            let actual = ExactSourceResidual::new(p, &rows, &center, &coefficients, 2).unwrap();
            let exact_dyadic = |a: &C| Gaussian::new(a.re.to_rational(), a.im.to_rational());
            let (exact, _) = construct(
                FloatField::from_rep(gaussian(0, 1)),
                &rows,
                exact_dyadic(&center),
                &coefficients,
                2,
                exact_dyadic,
                Clone::clone,
            );
            for (actual, exact) in actual.residuals.iter().zip(exact) {
                assert_eq!(actual.coefficients().len(), exact.coefficients().len());
                for (ball, value) in actual.coefficients().iter().zip(exact.coefficients()) {
                    assert!(ball.re.lower_bound().to_rational() <= value.re);
                    assert!(ball.re.upper_bound().to_rational() >= value.re);
                    assert!(ball.im.is_zero() && value.im.is_zero());
                }
            }
            let mut complex_center = center.clone();
            complex_center.im = p.tolerance(1000);
            assert!(!all_real(&rows, &complex_center, &coefficients));
            let mut complex_coefficients = coefficients;
            complex_coefficients[1][0].im = p.tolerance(1000);
            assert!(!all_real(&rows, &center, &complex_coefficients));
            let mut complex_rows = rows;
            complex_rows[0].entries[1].1[0].im = Rational::from((1, 97));
            assert!(!all_real(&complex_rows, &center, &complex_coefficients));
        }
    }
}
