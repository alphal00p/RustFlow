//! Root-aware whole-disk differential-defect estimates.
//! Polynomial products are owned by Symbolica; root errors are controlled via
//! their existing rational scalar ODEs and propagated through products without
//! subtracting nearly equal positive bounds. These are finite-MPFR estimates,
//! not directed intervals or a global accumulated solution-error certificate.
use super::{CompiledAlgebraicSystem, NumericRational};
use crate::fixed_series::FixedComplexRing;
use crate::ode::residual::RationalResidualChart;
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::collections::BTreeMap;
use std::sync::Arc;
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

type Polynomial = UnivariatePolynomial<FixedComplexRing>;
fn polynomial(p: Precision, coefficients: Vec<C>) -> Polynomial {
    Polynomial::from_coefficients(
        &FixedComplexRing::new(p),
        coefficients,
        Arc::new(PolyVariable::Temporary(0)),
    )
}
fn rational_polynomial(p: Precision, polynomial: &Polynomial) -> NumericRational {
    NumericRational {
        numerator: if polynomial.is_zero() {
            vec![p.zero()]
        } else {
            polynomial.coefficients().to_vec()
        },
        denominator: vec![p.i(1)],
    }
}
fn bound(p: Precision, value: &NumericRational, center: &C, radius: &Float) -> Result<Float> {
    crate::diffexp::rational_disk_bound(p, value, center, radius)?
        .ok_or_else(|| Error::Accuracy("cannot bound registered-root chart on trial disk".into()))
}

pub(super) struct AlgebraicResidualChart {
    center: C,
    main: Option<RationalResidualChart>,
    root: Option<RationalResidualChart>,
    root_polynomials: Vec<NumericRational>,
    solution_polynomials: Vec<NumericRational>,
}
impl AlgebraicResidualChart {
    /// Analytic-origin initialization checks root charts independently of its
    /// regular-singular solution recurrence. No solution majorant is asserted.
    pub(super) fn roots_only(
        c: &CompiledAlgebraicSystem,
        center: &C,
        coefficients: &[Vec<C>],
    ) -> Result<Self> {
        if coefficients.is_empty() || coefficients.iter().any(|r| r.len() != c.roots.len()) {
            return Err(Error::InvalidInput("root-only residual dimensions".into()));
        }
        let p = c.p;
        let root = c
            .root_system
            .as_ref()
            .map(|system| {
                RationalResidualChart::new(p, &system.polynomial_rows, center, coefficients, 1)
            })
            .transpose()?;
        let root_polynomials = (0..c.roots.len())
            .map(|i| NumericRational {
                numerator: coefficients.iter().map(|r| r[i].clone()).collect(),
                denominator: vec![p.i(1)],
            })
            .collect();
        Ok(Self {
            center: center.clone(),
            main: None,
            root,
            root_polynomials,
            solution_polynomials: Vec::new(),
        })
    }
    pub(super) fn new(
        c: &CompiledAlgebraicSystem,
        center: &C,
        solutions: &[Vec<C>],
        root_coefficients: &[Vec<C>],
    ) -> Result<Self> {
        let p = c.p;
        let size = c
            .size
            .checked_mul(c.count)
            .ok_or_else(|| Error::Limit("root residual dimension overflow".into()))?;
        if solutions.is_empty()
            || solutions.iter().any(|r| r.len() != size)
            || root_coefficients.is_empty()
            || root_coefficients.iter().any(|r| r.len() != c.roots.len())
            || c.residual_rows.len() != c.size
        {
            return Err(Error::InvalidInput("root residual chart dimensions".into()));
        }
        let solutions = (0..size)
            .map(|i| polynomial(p, solutions.iter().map(|r| r[i].clone()).collect()))
            .collect::<Vec<_>>();
        let roots = (0..c.roots.len())
            .map(|i| polynomial(p, root_coefficients.iter().map(|r| r[i].clone()).collect()))
            .collect::<Vec<_>>();
        let root = c
            .root_system
            .as_ref()
            .map(|system| {
                RationalResidualChart::new(p, &system.polynomial_rows, center, root_coefficients, 1)
            })
            .transpose()?;
        let mut root_products = BTreeMap::<Vec<usize>, Polynomial>::new();
        root_products.insert(Vec::new(), polynomial(p, vec![p.i(1)]));
        for term in &c.terms {
            if !root_products.contains_key(&term.roots) {
                let mut product = polynomial(p, vec![p.i(1)]);
                for &index in &term.roots {
                    product = &product * &roots[index];
                }
                root_products.insert(term.roots.clone(), product);
            }
        }
        let mut residuals = vec![
            NumericRational {
                numerator: vec![p.zero()],
                denominator: vec![p.i(1)]
            };
            size
        ];
        for (i, row) in c.residual_rows.iter().enumerate() {
            let denominator = polynomial(p, row.denominator.clone()).shift_var(center);
            // First collect q*P by root monomial. This shares root polynomial
            // products across sparse terms and keeps exact structural zeros.
            let mut grouped = BTreeMap::<(usize, Vec<usize>), Polynomial>::new();
            for (index, coefficients) in &row.entries {
                let term = &c.terms[*index];
                if term.row != i {
                    return Err(Error::InvalidInput("root residual term row changed".into()));
                }
                let numerator = polynomial(p, coefficients.clone()).shift_var(center);
                for epsilon in term.shift..c.count {
                    let source = (epsilon - term.shift) * c.size + term.column;
                    let key = (epsilon, term.roots.clone());
                    let contribution = &numerator * &solutions[source];
                    let previous = grouped.entry(key).or_insert_with(|| contribution.zero());
                    *previous = &*previous + &contribution;
                }
            }
            for epsilon in 0..c.count {
                let output = epsilon * c.size + i;
                let mut numerator = &denominator * &solutions[output].derivative();
                for ((order, roots), contribution) in &grouped {
                    if *order == epsilon {
                        numerator = &numerator - &(contribution * &root_products[roots]);
                    }
                }
                residuals[output] = NumericRational {
                    numerator: if numerator.is_zero() {
                        vec![p.zero()]
                    } else {
                        numerator.coefficients().to_vec()
                    },
                    denominator: if denominator.is_zero() {
                        vec![p.zero()]
                    } else {
                        denominator.coefficients().to_vec()
                    },
                };
            }
        }
        Ok(Self {
            center: center.clone(),
            main: Some(RationalResidualChart::from_residuals(residuals)?),
            root,
            root_polynomials: roots.iter().map(|q| rational_polynomial(p, q)).collect(),
            solution_polynomials: solutions
                .iter()
                .map(|q| rational_polynomial(p, q))
                .collect(),
        })
    }

    pub(super) fn root_error_bounds(
        &self,
        c: &CompiledAlgebraicSystem,
        step: &C,
    ) -> Result<Vec<Float>> {
        let p = c.p;
        let radius = p.norm(step);
        let Some(root) = &self.root else {
            return Ok(Vec::new());
        };
        let defects = root.defect_bounds(p, step)?;
        let available = (u64::from(p.bits) * 1000 / 3322) as u32;
        let rounding = p.tolerance(available.saturating_sub(3));
        c.roots
            .iter()
            .enumerate()
            .zip(defects)
            .map(|((i, root), defect)| {
                let norm = bound(p, &root.logarithmic_derivative, &self.center, &radius)?;
                let integral = norm * &radius;
                let growth = crate::diffexp::amplification_from_integral(p, &integral, None)?;
                // The starting root is reprojected from the exact radicand at each
                // chart. Retain an arithmetic rounding allowance rather than
                // silently treating that working-precision value as exact.
                let initial = rounding.clone() * p.norm(&self.root_polynomials[i].numerator[0]);
                let error = (defect + initial) * growth;
                if error.is_finite() {
                    Ok(error)
                } else {
                    Err(Error::Accuracy("nonfinite root approximation error".into()))
                }
            })
            .collect()
    }

    pub(super) fn defect_bounds(
        &self,
        c: &CompiledAlgebraicSystem,
        step: &C,
    ) -> Result<Vec<Float>> {
        let p = c.p;
        let radius = p.norm(step);
        let mut result = self
            .main
            .as_ref()
            .ok_or_else(|| {
                Error::InvalidInput(
                    "root-only chart cannot accept a solution transport step".into(),
                )
            })?
            .defect_bounds(p, step)?;
        let errors = self.root_error_bounds(c, step)?;
        if errors.is_empty() {
            return Ok(result);
        }
        let root_bounds = self
            .root_polynomials
            .iter()
            .map(|q| bound(p, q, &p.zero(), &radius))
            .collect::<Result<Vec<_>>>()?;
        let solution_bounds = self
            .solution_polynomials
            .iter()
            .map(|q| bound(p, q, &p.zero(), &radius))
            .collect::<Result<Vec<_>>>()?;
        let mut product_errors = BTreeMap::<Vec<usize>, Float>::new();
        product_errors.insert(Vec::new(), p.real(0));
        for term in &c.terms {
            let error = product_errors.entry(term.roots.clone()).or_insert_with(|| {
                let mut approximate = p.real(1);
                let mut error = p.real(0);
                for &index in &term.roots {
                    error = error * (root_bounds[index].clone() + &errors[index])
                        + approximate.clone() * &errors[index];
                    approximate *= &root_bounds[index];
                }
                error
            });
            if *error == p.real(0) {
                continue;
            }
            let coefficient = bound(p, &term.coefficient, &self.center, &radius)?;
            let kernel_error = coefficient * &*error * &radius;
            for epsilon in term.shift..c.count {
                let source = (epsilon - term.shift) * c.size + term.column;
                result[epsilon * c.size + term.row] +=
                    kernel_error.clone() * &solution_bounds[source];
            }
        }
        if result.iter().any(|x| !x.is_finite()) {
            return Err(Error::Accuracy(
                "nonfinite algebraic whole-segment defect".into(),
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebraic::{AlgebraicRun, AlgebraicSystem, RootSeed, SquareRoot};
    use crate::diffexp::{EpsilonBoundary, EpsilonSystem};
    use crate::ode::{BoundaryData, SeriesSystem, shift_polynomial};
    use crate::{FlowOptions, RunContext};
    fn atom(s: &str) -> Atom {
        Atom::parse(s, "root_residual_tests", Default::default()).unwrap()
    }
    fn symbol(s: &str) -> Symbol {
        match atom(s).as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        }
    }
    #[test]
    fn root_scalar_residual_bounds_actual_polynomial_root_error() {
        let p = Precision::decimal(60).unwrap();
        let x = symbol("x");
        let r = symbol("r");
        let c = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![vec![vec![atom("r")]]],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: atom("1+x"),
            }],
            nonzero_conditions: vec![],
        }
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
        let run = AlgebraicRun {
            compiled: &c,
            seeds: &seeds,
        };
        let boundary = BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        };
        let state = run.initial_state(&boundary).unwrap();
        let (_, chart) = run
            .local_chart(&boundary.point, &boundary.values, 8, &state)
            .unwrap();
        let step = p.rational(&Rational::from((1, 3)));
        let errors = chart
            .residual
            .as_ref()
            .unwrap()
            .root_error_bounds(&c, &step)
            .unwrap();
        let approximate = crate::ode::evaluate_taylor(p, &chart.coefficients, &step).0[0].clone();
        let exact = C::new(p.rational(&Rational::from((4, 3))).re.sqrt(), p.real(0));
        let actual = p.norm(&p.sub(&exact, &approximate));
        assert!(actual > p.real(0));
        assert!(errors[0] >= actual);
        assert!(
            run.accepted_state(&chart, &step, &p.tolerance(20))
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn two_root_kernel_full_residual_bounds_physical_defect() {
        let p = Precision::decimal(60).unwrap();
        let x = symbol("x");
        let r = symbol("r");
        let s = symbol("s");
        let c = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![vec![vec![atom("0")]], vec![vec![atom("r*s")]]],
            },
            roots: vec![
                SquareRoot {
                    symbol: r,
                    radicand: atom("1+x"),
                },
                SquareRoot {
                    symbol: s,
                    radicand: atom("1-x"),
                },
            ],
            nonzero_conditions: vec![],
        }
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::from([(r, RootSeed::Principal), (s, RootSeed::Principal)]);
        let run = AlgebraicRun {
            compiled: &c,
            seeds: &seeds,
        };
        let boundary = BoundaryData {
            point: p.zero(),
            values: vec![p.i(1), p.zero()],
        };
        let state = run.initial_state(&boundary).unwrap();
        let (coefficients, chart) = run
            .local_chart(&boundary.point, &boundary.values, 8, &state)
            .unwrap();
        let h = p.rational(&Rational::from((1, 4)));
        let bounds = run.whole_segment_residual(&chart, &h).unwrap().unwrap();
        let mut derivative = p.zero();
        for k in (1..coefficients.len()).rev() {
            derivative = p.add(
                &p.mul(&derivative, &h),
                &p.scale(&coefficients[k][1], k as i64, 1),
            );
        }
        let exact = C::new(p.rational(&Rational::from((15, 16))).re.sqrt(), p.real(0));
        let actual = p.norm(&p.mul(&h, &p.sub(&derivative, &exact)));
        assert!(actual > p.real(0));
        assert!(bounds[1] >= actual);
    }
    #[test]
    fn rowwise_exact_common_denominators_preserve_term_mapping() {
        let p = Precision::decimal(60).unwrap();
        let x = symbol("x");
        let r = symbol("r");
        let c = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![
                    vec![vec![atom("r/(1+x)+1/(1-x)")]],
                    vec![vec![atom("r*x^1000/(1+x)^2")]],
                ],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: atom("1+x"),
            }],
            nonzero_conditions: vec![],
        }
        .compile(p)
        .unwrap();
        assert_eq!(c.residual_rows.len(), 1);
        assert_eq!(c.residual_rows[0].entries.len(), 3);
        assert!(
            c.residual_rows[0]
                .entries
                .iter()
                .any(|(_, a)| a.len() > 1000)
        );
        let center = p.parse("0.2", "0.1").unwrap();
        let row = &c.residual_rows[0];
        let denominator = shift_polynomial(p, &row.denominator, &center, 0)[0].clone();
        for (index, numerator) in &row.entries {
            let actual = p.div(&shift_polynomial(p, numerator, &center, 0)[0], &denominator);
            let expected = c.terms[*index].coefficient.series(p, &center, 0).unwrap()[0].clone();
            assert!(p.close(&actual, &expected, 50));
        }
    }
    #[test]
    fn constant_registered_root_does_not_hide_high_degree_main_defect() {
        let p = Precision::decimal(60).unwrap();
        let x = symbol("x");
        let r = symbol("r");
        let c = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![
                    vec![vec![atom("0")]],
                    vec![vec![atom("r*2*1001*1002*1003*x^1000*(x-1/2)*(x-1)")]],
                ],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: atom("1"),
            }],
            nonzero_conditions: vec![],
        }
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
        let options = FlowOptions {
            series_order: 80,
            max_steps: 1000,
            ..Default::default()
        };
        let answer = c
            .transport(
                &EpsilonBoundary {
                    point: p.zero(),
                    leading: 0,
                    coefficients: vec![vec![p.i(1)], vec![p.zero()]],
                },
                &[p.i(1)],
                &seeds,
                &options,
                &RunContext::default(),
                true,
            )
            .unwrap();
        assert!(answer.solution.diagnostics.steps > 1);
        assert!(p.close(&answer.solution.coefficients[1][0], &p.i(-999), 20));
    }
    #[test]
    fn no_root_algebraic_owner_still_requires_rational_majorant() {
        let p = Precision::decimal(60).unwrap();
        let x = symbol("x");
        let c = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![vec![vec![atom("0")]], vec![vec![atom("1/(1-x)")]]],
            },
            roots: vec![],
            nonzero_conditions: vec![],
        }
        .compile(p)
        .unwrap();
        let seeds = BTreeMap::new();
        let run = AlgebraicRun {
            compiled: &c,
            seeds: &seeds,
        };
        let boundary = BoundaryData {
            point: p.zero(),
            values: vec![p.i(1), p.zero()],
        };
        let state = run.initial_state(&boundary).unwrap();
        let (_, chart) = run
            .local_chart(&boundary.point, &boundary.values, 8, &state)
            .unwrap();
        assert_eq!(chart.coefficients.len(), 9);
        assert!(chart.coefficients.iter().all(Vec::is_empty));
        let bounds = run
            .whole_segment_residual(&chart, &p.rational(&Rational::from((1, 4))))
            .unwrap()
            .unwrap();
        assert_eq!(bounds.len(), 2);
        assert_eq!(bounds[0], p.real(0));
        assert!(bounds[1] > p.real(0));
    }
}
