//! Explicit benchmark harness: same shared transport controller, only the
//! source-defect owner changes. This is never part of a production selector.
use super::*;
use crate::ode::{
    BoundaryData, CompiledSystem, DifferentialSystem, SeriesSystem, TaylorCoordinate,
};
use crate::{FlowOptions, RunContext};
use std::cell::Cell;
use std::time::Instant;

struct Controlled<'a> {
    system: &'a CompiledSystem,
    adaptive: bool,
    integer_charts: Cell<usize>,
    ball_charts: Cell<usize>,
}
impl Controlled<'_> {
    fn make(
        &self,
        rows: &[ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
    ) -> Result<AdaptiveResidual> {
        let p = self.system.p;
        let chart = if self.adaptive {
            AdaptiveResidual::new(p, rows, center, coefficients, 1, &Default::default())?
        } else {
            AdaptiveResidual::Ball {
                source: ExactSourceResidual::new(p, rows, center, coefficients, 1)?,
                cancellation: Default::default(),
            }
        };
        match chart {
            AdaptiveResidual::Integer { .. } => {
                self.integer_charts.set(self.integer_charts.get() + 1)
            }
            AdaptiveResidual::Ball { .. } => self.ball_charts.set(self.ball_charts.get() + 1),
        }
        Ok(chart)
    }
}
impl SeriesSystem for Controlled<'_> {
    type State = ();
    type Chart = AdaptiveResidual;
    fn initial_state(&self, _: &BoundaryData) -> Result<()> {
        Ok(())
    }
    fn accepted_state(&self, _: &Self::Chart, _: &C, _: &Float) -> Result<Option<()>> {
        Ok(Some(()))
    }
    fn precision(&self) -> Precision {
        self.system.p
    }
    fn dimension(&self) -> usize {
        self.system.dimension()
    }
    fn poles(&self) -> &[C] {
        &self.system.poles
    }
    fn local_chart(
        &self,
        center: &C,
        values: &[C],
        order: usize,
        _: &(),
    ) -> Result<(Vec<Vec<C>>, Self::Chart)> {
        let coefficients = self.system.taylor(center, values, order)?;
        let chart = self.make(&self.system.exact_source_rows, center, &coefficients)?;
        Ok((coefficients, chart))
    }
    fn mapped_chart(
        &self,
        coordinate: &TaylorCoordinate,
        _: &C,
        values: &[C],
        order: usize,
        _: &(),
    ) -> Result<crate::ode::MappedTaylorChart<Self::Chart>> {
        let mapped = self.system.in_coordinate(coordinate)?;
        let center = self.system.p.zero();
        let coefficients = mapped.taylor(&center, values, order)?;
        let chart = self.make(&mapped.exact_source_rows, &center, &coefficients)?;
        Ok((coefficients, chart, mapped.poles))
    }
    fn whole_segment_residual(&self, chart: &Self::Chart, step: &C) -> Result<Option<Vec<Float>>> {
        chart.defect_bounds(self.system.p, step, None).map(Some)
    }
    fn whole_segment_residual_with_budget(
        &self,
        chart: &Self::Chart,
        step: &C,
        values: &[C],
        tolerance: &Float,
    ) -> Result<Option<Vec<Float>>> {
        chart
            .defect_bounds(self.system.p, step, Some((values, tolerance)))
            .map(Some)
    }
    fn rhs(&self, point: &C, values: &[C], _: &Self::Chart) -> Result<Vec<C>> {
        let p = self.system.p;
        self.system
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .zip(values)
                    .try_fold(p.zero(), |sum, (entry, value)| {
                        Ok(p.add(&sum, &p.mul(&entry.series(p, point, 0)?[0], value)))
                    })
            })
            .collect()
    }
}

#[test]
#[ignore = "explicit constructor and complete transport benchmark; ADAPTIVE_RESIDUAL_BENCH_DIR contains initial/endpoint fixture directories"]
fn adaptive_residual_performance_and_full_transport() {
    let directory = std::path::PathBuf::from(std::env::var("ADAPTIVE_RESIDUAL_BENCH_DIR").unwrap());
    let mut result = Vec::new();
    for chart_kind in ["initial", "endpoint"] {
        for name in ["upstream12", "paper27"] {
            let input: serde_json::Value = serde_json::from_slice(
                &std::fs::read(directory.join(chart_kind).join(format!("{name}.json"))).unwrap(),
            )
            .unwrap();
            for bits in [201, 333] {
                let p = Precision { bits };
                let parse = |a: &str| {
                    Atom::parse(a, "adaptive_residual_bench", Default::default()).unwrap()
                };
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
                let system = DifferentialSystem {
                    variable: symbol!("adaptive_residual_bench::eta"),
                    matrix,
                }
                .compile(p, &Default::default())
                .unwrap();
                let pair = |a: &serde_json::Value| {
                    p.parse(a[0].as_str().unwrap(), a[1].as_str().unwrap())
                        .unwrap()
                };
                let center = pair(&input["start"]);
                let end = pair(&input["end"]);
                let values = input["boundary"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(pair)
                    .collect::<Vec<_>>();
                let coefficients = system.taylor(&center, &values, 80).unwrap();
                let mut constructor = Vec::new();
                let mut transport = Vec::new();
                let mut comparisons = Vec::new();
                for iteration in 0..3 {
                    for adaptive in if iteration % 2 == 0 {
                        [false, true]
                    } else {
                        [true, false]
                    } {
                        let owner = Controlled {
                            system: &system,
                            adaptive,
                            integer_charts: Cell::new(0),
                            ball_charts: Cell::new(0),
                        };
                        let start = Instant::now();
                        let chart = owner
                            .make(&system.exact_source_rows, &center, &coefficients)
                            .unwrap();
                        let elapsed = start.elapsed().as_nanos();
                        let mode = if matches!(chart, AdaptiveResidual::Integer { .. }) {
                            "integer"
                        } else {
                            "ball"
                        };
                        constructor.push(
                            serde_json::json!({"adaptive":adaptive,"ns":elapsed,"mode":mode}),
                        );
                        if chart_kind == "initial" {
                            owner.integer_charts.set(0);
                            owner.ball_charts.set(0);
                            let mut fallback_owners = 0;
                            let options = FlowOptions {
                                digits: 20,
                                series_order: 80,
                                ..Default::default()
                            };
                            let boundary = BoundaryData {
                                point: center.clone(),
                                values: values.clone(),
                            };
                            let start = Instant::now();
                            let (flow, _) = crate::ode::transport_series_observed(
                                &owner,
                                &boundary,
                                std::slice::from_ref(&end),
                                &options,
                                &RunContext::default(),
                                None,
                                |chart, _| {
                                    if let AdaptiveResidual::Integer { fallback, .. } = chart {
                                        fallback_owners += usize::from(
                                            fallback.residual.lock().unwrap().is_some(),
                                        );
                                    }
                                },
                            )
                            .unwrap();
                            let elapsed = start.elapsed().as_nanos();
                            transport.push(serde_json::json!({"adaptive":adaptive,"ns":elapsed,"steps":flow.diagnostics.steps,"rejected":flow.diagnostics.rejected_steps,"integer_charts":owner.integer_charts.get(),"ball_charts":owner.ball_charts.get(),"charts_with_lazy_ball":fallback_owners}));
                            comparisons.push((adaptive, flow.values));
                        }
                    }
                }
                if let Some((_, first)) = comparisons.first() {
                    for (_, values) in &comparisons {
                        for (a, b) in first.iter().zip(values) {
                            assert!(p.close(a, b, 20), "transport mismatch {name} {bits}");
                        }
                    }
                }
                result.push(serde_json::json!({"name":name,"chart":chart_kind,"bits":bits,"constructor":constructor,"transport":transport}));
                std::fs::write(
                    directory.join("adaptive-comparison.json"),
                    serde_json::to_vec_pretty(&result).unwrap(),
                )
                .unwrap();
            }
        }
    }
}

#[test]
#[ignore = "explicit fresh precision/order and adverse-selector experiment; uses ADAPTIVE_RESIDUAL_BENCH_DIR"]
fn adaptive_residual_refinement_and_adverse_cases() {
    let directory = std::path::PathBuf::from(std::env::var("ADAPTIVE_RESIDUAL_BENCH_DIR").unwrap());
    let mut report = Vec::new();
    for name in ["upstream12", "paper27"] {
        let input: serde_json::Value = serde_json::from_slice(
            &std::fs::read(directory.join("initial").join(format!("{name}.json"))).unwrap(),
        )
        .unwrap();
        let mut previous: Option<Vec<C>> = None;
        for (bits, order) in [(201, 80), (397, 112)] {
            let p = Precision { bits };
            let parse =
                |a: &str| Atom::parse(a, "adaptive_refinement", Default::default()).unwrap();
            let eps = parse("eps");
            let epsilon = parse(input["epsilon"].as_str().unwrap());
            let matrix = input["matrix"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| {
                    r.as_array()
                        .unwrap()
                        .iter()
                        .map(|a| {
                            parse(a.as_str().unwrap())
                                .replace(eps.clone())
                                .with(epsilon.clone())
                                .together()
                                .cancel()
                        })
                        .collect()
                })
                .collect();
            // Compile and rebuild all source/parameter/Taylor values at each
            // precision; never raise the precision of a previous numerical DE.
            let system = DifferentialSystem {
                variable: symbol!("adaptive_refinement::eta"),
                matrix,
            }
            .compile(p, &Default::default())
            .unwrap();
            let pair = |a: &serde_json::Value| {
                p.parse(a[0].as_str().unwrap(), a[1].as_str().unwrap())
                    .unwrap()
            };
            let boundary = BoundaryData {
                point: pair(&input["start"]),
                values: input["boundary"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(pair)
                    .collect(),
            };
            let end = pair(&input["end"]);
            let options = FlowOptions {
                digits: 20,
                series_order: order,
                ..Default::default()
            };
            for adaptive in [false, true] {
                let owner = Controlled {
                    system: &system,
                    adaptive,
                    integer_charts: Cell::new(0),
                    ball_charts: Cell::new(0),
                };
                let start = Instant::now();
                let flow = crate::ode::transport_series(
                    &owner,
                    &boundary,
                    std::slice::from_ref(&end),
                    &options,
                    &RunContext::default(),
                    None,
                )
                .unwrap();
                let elapsed = start.elapsed().as_nanos();
                if let Some(reference) = &previous {
                    for (a, b) in reference.iter().zip(&flow.values) {
                        assert!(
                            p.close(a, b, 20),
                            "independent precision/order disagreement {name} {bits} {adaptive}"
                        );
                    }
                } else {
                    previous = Some(flow.values.clone());
                }
                report.push(serde_json::json!({"kind":"refinement","name":name,"bits":bits,"order":order,"adaptive":adaptive,"ns":elapsed,"steps":flow.diagnostics.steps,"rejected":flow.diagnostics.rejected_steps,"integer_charts":owner.integer_charts.get(),"ball_charts":owner.ball_charts.get(),"values":flow.values.iter().map(|a|[a.re.as_raw().to_string(),a.im.as_raw().to_string()]).collect::<Vec<_>>() }));
            }
            std::fs::write(
                directory.join("refinement-adverse.json"),
                serde_json::to_vec_pretty(&report).unwrap(),
            )
            .unwrap();
        }
    }
    for (name, text, zero_center) in [
        ("nondyadic_complex", "(1/3+2*i/7+x/11)/(1-x/5)", false),
        ("sparse_degree257", "x^257", true),
        ("repeated_degree64", "1/(1-x)^64", false),
        (
            "rational_growth_degree32",
            "1/(1-(2^512+1)/(2^512+3)*x)^32",
            false,
        ),
    ] {
        for bits in [201, 397] {
            let p = Precision { bits };
            let parse = |a: &str| Atom::parse(a, "adaptive_adverse", Default::default()).unwrap();
            let matrix =
                vec![vec![parse(text).replace(parse("i")).with(Atom::num(
                    Gaussian::new(Rational::zero(), Rational::one()),
                ))]];
            let system = DifferentialSystem {
                variable: symbol!("adaptive_adverse::x"),
                matrix,
            }
            .compile(p, &Default::default())
            .unwrap();
            let center = if zero_center {
                p.zero()
            } else {
                p.parse("0.3", "0.125").unwrap()
            };
            let end = if zero_center {
                p.rational(&Rational::from((1, 2)))
            } else {
                p.add(&center, &p.parse("0.0000000001", "0").unwrap())
            };
            let boundary = BoundaryData {
                point: center.clone(),
                values: vec![p.i(1)],
            };
            let coefficients = system.taylor(&center, &boundary.values, 80).unwrap();
            let mut reference: Option<Vec<C>> = None;
            let mut reference_error: Option<String> = None;
            for iteration in 0..3 {
                for adaptive in if iteration % 2 == 0 {
                    [false, true]
                } else {
                    [true, false]
                } {
                    let owner = Controlled {
                        system: &system,
                        adaptive,
                        integer_charts: Cell::new(0),
                        ball_charts: Cell::new(0),
                    };
                    let start = Instant::now();
                    let chart = owner
                        .make(&system.exact_source_rows, &center, &coefficients)
                        .unwrap();
                    let constructor_ns = start.elapsed().as_nanos();
                    let mode = if matches!(chart, AdaptiveResidual::Integer { .. }) {
                        "integer"
                    } else {
                        "ball"
                    };
                    owner.integer_charts.set(0);
                    owner.ball_charts.set(0);
                    let start = Instant::now();
                    let flow = crate::ode::transport_series(
                        &owner,
                        &boundary,
                        std::slice::from_ref(&end),
                        &FlowOptions::default(),
                        &RunContext::default(),
                        None,
                    );
                    let transport_ns = start.elapsed().as_nanos();
                    let (steps, outcome) = match flow {
                        Ok(flow) => {
                            assert!(
                                reference_error.is_none(),
                                "changed adverse outcome {name} {bits}"
                            );
                            if let Some(reference) = &reference {
                                for (a, b) in reference.iter().zip(&flow.values) {
                                    assert!(p.close(a, b, 20));
                                }
                            } else {
                                reference = Some(flow.values.clone());
                            }
                            (Some(flow.diagnostics.steps), "success".to_owned())
                        }
                        Err(error) => {
                            assert!(reference.is_none(), "changed adverse outcome {name} {bits}");
                            assert!(
                                matches!(error, Error::InsufficientPrecision { .. }),
                                "unexpected {error}"
                            );
                            let error = format!("{error:?}");
                            if let Some(expected) = &reference_error {
                                assert_eq!(&error, expected);
                            } else {
                                reference_error = Some(error.clone());
                            }
                            (None, error)
                        }
                    };
                    report.push(serde_json::json!({"kind":"adverse","name":name,"bits":p.bits,"adaptive":adaptive,"constructor_ns":constructor_ns,"mode":mode,"transport_ns":transport_ns,"steps":steps,"outcome":outcome,"integer_charts":owner.integer_charts.get(),"ball_charts":owner.ball_charts.get()}));
                }
            }
            std::fs::write(
                directory.join("refinement-adverse.json"),
                serde_json::to_vec_pretty(&report).unwrap(),
            )
            .unwrap();
        }
    }
}
