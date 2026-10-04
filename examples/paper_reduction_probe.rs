//! Inspect one native IBP search without recursively expanding its raw RHS.
//!
//! cargo run --release --example paper_reduction_probe -- [numeric|domains]
//!   [three comma-separated one-based lines] [source-order: input|terms]
//!
//! Fixed paper kinematics, epsilon=1/2700 and symbolic eta on propagator 5.
//! Output is JSON. Residuals are search residuals, never certified masters.
//! The optional domains lane fixes numerator powers and leaves denominator
//! powers symbolic; its explicit limits may return incomplete search.
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::time::Instant;

use rustred::algebra::{Coefficient, CoefficientContext};
use rustred::sector::{
    Mask,
    zero::{Analyzer, Decision},
};
use rustred::solver::{
    CoordinateCase, NumericalExactBackend, SearchOptions, SectorConfig, SectorSolveOptions,
    SectorSolver, SourceDiscoveryStrategy, SourceRowFeature, SourceRowPriority, SourceSystem,
};
use serde_json::json;
use symbolica::prelude::*;
use symbolica_amflow::{MassMode, benchmarks};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

// This diagnostic has just one symbolic coefficient, eta. Constructing the
// native input here keeps experimental search controls out of the library API.
fn native_family() -> Result<rustred::family::IntegralFamily> {
    let (family, _) = benchmarks::paper_two_loop()?;
    let eta = symbol!("paper_reduction_probe::eta");
    let (family, _) = family.deform(eta, &MassMode::Propagators(vec![4]))?;
    let context = CoefficientContext::try_new(["eta"])?;
    let replacement = BTreeMap::from([(
        Atom::var(eta),
        context.parameter("eta").unwrap().to_expression(),
    )]);
    let coefficient = |a: &Atom| -> Result<Coefficient> {
        Ok(
            symbolica_amflow::family::substitute(a, &replacement).try_to_rational_polynomial(
                &Q,
                &Z,
                Some(context.one().get_variables().clone()),
            )?,
        )
    };
    let propagators = family
        .propagators
        .iter()
        .map(|p| {
            Ok(rustred::family::AffineDenominator::new(
                coefficient(&p.constant)?,
                p.scalar_products
                    .iter()
                    .map(&coefficient)
                    .collect::<Result<_>>()?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let gram = family
        .external_gram
        .iter()
        .map(|r| r.iter().map(&coefficient).collect::<Result<Vec<_>>>())
        .collect::<Result<Vec<_>>>()?;
    let dimension = coefficient(&Atom::num(Rational::from(4) - Rational::from((1, 1350))))?;
    let shifts = vec![context.zero(); propagators.len()];
    Ok(rustred::family::IntegralFamily::new(
        family.name,
        family.loops,
        family.external,
        context,
        dimension,
        propagators,
        gram,
        shifts,
    )?)
}

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let mode = args.first().map_or("numeric", String::as_str);
    if !["numeric", "domains"].contains(&mode) {
        return Err("mode must be numeric or domains".into());
    }
    let lines = args
        .get(1)
        .map_or("2,5,7", String::as_str)
        .split(',')
        .map(str::parse::<usize>)
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if lines.len() != 3
        || lines.iter().any(|&i| !(1..=7).contains(&i))
        || lines.iter().collect::<BTreeSet<_>>().len() != 3
    {
        return Err("three distinct physical line indices in 1..=7 are required".into());
    }
    let source_order = args.get(2).map_or("input", String::as_str);
    let strategy = match source_order {
        "input" => SourceDiscoveryStrategy::InputOrder,
        "terms" => SourceDiscoveryStrategy::Features(vec![
            SourceRowPriority {
                feature: SourceRowFeature::Terms,
                descending: false,
            },
            SourceRowPriority {
                feature: SourceRowFeature::CoefficientMonomials,
                descending: false,
            },
        ]),
        _ => return Err("source order must be input or terms".into()),
    };
    let started = Instant::now();
    let family = native_family()?;
    let sources = SourceSystem::<9>::from_family_with_lorentz(&family, false)?;
    let analyzer = Analyzer::try_unrestricted(&family)?;
    let sector: [bool; 9] = std::array::from_fn(|i| lines.contains(&(i + 1)));
    let mut zero = Vec::new();
    for bits in 0..512 {
        let child: [bool; 9] = std::array::from_fn(|i| bits & (1 << i) != 0);
        if child.iter().zip(sector).all(|(&a, b)| !a || b)
            && matches!(
                analyzer.analyze(&Mask::try_new(child)?)?,
                Decision::ProvedZero(_)
            )
        {
            zero.push(child);
        }
    }
    let solver = SectorSolver::new(
        &sources,
        sector,
        SectorConfig {
            zero_sectors: Arc::from(zero),
            numerical_exact_backend: NumericalExactBackend::SparseFactorized,
            source_discovery: strategy,
            ..Default::default()
        },
    )?;
    let preparation_ms = started.elapsed().as_millis();
    let inactive = (0..9).filter(|&i| !sector[i]).collect::<Vec<_>>();
    let corner: [i16; 9] = sector.map(i16::from);
    let mut targets = Vec::new();
    for (rank, dots) in [(0, 0), (3, 0), (8, 0), (3, 4)] {
        let mut target = corner;
        target[inactive[0]] = -rank;
        target[lines[0] - 1] += dots;
        targets.push(target);
    }
    let search = SearchOptions {
        max_depth: Some(3),
        ..Default::default()
    };
    if mode == "numeric" {
        let result = solver.solve_numeric_cases(
            targets
                .iter()
                .map(|v| CoordinateCase::new(v.map(Some)))
                .collect::<std::result::Result<_, _>>()?,
            search,
        )?;
        let mut rhs = BTreeSet::new();
        let mut rule_widths = Vec::new();
        let mut max_rank = 0;
        let mut max_dots = 0;
        for rule in &result.rules {
            rule_widths.push(rule.rhs.len());
            for term in &rule.rhs {
                let p: [i16; 9] = std::array::from_fn(|i| term.integral.powers()[i].value());
                max_rank = max_rank.max(p.iter().map(|&n| i32::from(-n.min(0))).sum::<i32>());
                max_dots = max_dots.max(p.iter().map(|&n| i32::from((n - 1).max(0))).sum::<i32>());
                rhs.insert(p);
            }
        }
        let by_lines = rhs
            .iter()
            .fold(BTreeMap::<usize, usize>::new(), |mut counts, p| {
                *counts
                    .entry(p.iter().filter(|&&n| n > 0).count())
                    .or_default() += 1;
                counts
            });
        println!(
            "{}",
            json!({
                "mode": mode, "lines": lines, "source_order": source_order,
                "epsilon": "1/2700", "targets": targets, "preparation_ms": preparation_ms,
                "elapsed_ms": result.stats.elapsed.as_millis(),
                "exact_ms": result.stats.exact_materialization.as_millis(),
                "seeds": result.stats.seeds, "rows": result.stats.rows,
                "exact_trace_rows": result.stats.exact_trace_rows,
                "direct_rules": result.stats.direct_rules, "modular_rules": result.stats.modular_rules,
                "residuals": result.residuals.iter().map(|c| c.integral().powers().iter().map(|p| p.value()).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "rule_widths": rule_widths, "unique_rhs": rhs.len(), "rhs_by_lines": by_lines,
                "max_rhs_rank": max_rank, "max_rhs_dots": max_dots,
                "complete_reduction": false,
            })
        );
    } else {
        // A single rank-three numerator assignment with three symbolic positive
        // powers, covering its exceptional branches under explicit work limits.
        let target = targets[1];
        let case = CoordinateCase::new(std::array::from_fn(|i| (!sector[i]).then_some(target[i])))?;
        let started = Instant::now();
        let result = solver.solve_domains(
            vec![case.into()],
            SectorSolveOptions {
                symbolic: search,
                numerical_depth: 3,
                max_symbolic_cases: Some(32),
                ..Default::default()
            },
        );
        match result {
            Ok(result) => println!(
                "{}",
                json!({
                    "mode": mode, "lines": lines, "source_order": source_order,
                    "epsilon": "1/2700", "preparation_ms": preparation_ms,
                    "elapsed_ms": started.elapsed().as_millis(), "rules": result.rules.len(),
                    "symbolic_cases": result.stats.symbolic_cases,
                    "symbolic_rows": result.stats.symbolic_rows,
                    "finite_residuals": result.finite_residuals.iter().map(|i| i.powers().iter().map(|p| p.value()).collect::<Vec<_>>()).collect::<Vec<_>>(),
                    "complete_reduction": false,
                })
            ),
            Err(error) => {
                println!(
                    "{}",
                    json!({
                        "mode": mode, "lines": lines, "source_order": source_order,
                        "epsilon": "1/2700", "preparation_ms": preparation_ms,
                        "elapsed_ms": started.elapsed().as_millis(), "error": error.to_string(),
                        "complete_reduction": false,
                    })
                );
                return Err(error.into());
            }
        }
    }
    Ok(())
}
