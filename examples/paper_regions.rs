//! Evaluate the first large-mass region coefficients of the four paper targets.
//! This checks boundary construction without reducing the full target family.
use symbolica::prelude::*;
use symbolica_amflow::*;

fn main() -> Result<()> {
    let mut parametric_rules = false;
    let mut symmetry_rules = false;
    let mut max_exact_frontier = 0;
    let mut max_backward_frontier = RustRedBackend::default().max_backward_frontier;
    let mut max_targets = 32768;
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--parametric-rules" => parametric_rules = true,
            "--symmetry-rules" => symmetry_rules = true,
            "--max-exact-frontier" | "--max-backward-frontier" | "--max-targets" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| Error::InvalidInput(format!("missing value for {argument}")))?
                    .parse::<usize>()
                    .map_err(|e| Error::InvalidInput(format!("{argument}: {e}")))?;
                match argument.as_str() {
                    "--max-exact-frontier" => max_exact_frontier = value,
                    "--max-backward-frontier" => max_backward_frontier = value,
                    _ => max_targets = value,
                }
            }
            "--help" => {
                println!(
                    "paper_regions [--parametric-rules] [--symmetry-rules] [--max-exact-frontier N] [--max-backward-frontier N] [--max-targets N]"
                );
                return Ok(());
            }
            _ => return Err(Error::InvalidInput(format!("unknown option {argument}"))),
        }
    }
    let (family, targets) = benchmarks::paper_two_loop()?;
    let options = FlowOptions {
        cache_directory: Some(".amflow-cache/regions".into()),
        ..Default::default()
    };
    let backend = RustRedBackend {
        parametric_rules,
        symmetry_rules,
        max_depth: 3,
        max_targets,
        max_sector_batch: 32,
        max_exact_frontier,
        max_backward_frontier,
        checkpoints: Some(".amflow-cache/regions".into()),
        ..Default::default()
    };
    let context = RunContext {
        progress: Some(std::sync::Arc::new(|event| {
            if !matches!(event, Progress::Step { .. }) {
                eprintln!("{event:?}");
            }
        })),
        ..Default::default()
    };
    let evaluator = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let (_, shifted) = family.deform(symbol!("paper_regions_eta"), &options.mass_mode)?;
    let regions = regions::enumerate_regions(&family, 10000)?;
    let epsilon = Rational::from((1, 2700));
    for (index, target) in targets.iter().enumerate() {
        for region in &regions {
            let expansion = regions::expand_region(&family, target, &shifted, region, 2)?;
            for (order, coefficient) in expansion.coefficients.iter().enumerate() {
                let factors = integrand::factor_region(
                    coefficient,
                    &expansion.coordinates,
                    &family,
                    &region.hard,
                    10000,
                )?;
                let mut values = Vec::new();
                for digits in [60, 80] {
                    let p = Precision::decimal(digits)?;
                    let parameters = ahash::HashMap::from_iter([(
                        Atom::var(family.epsilon),
                        p.rational(&epsilon),
                    )]);
                    let mut value = p.zero();
                    for term in &factors {
                        if term
                            .factors
                            .iter()
                            .map(|f| recursive::scaleless(&f.family, &f.integral))
                            .collect::<Result<Vec<_>>>()?
                            .contains(&true)
                        {
                            continue;
                        }
                        let mut v = p.eval(&term.coefficient, &parameters)?;
                        for f in &term.factors {
                            v = p.mul(
                                &v,
                                &evaluator.evaluate(&f.family, &f.integral, &epsilon, p)?,
                            );
                        }
                        value = p.add(&value, &v);
                    }
                    let determinant = p.eval(&expansion.jacobian_determinant, &parameters)?;
                    let jacobian = p.pow(
                        &ComplexFloat::new(p.norm(&determinant), p.real(0)),
                        &p.rational(
                            &(Rational::from(family.dimension) - &epsilon * &Rational::from(2)),
                        ),
                    );
                    values.push(p.mul(&value, &jacobian));
                }
                let p = Precision::decimal(80)?;
                if !p.close(&values[0], &values[1], 20) {
                    return Err(Error::Accuracy(
                        "region coefficient precision refinement".into(),
                    ));
                }
                if values[1] != p.zero() {
                    println!(
                        "target {index}, hard {:?}, eta^({}), t^{order}: {}",
                        region.hard_branches, expansion.eta_power, values[1]
                    );
                }
            }
        }
    }
    println!("all region coefficients stable to 20 digits at increased precision/order");
    Ok(())
}
