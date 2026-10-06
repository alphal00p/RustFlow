//! Read-only experiment: current finite-epsilon kernel and physical error owner.
use serde_json::json;
use std::collections::BTreeMap;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::{
    diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon},
    kinematics::KinematicSystem,
    transport_cache::*,
    *,
};

fn main() -> Result<()> {
    let x = symbol!("epsilon_singular_probe::x");
    let epsilon = symbol!("epsilon_singular_probe::epsilon");
    let e = Atom::var(epsilon);
    let original = DifferentialSystem {
        variable: x,
        matrix: vec![
            vec![Atom::one() / &e, Atom::one() / &e],
            vec![-Atom::one() / &e, -Atom::one() / &e],
        ],
    };
    let gauge = vec![
        vec![Atom::one() / &e, Atom::one()],
        vec![-Atom::one() / &e, Atom::new()],
    ];
    let regular = original.change_basis(&gauge)?;
    assert_eq!(
        regular.matrix,
        vec![
            vec![Atom::new(), Atom::one()],
            vec![Atom::new(), Atom::new()]
        ]
    );
    let unregularized = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(x, original.matrix)]),
    };
    let diagonal_error = match EpsilonShearing::regularize(&unregularized, &RunContext::default()) {
        Err(Error::Unsupported(message)) => message,
        _ => panic!("nilpotent diagonal epsilon poles must not admit diagonal shearing"),
    };
    let options = FlowOptions::default();
    let p = Precision::decimal(90)?;
    let basis = [
        symbol!("epsilon_singular_probe::I").call(1),
        symbol!("epsilon_singular_probe::I").call(2),
    ];
    let mut cases = Vec::new();
    for denominator in [100_i64, 2500, 5800] {
        let sample = Rational::from((1, denominator));
        for distance in [Rational::from(1), Rational::from((1, 1000))] {
            let matrix = vec![
                vec![Atom::num(denominator), Atom::num(denominator)],
                vec![Atom::num(-denominator), Atom::num(-denominator)],
            ];
            let system = EpsilonSystem {
                variable: x,
                matrices: vec![matrix.clone()],
            };
            let endpoint = Atom::num(distance.clone());
            let step = p.rational(&distance);
            let expected = [
                p.add(&p.i(1), &p.scale(&step, denominator, 1)),
                p.scale(&step, -denominator, 1),
            ];
            let start = Instant::now();
            let low = transport_epsilon(
                &system,
                |p| {
                    Ok(EpsilonBoundary {
                        point: p.zero(),
                        leading: 0,
                        coefficients: vec![vec![p.i(1), p.zero()]],
                    })
                },
                std::slice::from_ref(&endpoint),
                &options,
                &RunContext::default(),
                true,
            )?;
            let low_seconds = start.elapsed().as_secs_f64();
            for (value, expected) in low.coefficients[0].iter().zip(&expected) {
                assert!(p.close(value, expected, 50));
            }
            let amplification = system
                .compile(p, &Default::default())?
                .error_amplification(&p.zero(), &step)?;
            let log10_amplification =
                p.log(&Complex::new(amplification.clone(), p.real(0))).re / p.log(&p.i(10)).re;
            let exact_transfer_norm = Rational::one() + Rational::from(2 * denominator) * &distance;
            let exact_input_error = exact_transfer_norm.to_multi_prec_float(p.bits)
                * (p.tolerance(30) + p.tolerance(35));
            let flow = RustFlow::new(
                KinematicSystem {
                    epsilon,
                    derivatives: BTreeMap::from([(x, matrix)]),
                },
                &basis,
                &Atom::one(),
                Prescription::PlusI0,
                &format!("real x; exact epsilon={sample}"),
            )?;
            let mut cache = RustFlowCache::default();
            cache.insert(CachedBoundary { identity:flow.identity().clone(),point:CachedPoint::Exact(BTreeMap::from([(x,Atom::new())])),kind:PointKind::Physical,range:EpsilonRange::new(0,0)?,coefficients:vec![vec![p.i(1),p.zero()]],accuracy:BoundaryAccuracy::supplied(30,p.bits,vec![vec![p.tolerance(35);2]],"Probe source: 30-digit evidence, component errors 1e-35, 90-digit storage; no exact-boundary shortcut")? })?;
            let policy = ScaledDistance {
                scales: BTreeMap::new(),
                admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
            };
            let start = Instant::now();
            let physical = match flow.evaluate_to(
                &mut cache,
                &BTreeMap::from([(x, endpoint)]),
                EpsilonRange::new(0, 0)?,
                &options,
                &RunContext::default(),
                &policy,
            ) {
                Ok(result) => {
                    for (value, expected) in result.boundary.coefficients[0].iter().zip(&expected) {
                        assert!(p.close(value, expected, 50));
                    }
                    json!({"status":"accepted","verified_digits":result.boundary.accuracy.verified_digits(),"component_errors":result.boundary.accuracy.comparison_errors()[0].iter().map(ToString::to_string).collect::<Vec<_>>(),"steps":result.transport.as_ref().map(|s|s.diagnostics.steps),"inserted_points":result.inserted_points,"cache_entries":cache.len()})
                }
                Err(error) => {
                    assert_eq!(cache.len(), 1, "failed route must not mutate sample bank");
                    json!({"status":"rejected","error":error.to_string(),"error_debug":format!("{error:?}"),"cache_entries":cache.len()})
                }
            };
            let physical_seconds = start.elapsed().as_secs_f64();
            cases.push(json!({"epsilon":sample.to_string(),"distance":distance.to_string(),"expected":expected.iter().map(ToString::to_string).collect::<Vec<_>>(),"low_level":{"verified_digits":low.verified_digits,"steps":low.diagnostics.steps,"rejected_steps":low.diagnostics.rejected_steps,"working_bits":low.diagnostics.working_bits,"order":low.diagnostics.expansion_order,"seconds":low_seconds,"does_not_include_boundary_uncertainty":true},"unchanged_error_amplification":amplification.to_string(),"log10_error_amplification":log10_amplification.to_string(),"exact_nilpotent_transfer_infinity_norm":exact_transfer_norm.to_string(),"exact_transfer_propagated_source_error_for_comparison_only":exact_input_error.to_string(),"physical":physical,"physical_seconds":physical_seconds}));
        }
    }
    println!("{}",serde_json::to_string_pretty(&json!({"base_commit":"6a774c0","scope":"Unchanged engine probe; exact rational specialization at epsilon, full existing physical source uncertainty retained; closed-form transfer bound is comparison only and not supplied to controller.","source_verified_digits":30,"source_error_per_component":"1e-35","source_storage_decimal_digits":90,"options":format!("{options:?}"),"diagonal_shearing_error":diagonal_error,"exact_nondiagonal_gauge_verified":true,"cases":cases})).unwrap());
    Ok(())
}
