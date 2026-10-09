//! Integration checks starting from generated occupied sources and regions.
use super::*;
use crate::finite_density::guarded::GuardedDiscoveryOptions;
use std::collections::BTreeMap;

#[test]
fn native_occupied_endpoint_reconstructs_singular_weights_and_rejects_bad_limits() {
    let input: crate::finite_density::DensityInput = serde_json::from_str(include_str!(
        "../../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    let options = FlowOptions {
        digits: 20,
        guard_digits: 24,
        series_order: 60,
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let context = RunContext::default();
    let flow = PreparedOccupiedFlow::<7>::prepare(
        &input.prepare().unwrap(),
        &[0],
        &options,
        WeightedClosureOptions {
            max_rounds: 12,
            discovery: GuardedDiscoveryOptions {
                max_depth: 3,
                max_domains: 8192,
                sample_seed: 0,
            },
            ..Default::default()
        },
        &context,
    )
    .unwrap();
    let epsilon = Rational::from((4, 5));
    let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
    let expected = flow.evaluate(&epsilon, &options, &context).unwrap();
    assert!(p.norm(&expected[0]) > p.tolerance(10));
    let original = &flow.closed.reduced;
    let n = original.basis.len();
    assert!(n >= 2);
    let j = original
        .basis
        .iter()
        .position(|i| {
            original.targets[0].get(i).is_some_and(|weight| {
                crate::frobenius::valuation(weight, flow.closed.variable).unwrap() <= 0
            })
        })
        .unwrap();
    let k = usize::from(j == 0);
    let eta = Atom::var(flow.closed.variable);
    // Exact invertible change of basis on eta>0:
    // new[j]=old[k], new[k]=old[k]-eta*old[j].
    // Thus old[j]=(new[j]-new[k])/eta. The two separate endpoint
    // components agree at zero, but their difference carries the target.
    let mut backward = (0..n)
        .map(|i| {
            (0..n)
                .map(|l| Atom::num(i64::from(i == l)))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    backward[j][j] = eta.clone().pow(-1);
    backward[j][k] = -eta.clone().pow(-1);
    backward[k][j] = Atom::num(1);
    backward[k][k] = Atom::num(0);
    let system = flow
        .system
        .as_ref()
        .unwrap()
        .change_basis(&backward)
        .unwrap();
    let mut reduced = original.clone();
    reduced.matrix = system.matrix.clone();
    reduced.nonzero_conditions.push(eta.clone());
    reduced.targets = original
        .targets
        .iter()
        .map(|target| {
            (0..n)
                .filter_map(|column| {
                    let weight = (0..n)
                        .fold(Atom::num(0), |sum, row| {
                            sum + target
                                .get(&original.basis[row])
                                .cloned()
                                .unwrap_or_default()
                                * &backward[row][column]
                        })
                        .together()
                        .cancel();
                    (!weight.is_zero()).then(|| (original.basis[column].clone(), weight))
                })
                .collect::<BTreeMap<_, _>>()
        })
        .collect();
    assert!(
        reduced.targets[0]
            .values()
            .any(|weight| crate::frobenius::valuation(weight, flow.closed.variable).unwrap() < 0)
    );
    let transport = ConnectionTransport::default();
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let boundary = OccupiedFlowBoundary::new(
        &backend,
        &options,
        &context,
        flow.epsilon,
        60,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let evaluate = |reduced: &crate::reduction::ReducedSystem, start_scale| {
        transport.evaluate(
            ConnectionRequest {
                system: &system,
                reduced,
                epsilon: flow.epsilon,
                loops: flow.family.loops(),
                positive_mass_contour: true,
                restrict_to_right_half_plane: true,
                start_scale,
            },
            &epsilon,
            &options,
            &context,
            |_| Ok(()),
            |new_basis, precision| {
                // At infinity z=1/eta, old=T(z)*new is a polynomial map.
                // Map every known coefficient and every log power exactly. Keep
                // the common known depth; no unknown tail is filled with zeros.
                let mut old_basis = new_basis.clone();
                for (old, new) in old_basis.columns.iter_mut().zip(&new_basis.columns) {
                    for (order, logs) in old.coefficients.iter_mut().enumerate() {
                        let current = &new.coefficients[order];
                        let previous = order.checked_sub(1).map(|r| &new.coefficients[r]);
                        let log_count = current.len().max(previous.map_or(0, Vec::len));
                        logs.resize(log_count, vec![precision.zero(); n]);
                        for (log, row) in logs.iter_mut().enumerate() {
                            row[k] = current
                                .get(log)
                                .map_or_else(|| precision.zero(), |v| v[j].clone());
                            row[j] = previous
                                .and_then(|v| v.get(log))
                                .map_or_else(|| precision.zero(), |v| precision.sub(&v[j], &v[k]));
                        }
                    }
                }
                Ok(boundary
                    .constants(
                        &flow.family,
                        &original.basis,
                        &flow.shifted,
                        &epsilon,
                        precision,
                        &old_basis,
                    )?
                    .constants)
            },
        )
    };
    let d = Rational::from(options.dimension) - &epsilon * &Rational::from(2);
    let normalization = p
        .eval(
            &native_measure_to_euclidean(2, 1, &Atom::num(d)).unwrap(),
            &ahash::HashMap::default(),
        )
        .unwrap();
    for start_scale in [8, 12] {
        let actual = evaluate(&reduced, start_scale).unwrap();
        for (value, expected) in actual.iter().zip(&expected) {
            assert!(p.close(&p.mul(value, &normalization), expected, 16));
        }
    }
    // Same source-derived solution and integrated boundary, now requesting a
    // genuinely divergent observable. The nonzero scalar target proves this
    // failure is not an arbitrary component with an accidental zero limit.
    let mut divergent = reduced.clone();
    divergent.targets.truncate(1);
    for weight in divergent.targets[0].values_mut() {
        *weight = (&*weight / &eta).together().cancel();
    }
    assert!(
        matches!(evaluate(&divergent, 8), Err(Error::Numerical(message))
        if message.contains("uncancelled physical endpoint divergence"))
    );
    let mut shallow = reduced;
    shallow.targets.truncate(1);
    for weight in shallow.targets[0].values_mut() {
        *weight = (&*weight * eta.clone().pow(-200)).together().cancel();
    }
    assert!(
        matches!(evaluate(&shallow, 8), Err(Error::Accuracy(message))
        if message.contains("endpoint series does not reach"))
    );
}
