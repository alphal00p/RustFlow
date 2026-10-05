use super::*;
use crate::reduction::{LinearCombination, Reduction};
use std::sync::atomic::{AtomicUsize, Ordering};

struct NoReduction;
impl ReductionBackend for NoReduction {
    fn identity(&self) -> String {
        "supplied-connection-test".into()
    }
    fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
        Err(Error::Reduction(
            "unexpected IBP search in supplied test".into(),
        ))
    }
}

#[test]
fn supplied_target_maps_share_samples_and_cancel_before_fitting() -> Result<()> {
    let epsilon = symbol!("supplied_projection_eps");
    let mass = symbol!("supplied_projection_mass");
    let eta = symbol!("supplied_projection_eta");
    let family = IntegralFamily {
        name: "supplied_projection_tadpole".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -Atom::var(mass),
            scalar_products: vec![Atom::one()],
        }],
        physical_propagators: 1,
        epsilon,
        dimension: 4,
    };
    let first = Integral(vec![1]);
    let second = Integral(vec![2]);
    let preparations = AtomicUsize::new(0);
    let prepare = |family: &IntegralFamily,
                   targets: &[Integral],
                   options: &FlowOptions,
                   context: &RunContext| {
        preparations.fetch_add(1, Ordering::Relaxed);
        assert_eq!(family.propagators[0].constant, Atom::num(-2));
        let derivative = (Atom::one() - Atom::var(epsilon)) / (Atom::num(2) + Atom::var(eta));
        let maps = BTreeMap::from([
            (
                first.clone(),
                LinearCombination::from([(first.clone(), Atom::one())]),
            ),
            (
                second.clone(),
                LinearCombination::from([(first.clone(), derivative.clone())]),
            ),
        ]);
        PreparedFlow::from_supplied(
            family,
            targets,
            &KinematicPoint::default(),
            SuppliedAuxiliarySystem {
                variable: eta,
                reduced: ReducedSystem {
                    basis: vec![first.clone()],
                    matrix: vec![vec![derivative]],
                    targets: targets.iter().map(|target| maps[target].clone()).collect(),
                    nonzero_conditions: vec![Atom::num(2) + Atom::var(eta)],
                    candidates: maps,
                    transformations: vec![],
                },
                deformation_mask: vec![true],
                provenance: "analytic tadpole scaling and exact power recurrence".into(),
            },
            options,
            context,
        )
    };
    let rows = vec![
        LinearCombination::from([(first.clone(), Atom::one())]),
        LinearCombination::from([
            (first.clone(), Atom::one()),
            (second.clone(), Atom::num(-2)),
        ]),
    ];
    let fits = solve_integral_projections_with_preparer(
        &[(family, rows)],
        &KinematicPoint(BTreeMap::from([(Atom::var(mass), Atom::num(2))])),
        0,
        &FlowOptions::default(),
        &NoReduction,
        &RunContext::default(),
        &Default::default(),
        Some(&prepare),
    )?;
    assert_eq!(preparations.load(Ordering::Relaxed), 1);
    let p = Precision::decimal(80)?;
    assert_eq!(fits[0].verified_digits, Some(20));
    assert_eq!(fits[1].verified_digits, Some(20));
    assert!(p.close(&fits[0].coefficients[&-1], &p.i(2), 20));
    assert!(p.close(&fits[1].coefficients[&-1], &p.zero(), 20));
    // I1 - 2 I2 = epsilon I1 at m²=2, so its finite term is exactly 2.
    assert!(p.close(&fits[1].coefficients[&0], &p.i(2), 20));
    Ok(())
}
