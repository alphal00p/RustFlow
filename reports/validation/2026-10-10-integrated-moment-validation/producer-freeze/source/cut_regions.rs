//! Region coordinates compatible with an oriented compact final-state measure.
use crate::algebra::{determinant, rref};
use crate::cuts::{CutDefinition, CutFamily, LoopPrescription};
use crate::regions::LoopRegion;
use crate::*;
use symbolica::prelude::*;

/// Keep the native branch enumeration, then choose an equivalent basis for
/// each admitted hard space which fixes every original real loop coordinate.
/// Real momenta have compact support and cannot acquire a hard component.
pub(crate) fn enumerate(family: &CutFamily, context: &RunContext) -> Result<Vec<LoopRegion>> {
    let ordinary = family.family();
    let count = ordinary.loops.len();
    let virtual_slots = family
        .cuts()
        .loop_prescriptions()
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (*p == LoopPrescription::PlusI0).then_some(i))
        .collect::<Vec<_>>();
    let mut out = Vec::new();
    for region in regions::enumerate_regions(ordinary, 10000)? {
        context.cancellation.check()?;
        let hard_columns = (0..count).filter(|&i| region.hard[i]).collect::<Vec<_>>();
        let mut admitted = true;
        for definition in family.cuts().lines().values() {
            let CutDefinition::PositiveEnergy { momentum } = definition else {
                return Err(Error::Unsupported(
                    "region analysis needs positive-energy cuts".into(),
                ));
            };
            for &column in &hard_columns {
                let component = (0..count)
                    .fold(Atom::zero(), |sum, i| {
                        sum + Atom::num(momentum.loops[i].clone())
                            * &region.transformation[i][column]
                    })
                    .together()
                    .cancel();
                admitted &= component.is_zero();
            }
        }
        if !admitted {
            continue;
        }
        // The complete cut measure spans the real subspace. Check the stronger
        // coordinate statement as well; do not infer it numerically.
        if (0..count).any(|i| {
            !virtual_slots.contains(&i)
                && hard_columns
                    .iter()
                    .any(|&j| !region.transformation[i][j].is_zero())
        }) {
            return Err(Error::Unsupported(
                "hard cut region changes a real integration direction".into(),
            ));
        }
        out.push(adapt_compact_region(region, count, &virtual_slots)?);
    }
    if out.is_empty() {
        return Err(Error::Unsupported(
            "no regions preserve the complete final-state measure".into(),
        ));
    }
    context.cancellation.check()?;
    Ok(out)
}

/// Enumerate correlated native branch regions with every occupied coordinate
/// fixed and soft. This shares the final-state owner's exact virtual-space
/// projection; it assumes no total final-state momentum constraint.
pub(crate) fn enumerate_compact(
    ordinary: &IntegralFamily,
    compact: &[bool],
    context: &RunContext,
) -> Result<Vec<LoopRegion>> {
    let count = ordinary.loops.len();
    if compact.len() != count {
        return Err(Error::InvalidInput("compact loop mask dimensions".into()));
    }
    let virtual_slots = (0..count).filter(|&i| !compact[i]).collect::<Vec<_>>();
    let mut out = Vec::new();
    for region in regions::enumerate_regions(ordinary, 10000)? {
        context.cancellation.check()?;
        if (0..count).any(|i| {
            compact[i]
                && (0..count).any(|j| region.hard[j] && !region.transformation[i][j].is_zero())
        }) {
            continue;
        }
        out.push(adapt_compact_region(region, count, &virtual_slots)?);
    }
    if out.is_empty() {
        return Err(Error::Unsupported(
            "no regions preserve the occupied coordinates".into(),
        ));
    }
    Ok(out)
}

fn adapt_compact_region(
    region: LoopRegion,
    count: usize,
    virtual_slots: &[usize],
) -> Result<LoopRegion> {
    let hard_columns = (0..count).filter(|&i| region.hard[i]).collect::<Vec<_>>();
    let mut basis = hard_columns
        .iter()
        .map(|&j| {
            virtual_slots
                .iter()
                .map(|&i| region.transformation[i][j].clone())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let hard_count = basis.len();
    if rref(basis.clone()).1.len() != hard_count {
        return Err(Error::InvalidInput(
            "dependent virtual hard directions".into(),
        ));
    }
    // Project the enumerator's remaining columns onto the virtual space.
    // Complete the hard basis without mixing real and virtual measures.
    for column in 0..count {
        let vector = virtual_slots
            .iter()
            .map(|&i| region.transformation[i][column].clone())
            .collect::<Vec<_>>();
        let mut candidate = basis.clone();
        candidate.push(vector.clone());
        if rref(candidate).1.len() > basis.len() {
            basis.push(vector);
        }
    }
    if basis.len() != virtual_slots.len() {
        return Err(Error::InvalidInput(
            "virtual region coordinates do not span their measure".into(),
        ));
    }
    let mut transformation = vec![vec![Atom::zero(); count]; count];
    let mut hard = vec![false; count];
    for (i, row) in transformation.iter_mut().enumerate() {
        if !virtual_slots.contains(&i) {
            row[i] = Atom::one();
        }
    }
    for (column, vector) in basis.iter().enumerate() {
        let slot = virtual_slots[column];
        hard[slot] = column < hard_count;
        for (row, &original) in virtual_slots.iter().enumerate() {
            transformation[original][slot] = vector[row].clone();
        }
    }
    let jacobian_determinant = determinant(transformation.clone());
    if jacobian_determinant.is_zero() {
        return Err(Error::InvalidInput("singular adapted cut region".into()));
    }
    Ok(LoopRegion {
        transformation,
        hard,
        hard_branches: region.hard_branches,
        jacobian_determinant,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cuts::{CutLine, CutMetadata, MomentumRouting};

    #[test]
    fn virtual_difference_can_stay_soft_with_exact_real_coordinates() -> Result<()> {
        let gram = vec![vec![Atom::num(25)]];
        let specifications = [
            ([1, 0, 0], 0, 1),
            ([-1, 0, 0], 1, 1),
            ([0, 1, 0], 0, 2),
            ([0, 0, 1], 0, 3),
            ([0, 1, -1], 0, 4),
            ([-1, 1, 0], 0, 0),
            ([-1, 0, 1], 0, 0),
            ([0, -1, 0], 1, 0),
            ([0, 0, -1], 1, 0),
        ];
        let ordinary = IntegralFamily {
            name: "cut_region_common_virtual_hard".into(),
            loops: vec!["r".into(), "k".into(), "l".into()],
            external: vec!["P".into()],
            external_gram: gram.clone(),
            propagators: specifications
                .iter()
                .map(|(q, e, m)| Propagator::quadratic(q, &[*e], Atom::num(*m), &gram))
                .collect::<Result<_>>()?,
            physical_propagators: 5,
            epsilon: symbol!("cut_region::eps"),
            dimension: 4,
        };
        let family = CutFamily::new(
            ordinary,
            CutMetadata::new(
                [([1, 0, 0], 0), ([-1, 0, 0], 1)]
                    .into_iter()
                    .enumerate()
                    .map(|(propagator, (q, e))| CutLine {
                        propagator,
                        definition: CutDefinition::PositiveEnergy {
                            momentum: MomentumRouting {
                                loops: q.into_iter().map(Rational::from).collect(),
                                external: vec![e.into()],
                            },
                        },
                    }),
                vec![
                    LoopPrescription::Insensitive,
                    LoopPrescription::PlusI0,
                    LoopPrescription::PlusI0,
                ],
            )?,
        )?;
        let mut common_hard = false;
        for region in enumerate(&family, &RunContext::default())? {
            assert!(!region.hard[0]);
            assert_eq!(
                region.transformation[0],
                vec![Atom::one(), Atom::zero(), Atom::zero()]
            );
            assert_eq!(region.transformation[1][0], Atom::zero());
            assert_eq!(region.transformation[2][0], Atom::zero());
            assert_eq!(
                region.jacobian_determinant,
                determinant(region.transformation.clone())
            );
            let transformed = family.family().transform_loops(&region.transformation)?;
            let hard_branch = |slot| -> Result<bool> {
                Ok(regions::branch(&transformed.propagators[slot], 3)?
                    .iter()
                    .zip(&region.hard)
                    .any(|(a, h)| *h && !a.is_zero()))
            };
            common_hard |= hard_branch(2)? && hard_branch(3)? && !hard_branch(4)?;
        }
        assert!(common_hard, "missed the hard k,l with soft k-l region");
        Ok(())
    }
}
