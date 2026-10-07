//! Bounded discovery and exact transport of routing automorphisms.
//!
//! Branch displacements suggest loop translations and external reflections;
//! they do not prove an integral identity. RustRed verifies each momentum map
//! and compiles its action on the complete physical root. Concrete transport
//! expands only finite numerator powers, and every retained term must strictly
//! descend under the caller's actual reduction order. Any unsupported map,
//! exhausted search, or transport limit leaves the ordinary reducer available.

use std::cmp::Ordering;
use std::sync::Arc;

use rustred::algebra::{Coefficient, CoefficientContext};
use rustred::family::{IntegralFamily, IntegralKey, ScalarProductCoordinate};
use rustred::sector::Mask;
use rustred::sector::symmetry::integral_transport::{self, ExpansionLimits, Prepared};
use rustred::sector::symmetry::{self, CoefficientMatrix, DenominatorAction, MomentumMap};
use rustred::solver::{Integral, IntegralOrder, Term};
use symbolica::domains::rational::RationalField;
use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
use symbolica::prelude::*;
use symbolica::tensors::matrix::Matrix;

/// Bounds for an optional accelerator, not limits on integral evaluation.
#[derive(Clone, Copy)]
pub(super) struct Limits {
    pub max_candidates: usize,
    pub max_external_maps: usize,
    pub max_denominators: usize,
    pub max_momenta: usize,
    pub verification: symmetry::Limits,
    pub expansion: ExpansionLimits,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_candidates: 256,
            max_external_maps: 32,
            max_denominators: 64,
            max_momenta: 16,
            verification: symmetry::Limits {
                max_matrix_entries: 16_384,
                max_exact_operations: 1_000_000,
                max_symbolica_single_matrix_entries: 16_384,
                max_symbolica_live_matrix_entries: 65_536,
                max_nonzero_conditions: 4096,
                max_condition_sources: 16_384,
                ..Default::default()
            },
            expansion: ExpansionLimits {
                max_total_power: 64,
                max_endpoints: 512,
                max_native_polynomial_terms: 4096,
                ..Default::default()
            },
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Statistics {
    pub candidates: usize,
    pub automorphisms: usize,
    pub applied: usize,
    pub missed: usize,
    pub transport_failures: usize,
    pub search_limited: bool,
}

pub(super) struct Bank {
    fingerprint: String,
    arity: usize,
    maps: Vec<Prepared>,
    expansion: ExpansionLimits,
    pub statistics: Statistics,
}

impl Bank {
    /// Discovery is deliberately incomplete. Symbolic Gram entries, non-square
    /// branches, and missing independent branch anchors simply give no maps.
    pub fn discover(family: &IntegralFamily, physical: usize, limits: Limits) -> Self {
        let mut bank = Self {
            fingerprint: family.fingerprint().to_owned(),
            arity: family.denominator_count(),
            maps: Vec::new(),
            expansion: limits.expansion,
            statistics: Statistics::default(),
        };
        if physical > bank.arity
            || bank.arity > limits.max_denominators
            || family.loop_count() > limits.max_momenta
            || family.external_count() > limits.max_momenta
        {
            bank.statistics.search_limited = true;
            return bank;
        }
        if family.loop_count() == 0 || family.external_count() == 0 {
            return bank;
        }
        let Some(candidates) = candidates(family, physical, limits) else {
            return bank;
        };
        bank.statistics.search_limited = candidates.limited;
        // IntegralFamily intentionally is not Clone. Rebuild through its public
        // validator while retaining the exact coefficient context and ordering.
        let Ok(target) = IntegralFamily::new(
            family.name().to_owned(),
            family.loop_momenta().to_vec(),
            family.external_momenta().to_vec(),
            family.coefficient_context().clone(),
            family.dimension().clone(),
            family.denominators().to_vec(),
            family.external_gram().to_vec(),
            family.power_shifts().to_vec(),
        ) else {
            return bank;
        };
        let target = Arc::new(target);
        let Ok(root) = Mask::try_new((0..bank.arity).map(|i| i < physical)) else {
            return bank;
        };
        for candidate in candidates.maps {
            bank.statistics.candidates += 1;
            let Ok(map) = symmetry::verify(family, &target, candidate, limits.verification) else {
                continue;
            };
            // A momentum symmetry that fixes every denominator cannot yield a
            // strictly descending integral. Discard it before concrete calls.
            if map.row_actions().iter().enumerate().all(|(axis, action)| {
                matches!(action, DenominatorAction::Monomial { target, scale }
                    if *target == axis && *scale == family.coefficient_context().one())
            }) {
                continue;
            }
            // This checks masses and eta, a unit Jacobian, a unit bijection of
            // ALL physical denominators, rational affine numerator rows, and
            // every nonzero condition. Unresolved symbolic guards are rejected;
            // no condition is discarded when a transport is returned below.
            if let Ok(prepared) = integral_transport::compile(
                family,
                Arc::clone(&target),
                Arc::new(map),
                root.clone(),
                root.clone(),
                limits.expansion,
            ) {
                bank.maps.push(prepared);
            }
        }
        bank.statistics.automorphisms = bank.maps.len();
        bank
    }

    /// Return a complete exact identity or no rule. Partial transport results,
    /// overflowing compact keys, and nondecreasing terms are never installed.
    pub fn apply<const N: usize>(
        &mut self,
        family: &IntegralFamily,
        powers: [i16; N],
        order: &IntegralOrder<N>,
    ) -> Option<Vec<Term<N, Coefficient>>> {
        let result = self.apply_inner(family, powers, order);
        if result.is_some() {
            self.statistics.applied += 1;
        } else {
            self.statistics.missed += 1;
        }
        result
    }

    fn apply_inner<const N: usize>(
        &mut self,
        family: &IntegralFamily,
        powers: [i16; N],
        order: &IntegralOrder<N>,
    ) -> Option<Vec<Term<N, Coefficient>>> {
        if self.arity > N
            || powers[self.arity..].iter().any(|&value| value != 0)
            || self.fingerprint != family.fingerprint()
        {
            return None;
        }
        let lhs = Integral::numeric(powers).ok()?;
        let key = IntegralKey::try_new(powers[..self.arity].iter().copied().map(i64::from)).ok()?;
        for map in &self.maps {
            let transported = match map.transport(&key, self.expansion) {
                Ok(value) => value,
                Err(_) => {
                    self.statistics.transport_failures += 1;
                    continue;
                }
            };
            let terms = transported
                .terms()
                .iter()
                .map(|term| {
                    let mut values = term
                        .key()
                        .powers()
                        .iter()
                        .copied()
                        .map(i16::try_from)
                        .collect::<Result<Vec<_>, _>>()
                        .ok()?;
                    values.resize(N, 0);
                    let integral = Integral::numeric(values.try_into().ok()?).ok()?;
                    (order.compare(&lhs, &integral) == Ordering::Less).then(|| Term {
                        integral,
                        coefficient: term.coefficient().clone(),
                    })
                })
                .collect::<Option<Vec<_>>>();
            if let Some(terms) = terms {
                return Some(terms);
            }
        }
        None
    }
}

type RationalMatrix = Matrix<RationalField>;

struct Branch {
    direction: Vec<Rational>,
    scale: Rational,
    displacements: Vec<Vec<Rational>>,
}

struct Candidates {
    maps: Vec<MomentumMap>,
    limited: bool,
}

fn constant(value: &Coefficient) -> Option<Rational> {
    value.is_constant().then(|| {
        Rational::from_unchecked(
            value.numerator.get_constant(),
            value.denominator.get_constant(),
        )
    })
}

fn native_matrix(
    context: &CoefficientContext,
    matrix: RationalMatrix,
) -> Option<CoefficientMatrix> {
    let rows = matrix.nrows();
    let columns = matrix.ncols();
    let template = context.one();
    CoefficientMatrix::try_new(
        rows,
        columns,
        matrix.into_vec().into_iter().map(|value| {
            Coefficient::from_num_den(
                template.numerator.constant(value.numerator()),
                template.numerator.constant(value.denominator()),
                &Z,
                false,
            )
        }),
    )
    .ok()
}

/// For D = scale * (v.l + w.p)^2 + constant, the sign convention is
/// v[pivot] = 1. Coefficients outside a rational rank-one branch are unsupported.
fn branches(family: &IntegralFamily, physical: usize) -> Option<Vec<Branch>> {
    let loops = family.loop_count();
    let external = family.external_count();
    let mut groups: Vec<Branch> = Vec::new();
    for denominator in &family.denominators()[..physical] {
        let mut square = vec![vec![Rational::zero(); loops]; loops];
        let mut mixed = vec![vec![Rational::zero(); external]; loops];
        for (coordinate, value) in family.coordinates().iter().zip(denominator.coefficients()) {
            let value = constant(value)?;
            match *coordinate {
                ScalarProductCoordinate::LoopLoop { left, right } => {
                    let value = if left == right {
                        value
                    } else {
                        value / Rational::from(2)
                    };
                    square[left][right] = value.clone();
                    square[right][left] = value;
                }
                ScalarProductCoordinate::LoopExternal {
                    loop_index,
                    external_index,
                } => {
                    mixed[loop_index][external_index] = value;
                }
            }
        }
        let pivot = (0..loops).find(|&i| !square[i][i].is_zero())?;
        let scale = square[pivot][pivot].clone();
        let direction = square[pivot].iter().map(|q| q / &scale).collect::<Vec<_>>();
        let twice_scale = &scale * &Rational::from(2);
        let displacement = mixed[pivot]
            .iter()
            .map(|q| q / &twice_scale)
            .collect::<Vec<_>>();
        if (0..loops).any(|i| {
            (0..loops).any(|j| square[i][j] != &scale * &direction[i] * &direction[j])
                || (0..external).any(|j| {
                    mixed[i][j] != Rational::from(2) * &scale * &direction[i] * &displacement[j]
                })
        }) {
            return None;
        }
        if let Some(group) = groups
            .iter_mut()
            .find(|g| g.scale == scale && g.direction == direction)
        {
            if !group.displacements.contains(&displacement) {
                group.displacements.push(displacement);
            }
        } else {
            groups.push(Branch {
                direction,
                scale,
                displacements: vec![displacement],
            });
        }
    }
    Some(groups)
}

fn candidates(family: &IntegralFamily, physical: usize, limits: Limits) -> Option<Candidates> {
    let mut result = Candidates {
        maps: Vec::new(),
        limited: false,
    };
    if limits.max_candidates == 0 || limits.max_external_maps == 0 {
        result.limited = true;
        return Some(result);
    }
    let loops = family.loop_count();
    let external = family.external_count();
    let groups = branches(family, physical)?;
    let gram = family
        .external_gram()
        .iter()
        .map(|row| row.iter().map(constant).collect::<Option<Vec<_>>>())
        .collect::<Option<Vec<_>>>()?;
    let mut external_maps = vec![Matrix::identity(external as u32, Q)];
    for group in &groups {
        for (i, left) in group.displacements.iter().enumerate() {
            for right in group.displacements.iter().skip(i + 1) {
                let difference = left
                    .iter()
                    .zip(right)
                    .map(|(l, r)| l - r)
                    .collect::<Vec<_>>();
                let gram_difference = gram
                    .iter()
                    .map(|row| {
                        row.iter()
                            .zip(&difference)
                            .fold(Rational::zero(), |sum, (g, d)| sum + g * d)
                    })
                    .collect::<Vec<_>>();
                let norm = difference
                    .iter()
                    .zip(&gram_difference)
                    .fold(Rational::zero(), |sum, (d, g)| sum + d * g);
                if norm.is_zero() {
                    continue;
                }
                let reflected = Matrix::from_nested_vec(
                    (0..external)
                        .map(|i| {
                            (0..external)
                                .map(|j| {
                                    Rational::from(i64::from(i == j))
                                        - Rational::from(2) * &gram_difference[i] * &difference[j]
                                            / &norm
                                })
                                .collect()
                        })
                        .collect(),
                    Q,
                )
                .ok()?;
                if !external_maps.contains(&reflected) {
                    if external_maps.len() == limits.max_external_maps {
                        result.limited = true;
                    } else {
                        external_maps.push(reflected);
                    }
                }
            }
        }
    }
    let mut anchors: Vec<&Branch> = Vec::new();
    for group in &groups {
        let rows = anchors
            .iter()
            .map(|b| b.direction.clone())
            .chain(std::iter::once(group.direction.clone()))
            .collect::<Vec<_>>();
        if Matrix::from_nested_vec(rows, Q).ok()?.rank() > anchors.len() {
            anchors.push(group);
        }
        if anchors.len() == loops {
            break;
        }
    }
    if anchors.len() != loops {
        return None;
    }
    let basis =
        Matrix::from_nested_vec(anchors.iter().map(|b| b.direction.clone()).collect(), Q).ok()?;
    let inverse = basis.inv().ok()?;
    let source = Matrix::from_nested_vec(
        anchors.iter().map(|b| b.displacements[0].clone()).collect(),
        Q,
    )
    .ok()?;
    let context = family.coefficient_context();
    for external_map in external_maps {
        let mapped = &source * &external_map;
        let mut destinations = vec![0; loops];
        loop {
            if result.maps.len() == limits.max_candidates {
                result.limited = true;
                return Some(result);
            }
            let target = Matrix::from_nested_vec(
                anchors
                    .iter()
                    .zip(&destinations)
                    .map(|(b, &i)| b.displacements[i].clone())
                    .collect(),
                Q,
            )
            .ok()?;
            let translation = &inverse * &(&target - &mapped);
            let candidate = MomentumMap::new(
                native_matrix(context, Matrix::identity(loops as u32, Q))?,
                native_matrix(context, translation)?,
                native_matrix(context, external_map.clone())?,
            );
            if !result.maps.contains(&candidate) {
                result.maps.push(candidate);
            }
            let mut axis = 0;
            while axis < loops {
                destinations[axis] += 1;
                if destinations[axis] < anchors[axis].displacements.len() {
                    break;
                }
                destinations[axis] = 0;
                axis += 1;
            }
            if axis == loops {
                break;
            }
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn paper(shifted: Vec<usize>) -> IntegralFamily {
        let (family, _) = crate::benchmarks::paper_two_loop().unwrap();
        let (family, _) = family
            .deform(
                symbol!("symmetry_tests::eta"),
                &crate::MassMode::Propagators(shifted),
            )
            .unwrap();
        family
            .convert_at_epsilon(Some(&Rational::from((1, 2700))))
            .unwrap()
            .family
    }

    fn order<const N: usize>(powers: [i16; N]) -> IntegralOrder<N> {
        IntegralOrder::new(powers.map(|n| n > 0), [false; N])
    }

    #[test]
    fn generic_equal_mass_bubble_uses_the_same_discovery() {
        let context = CoefficientContext::try_new(["symmetry_bubble_eta"]).unwrap();
        let mass = &context.integer(1) + &context.parameter("symmetry_bubble_eta").unwrap();
        let family = IntegralFamily::new(
            "equal_mass_bubble".to_owned(),
            vec!["l".to_owned()],
            vec!["p".to_owned()],
            context.clone(),
            context.integer(4),
            vec![
                rustred::family::AffineDenominator::new(
                    -mass.clone(),
                    vec![context.one(), context.zero()],
                ),
                rustred::family::AffineDenominator::new(
                    context.integer(3) - mass,
                    vec![context.one(), context.integer(2)],
                ),
            ],
            vec![vec![context.integer(3)]],
            vec![context.zero(); 2],
        )
        .unwrap();
        let mut bank = Bank::discover(&family, 2, Limits::default());
        assert_eq!(bank.statistics.automorphisms, 1);
        let terms = bank.apply(&family, [2, 0], &order([2, 0])).unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].integral, Integral::numeric([0, 2]).unwrap());
        assert_eq!(terms[0].coefficient, context.one());
        assert!(bank.apply(&family, [0, 2], &order([0, 2])).is_none());
    }

    #[test]
    fn discovers_routing_and_transports_the_exact_numerator_square() {
        let family = paper(vec![4]);
        let mut bank = Bank::discover(&family, 7, Limits::default());
        assert_eq!(bank.statistics.candidates, 36);
        assert_eq!(bank.statistics.automorphisms, 1);
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, -2];
        let terms = bank.apply(&family, powers, &order(powers)).unwrap();
        let found = terms
            .into_iter()
            .map(|t| {
                (
                    t.integral.powers().map(|p| p.value()),
                    constant(&t.coefficient).unwrap(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        // D9 -> D9 + 2 D6 - 2 D4. No denominator was expanded.
        let expected = BTreeMap::from([
            ([0, 0, 1, 0, 1, 0, 1, 0, -2], Rational::from(1)),
            ([0, 0, 1, 0, 1, -2, 1, 0, 0], Rational::from(4)),
            ([0, 0, 1, -2, 1, 0, 1, 0, 0], Rational::from(4)),
            ([0, 0, 1, 0, 1, -1, 1, 0, -1], Rational::from(4)),
            ([0, 0, 1, -1, 1, 0, 1, 0, -1], Rational::from(-4)),
            ([0, 0, 1, -1, 1, -1, 1, 0, 0], Rational::from(-8)),
        ]);
        assert_eq!(found, expected);
    }

    #[test]
    fn exact_transport_composed_twice_is_the_identity() {
        let family = paper(vec![4]);
        let bank = Bank::discover(&family, 7, Limits::default());
        let key = IntegralKey::try_new([1, -1, 0, -1, 2, 0, 1, -1, -3]).unwrap();
        let map = &bank.maps[0];
        let first = map.transport(&key, bank.expansion).unwrap();
        let mut composition = BTreeMap::<Vec<i64>, Coefficient>::new();
        for term in first.terms() {
            for output in map.transport(term.key(), bank.expansion).unwrap().terms() {
                let coefficient = term.coefficient() * output.coefficient();
                let entry = composition
                    .entry(output.key().powers().to_vec())
                    .or_insert_with(|| family.coefficient_context().zero());
                *entry = &*entry + &coefficient;
            }
        }
        composition.retain(|_, c| !c.is_zero());
        assert_eq!(
            composition,
            BTreeMap::from([(key.powers().to_vec(), family.coefficient_context().one())])
        );
    }

    #[test]
    fn noninvariant_eta_and_gram_and_foreign_family_are_rejected() {
        let family = paper(vec![4]);
        let mut bank = Bank::discover(&family, 7, Limits::default());
        let changed_eta = paper(vec![0]);
        assert_eq!(
            Bank::discover(&changed_eta, 7, Limits::default())
                .statistics
                .automorphisms,
            0
        );
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, 0];
        assert!(bank.apply(&changed_eta, powers, &order(powers)).is_none());
        let mut gram = family.external_gram().to_vec();
        gram[0][0] = family.coefficient_context().integer(2);
        let changed_gram = IntegralFamily::new(
            family.name().to_owned(),
            family.loop_momenta().to_vec(),
            family.external_momenta().to_vec(),
            family.coefficient_context().clone(),
            family.dimension().clone(),
            family.denominators().to_vec(),
            gram,
            family.power_shifts().to_vec(),
        )
        .unwrap();
        assert!(
            symmetry::verify(
                &changed_gram,
                &changed_gram,
                bank.maps[0].verified_map().momentum().clone(),
                Limits::default().verification
            )
            .is_err()
        );
    }

    #[test]
    fn positive_affine_denominators_and_nondescending_rules_are_never_used() {
        let family = paper(vec![4]);
        let mut bank = Bank::discover(&family, 7, Limits::default());
        for powers in [
            [0, 0, 1, 0, 1, 0, 1, 0, 0], // reverse sector direction
            [0, 1, 0, 0, 1, 0, 1, 0, 0], // unchanged integral
            [1, 0, 0, 0, 1, 0, 1, 0, 1], // positive affine ISP
        ] {
            assert!(bank.apply(&family, powers, &order(powers)).is_none());
        }
        assert_eq!(bank.statistics.transport_failures, 1);
        assert_eq!(bank.statistics.applied, 0);
    }

    #[test]
    fn bounded_search_and_expansion_leave_no_partial_identity() {
        let family = paper(vec![4]);
        let none = Bank::discover(
            &family,
            7,
            Limits {
                max_candidates: 0,
                ..Limits::default()
            },
        );
        assert!(none.statistics.search_limited);
        assert!(none.maps.is_empty());
        let limits = Limits {
            expansion: ExpansionLimits {
                max_endpoints: 1,
                ..Limits::default().expansion
            },
            ..Limits::default()
        };
        let mut bank = Bank::discover(&family, 7, limits);
        assert_eq!(bank.statistics.automorphisms, 1);
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, -2];
        assert!(bank.apply(&family, powers, &order(powers)).is_none());
        assert_eq!(bank.statistics.transport_failures, 1);
        bank.expansion = ExpansionLimits {
            max_total_power: 128,
            ..Limits::default().expansion
        };
        let overflow = [1, 0, 0, -64, 1, 0, 1, 0, -1];
        assert!(bank.apply(&family, overflow, &order(overflow)).is_none());
    }
}
