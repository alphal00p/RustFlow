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
use symbolica::utils::resources;

fn verification_failure(error: symmetry::Error) -> crate::Result<()> {
    if let Some(resource) = crate::error::preserve_resource(&error) {
        return Err(resource);
    }
    match error {
        symmetry::Error::SingularLoopMap
        | symmetry::Error::SingularExternalMap
        | symmetry::Error::ExternalGramMismatch { .. }
        | symmetry::Error::DenominatorReplayMismatch { .. }
        | symmetry::Error::ResourceLimit { .. } => Ok(()),
        error => Err(super::native_error(error)),
    }
}

fn transport_failure(error: integral_transport::Error) -> crate::Result<()> {
    if let Some(resource) = crate::error::preserve_resource(&error) {
        return Err(resource);
    }
    use integral_transport::{Error as Transport, ExpansionError as Expansion};
    use rustred::algebra::ExactAlgebraError as Algebra;
    match error {
        Transport::NonUnitJacobian
        | Transport::AnalyticPowerShift { .. }
        | Transport::UnresolvedCondition { .. }
        | Transport::NonconstantCoefficient { .. }
        | Transport::NonUnitActiveRow { .. }
        | Transport::ActiveBijection
        | Transport::PositiveOutsideRoot { .. }
        | Transport::Expansion(
            Expansion::ResourceLimit { .. }
            | Expansion::NativeExponentLimit { .. }
            | Expansion::NativePolynomialSupportExceeded { .. }
            | Expansion::PowerShiftUnderflow { .. }
            | Expansion::NativeExponentDegreeOverflow { .. }
            | Expansion::ExactAlgebra(
                Algebra::ResourceLimit { .. }
                | Algebra::ExponentLimit { .. }
                | Algebra::ExponentArithmeticOverflow { .. },
            ),
        ) => Ok(()),
        error => Err(super::native_error(error)),
    }
}

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
    pub fn discover(
        family: &IntegralFamily,
        physical: usize,
        limits: Limits,
    ) -> crate::Result<Self> {
        resources::checkpoint()?;
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
            return Ok(bank);
        }
        if family.loop_count() == 0 || family.external_count() == 0 {
            return Ok(bank);
        }
        let Some(candidates) = candidates(family, physical, limits)? else {
            return Ok(bank);
        };
        bank.statistics.search_limited = candidates.limited;
        // IntegralFamily intentionally is not Clone. Rebuild through its public
        // validator while retaining the exact coefficient context and ordering.
        let target = IntegralFamily::new(
            family.name().to_owned(),
            family.loop_momenta().to_vec(),
            family.external_momenta().to_vec(),
            family.coefficient_context().clone(),
            family.dimension().clone(),
            family.denominators().to_vec(),
            family.external_gram().to_vec(),
            family.power_shifts().to_vec(),
        )
        .map_err(super::native_error)?;
        let target = Arc::new(target);
        let root =
            Mask::try_new((0..bank.arity).map(|i| i < physical)).map_err(super::native_error)?;
        for candidate in candidates.maps {
            resources::checkpoint()?;
            bank.statistics.candidates += 1;
            let map = match symmetry::verify(family, &target, candidate, limits.verification) {
                Ok(map) => map,
                Err(error) => {
                    verification_failure(error)?;
                    continue;
                }
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
            match integral_transport::compile(
                family,
                Arc::clone(&target),
                Arc::new(map),
                root.clone(),
                root.clone(),
                limits.expansion,
            ) {
                Ok(prepared) => bank.maps.push(prepared),
                Err(error) => transport_failure(error)?,
            }
        }
        bank.statistics.automorphisms = bank.maps.len();
        Ok(bank)
    }

    /// Return a complete exact identity or no rule. Partial transport results,
    /// overflowing compact keys, and nondecreasing terms are never installed.
    pub fn apply<const N: usize>(
        &mut self,
        family: &IntegralFamily,
        powers: [i16; N],
        order: &IntegralOrder<N>,
    ) -> crate::Result<Option<Vec<Term<N, Coefficient>>>> {
        let result = self.apply_inner(family, powers, order)?;
        if result.is_some() {
            self.statistics.applied += 1;
        } else {
            self.statistics.missed += 1;
        }
        Ok(result)
    }

    fn apply_inner<const N: usize>(
        &mut self,
        family: &IntegralFamily,
        powers: [i16; N],
        order: &IntegralOrder<N>,
    ) -> crate::Result<Option<Vec<Term<N, Coefficient>>>> {
        resources::checkpoint()?;
        if self.arity > N
            || powers[self.arity..].iter().any(|&value| value != 0)
            || self.fingerprint != family.fingerprint()
        {
            return Ok(None);
        }
        // Compact-key admission is representation-only (PowerError), not a
        // fallible resource operation. Unsupported powers retain the fallback.
        let Ok(lhs) = Integral::numeric(powers) else {
            return Ok(None);
        };
        let key = IntegralKey::try_new(powers[..self.arity].iter().copied().map(i64::from))
            .map_err(super::native_error)?;
        for map in &self.maps {
            resources::checkpoint()?;
            let transported = match map.transport(&key, self.expansion) {
                Ok(value) => value,
                Err(error) => {
                    transport_failure(error)?;
                    self.statistics.transport_failures += 1;
                    continue;
                }
            };
            let mut terms = Vec::new();
            let mut admitted = true;
            for term in transported.terms() {
                let Some(mut values) = term
                    .key()
                    .powers()
                    .iter()
                    .copied()
                    .map(i16::try_from)
                    .collect::<Result<Vec<_>, _>>()
                    .ok()
                else {
                    admitted = false;
                    break;
                };
                values.resize(N, 0);
                let Ok(values) = values.try_into() else {
                    admitted = false;
                    break;
                };
                let Ok(integral) = Integral::numeric(values) else {
                    admitted = false;
                    break;
                };
                if order.compare(&lhs, &integral) != Ordering::Less {
                    admitted = false;
                    break;
                }
                terms.push(Term {
                    integral,
                    coefficient: term.coefficient().clone(),
                });
            }
            if admitted {
                return Ok(Some(terms));
            }
        }
        Ok(None)
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
) -> crate::Result<CoefficientMatrix> {
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
    .map_err(super::native_error)
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

fn candidates(
    family: &IntegralFamily,
    physical: usize,
    limits: Limits,
) -> crate::Result<Option<Candidates>> {
    resources::checkpoint()?;
    let mut result = Candidates {
        maps: Vec::new(),
        limited: false,
    };
    if limits.max_candidates == 0 || limits.max_external_maps == 0 {
        result.limited = true;
        return Ok(Some(result));
    }
    let loops = family.loop_count();
    let external = family.external_count();
    let Some(groups) = branches(family, physical) else {
        return Ok(None);
    };
    let Some(gram) = family
        .external_gram()
        .iter()
        .map(|row| row.iter().map(constant).collect::<Option<Vec<_>>>())
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let mut external_maps = vec![Matrix::identity(external as u32, Q)];
    for group in &groups {
        for (i, left) in group.displacements.iter().enumerate() {
            for right in group.displacements.iter().skip(i + 1) {
                resources::checkpoint()?;
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
                .map_err(crate::Error::Reduction)?;
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
        resources::checkpoint()?;
        let rows = anchors
            .iter()
            .map(|b| b.direction.clone())
            .chain(std::iter::once(group.direction.clone()))
            .collect::<Vec<_>>();
        if Matrix::from_nested_vec(rows, Q)
            .map_err(crate::Error::Reduction)?
            .rank()
            > anchors.len()
        {
            anchors.push(group);
        }
        if anchors.len() == loops {
            break;
        }
    }
    if anchors.len() != loops {
        return Ok(None);
    }
    let basis = Matrix::from_nested_vec(anchors.iter().map(|b| b.direction.clone()).collect(), Q)
        .map_err(crate::Error::Reduction)?;
    let inverse = basis.inv().map_err(super::native_error)?;
    let source = Matrix::from_nested_vec(
        anchors.iter().map(|b| b.displacements[0].clone()).collect(),
        Q,
    )
    .map_err(crate::Error::Reduction)?;
    let context = family.coefficient_context();
    for external_map in external_maps {
        let mapped = &source * &external_map;
        let mut destinations = vec![0; loops];
        loop {
            resources::checkpoint()?;
            if result.maps.len() == limits.max_candidates {
                result.limited = true;
                return Ok(Some(result));
            }
            let target = Matrix::from_nested_vec(
                anchors
                    .iter()
                    .zip(&destinations)
                    .map(|(b, &i)| b.displacements[i].clone())
                    .collect(),
                Q,
            )
            .map_err(crate::Error::Reduction)?;
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
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn shared_resource_causes_survive_discovery_and_transport_failures() {
        let expected = crate::ResourceLimitError::Memory {
            rss_bytes: 2,
            additional_bytes: 3,
            limit_bytes: 4,
        };
        let verification = symmetry::Error::ExactAlgebra(
            rustred::algebra::ExactAlgebraError::Resource(expected.clone()),
        );
        assert!(matches!(
            verification_failure(verification),
            Err(crate::Error::Resource(error)) if error == expected
        ));
        let transport =
            integral_transport::Error::Expansion(integral_transport::ExpansionError::ExactAlgebra(
                rustred::algebra::ExactAlgebraError::Resource(expected.clone()),
            ));
        assert!(matches!(
            transport_failure(transport),
            Err(crate::Error::Resource(error)) if error == expected
        ));
    }

    #[test]
    fn invalid_exact_algebra_is_not_an_optional_accelerator_miss() {
        assert!(matches!(
            verification_failure(symmetry::Error::ExactAlgebra(
                rustred::algebra::ExactAlgebraError::ZeroDenominator,
            )),
            Err(crate::Error::Reduction(_))
        ));
        assert!(matches!(
            transport_failure(integral_transport::Error::Expansion(
                integral_transport::ExpansionError::ExactAlgebra(
                    rustred::algebra::ExactAlgebraError::ZeroDenominator,
                ),
            )),
            Err(crate::Error::Reduction(_))
        ));
    }

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
        let mut bank = Bank::discover(&family, 2, Limits::default()).unwrap();
        assert_eq!(bank.statistics.automorphisms, 1);
        let terms = bank
            .apply(&family, [2, 0], &order([2, 0]))
            .unwrap()
            .unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].integral, Integral::numeric([0, 2]).unwrap());
        assert_eq!(terms[0].coefficient, context.one());
        assert!(
            bank.apply(&family, [0, 2], &order([0, 2]))
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn discovers_routing_and_transports_the_exact_numerator_square() {
        let family = paper(vec![4]);
        let mut bank = Bank::discover(&family, 7, Limits::default()).unwrap();
        assert_eq!(bank.statistics.candidates, 36);
        assert_eq!(bank.statistics.automorphisms, 1);
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, -2];
        let terms = bank
            .apply(&family, powers, &order(powers))
            .unwrap()
            .unwrap();
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
        let bank = Bank::discover(&family, 7, Limits::default()).unwrap();
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
        let mut bank = Bank::discover(&family, 7, Limits::default()).unwrap();
        let changed_eta = paper(vec![0]);
        assert_eq!(
            Bank::discover(&changed_eta, 7, Limits::default())
                .unwrap()
                .statistics
                .automorphisms,
            0
        );
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, 0];
        assert!(
            bank.apply(&changed_eta, powers, &order(powers))
                .unwrap()
                .is_none()
        );
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
        let mut bank = Bank::discover(&family, 7, Limits::default()).unwrap();
        for powers in [
            [0, 0, 1, 0, 1, 0, 1, 0, 0], // reverse sector direction
            [0, 1, 0, 0, 1, 0, 1, 0, 0], // unchanged integral
            [1, 0, 0, 0, 1, 0, 1, 0, 1], // positive affine ISP
        ] {
            assert!(
                bank.apply(&family, powers, &order(powers))
                    .unwrap()
                    .is_none()
            );
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
        )
        .unwrap();
        assert!(none.statistics.search_limited);
        assert!(none.maps.is_empty());
        let limits = Limits {
            expansion: ExpansionLimits {
                diagnostics_enabled: true,
                max_endpoints: 1,
                ..Limits::default().expansion
            },
            ..Limits::default()
        };
        let mut bank = Bank::discover(&family, 7, limits).unwrap();
        assert_eq!(bank.statistics.automorphisms, 1);
        let powers = [1, 0, 0, 0, 1, 0, 1, 0, -2];
        assert!(
            bank.apply(&family, powers, &order(powers))
                .unwrap()
                .is_none()
        );
        assert_eq!(bank.statistics.transport_failures, 1);
        bank.expansion = ExpansionLimits {
            diagnostics_enabled: true,
            max_total_power: 128,
            ..Limits::default().expansion
        };
        let overflow = [1, 0, 0, -64, 1, 0, 1, 0, -1];
        assert!(
            bank.apply(&family, overflow, &order(overflow))
                .unwrap()
                .is_none()
        );
    }
}
