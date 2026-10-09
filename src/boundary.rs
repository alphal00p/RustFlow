//! Boundary providers and matching of asymptotic data to a Frobenius basis.
use crate::frobenius::FrobeniusBasis;
use crate::numeric::solve_constraints;
use crate::{ComplexFloat as C, Error, Integral, IntegralFamily, Precision, Result};
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct LeadingBoundary {
    pub exponent: Atom,
    pub coefficient: C,
}

pub trait BoundaryProvider: Send + Sync {
    /// Match the explicitly prepared deformation. Custom providers which use
    /// `constants` retain their supplied-boundary contract; automatic providers
    /// must use this mask rather than reconstructing it from evaluation options.
    #[allow(clippy::too_many_arguments)]
    fn constants_with_deformation(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        shifted: &[bool],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<C>> {
        validate_deformation_mask(family, shifted)?;
        self.constants(family, basis, epsilon, p, solutions)
    }

    fn constants(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<C>> {
        solutions.match_leading(&self.leading(family, basis, epsilon, p)?)
    }
    /// Data at z=1/eta -> 0; one entry per integral in the supplied basis.
    fn leading(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<LeadingBoundary>>;
}

/// All-hard boundary for one-loop families with unit quadratic coefficients.
/// Its coefficients are analytic massive-vacuum integrals, not fitted values.
#[derive(Clone, Copy, Debug, Default)]
pub struct OneLoopBoundary;
impl BoundaryProvider for OneLoopBoundary {
    fn constants_with_deformation(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        shifted: &[bool],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<C>> {
        validate_deformation_mask(family, shifted)?;
        if !all_physical_lines_shifted(family, shifted) {
            return Err(Error::Unsupported(
                "OneLoopBoundary requires all physical propagators to be shifted; use RecursiveBoundary for partial deformations".into(),
            ));
        }
        self.constants(family, basis, epsilon, p, solutions)
    }

    fn leading(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<LeadingBoundary>> {
        if family.loops.len() != 1 {
            return Err(Error::Unsupported(
                "OneLoopBoundary supports one-loop families; use RecursiveBoundary for multiloop families".into(),
            ));
        }
        for d in &family.propagators {
            if d.scalar_products.first() != Some(&Atom::num(1)) {
                return Err(Error::Unsupported(
                    "one-loop boundary requires unit coefficient of l^2".into(),
                ));
            }
        }
        let half_dimension = Rational::from((family.dimension, 2)) - epsilon;
        basis
            .iter()
            .map(|integral| {
                let power = integral.0[..family.physical_propagators]
                    .iter()
                    .map(|&v| v as i64)
                    .sum::<i64>();
                let rank = -integral.0[family.physical_propagators..]
                    .iter()
                    .map(|&v| v as i64)
                    .sum::<i64>();
                if power <= 0 {
                    return Ok(LeadingBoundary {
                        exponent: Atom::num(power - rank) - Atom::num((family.dimension, 2))
                            + Atom::var(family.epsilon),
                        coefficient: p.zero(),
                    });
                }
                let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re);
                let coefficient = p.div(
                    &p.mul(
                        &gamma(&half_dimension + &Rational::from(rank))?,
                        &gamma(Rational::from(power - rank) - &half_dimension)?,
                    ),
                    &p.mul(
                        &gamma(half_dimension.clone())?,
                        &gamma(Rational::from(power))?,
                    ),
                );
                let coefficient = if (power + rank) % 2 == 0 {
                    coefficient
                } else {
                    p.neg(&coefficient)
                };
                Ok(LeadingBoundary {
                    exponent: Atom::num(power - rank) - Atom::num((family.dimension, 2))
                        + Atom::var(family.epsilon),
                    coefficient,
                })
            })
            .collect()
    }
}

pub(crate) fn validate_deformation_mask(family: &IntegralFamily, shifted: &[bool]) -> Result<()> {
    if shifted.len() != family.propagators.len()
        || family.physical_propagators > shifted.len()
        || !shifted.iter().any(|&shift| shift)
        || shifted[family.physical_propagators..]
            .iter()
            .any(|&shift| shift)
    {
        return Err(Error::InvalidInput(
            "auxiliary deformation must select physical propagators with one mask slot per family propagator".into(),
        ));
    }
    Ok(())
}

pub(crate) fn all_physical_lines_shifted(family: &IntegralFamily, shifted: &[bool]) -> bool {
    shifted[..family.physical_propagators]
        .iter()
        .all(|&shift| shift)
}

impl FrobeniusBasis {
    /// Fix constants from asymptotic leading coefficients and absence of other
    /// power sectors. Only valid when the provider supplies all contributing regions.
    pub fn match_leading(&self, leading: &[LeadingBoundary]) -> Result<Vec<C>> {
        let p = self.precision;
        let n = self.columns.len();
        if leading.len() != n {
            return Err(Error::InvalidInput("boundary basis dimension".into()));
        }
        let mut equations = Vec::new();
        let mut rhs = Vec::new();
        for (j, column) in self.columns.iter().enumerate() {
            if !leading.iter().any(|b| {
                (&b.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                    .is_ok()
            }) {
                let mut row = vec![p.zero(); n];
                row[j] = p.i(1);
                equations.push(row);
                rhs.push(p.zero());
            }
        }
        for (i, boundary) in leading.iter().enumerate() {
            let mut max_preceding = 0;
            for column in &self.columns {
                if let Ok(k) = (&boundary.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                {
                    max_preceding = max_preceding.max(k);
                }
            }
            if max_preceding > self.columns[0].coefficients.len() as i64 {
                return Err(Error::Accuracy("insufficient boundary series order".into()));
            }
            for before in 0..=max_preceding.max(0) {
                let mut row = vec![p.zero(); n];
                for (j, column) in self.columns.iter().enumerate() {
                    if let Ok(k) = (&boundary.exponent - &column.exponent - Atom::num(before))
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<i64>()
                        && k >= 0
                    {
                        let coeff = column.coefficients.get(k as usize).ok_or_else(|| {
                            Error::Accuracy("boundary expansion too short".into())
                        })?;
                        row[j] = coeff[0][i].clone();
                    }
                }
                equations.push(row);
                rhs.push(if before == 0 {
                    boundary.coefficient.clone()
                } else {
                    p.zero()
                });
            }
        }
        solve_constraints(p, equations, rhs, n)
    }
}

/// One Taylor series within a dimensional region at z=1/eta=0.
/// This type certifies vanishing logarithmic terms; use [`LogRegionBoundary`]
/// when logarithmic coefficients may be nonzero or unknown.
#[derive(Clone, Debug)]
pub struct RegionBoundary {
    pub exponent: Atom,
    pub coefficients: Vec<C>,
}

/// Power/log coefficients within one symbolic dimensional region.
///
/// `coefficients[k][l]` multiplies `z^(exponent+k) log(z)^l`. `None`, a
/// missing logarithm, and a missing Taylor order are unknown coefficients,
/// never inferred zeros. Powers below `exponent` are known to vanish.
/// `max_log_power = Some(m)` additionally certifies that every logarithmic
/// power greater than `m` vanishes, including at uncomputed Taylor orders.
/// Regions absent from a component are known to vanish, so callers must supply
/// the complete list of possible exponent classes, even for unknown regions.
#[derive(Clone, Debug)]
pub struct LogRegionBoundary {
    pub exponent: Atom,
    pub coefficients: Vec<Vec<Option<C>>>,
    pub max_log_power: Option<usize>,
}

impl From<&RegionBoundary> for LogRegionBoundary {
    fn from(region: &RegionBoundary) -> Self {
        Self {
            exponent: region.exponent.clone(),
            coefficients: region
                .coefficients
                .iter()
                .map(|value| vec![Some(value.clone())])
                .collect(),
            max_log_power: Some(0),
        }
    }
}

impl FrobeniusBasis {
    /// Match known power/log coefficients without specializing symbolic
    /// epsilon-dependent exponents. Overlapping regions are summed only when
    /// every contribution to that coefficient is known. Missing rank and
    /// insufficient Frobenius depth are explicit errors.
    pub fn match_log_regions(&self, data: &[Vec<LogRegionBoundary>]) -> Result<Vec<C>> {
        let p = self.precision;
        let n = self.columns.len();
        if data.len() != n
            || self.columns.iter().any(|column| {
                column
                    .coefficients
                    .iter()
                    .flatten()
                    .any(|row| row.len() != n || row.iter().any(|value| !p.finite(value)))
            })
        {
            return Err(Error::InvalidInput(
                "region boundary dimension or nonfinite basis".into(),
            ));
        }
        for region in data.iter().flatten() {
            if region.coefficients.iter().any(|logs| {
                region
                    .max_log_power
                    .is_some_and(|max| logs.len().saturating_sub(1) > max)
                    || logs.iter().flatten().any(|value| !p.finite(value))
            }) {
                return Err(Error::InvalidInput(
                    "region coefficients exceed the declared logarithmic degree or are nonfinite"
                        .into(),
                ));
            }
        }
        let mut rows = Vec::new();
        let mut rhs = Vec::new();
        for (i, regions) in data.iter().enumerate() {
            let mut powers = std::collections::BTreeSet::new();
            // Available solution powers supply known-zero and logarithmic
            // constraints, including below an uncomputed region's leading power.
            for column in &self.columns {
                for k in 0..column.coefficients.len() {
                    powers.insert((&column.exponent + Atom::num(k as i64)).together().cancel());
                }
            }
            // Supplied known coefficients outside that range require more
            // Frobenius depth; entirely unknown data must not demand it.
            for region in regions {
                for (k, logs) in region.coefficients.iter().enumerate() {
                    if logs.iter().any(Option::is_some) {
                        powers.insert((&region.exponent + Atom::num(k as i64)).together().cancel());
                    }
                }
            }
            for exponent in powers {
                let offsets = regions
                    .iter()
                    .map(|region| integer_offset(&exponent, &region.exponent))
                    .collect::<Vec<_>>();
                let absent_class = offsets.iter().all(Option::is_none);
                let mut logarithms = vec![vec![p.zero(); n]];
                let mut complete = true;
                for (j, column) in self.columns.iter().enumerate() {
                    let Some(k) = integer_offset(&exponent, &column.exponent)
                        .and_then(|k| usize::try_from(k).ok())
                    else {
                        continue;
                    };
                    let Some(terms) = column.coefficients.get(k) else {
                        complete = false;
                        continue;
                    };
                    logarithms.resize(logarithms.len().max(terms.len()), vec![p.zero(); n]);
                    for (l, row) in terms.iter().enumerate() {
                        logarithms[l][j] = row[i].clone();
                    }
                }
                for (region, offset) in regions.iter().zip(&offsets) {
                    if let Some(logs) = offset
                        .and_then(|k| usize::try_from(k).ok())
                        .and_then(|k| region.coefficients.get(k))
                    {
                        logarithms.resize(logarithms.len().max(logs.len()), vec![p.zero(); n]);
                    }
                }
                for (l, row) in logarithms.into_iter().enumerate() {
                    let mut desired = p.zero();
                    let mut known = true;
                    for (region, offset) in regions.iter().zip(&offsets) {
                        let Some(k) = offset.and_then(|k| usize::try_from(k).ok()) else {
                            // Distinct exponent classes and powers below the
                            // region's leading power contribute exactly zero.
                            continue;
                        };
                        if region.max_log_power.is_some_and(|max| l > max) {
                            continue;
                        }
                        let Some(value) = region
                            .coefficients
                            .get(k)
                            .and_then(|logs| logs.get(l))
                            .and_then(Option::as_ref)
                        else {
                            known = false;
                            break;
                        };
                        desired = p.add(&desired, value);
                    }
                    if !known || (!complete && absent_class) {
                        continue;
                    }
                    if !complete {
                        return Err(Error::Accuracy(
                            "insufficient Frobenius boundary order for a known power/log coefficient"
                                .into(),
                        ));
                    }
                    rows.push(row);
                    rhs.push(desired);
                }
            }
        }
        for (j, column) in self.columns.iter().enumerate() {
            if !data
                .iter()
                .flatten()
                .any(|region| integer_offset(&region.exponent, &column.exponent).is_some())
            {
                let mut row = vec![p.zero(); n];
                row[j] = p.i(1);
                rows.push(row);
                rhs.push(p.zero());
            }
        }
        solve_constraints(p, rows, rhs, n)
    }
}

fn integer_offset(a: &Atom, b: &Atom) -> Option<i64> {
    (a - b).together().cancel().to_string().parse().ok()
}

#[cfg(test)]
mod logarithmic_tests {
    use super::*;
    use crate::frobenius::FrobeniusColumn;

    fn jordan(p: Precision) -> FrobeniusBasis {
        FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: parse!("boundary_log_eps"),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: parse!("boundary_log_eps"),
                    coefficients: vec![vec![vec![p.zero(), p.i(1)], vec![p.i(1), p.zero()]]],
                },
            ],
            precision: p,
        }
    }

    #[test]
    fn known_nonzero_logarithm_fixes_a_jordan_solution() {
        let p = Precision::decimal(60).unwrap();
        let basis = jordan(p);
        let constants = basis
            .match_log_regions(&[
                vec![LogRegionBoundary {
                    exponent: parse!("boundary_log_eps"),
                    coefficients: vec![vec![Some(p.i(7)), Some(p.i(3))]],
                    max_log_power: Some(1),
                }],
                vec![LogRegionBoundary {
                    exponent: parse!("boundary_log_eps"),
                    coefficients: vec![],
                    max_log_power: None,
                }],
            ])
            .unwrap();
        assert!(p.close(&constants[0], &p.i(7), 40));
        assert!(p.close(&constants[1], &p.i(3), 40));
    }

    #[test]
    fn unknown_constant_is_not_zero_when_logarithm_is_known() {
        let p = Precision::decimal(60).unwrap();
        let basis = jordan(p);
        let result = basis.match_log_regions(&[
            vec![LogRegionBoundary {
                exponent: parse!("boundary_log_eps"),
                coefficients: vec![vec![None, Some(p.i(3))]],
                max_log_power: Some(1),
            }],
            vec![LogRegionBoundary {
                exponent: parse!("boundary_log_eps"),
                coefficients: vec![vec![Some(p.i(3))]],
                max_log_power: Some(0),
            }],
        ]);
        assert!(matches!(result, Err(Error::IncompleteReduction(_))));
    }

    #[test]
    fn overlapping_uncomputed_region_cannot_supply_a_zero() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![FrobeniusColumn {
                exponent: parse!("boundary_overlap_eps+1"),
                coefficients: vec![vec![vec![p.i(1)]]],
            }],
            precision: p,
        };
        let mut data = vec![vec![
            LogRegionBoundary {
                exponent: parse!("boundary_overlap_eps"),
                coefficients: vec![vec![Some(p.zero())]],
                max_log_power: Some(0),
            },
            LogRegionBoundary {
                exponent: parse!("boundary_overlap_eps+1"),
                coefficients: vec![vec![Some(p.i(3))]],
                max_log_power: Some(0),
            },
        ]];
        assert!(matches!(
            basis.match_log_regions(&data),
            Err(Error::IncompleteReduction(_))
        ));
        data[0][0].coefficients.push(vec![Some(p.i(2))]);
        assert!(p.close(&basis.match_log_regions(&data).unwrap()[0], &p.i(5), 40));
    }

    #[test]
    fn known_data_require_depth_but_unknown_orders_do_not() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![FrobeniusColumn {
                exponent: Atom::new(),
                coefficients: vec![vec![vec![p.i(1)]]],
            }],
            precision: p,
        };
        let mut data = vec![vec![LogRegionBoundary {
            exponent: Atom::new(),
            coefficients: vec![vec![Some(p.i(7))], vec![None]],
            max_log_power: None,
        }]];
        assert!(p.close(&basis.match_log_regions(&data).unwrap()[0], &p.i(7), 40));
        data[0][0].coefficients[1][0] = Some(p.i(8));
        assert!(matches!(
            basis.match_log_regions(&data),
            Err(Error::Accuracy(_))
        ));
    }

    #[test]
    fn dimensional_classes_remain_distinct() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: parse!("boundary_classes_eps"),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: parse!("2*boundary_classes_eps"),
                    coefficients: vec![vec![vec![p.i(1), p.i(1)]]],
                },
            ],
            precision: p,
        };
        let region = |exponent, value| LogRegionBoundary {
            exponent,
            coefficients: vec![vec![Some(value)]],
            max_log_power: Some(0),
        };
        let data = vec![
            vec![
                region(parse!("boundary_classes_eps"), p.i(2)),
                region(parse!("2*boundary_classes_eps"), p.i(3)),
            ],
            vec![region(parse!("2*boundary_classes_eps"), p.i(3))],
        ];
        let constants = basis.match_log_regions(&data).unwrap();
        assert!(p.close(&constants[0], &p.i(2), 40));
        assert!(p.close(&constants[1], &p.i(3), 40));
    }

    #[test]
    fn missing_logarithm_and_declared_zero_logarithm_differ() {
        let p = Precision::decimal(60).unwrap();
        let basis = jordan(p);
        let mut data = vec![
            vec![LogRegionBoundary {
                exponent: parse!("boundary_log_eps"),
                coefficients: vec![vec![Some(p.i(7))]],
                max_log_power: None,
            }],
            vec![LogRegionBoundary {
                exponent: parse!("boundary_log_eps"),
                coefficients: vec![],
                max_log_power: None,
            }],
        ];
        assert!(matches!(
            basis.match_log_regions(&data),
            Err(Error::IncompleteReduction(_))
        ));
        data[0][0].max_log_power = Some(0);
        let constants = basis.match_log_regions(&data).unwrap();
        assert!(p.close(&constants[0], &p.i(7), 40));
        assert!(p.close(&constants[1], &p.zero(), 40));
    }
}

impl FrobeniusBasis {
    /// Match every supplied power and require logarithmic coefficients to vanish.
    /// Distinct epsilon-dependent exponent classes are kept separate.
    /// For nonzero or unknown logarithms use [`Self::match_log_regions`].
    pub fn match_regions(&self, data: &[Vec<RegionBoundary>]) -> Result<Vec<C>> {
        let p = self.precision;
        let n = self.columns.len();
        if data.len() != n {
            return Err(Error::InvalidInput("region boundary dimension".into()));
        }
        let mut rows = Vec::new();
        let mut rhs = Vec::new();
        for (i, regions) in data.iter().enumerate() {
            let mut powers = std::collections::BTreeSet::new();
            for region in regions {
                for k in 0..region.coefficients.len() {
                    powers.insert((&region.exponent + Atom::num(k as i64)).together().cancel());
                }
            }
            // A region starts at its stated leading power. Lower powers in
            // the same indicial class are known zeros, not missing data.
            for column in &self.columns {
                let end = regions
                    .iter()
                    .filter_map(|r| {
                        let end = (&r.exponent + Atom::num(r.coefficients.len() as i64 - 1)
                            - &column.exponent)
                            .together()
                            .cancel()
                            .to_string()
                            .parse::<usize>()
                            .ok()?;
                        // An uncomputed region states its leading power only.
                        // Use available lower known-zero coefficients without
                        // demanding a longer expansion for an unused component.
                        Some(if r.coefficients.is_empty() {
                            end.min(column.coefficients.len().saturating_sub(1))
                        } else {
                            end
                        })
                    })
                    .max();
                if let Some(end) = end {
                    if end >= column.coefficients.len() {
                        return Err(Error::Accuracy("insufficient boundary series".into()));
                    }
                    for k in 0..=end {
                        powers.insert((&column.exponent + Atom::num(k as i64)).together().cancel());
                    }
                }
            }
            // A class absent in this component is identically zero even when
            // it is allowed in another component. These constraints can fix
            // constants without evaluating any region integral here.
            for column in &self.columns {
                if !regions.iter().any(|r| {
                    (&r.exponent - &column.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<i64>()
                        .is_ok()
                }) {
                    for (k, logs) in column.coefficients.iter().enumerate() {
                        if logs.iter().any(|row| row[i] != p.zero()) {
                            powers.insert(
                                (&column.exponent + Atom::num(k as i64)).together().cancel(),
                            );
                        }
                    }
                }
            }
            'powers: for exponent in powers {
                let absent_class = !regions.iter().any(|r| {
                    (&exponent - &r.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<i64>()
                        .is_ok()
                });
                // Do not mistake an uncomputed subleading term of another
                // region in the same class for a known zero.
                if regions.iter().any(|r| {
                    (&exponent - &r.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                        .is_ok_and(|k| k >= r.coefficients.len())
                }) {
                    continue;
                }
                let mut desired = p.zero();
                for region in regions {
                    if let Ok(k) = (&exponent - &region.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                        && let Some(value) = region.coefficients.get(k)
                    {
                        desired = p.add(&desired, value);
                    }
                }
                let mut logarithms = vec![vec![p.zero(); n]];
                for (j, column) in self.columns.iter().enumerate() {
                    if let Ok(k) = (&exponent - &column.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                    {
                        let Some(terms) = column.coefficients.get(k) else {
                            if absent_class {
                                continue 'powers;
                            }
                            return Err(Error::Accuracy(
                                "insufficient Frobenius boundary order".into(),
                            ));
                        };
                        logarithms.resize(logarithms.len().max(terms.len()), vec![p.zero(); n]);
                        for (l, row) in terms.iter().enumerate() {
                            logarithms[l][j] = row[i].clone();
                        }
                    }
                }
                for (l, row) in logarithms.into_iter().enumerate() {
                    rows.push(row);
                    rhs.push(if l == 0 { desired.clone() } else { p.zero() });
                }
            }
        }
        for (j, column) in self.columns.iter().enumerate() {
            if !data.iter().flatten().any(|r| {
                (&r.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                    .is_ok()
            }) {
                let mut row = vec![p.zero(); n];
                row[j] = p.i(1);
                rows.push(row);
                rhs.push(p.zero());
            }
        }
        solve_constraints(p, rows, rhs, n)
    }
}

// Greedy row-rank analysis of the Frobenius solution determines which region
// coefficients are actually needed. Most components need no boundary integral:
// the differential equations already determine them from a smaller set.
impl FrobeniusBasis {
    pub(crate) fn region_orders(
        &self,
        powers: &[Vec<Atom>],
        max_half_order: usize,
    ) -> Result<Vec<Vec<Option<usize>>>> {
        use std::collections::BTreeMap;
        let p = self.precision;
        let n = self.columns.len();
        if powers.len() != n
            || self.columns.iter().any(|column| {
                column
                    .coefficients
                    .iter()
                    .flatten()
                    .any(|row| row.len() != n)
            })
        {
            return Err(Error::InvalidInput("region planning dimension".into()));
        }
        let zero = p.zero();
        let mut orders: Vec<Vec<Option<usize>>> = powers
            .iter()
            .map(|regions| vec![None; regions.len()])
            .collect();
        let mut independent = IndependentBoundaryRows::new(p);
        // Classes absent from every region have zero integration constants.
        for (j, column) in self.columns.iter().enumerate() {
            if !powers
                .iter()
                .flatten()
                .any(|power| half_offset(power, &column.exponent).is_some())
            {
                let mut row = vec![zero.clone(); n];
                row[j] = p.i(1);
                independent.insert(row);
            }
        }
        // A candidate is a coefficient in one component of the full solution.
        // Store only its identity; dense rows are built lazily after sorting.
        let mut candidates = BTreeMap::<(usize, Atom), ()>::new();
        for column in &self.columns {
            for (k, logs) in column.coefficients.iter().enumerate() {
                let exponent = (&column.exponent + Atom::num(k as i64)).together().cancel();
                for i in 0..n {
                    if logs.iter().any(|row| row[i] != zero) {
                        candidates.insert((i, exponent.clone()), ());
                    }
                }
            }
        }
        let mut candidates = candidates
            .into_keys()
            .filter_map(|(i, exponent)| {
                let required = powers[i]
                    .iter()
                    .map(|power| {
                        half_offset(&exponent, power).and_then(|k| usize::try_from(k).ok())
                    })
                    .collect::<Vec<_>>();
                let cost = required.iter().flatten().copied().max().unwrap_or(0);
                (cost <= max_half_order).then_some((cost, i, exponent, required))
            })
            .collect::<Vec<_>>();
        candidates.sort_by(|a, b| (a.0, a.1, &a.2).cmp(&(b.0, b.1, &b.2)));
        for (_, i, exponent, required) in candidates {
            if independent.rows.len() == n {
                return Ok(orders);
            }
            let mut rows = vec![vec![zero.clone(); n]];
            let mut complete = true;
            for (j, column) in self.columns.iter().enumerate() {
                let Some(k) = (&exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<usize>()
                    .ok()
                else {
                    continue;
                };
                let Some(logs) = column.coefficients.get(k) else {
                    complete = false;
                    break;
                };
                rows.resize(rows.len().max(logs.len()), vec![zero.clone(); n]);
                for (l, values) in logs.iter().enumerate() {
                    rows[l][j] = values[i].clone();
                }
            }
            if !complete {
                continue;
            }
            let mut useful = false;
            for row in rows {
                useful |= independent.insert(row);
            }
            if useful {
                for (order, needed) in orders[i].iter_mut().zip(required) {
                    if let Some(needed) = needed {
                        *order = Some(order.unwrap_or(0).max(needed));
                    }
                }
            }
        }
        if independent.rows.len() == n {
            Ok(orders)
        } else {
            Err(Error::IncompleteReduction(format!(
                "boundary coefficients through half-order {max_half_order} fix {} of {n} constants; increase the Frobenius order or boundary limit",
                independent.rows.len()
            )))
        }
    }
}

fn half_offset(a: &Atom, b: &Atom) -> Option<i64> {
    (Atom::num(2) * (a - b))
        .together()
        .cancel()
        .to_string()
        .parse()
        .ok()
}

struct IndependentBoundaryRows {
    p: Precision,
    rows: Vec<(usize, Vec<C>)>,
}
impl IndependentBoundaryRows {
    fn new(p: Precision) -> Self {
        Self {
            p,
            rows: Vec::new(),
        }
    }
    fn insert(&mut self, mut row: Vec<C>) -> bool {
        let p = self.p;
        for (pivot, previous) in &self.rows {
            let factor = row[*pivot].clone();
            for j in *pivot..row.len() {
                row[j] = p.sub(&row[j], &p.mul(&factor, &previous[j]));
            }
        }
        let Some(pivot) = row.iter().position(|v| p.norm(v) > p.tolerance(p.bits / 5)) else {
            return false;
        };
        row[..pivot].fill(p.zero());
        let divisor = row[pivot].clone();
        for value in &mut row[pivot..] {
            *value = p.div(value, &divisor);
        }
        // Insertion can discover a pivot left of an older pivot. Keep echelon
        // order so later rows eliminate that new pivot first.
        let index = self.rows.partition_point(|(column, _)| *column < pivot);
        self.rows.insert(index, (pivot, row));
        true
    }
}

#[cfg(test)]
mod planning_tests {
    use super::*;
    use crate::frobenius::FrobeniusColumn;

    #[test]
    fn avoids_high_order_components_already_fixed_by_the_system() {
        let p = Precision::decimal(60).unwrap();
        // Fundamental matrix [[1,0],[x^12,1]]. Its two constants are
        // determined at order zero; expanding the second component to x^12
        // would ask for a redundant, potentially expensive boundary integral.
        let mut first = vec![vec![vec![p.zero(); 2]]; 13];
        first[0][0][0] = p.i(1);
        first[12][0][1] = p.i(1);
        let mut second = vec![vec![vec![p.zero(); 2]]; 13];
        second[0][0][1] = p.i(1);
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: first,
                },
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: second,
                },
            ],
            precision: p,
        };
        let powers = vec![vec![Atom::new()]; 2];
        assert_eq!(
            basis.region_orders(&powers, 0).unwrap(),
            vec![vec![Some(0)]; 2]
        );
        let values = basis
            .match_regions(&[
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![p.i(2)],
                }],
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![p.i(3)],
                }],
            ])
            .unwrap();
        assert!(p.close(&values[0], &p.i(2), 40));
        assert!(p.close(&values[1], &p.i(3), 40));
    }

    #[test]
    fn includes_every_overlapping_region_at_the_selected_power() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![FrobeniusColumn {
                exponent: Atom::num(1),
                coefficients: vec![vec![vec![p.i(1)]]],
            }],
            precision: p,
        };
        let powers = vec![vec![Atom::new(), Atom::num(1)]];
        assert!(matches!(
            basis.region_orders(&powers, 1),
            Err(Error::IncompleteReduction(_))
        ));
        assert_eq!(
            basis.region_orders(&powers, 2).unwrap(),
            vec![vec![Some(2), Some(0)]]
        );
        let values = basis
            .match_regions(&[vec![
                RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![p.zero(), p.i(2)],
                },
                RegionBoundary {
                    exponent: Atom::num(1),
                    coefficients: vec![p.i(3)],
                },
            ]])
            .unwrap();
        assert!(p.close(&values[0], &p.i(5), 40));
        // Leaving the first region uncomputed must not pretend its x^1 term
        // vanishes and return the second region's coefficient alone.
        assert!(matches!(
            basis.match_regions(&[vec![
                RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![]
                },
                RegionBoundary {
                    exponent: Atom::num(1),
                    coefficients: vec![p.i(3)]
                },
            ]]),
            Err(Error::IncompleteReduction(_))
        ));
    }

    #[test]
    fn uses_log_constraints_without_evaluating_redundant_components() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: vec![vec![vec![p.zero(), p.i(1)], vec![p.i(1), p.zero()]]],
                },
            ],
            precision: p,
        };
        let powers = vec![vec![Atom::new()]; 2];
        assert_eq!(
            basis.region_orders(&powers, 0).unwrap(),
            vec![vec![Some(0)], vec![None]]
        );
        let values = basis
            .match_regions(&[
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![p.i(7)],
                }],
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![],
                }],
            ])
            .unwrap();
        assert!(p.close(&values[0], &p.i(7), 40));
        assert!(p.close(&values[1], &p.zero(), 40));
    }

    #[test]
    fn per_component_absent_classes_are_known_zero_constraints() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: parse!("component_eps"),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: parse!("component_eps"),
                    coefficients: vec![vec![vec![p.i(1), p.i(1)]]],
                },
            ],
            precision: p,
        };
        let powers = vec![vec![Atom::new()], vec![parse!("component_eps")]];
        assert_eq!(
            basis.region_orders(&powers, 0).unwrap(),
            vec![vec![None], vec![Some(0)]]
        );
        let values = basis
            .match_regions(&[
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![],
                }],
                vec![RegionBoundary {
                    exponent: parse!("component_eps"),
                    coefficients: vec![p.i(5)],
                }],
            ])
            .unwrap();
        assert!(p.close(&values[0], &p.i(-5), 40));
        assert!(p.close(&values[1], &p.i(5), 40));
    }

    #[test]
    fn unused_late_region_does_not_require_a_longer_frobenius_series() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: Atom::new(),
                    coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
                },
            ],
            precision: p,
        };
        let powers = vec![vec![Atom::new()], vec![Atom::num(100)]];
        assert_eq!(
            basis.region_orders(&powers, 0).unwrap(),
            vec![vec![Some(0)], vec![None]]
        );
        let values = basis
            .match_regions(&[
                vec![RegionBoundary {
                    exponent: Atom::new(),
                    coefficients: vec![p.i(5)],
                }],
                vec![RegionBoundary {
                    exponent: Atom::num(100),
                    coefficients: vec![],
                }],
            ])
            .unwrap();
        assert!(p.close(&values[0], &p.i(5), 40));
        assert!(p.close(&values[1], &p.zero(), 40));
    }

    #[test]
    fn retains_half_integer_and_absent_dimensional_classes() {
        let p = Precision::decimal(60).unwrap();
        let basis = FrobeniusBasis {
            columns: vec![
                FrobeniusColumn {
                    exponent: Atom::num((1, 2)),
                    coefficients: vec![vec![vec![p.i(1), p.zero()]]],
                },
                FrobeniusColumn {
                    exponent: parse!("planning_eps"),
                    coefficients: vec![vec![vec![p.zero(), p.i(1)]]],
                },
            ],
            precision: p,
        };
        let powers = vec![vec![Atom::new()]; 2];
        assert_eq!(
            basis.region_orders(&powers, 1).unwrap(),
            vec![vec![Some(1)], vec![None]]
        );
        let values = basis
            .match_regions(&[
                vec![
                    RegionBoundary {
                        exponent: Atom::new(),
                        coefficients: vec![p.zero()],
                    },
                    RegionBoundary {
                        exponent: Atom::num((1, 2)),
                        coefficients: vec![p.i(4)],
                    },
                ],
                vec![
                    RegionBoundary {
                        exponent: Atom::new(),
                        coefficients: vec![],
                    },
                    RegionBoundary {
                        exponent: Atom::num((1, 2)),
                        coefficients: vec![],
                    },
                ],
            ])
            .unwrap();
        assert!(p.close(&values[0], &p.i(4), 40));
        assert!(p.close(&values[1], &p.zero(), 40));
    }
}
