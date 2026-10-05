//! Exact target reductions applied to cached Laurent master coefficients.
use crate::transport_cache::{CachedBoundary, CachedPoint, EpsilonRange};
use crate::{ComplexFloat, Error, Precision, PreparedPhysicalFamily, Result};
use std::collections::BTreeMap;
use symbolica::poly::series::{Series, SeriesDepth};
use symbolica::prelude::*;

/// A target value projected from an independently verified master boundary.
/// Errors include inherited master uncertainty and precision refinement of the
/// projection. No new epsilon fit or independent master calculation is implied.
#[derive(Clone, Debug)]
pub struct ProjectedLaurentExpansion {
    pub coefficients: BTreeMap<i32, ComplexFloat>,
    pub absolute_errors: BTreeMap<i32, Float>,
    pub verified_digits: u32,
    pub working_bits: u32,
}

type ExactSeries = Series<symbolica::domains::atom::AtomField>;

fn expansion(coefficient: &Atom, epsilon: Symbol, depth: SeriesDepth) -> Result<ExactSeries> {
    crate::family::scalar_symbols(coefficient.as_view(), &mut Default::default())?;
    let coefficient = crate::family::encode_complex(coefficient);
    let _: RationalPolynomial<IntegerRing, u16> = coefficient
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::Unsupported(format!("rational target coefficient required: {e}")))?;
    coefficient
        .series(epsilon, 0, depth)
        .map_err(|e| Error::Unsupported(format!("target coefficient series: {e}")))
}

fn valuation(coefficient: &Atom, epsilon: Symbol) -> Result<Option<i32>> {
    if coefficient.is_zero() {
        return Ok(None);
    }
    let series = expansion(coefficient, epsilon, SeriesDepth::relative(1))?;
    if series.is_zero() {
        return Err(Error::Numerical(
            "nonzero target weight has no leading series coefficient".into(),
        ));
    }
    Ok(Some(
        series
            .get_trailing_exponent()
            .to_string()
            .parse()
            .map_err(|_| Error::Limit("target weight epsilon valuation exceeds i32".into()))?,
    ))
}

impl PreparedPhysicalFamily {
    pub(crate) fn ordinary_master_leading(&self) -> Result<i32> {
        let family = self.family();
        let mut symbols = std::collections::BTreeSet::new();
        for coefficient in family.external_gram.iter().flatten().chain(
            family
                .propagators
                .iter()
                .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
        ) {
            crate::family::scalar_symbols(coefficient.as_view(), &mut symbols)?;
        }
        if symbols.contains(&Atom::var(family.epsilon)) {
            return Err(Error::Unsupported("automatic master pole bounds require epsilon-independent propagators and kinematics".into()));
        }
        i32::try_from(family.loops.len())
            .ok()
            .and_then(|n| n.checked_mul(-2))
            .ok_or_else(|| Error::Limit("master Laurent bound overflow".into()))
    }

    /// Determine how far masters must be known before reducing targets. A pole
    /// in a reduction coefficient requires additional positive epsilon orders;
    /// a cache truncated at the target's last order is not sufficient in general.
    pub fn required_master_range(
        &self,
        point: &BTreeMap<Symbol, Atom>,
        targets: EpsilonRange,
    ) -> Result<EpsilonRange> {
        EpsilonRange::new(targets.leading, targets.last)?;
        self.flow()
            .identity()
            .validate_point(&CachedPoint::Exact(point.clone()))?;
        let substitutions = point
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect();
        let leading = self.ordinary_master_leading()?;
        let mut last = leading;
        let columns = self
            .basis()
            .iter()
            .enumerate()
            .map(|(column, integral)| (integral, column))
            .collect::<BTreeMap<_, _>>();
        for (integral, weight) in self.target_reductions().iter().flat_map(|row| row.iter()) {
            let column = *columns.get(integral).ok_or_else(|| {
                Error::IncompleteReduction(
                    "target reduction leaves the physical master basis".into(),
                )
            })?;
            let weight = self.cache_target_weight(column, weight);
            let weight = crate::family::substitute(&weight, &substitutions)
                .together()
                .cancel();
            if let Some(power) = valuation(&weight, self.family().epsilon)? {
                last = last.max(
                    targets
                        .last
                        .checked_sub(power)
                        .ok_or_else(|| Error::Limit("target Laurent order overflow".into()))?,
                );
            }
        }
        EpsilonRange::new(leading, last)
    }

    /// Apply the retained exact target reduction at this boundary's coordinates.
    /// Insufficient master orders and excessive propagated uncertainty are typed
    /// errors; missing positive orders are never silently replaced by zero.
    pub fn project_targets(
        &self,
        boundary: &CachedBoundary,
        range: EpsilonRange,
        digits: u32,
    ) -> Result<Vec<ProjectedLaurentExpansion>> {
        EpsilonRange::new(range.leading, range.last)?;
        boundary.validate()?;
        if boundary.identity.key() != self.flow().identity().key() {
            return Err(Error::InvalidInput(
                "target projection has a different physical basis identity".into(),
            ));
        }
        if digits == 0 || digits > boundary.accuracy.verified_digits() {
            return Err(Error::Accuracy(
                "target projection requires sufficient inherited accuracy".into(),
            ));
        }
        let master_leading = self.ordinary_master_leading()?;
        if boundary.range.leading > master_leading {
            return Err(Error::InvalidInput(format!(
                "target projection requires master coefficients from epsilon order {master_leading}; omitted pole coefficients cannot be assumed zero"
            )));
        }
        let coordinates = boundary.point.restart_coordinates()?;
        let substitutions = coordinates
            .into_iter()
            .map(|(s, a)| (Atom::var(s), a))
            .collect();
        let p = Precision {
            bits: boundary.accuracy.working_bits(),
        };
        let refined = Precision {
            bits: p
                .bits
                .checked_add(64)
                .ok_or_else(|| Error::Limit("target projection precision overflow".into()))?,
        };
        let columns = self
            .basis()
            .iter()
            .enumerate()
            .map(|(column, integral)| (integral, column))
            .collect::<BTreeMap<_, _>>();
        let mut result = Vec::with_capacity(self.target_reductions().len());
        for target in self.target_reductions() {
            let mut rows = Vec::new();
            for (integral, weight) in target {
                let column = *columns.get(integral).ok_or_else(|| {
                    Error::IncompleteReduction(
                        "target reduction leaves the physical master basis".into(),
                    )
                })?;
                let weight = self.cache_target_weight(column, weight);
                let weight = crate::family::substitute(&weight, &substitutions)
                    .together()
                    .cancel();
                let Some(leading) = valuation(&weight, self.family().epsilon)? else {
                    continue;
                };
                let required_last = range
                    .last
                    .checked_sub(leading)
                    .ok_or_else(|| Error::Limit("target Laurent order overflow".into()))?;
                if required_last > boundary.range.last {
                    return Err(Error::InvalidInput(format!(
                        "target reduction requires master epsilon order {required_last}, but boundary ends at {}",
                        boundary.range.last
                    )));
                }
                let last_weight = range
                    .last
                    .checked_sub(boundary.range.leading)
                    .ok_or_else(|| Error::Limit("target weight order overflow".into()))?;
                if leading > last_weight {
                    continue;
                }
                let depth = i64::from(last_weight) + 1;
                let series =
                    expansion(&weight, self.family().epsilon, SeriesDepth::absolute(depth))?;
                rows.push((column, leading, series));
            }
            let count = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
            let mut values = vec![p.zero(); count];
            let mut high_values = vec![refined.zero(); count];
            let mut input_errors = vec![refined.zero(); count];
            for (column, leading, series) in &rows {
                // Shift both Laurent series to ordinary series, then delegate
                // the products to Symbolica. Keep the absolute truncation length:
                // trailing zeros must not turn unknown higher terms into zeros.
                let offset = i64::from(boundary.range.leading) + i64::from(*leading);
                let length = i64::from(range.last) - offset + 1;
                if length <= 0 {
                    continue;
                }
                let length = usize::try_from(length)
                    .map_err(|_| Error::Limit("target product length overflow".into()))?;
                if length > boundary.coefficients.len() {
                    return Err(Error::InvalidInput(
                        "target product needs additional master coefficients".into(),
                    ));
                }
                let mut weights_low = Vec::with_capacity(length);
                let mut weights_high = Vec::with_capacity(length);
                let mut weights_norm = Vec::with_capacity(length);
                let mut masters = Vec::with_capacity(length);
                let mut errors = Vec::with_capacity(length);
                for index in 0..length {
                    let power = i64::from(*leading) + index as i64;
                    let weight = series.coefficient(Rational::from(power)).ok_or_else(|| {
                        Error::Accuracy(
                            "target weight coefficient lies in the unknown series remainder".into(),
                        )
                    })?;
                    let low = p.eval(&weight, &Default::default())?;
                    let high = refined.eval(&weight, &Default::default())?;
                    weights_norm.push(ComplexFloat::new(refined.norm(&high), refined.real(0)));
                    weights_low.push(low);
                    weights_high.push(high);
                    let master = &boundary.coefficients[index][*column];
                    let scale = maximum(refined.norm(master), refined.real(1));
                    let floor = refined.tolerance(boundary.accuracy.verified_digits()) * scale;
                    let error = maximum(
                        boundary.accuracy.comparison_errors()[index][*column].clone(),
                        floor,
                    );
                    errors.push(ComplexFloat::new(error, refined.real(0)));
                    masters.push(master.clone());
                }
                let multiply = |p, left: &[ComplexFloat], right: &[ComplexFloat]| {
                    let a = crate::fixed_series::fixed_series(p, self.family().epsilon, left)?;
                    let b = crate::fixed_series::fixed_series(p, self.family().epsilon, right)?;
                    crate::fixed_series::coefficients(&(&a * &b), length)
                };
                let low = multiply(p, &weights_low, &masters)?;
                let high = multiply(refined, &weights_high, &masters)?;
                let propagated = multiply(refined, &weights_norm, &errors)?;
                for (index, power) in (range.leading..=range.last).enumerate() {
                    let source = i64::from(power) - offset;
                    if source < 0 {
                        continue;
                    }
                    let source = usize::try_from(source)
                        .map_err(|_| Error::Limit("target product index overflow".into()))?;
                    values[index] = p.add(&values[index], &low[source]);
                    high_values[index] = refined.add(&high_values[index], &high[source]);
                    input_errors[index] = refined.add(&input_errors[index], &propagated[source]);
                }
            }
            let mut coefficients = BTreeMap::new();
            let mut absolute_errors = BTreeMap::new();
            for (index, power) in (range.leading..=range.last).enumerate() {
                let value = &values[index];
                let high = &high_values[index];
                if !p.finite(value) || !refined.finite(high) {
                    return Err(Error::Numerical("nonfinite projected target".into()));
                }
                let error = Float::with_val(
                    p.bits,
                    input_errors[index].re.as_raw()
                        + refined.norm(&refined.sub(high, value)).as_raw(),
                );
                let scale = maximum(p.norm(value), p.real(1));
                if !error.is_finite() || error > p.tolerance(digits) * scale {
                    return Err(Error::Accuracy("target reduction amplifies cached master uncertainty beyond the requested accuracy".into()));
                }
                coefficients.insert(power, value.clone());
                absolute_errors.insert(power, error);
            }
            result.push(ProjectedLaurentExpansion {
                coefficients,
                absolute_errors,
                verified_digits: digits,
                working_bits: p.bits,
            });
        }
        Ok(result)
    }
}

fn maximum(a: Float, b: Float) -> Float {
    if a > b { a } else { b }
}
