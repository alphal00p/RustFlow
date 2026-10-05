//! Local Taylor coordinates. Charts do not change the physical continuation path.
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

/// Choice of regular, boundary-centered Taylor coordinates. Off-center
/// predivision and singular-center matching are separate operations.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LocalCoordinate {
    #[default]
    Identity,
    BalancedMobius,
}

/// A Taylor polynomial's coordinate, distinct from its physical interval.
/// Möbius coefficients are exact binary-rationalized *chart choices*, not
/// additional precision evidence for the physical input or boundary values.
#[derive(Clone, Debug, Default)]
pub enum TaylorCoordinate {
    #[default]
    Identity,
    Mobius(MobiusCoordinate),
}

#[derive(Clone, Debug)]
pub struct MobiusCoordinate {
    center: Atom,
    direction: Atom,
    a: Rational,
    b: Rational,
}

fn exact(p: Precision, value: &C) -> Result<Atom> {
    if !p.finite(value) {
        return Err(Error::InvalidInput(
            "nonfinite local chart coordinate".into(),
        ));
    }
    Ok(Atom::num(Complex::new(
        value.re.to_rational(),
        value.im.to_rational(),
    )))
}

impl TaylorCoordinate {
    /// Choose a real fractional-linear map on the current straight physical
    /// leg. Projected singularities only choose the map; the transformed
    /// connection's actual complex poles determine its convergence radius.
    pub(crate) fn balanced(p: Precision, center: &C, target: &C, poles: &[C]) -> Result<Self> {
        let direction = p.sub(target, center);
        if direction == p.zero() || !p.finite(&direction) {
            return Err(Error::InvalidInput("degenerate local chart leg".into()));
        }
        let mut left: Option<Rational> = None;
        let mut right: Option<Rational> = None;
        for pole in poles {
            let z = p.div(&p.sub(pole, center), &direction);
            if !p.finite(&z) {
                return Err(Error::Accuracy("nonfinite projected chart pole".into()));
            }
            // Boundary-centered charts cannot stop at the real projection of
            // a nonreal pole and then recenter ahead of it (predivision does
            // that). Symmetric radial projections avoid an artificial barrier
            // at a perfectly regular physical point. These only choose the map;
            // all actual transformed poles still bound its convergence disk.
            let projections = if z.im == p.real(0) {
                vec![z.re.to_rational()]
            } else {
                let radius = p.norm(&z).to_rational();
                vec![-radius.clone(), radius]
            };
            for value in projections {
                if value < 0 {
                    if left.as_ref().is_none_or(|old| value > *old) {
                        left = Some(value);
                    }
                } else if value > 0 && right.as_ref().is_none_or(|old| value < *old) {
                    right = Some(value);
                }
            }
        }
        let (a, b) = match (left, right) {
            (Some(l), Some(r)) => (
                (&l * &r * Rational::from(2)) / (&l - &r),
                (&l + &r) / (&l - &r),
            ),
            (None, Some(r)) => (r * Rational::from(2), Rational::from(1)),
            (Some(l), None) => (-l * Rational::from(2), Rational::from(-1)),
            (None, None) => return Ok(Self::Identity),
        };
        if a <= 0 {
            return Err(Error::Accuracy("degenerate balanced local map".into()));
        }
        // Use the exact difference of the rounded endpoint coordinates, not a
        // separately rounded subtraction, so the declared chart leg is exact.
        let c = exact(p, center)?;
        let d = exact(p, target)? - &c;
        Ok(Self::Mobius(MobiusCoordinate {
            center: c,
            direction: d,
            a,
            b,
        }))
    }

    pub(crate) fn expression(&self, variable: Symbol) -> Result<Atom> {
        let Self::Mobius(m) = self else {
            return Err(Error::InvalidInput(
                "identity charts use physical offsets".into(),
            ));
        };
        let y = Atom::var(variable);
        Ok((&m.center
            + &m.direction * Atom::num(m.a.clone()) * &y
                / (Atom::num(1) + Atom::num(m.b.clone()) * y))
            .together()
            .cancel())
    }
    pub(crate) fn denominator(&self, variable: Symbol) -> Atom {
        match self {
            Self::Identity => Atom::num(1),
            Self::Mobius(m) => Atom::num(1) + Atom::num(m.b.clone()) * Atom::var(variable),
        }
    }
    /// Polynomial argument at a physical point. `center` remains the physical
    /// known boundary even when the polynomial uses a fractional-linear map.
    pub fn local_point(&self, p: Precision, center: &C, point: &C) -> Result<C> {
        match self {
            Self::Identity => Ok(p.sub(point, center)),
            Self::Mobius(m) => {
                if exact(p, center)? != m.center {
                    return Err(Error::InvalidInput(
                        "saved chart and physical center differ".into(),
                    ));
                }
                let t = ((exact(p, point)? - &m.center) / &m.direction)
                    .together()
                    .cancel();
                let denominator = (Atom::num(m.a.clone()) - Atom::num(m.b.clone()) * &t)
                    .together()
                    .cancel();
                if denominator.is_zero() {
                    return Err(Error::Accuracy(
                        "physical point is at local-map infinity".into(),
                    ));
                }
                let result = p.eval(&(t / denominator).together().cancel(), &Default::default())?;
                if !p.finite(&result) {
                    return Err(Error::Accuracy("nonfinite inverse local map".into()));
                }
                Ok(result)
            }
        }
    }
    pub(crate) fn physical_point(&self, p: Precision, center: &C, local: &C) -> Result<C> {
        match self {
            Self::Identity => Ok(p.add(center, local)),
            Self::Mobius(m) => {
                let y = exact(p, local)?;
                let denominator = Atom::num(1) + Atom::num(m.b.clone()) * &y;
                if denominator.is_zero() {
                    return Err(Error::Accuracy("local map reached infinity".into()));
                }
                let result = p.eval(
                    &(&m.center + &m.direction * Atom::num(m.a.clone()) * y / denominator)
                        .together()
                        .cancel(),
                    &Default::default(),
                )?;
                if !p.finite(&result) {
                    return Err(Error::Accuracy("nonfinite forward local map".into()));
                }
                Ok(result)
            }
        }
    }
    /// dx/dy at a physical point, used to compare the polynomial derivative
    /// against the original physical differential equation.
    pub(crate) fn jacobian_at(&self, p: Precision, center: &C, point: &C) -> Result<C> {
        match self {
            Self::Identity => Ok(p.complex(1, 0)),
            Self::Mobius(m) => {
                let y = exact(p, &self.local_point(p, center, point)?)?;
                let denominator = Atom::num(1) + Atom::num(m.b.clone()) * y;
                let result = p.eval(
                    &(&m.direction * Atom::num(m.a.clone()) / denominator.pow(2))
                        .together()
                        .cancel(),
                    &Default::default(),
                )?;
                if !p.finite(&result) || result == p.zero() {
                    return Err(Error::Accuracy("singular local-map Jacobian".into()));
                }
                Ok(result)
            }
        }
    }

    /// Substitute every original domain independently, before matrix/Jacobian
    /// cancellation. The artificial map pole is always retained as a guard.
    pub(crate) fn pullback_rows(
        &self,
        variable: Symbol,
        rows: &[Vec<Atom>],
    ) -> Result<(Vec<Vec<Atom>>, Vec<Atom>)> {
        let expression = self.expression(variable)?;
        let jacobian = expression.derivative(variable).together().cancel();
        let replacements = BTreeMap::from([(Atom::var(variable), expression)]);
        let mut symbols = BTreeSet::new();
        for a in rows.iter().flatten() {
            crate::family::scalar_symbols(a.as_view(), &mut symbols)?;
        }
        let variables = symbols
            .into_iter()
            .filter_map(|a| match a.as_view() {
                AtomView::Var(v) => Some(v.get_symbol()),
                _ => None,
            })
            .collect();
        let original = crate::physical_conditions::rational_denominator_conditions(
            &rows.iter().flatten().cloned().collect::<Vec<_>>(),
            &variables,
        )?;
        let mut guards = original
            .iter()
            .map(|g| {
                crate::family::substitute(g, &replacements)
                    .together()
                    .cancel()
            })
            .collect::<Vec<_>>();
        guards.push(self.denominator(variable));
        if guards.iter().any(|guard| guard.is_zero()) {
            return Err(Error::InvalidInput(
                "local map lies on an original source domain exclusion".into(),
            ));
        }
        let rows = rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|a| {
                        (&jacobian * crate::family::substitute(a, &replacements))
                            .together()
                            .cancel()
                    })
                    .collect()
            })
            .collect();
        Ok((rows, guards))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_balanced_map_images_and_inverse() -> Result<()> {
        let p = Precision::decimal(60)?;
        let map = TaylorCoordinate::balanced(p, &p.i(1), &p.i(5), &[p.i(-3), p.i(13)])?;
        let y = symbol!("mobius_map_unit::y");
        let expression = map.expression(y)?;
        for (local, physical) in [(-1, -3), (0, 1), (1, 13)] {
            let exact = crate::family::substitute(
                &expression,
                &BTreeMap::from([(Atom::var(y), Atom::num(local))]),
            );
            assert!((exact - Atom::num(physical)).together().cancel().is_zero());
            assert!(p.close(
                &map.local_point(p, &p.i(1), &p.i(physical))?,
                &p.i(local),
                50
            ));
        }
        let local = p.rational(&Rational::from((1, 5)));
        let physical = map.physical_point(p, &p.i(1), &local)?;
        assert!(p.close(&map.local_point(p, &p.i(1), &physical)?, &local, 50));
        let derivative = p.eval(
            &crate::family::substitute(
                &expression.derivative(y),
                &BTreeMap::from([(Atom::var(y), exact(p, &local)?)]),
            ),
            &Default::default(),
        )?;
        assert!(p.close(&map.jacobian_at(p, &p.i(1), &physical)?, &derivative, 50));
        Ok(())
    }

    #[test]
    fn one_sided_limits_and_entire_identity_are_explicit() -> Result<()> {
        let p = Precision::decimal(50)?;
        for (pole, local) in [(-2, -1), (2, 1)] {
            let map = TaylorCoordinate::balanced(p, &p.zero(), &p.i(1), &[p.i(pole)])?;
            assert!(p.close(
                &map.physical_point(p, &p.zero(), &p.i(local))?,
                &p.i(pole),
                40
            ));
        }
        assert!(matches!(
            TaylorCoordinate::balanced(p, &p.zero(), &p.i(1), &[])?,
            TaylorCoordinate::Identity
        ));
        assert!(TaylorCoordinate::balanced(p, &p.zero(), &p.zero(), &[]).is_err());
        Ok(())
    }

    #[test]
    fn a_nonreal_pole_does_not_create_a_fake_real_barrier() -> Result<()> {
        let p = Precision::decimal(50)?;
        let pole = p.parse("1e-100", "1")?;
        let map = TaylorCoordinate::balanced(p, &p.zero(), &p.i(1), &[pole])?;
        let next = map.physical_point(p, &p.zero(), &p.rational(&Rational::from((1, 3))))?;
        assert!(next.re > p.rational(&Rational::from((1, 4))).re);
        assert_eq!(next.im, p.real(0));
        Ok(())
    }

    #[test]
    fn map_pole_survives_a_canceling_jacobian() -> Result<()> {
        let p = Precision::decimal(50)?;
        let x = symbol!("mobius_map_guard::x");
        let map = TaylorCoordinate::balanced(p, &p.zero(), &p.i(2), &[p.i(1)])?;
        let rows = vec![vec![Atom::one() / (Atom::var(x) - Atom::one()).pow(2)]];
        let (transformed, guards) = map.pullback_rows(x, &rows)?;
        let excluded = map.denominator(x);
        assert!(guards.contains(&excluded));
        // Here B(y)=2/(y-1)^2 is regular at the artificial map pole y=-1.
        let at_pole = crate::family::substitute(
            &transformed[0][0],
            &BTreeMap::from([(Atom::var(x), Atom::num(-1))]),
        );
        assert!((at_pole - Atom::num((1, 2))).together().cancel().is_zero());
        assert!(map.physical_point(p, &p.zero(), &p.i(-1)).is_err());
        Ok(())
    }
}
