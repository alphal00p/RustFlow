//! Real-axis contours with explicitly declared polynomial i0 prescriptions.
use crate::{ComplexFloat as C, Error, Precision, Prescription, Result, RunContext};
use std::{collections::BTreeSet, sync::Arc};
use symbolica::domains::float::{FloatField, RealBall};
use symbolica::domains::rational::RationalField;
use symbolica::poly::univariate::{IsolatedRoot, RootLocation, UnivariatePolynomial};
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct PolynomialPrescription {
    pub polynomial: Atom,
    pub prescription: Prescription,
}

#[derive(Clone, Debug)]
pub struct PrescribedContour {
    pub variable: Symbol,
    pub prescriptions: Vec<PolynomialPrescription>,
    /// Absolute imaginary side of an otherwise unprescribed real singularity.
    pub unprescribed_side: Prescription,
}

#[derive(Clone, Debug)]
pub struct ContourCrossing {
    pub point: C,
    pub side: Prescription,
    /// Empty for a singularity assigned the caller's default side.
    pub prescription_indices: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct PlannedContour {
    pub waypoints: Vec<C>,
    pub crossings: Vec<ContourCrossing>,
}

fn exact_polynomial(a: &Atom, variable: Symbol) -> Result<UnivariatePolynomial<RationalField>> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(a.as_view(), &mut symbols)?;
    if !symbols.is_subset(&BTreeSet::from([Atom::var(variable)])) {
        return Err(Error::Unsupported(
            "contour polynomials require exact real coefficients and only the declared variable"
                .into(),
        ));
    }
    let fraction: RationalPolynomial<IntegerRing, u16> =
        a.try_to_rational_polynomial(&Q, &Z, None).map_err(|e| {
            Error::Unsupported(format!(
                "prescription requires an exact real polynomial: {e}"
            ))
        })?;
    if !fraction.denominator.is_constant()
        || fraction
            .numerator
            .variables()
            .iter()
            .any(|v| *v != PolyVariable::Symbol(variable))
    {
        return Err(Error::Unsupported(
            "prescription requires an exact univariate polynomial".into(),
        ));
    }
    let denominator = Rational::from(fraction.denominator.get_constant());
    if fraction.numerator.is_constant() {
        return Ok(UnivariatePolynomial::from_coefficients(
            &Q,
            vec![Rational::from(fraction.numerator.get_constant()) / denominator],
            Arc::new(PolyVariable::Symbol(variable)),
        ));
    }
    Ok(fraction
        .numerator
        .to_univariate_from_univariate(0)
        .map_coeff(|c| Rational::from(c.clone()) / &denominator, Q))
}

fn side(sign: bool) -> Prescription {
    if sign {
        Prescription::PlusI0
    } else {
        Prescription::MinusI0
    }
}
// Refine real certificates only when their present width is insufficient.
// Geometric accuracy is determined by separation, not requested decimal digits.
fn narrow_real_root(
    root: &mut IsolatedRoot,
    target: &Rational,
    context: &RunContext,
) -> Result<()> {
    if target <= &Rational::zero() {
        return Err(Error::Accuracy("no positive root-separation bound".into()));
    }
    context.cancellation.check()?;
    if root.enclosure().radius() > target {
        *root = root.clone().refined(target);
        context.cancellation.check()?;
    }
    Ok(())
}

// P(endpoint+z) = c0 + sum(c_k z^k). For |z| <= delta <= 1,
// the nonconstant terms are bounded by delta * sum |c_k| <= |c0|/2.
// This exact disk contains no zero, so it bounds the required real-root
// enclosure width without depending on the initial isolation radius.
fn endpoint_separation(
    polynomial: &UnivariatePolynomial<RationalField>,
    endpoint: &Rational,
    context: &RunContext,
) -> Result<Rational> {
    context.cancellation.check()?;
    let shifted = polynomial.shift_var(endpoint);
    let coefficients = shifted.coefficients();
    let constant = coefficients
        .first()
        .cloned()
        .unwrap_or_else(Rational::zero)
        .abs();
    if constant.is_zero() {
        return Err(Error::InvalidInput(
            "a contour endpoint is a singular zero".into(),
        ));
    }
    let sum = coefficients
        .iter()
        .skip(1)
        .fold(Rational::zero(), |sum, c| sum + c.clone().abs());
    context.cancellation.check()?;
    Ok(if sum.is_zero() {
        Rational::one()
    } else {
        Rational::one().min(constant / (Rational::from(2) * sum))
    })
}

fn real_crossing(
    polynomial: &UnivariatePolynomial<RationalField>,
    root: &mut IsolatedRoot,
    lower: &Rational,
    upper: &Rational,
    context: &RunContext,
) -> Result<bool> {
    if !matches!(
        root.classify_location(),
        RootLocation::Real | RootLocation::Zero
    ) {
        return Ok(false);
    }
    for attempt in 0..2 {
        let disk = root.enclosure();
        let a = &disk.center().re - disk.radius();
        let b = &disk.center().re + disk.radius();
        if b < *lower || a > *upper {
            return Ok(false);
        }
        if a > *lower && b < *upper {
            return Ok(true);
        }
        if attempt == 1 || disk.radius().is_zero() {
            break;
        }
        let separation = endpoint_separation(polynomial, lower, context)?
            .min(endpoint_separation(polynomial, upper, context)?);
        let target = separation / Rational::from(4);
        narrow_real_root(root, &target, context)?;
    }
    Err(Error::Accuracy(
        "singular root cannot be separated from the contour endpoint".into(),
    ))
}

fn certified_derivative_positive(
    polynomial: &UnivariatePolynomial<RationalField>,
    root: &mut IsolatedRoot,
    p: Precision,
    context: &RunContext,
) -> Result<bool> {
    let derivative = polynomial.derivative();
    let mut bits = p.bits;
    for _ in 0..4 {
        context.cancellation.check()?;
        let value = derivative
            .map_coeff(
                |c| RealBall::from_rational_bounds(c, c, bits),
                FloatField::from_rep(RealBall::exact(Float::with_val(bits, 0))),
            )
            .evaluate(&root.enclosure().to_ball(bits).re);
        if value.is_strictly_positive() {
            return Ok(true);
        }
        if value.is_strictly_negative() {
            return Ok(false);
        }
        let target = root.enclosure().radius() / &Rational::from(4);
        if !target.is_zero() {
            narrow_real_root(root, &target, context)?;
        }
        bits = bits
            .checked_add(64)
            .ok_or_else(|| Error::Limit("contour sign precision overflow".into()))?;
    }
    Err(Error::Accuracy(
        "could not certify the sign of the prescribed polynomial derivative".into(),
    ))
}

fn distance_to_disk(
    center: &Rational,
    other_center: &Complex<Rational>,
    other_radius: &Rational,
) -> Rational {
    (&other_center.re - center)
        .abs()
        .max(other_center.im.clone().abs())
        - other_radius
}

fn rounded_triangle(
    p: Precision,
    center: &C,
    interval: &[Rational; 2],
    radius: &Rational,
    distance: &Rational,
    chosen: Prescription,
    direction: i64,
) -> Result<[C; 3]> {
    let exact_center = center.re.to_rational();
    let numeric_radius = p.rational(radius).re;
    let step = C::new(numeric_radius.clone(), p.real(0));
    let imaginary = if chosen == Prescription::PlusI0 {
        numeric_radius
    } else {
        -numeric_radius
    };
    let triangle = [
        p.sub(center, &p.scale(&step, direction, 1)),
        p.add(center, &C::new(p.real(0), imaginary)),
        p.add(center, &p.scale(&step, direction, 1)),
    ];
    // The convex hull of these actual rounded vertices lies inside an L1
    // disk excluding every other root. Its real entry/exit bracket the whole
    // certified interval of the selected real root.
    if triangle.iter().any(|vertex| {
        !p.finite(vertex)
            || (&vertex.re.to_rational() - &exact_center).abs() + vertex.im.to_rational().abs()
                >= *distance
    }) {
        return Err(Error::Accuracy(
            "rounded contour leaves the certified obstacle-free disk".into(),
        ));
    }
    let entry = triangle[0].re.to_rational();
    let exit = triangle[2].re.to_rational();
    if entry.clone().min(exit.clone()) >= interval[0]
        || entry.max(exit) <= interval[1]
        || triangle[1].im == p.real(0)
    {
        return Err(Error::Accuracy(
            "rounded contour fails to separate the real zero".into(),
        ));
    }
    Ok(triangle)
}

impl PrescribedContour {
    /// Construct disjoint triangular detours between regular real endpoints.
    /// At each simple real zero a prescription on P becomes the absolute x-side
    /// sign(i0) * sign(P'). Reverse transport uses the same geometric path.
    /// Multiple prescribed zeros and inconsistent prescriptions at a shared zero are errors.
    /// All roots of the exact singularity and prescription polynomials constrain
    /// detour radii, including nonreal roots. Numerical pole estimates cannot be
    /// substituted for these polynomials: proximity does not prove root identity.
    pub fn plan(
        &self,
        p: Precision,
        singularities: &[Atom],
        start: &C,
        end: &C,
    ) -> Result<PlannedContour> {
        self.plan_with_context(p, singularities, start, end, &RunContext::default())
    }

    /// Plan with cancellation checks between native polynomial/root operations.
    /// Native isolation and refinement are synchronous calls and cannot be
    /// preempted by this token while they are running.
    pub fn plan_with_context(
        &self,
        p: Precision,
        singularities: &[Atom],
        start: &C,
        end: &C,
        context: &RunContext,
    ) -> Result<PlannedContour> {
        context.cancellation.check()?;
        if p.bits < 32
            || !p.finite(start)
            || !p.finite(end)
            || start.im != p.real(0)
            || end.im != p.real(0)
        {
            return Err(Error::InvalidInput(
                "prescribed contours require finite real endpoints and at least 32 working bits"
                    .into(),
            ));
        }
        let (lower, upper) = if start.re <= end.re {
            (start.re.to_rational(), end.re.to_rational())
        } else {
            (end.re.to_rational(), start.re.to_rational())
        };
        struct Obstacle {
            root: IsolatedRoot,
            crossing: bool,
            side: Option<Prescription>,
            indices: Vec<usize>,
        }
        let mut obstacles = Vec::<Obstacle>::new();
        let all = singularities.iter().map(|a| (a, None)).chain(
            self.prescriptions
                .iter()
                .enumerate()
                .map(|(i, d)| (&d.polynomial, Some(i))),
        );
        for (expression, declaration) in all {
            context.cancellation.check()?;
            let polynomial = exact_polynomial(expression, self.variable)?;
            if polynomial.is_zero() {
                return Err(Error::InvalidInput(
                    "a contour polynomial is identically zero".into(),
                ));
            }
            if polynomial.evaluate(&lower).is_zero() || polynomial.evaluate(&upper).is_zero() {
                return Err(Error::InvalidInput(
                    "a contour endpoint is a singular zero".into(),
                ));
            }
            if polynomial.is_constant() {
                continue;
            }
            let roots = polynomial.isolate_roots();
            context.cancellation.check()?;
            for (mut root, multiplicity) in roots {
                context.cancellation.check()?;
                // Isolation already returns certified disks. Complex obstacles
                // retain these disks; neither rounded center equality nor an
                // unnecessarily tiny requested radius is used as a certificate.
                let inside = real_crossing(&polynomial, &mut root, &lower, &upper, context)?;
                let chosen = if let Some(index) = declaration.filter(|_| inside) {
                    if multiplicity != 1 {
                        return Err(Error::Unsupported("multiple prescribed zeros require a local deformation beyond a simple contour detour".into()));
                    }
                    let positive =
                        certified_derivative_positive(&polynomial, &mut root, p, context)?;
                    Some(side(
                        positive
                            == (self.prescriptions[index].prescription == Prescription::PlusI0),
                    ))
                } else {
                    None
                };
                if let Some(old) = obstacles.iter_mut().find(|old| old.root == root) {
                    if root.enclosure().radius() < old.root.enclosure().radius() {
                        old.root = root.clone();
                    }
                    if let Some(chosen) = chosen {
                        if old.side.is_some_and(|previous| previous != chosen) {
                            return Err(Error::Unsupported(
                                "conflicting polynomial prescriptions at a shared real zero".into(),
                            ));
                        }
                        old.side = Some(chosen);
                        old.indices.push(declaration.unwrap());
                    }
                } else {
                    obstacles.push(Obstacle {
                        root,
                        crossing: inside,
                        side: chosen,
                        indices: declaration
                            .filter(|_| chosen.is_some())
                            .into_iter()
                            .collect(),
                    });
                }
            }
        }
        let direction = if start.re <= end.re { 1 } else { -1 };
        let mut triangles = Vec::new();
        for index in 0..obstacles.len() {
            context.cancellation.check()?;
            if !obstacles[index].crossing {
                continue;
            }
            let mut accepted = None;
            for attempt in 0..5 {
                let disk = obstacles[index].root.enclosure();
                let center = C::new(p.rational(&disk.center().re).re, p.real(0));
                let exact_center = center.re.to_rational();
                let uncertainty = disk.radius() + &(&disk.center().re - &exact_center).abs();
                let margin = p.tolerance(p.bits / 5).to_rational()
                    * exact_center.clone().abs().max(Rational::one())
                    * Rational::from(16);
                let mut distance = (&exact_center - &lower)
                    .abs()
                    .min((&upper - &exact_center).abs());
                for (other_index, other) in obstacles.iter().enumerate() {
                    if other_index == index {
                        continue;
                    }
                    let other_disk = other.root.enclosure();
                    // max(|dx|, |dy|) is an exact lower bound on Euclidean
                    // center distance. Subtracting the actual certified radius
                    // gives a lower bound on distance to the other exact root.
                    let separation =
                        distance_to_disk(&exact_center, other_disk.center(), other_disk.radius());
                    distance = distance.min(separation);
                }
                if distance <= margin {
                    // Separate roots from different exact polynomials may
                    // initially have overlapping isolation disks. Refine only
                    // obstructing real certificates, then repeat the unchanged
                    // exact clearance proof. Complex disks remain conservative
                    // obstacles; no center matching or root merging is used.
                    let mut refined = false;
                    if attempt < 4 {
                        for (other_index, other) in obstacles.iter_mut().enumerate() {
                            if other_index == index {
                                continue;
                            }
                            let enclosure = other.root.enclosure();
                            if distance_to_disk(
                                &exact_center,
                                enclosure.center(),
                                enclosure.radius(),
                            ) > margin
                                || enclosure.radius().is_zero()
                            {
                                continue;
                            }
                            let previous = enclosure.radius().clone();
                            if matches!(
                                other.root.classify_location(),
                                RootLocation::Real | RootLocation::Zero
                            ) {
                                narrow_real_root(
                                    &mut other.root,
                                    &(&previous / &Rational::from(4)),
                                    context,
                                )?;
                                refined |= other.root.enclosure().radius() < &previous;
                            }
                        }
                    }
                    if refined {
                        continue;
                    }
                    return Err(Error::Accuracy(
                        "contour obstacle disks are unresolved at working precision".into(),
                    ));
                }
                let radius = (&distance - &margin) / Rational::from(8);
                if radius <= &uncertainty + &margin {
                    if attempt == 4 || disk.radius().is_zero() {
                        return Err(Error::Accuracy(
                            "contour detour is unresolved at working precision".into(),
                        ));
                    }
                    let target = (&radius - &margin) / Rational::from(4);
                    narrow_real_root(&mut obstacles[index].root, &target, context)?;
                    continue;
                }
                let chosen = obstacles[index].side.unwrap_or(self.unprescribed_side);
                let interval = [
                    &disk.center().re - disk.radius(),
                    &disk.center().re + disk.radius(),
                ];
                let triangle =
                    rounded_triangle(p, &center, &interval, &radius, &distance, chosen, direction)?;
                accepted = Some((
                    ContourCrossing {
                        point: center,
                        side: chosen,
                        prescription_indices: obstacles[index].indices.clone(),
                    },
                    triangle,
                ));
                break;
            }
            triangles.push(accepted.ok_or_else(|| {
                Error::Accuracy("could not resolve a real contour crossing".into())
            })?);
        }
        triangles.sort_by(|a, b| a.0.point.re.partial_cmp(&b.0.point.re).unwrap());
        for pair in triangles.windows(2) {
            let right = pair[0].1[0]
                .re
                .to_rational()
                .max(pair[0].1[2].re.to_rational());
            let left = pair[1].1[0]
                .re
                .to_rational()
                .min(pair[1].1[2].re.to_rational());
            if right >= left {
                return Err(Error::Accuracy("certified contour detours overlap".into()));
            }
        }
        if direction < 0 {
            triangles.reverse();
        }
        let mut waypoints = Vec::new();
        let mut crossings = Vec::new();
        for (crossing, triangle) in triangles {
            crossings.push(crossing);
            waypoints.extend(triangle);
        }
        context.cancellation.check()?;
        waypoints.push(end.clone());
        Ok(PlannedContour {
            waypoints,
            crossings,
        })
    }
}

#[cfg(test)]
mod disk_geometry_tests {
    use super::*;

    #[test]
    fn rounded_triangle_rejects_exact_margin_collapse_disks() {
        let p = Precision::decimal(60).unwrap();
        let c = Rational::from((1, 2));
        let margin = Rational::from((Integer::from(16), Integer::from(10).pow(43)));
        let delta = &margin + &Rational::from((Integer::one(), Integer::from(10).pow(68)));
        // These radius-zero disks are exact analytic certificates for the roots
        // 1/2 +/- i*delta of (x-1/2)^2 + delta^2.
        for sign in [-1, 1] {
            let imaginary = &delta * &Rational::from(sign);
            let other = Complex::new(c.clone(), imaginary.clone());
            assert!((&imaginary * &imaginary - &delta * &delta).is_zero());
            let distance = distance_to_disk(&c, &other, &Rational::zero());
            let radius = (&distance - &margin) / Rational::from(8);
            assert!(radius > Rational::zero());
            assert!(matches!(
                rounded_triangle(
                    p,
                    &p.rational(&c),
                    &[c.clone(), c.clone()],
                    &radius,
                    &distance,
                    Prescription::PlusI0,
                    1
                ),
                Err(Error::Accuracy(_))
            ));
        }
    }

    #[test]
    fn disk_radius_reduces_clearance_and_resolved_triangle_is_certified() {
        let p = Precision::decimal(60).unwrap();
        let c = Rational::from((1, 2));
        let other = Complex::new(c.clone(), Rational::from((1, 100)));
        let distance = distance_to_disk(&c, &other, &Rational::from((1, 1000)));
        assert_eq!(distance, Rational::from((9, 1000)));
        let radius = &distance / &Rational::from(8);
        let triangle = rounded_triangle(
            p,
            &p.rational(&c),
            &[c.clone(), c.clone()],
            &radius,
            &distance,
            Prescription::MinusI0,
            -1,
        )
        .unwrap();
        assert!(triangle[0].re > p.rational(&c).re);
        assert!(triangle[2].re < p.rational(&c).re);
        assert!(triangle[1].im < p.real(0));
    }
}
