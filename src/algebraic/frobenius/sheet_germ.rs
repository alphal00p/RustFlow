//! Exact connection between a Gaussian-leading formal germ and supplied sheets.
//! Native certified quadratic disks select signs; floating centers never do.
use super::*;
use crate::frobenius::ExactFrobeniusLimits;
use crate::frobenius::exact::{Gaussian, field, gaussian_bounded};
use symbolica::poly::univariate::{IsolatedRoot, RootLocation, UnivariatePolynomial};

fn rational_sqrt(value: &Rational) -> Option<Rational> {
    if value < &Rational::zero() {
        return None;
    }
    let n = value.numerator_ref().root(2);
    let d = value.denominator_ref().root(2);
    if &n * &n == *value.numerator_ref() && &d * &d == *value.denominator_ref() {
        Some(Rational::from((n, d)))
    } else {
        None
    }
}

fn gaussian_sqrt(value: &Gaussian) -> Result<Gaussian> {
    let unsupported = || {
        Error::Unsupported("constrained root endpoints currently require Gaussian-rational leading constants for each actual subgroup generator; native algebraic-extension coefficient domains remain open".into())
    };
    let norm =
        rational_sqrt(&(&value.re * &value.re + &value.im * &value.im)).ok_or_else(unsupported)?;
    let re = rational_sqrt(&((&norm + &value.re) / Rational::from(2))).ok_or_else(unsupported)?;
    let mut im =
        rational_sqrt(&((&norm - &value.re) / Rational::from(2))).ok_or_else(unsupported)?;
    if value.im < Rational::zero() {
        im = -im;
    }
    let result = Gaussian::new(re, im);
    if &result * &result != *value {
        return Err(unsupported());
    }
    Ok(result)
}

fn principal_root(
    value: &Gaussian,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<IsolatedRoot> {
    if value.is_zero() {
        return Err(Error::InvalidInput(
            "zero matching root in exact sheet proof".into(),
        ));
    }
    let bits = crate::frobenius::exact::height(value);
    if bits.saturating_add(1).saturating_mul(64) > limits.max_coefficient_bits {
        return Err(Error::Limit(
            "exact quadratic root isolation height estimate exceeds limit".into(),
        ));
    }
    context.cancellation.check()?;
    let polynomial = UnivariatePolynomial::from_coefficients(
        &field(),
        vec![
            -value.clone(),
            Gaussian::new(Rational::zero(), Rational::zero()),
            Gaussian::new(Rational::one(), Rational::zero()),
        ],
        std::sync::Arc::new(PolyVariable::Symbol(symbol!(
            "symbolica_amflow::sheet_quadratic"
        ))),
    );
    let roots = polynomial.isolate_roots();
    context.cancellation.check()?;
    for (mut root, _) in roots {
        let location = root.classify_location();
        for attempt in 0..12 {
            context.cancellation.check()?;
            let disk = root.enclosure();
            let coordinate = if location == RootLocation::Imaginary {
                &disk.center().im
            } else {
                &disk.center().re
            };
            if coordinate - disk.radius() > Rational::zero() {
                return Ok(root);
            }
            if coordinate + disk.radius() < Rational::zero() {
                break;
            }
            let bits = 32_u32.checked_shl(attempt).unwrap_or(u32::MAX);
            if u64::from(bits) > limits.max_coefficient_bits {
                break;
            }
            let tolerance = Rational::one() / Rational::from(Integer::from(2).pow(u64::from(bits)));
            root = root.refined(&tolerance);
        }
    }
    Err(Error::Accuracy(
        "native quadratic disks could not certify the principal sheet within resource limits"
            .into(),
    ))
}

fn excludes_zero(value: &ComplexBall) -> bool {
    value.re.center > value.re.radius
        || value.re.center < -value.re.radius.clone()
        || value.im.center > value.im.radius
        || value.im.center < -value.im.radius.clone()
}

fn integer_power(value: &ComplexBall, power: i64, p: Precision) -> ComplexBall {
    let field = FloatField::from_rep(stored_ball(&p.zero(), p));
    let base = if power < 0 {
        field.inv(value)
    } else {
        value.clone()
    };
    field.pow(&base, power.unsigned_abs())
}

pub(super) fn matching_root_balls(
    lifted: &RationalAlgebraicSystem,
    matching: &Atom,
    seeds: &BTreeMap<Symbol, RootSeed>,
    p: Precision,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Vec<ComplexBall>> {
    if lifted.source.roots.is_empty() {
        return Ok(Vec::new());
    }
    if u64::from(p.bits) > limits.max_coefficient_bits {
        return Err(Error::Limit(
            "matching root refinement exceeds exact work limit".into(),
        ));
    }
    let rules = BTreeMap::from([(Atom::var(lifted.rational.variable), matching.clone())]);
    let tolerance = Rational::one() / Rational::from(Integer::from(2).pow(u64::from(p.bits)));
    lifted
        .source
        .roots
        .iter()
        .map(|root| {
            context.cancellation.check()?;
            crate::frobenius::exact::preflight_operation(
                &[&root.radicand, matching],
                lifted.rational.variable,
                2,
                1,
                limits,
                context,
            )?;
            let value = gaussian_bounded(&substitute(&root.radicand, &rules), limits, context)?;
            let selected = principal_root(&value, limits, context)?.refined(&tolerance);
            context.cancellation.check()?;
            let ball = selected.enclosure().to_ball(p.bits);
            match seeds.get(&root.symbol).unwrap_or(&RootSeed::Principal) {
                RootSeed::Principal => Ok(ball),
                RootSeed::Opposite => Ok(-ball),
                _ => Err(Error::Unsupported(
                    "exact endpoint sheet proof requires discrete root germs".into(),
                )),
            }
        })
        .collect()
}

/// Exact valuation and leading coefficient use the existing native Gaussian
/// rational-polynomial converter, also used by branch-path certification.
pub(super) fn valuation(
    expression: &Atom,
    variable: Symbol,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<(i64, Gaussian)> {
    crate::frobenius::exact::preflight_operation(&[expression], variable, 1, 1, limits, context)?;
    let [numerator, denominator] = crate::root_path::rational_parts(expression, variable)?;
    context.cancellation.check()?;
    let n = numerator
        .coefficients()
        .iter()
        .position(|a| !a.is_zero())
        .ok_or_else(|| Error::InvalidInput("identically zero root monomial".into()))?;
    let d = denominator
        .coefficients()
        .iter()
        .position(|a| !a.is_zero())
        .ok_or_else(|| Error::InvalidInput("zero root denominator".into()))?;
    if n.max(d) > limits.max_order {
        return Err(Error::Limit(
            "exact root valuation exceeds order limit".into(),
        ));
    }
    Ok((
        n as i64 - d as i64,
        &numerator.coefficients()[n] / &denominator.coefficients()[d],
    ))
}

pub(super) fn monomial_radicand(
    lifted: &RationalAlgebraicSystem,
    mask: usize,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Atom> {
    let inputs = lifted
        .source
        .roots
        .iter()
        .enumerate()
        .filter(|(i, _)| mask & (1 << i) != 0)
        .map(|(_, root)| &root.radicand)
        .collect::<Vec<_>>();
    crate::frobenius::exact::preflight_operation(
        &inputs,
        lifted.rational.variable,
        2,
        1,
        limits,
        context,
    )?;
    Ok(super::intersection_factor(&lifted.source.roots, mask, mask))
}

/// Returns (valuation, selected leading coefficient, H), where the physical
/// monomial is a*z^(valuation/2)*sqrt(H), H(0)=1, on the shared log branch.
#[allow(clippy::too_many_arguments)]
pub(super) fn selected_germ(
    lifted: &RationalAlgebraicSystem,
    mask: usize,
    matching: &Atom,
    seeds: &BTreeMap<Symbol, RootSeed>,
    winding: i32,
    p: Precision,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<(i64, Gaussian, Atom)> {
    let z = lifted.rational.variable;
    let radicand = monomial_radicand(lifted, mask, limits, context)?;
    let (v, leading) = valuation(&radicand, z, limits, context)?;
    if crate::frobenius::exact::height(&leading)
        .saturating_add(1)
        .saturating_mul(16)
        > limits.max_coefficient_bits
    {
        return Err(Error::Limit(
            "exact leading square-root height estimate exceeds limit".into(),
        ));
    }
    crate::frobenius::exact::preflight_operation(&[&radicand, matching], z, 4, 1, limits, context)?;
    let a = gaussian_sqrt(&leading)?;
    let h = (&radicand / (Atom::num(leading.clone()) * Atom::var(z).pow(v)))
        .together()
        .cancel();
    let match_value = gaussian_bounded(matching, limits, context)?;
    let rules = BTreeMap::from([(Atom::var(z), matching.clone())]);
    let at_match = gaussian_bounded(&substitute(&radicand, &rules), limits, context)?;
    let t = symbol!("symbolica_amflow::sheet_radial_parameter");
    let path = (substitute(
        &radicand,
        &BTreeMap::from([(Atom::var(z), Atom::var(t) * matching)]),
    ) / Atom::var(t).pow(v))
    .together()
    .cancel();
    crate::frobenius::exact::preflight_operation(&[&path], t, 4, 1, limits, context)?;
    let flip = crate::root_path::principal_flip(&path, t, p, context)?.ok_or_else(|| {
        Error::InvalidInput("root germ radial path has a branch point or pole".into())
    })?;
    let at_zero = gaussian_bounded(
        &substitute(&path, &BTreeMap::from([(Atom::var(t), Atom::new())])),
        limits,
        context,
    )?;
    let start_root = principal_root(&at_zero, limits, context)?;
    let end_root = principal_root(&at_match, limits, context)?;
    let coordinate_root = principal_root(&match_value, limits, context)?;
    let mut physical_roots = Vec::new();
    let mut declared_flip = false;
    for (i, root) in lifted.source.roots.iter().enumerate() {
        if mask & (1 << i) == 0 {
            continue;
        }
        context.cancellation.check()?;
        match seeds.get(&root.symbol).unwrap_or(&RootSeed::Principal) {
            RootSeed::Principal => {},
            RootSeed::Opposite => declared_flip = !declared_flip,
            _ => return Err(Error::Unsupported("exact endpoint sheet proof needs discrete principal/opposite germs, not numerical sign hints".into())),
        }
        crate::frobenius::exact::preflight_operation(
            &[&root.radicand, matching],
            z,
            2,
            1,
            limits,
            context,
        )?;
        let value = gaussian_bounded(&substitute(&root.radicand, &rules), limits, context)?;
        physical_roots.push(principal_root(&value, limits, context)?);
    }
    for attempt in 0..10 {
        context.cancellation.check()?;
        let bits = 64_u32.checked_shl(attempt).unwrap_or(u32::MAX);
        if u64::from(bits) > limits.max_coefficient_bits {
            break;
        }
        let p = Precision { bits };
        let tolerance = Rational::one() / Rational::from(Integer::from(2).pow(u64::from(bits)));
        let ball = |r: &IsolatedRoot| r.clone().refined(&tolerance).enclosure().to_ball(bits);
        let coordinate = ball(&coordinate_root);
        if v < 0 && !excludes_zero(&coordinate) {
            continue;
        }
        let phase = &crate::ode::source::exact_ball(&a, p) * &integer_power(&coordinate, v, p);
        let start = ball(&start_root);
        let phase_flip = if excludes_zero(&(&phase + &start)) {
            Some(false)
        } else if excludes_zero(&(&phase - &start)) {
            Some(true)
        } else {
            None
        };
        let physical = physical_roots
            .iter()
            .fold(stored_ball(&p.i(1), p), |a, r| &a * &ball(r));
        let end = ball(&end_root);
        let product_flip = if excludes_zero(&(&physical + &end)) {
            Some(false)
        } else if excludes_zero(&(&physical - &end)) {
            Some(true)
        } else {
            None
        };
        if let (Some(phase_flip), Some(product_flip)) = (phase_flip, product_flip) {
            let winding_flip = v.rem_euclid(2) != 0 && winding.rem_euclid(2) != 0;
            return Ok((
                v,
                if phase_flip ^ flip ^ product_flip ^ declared_flip ^ winding_flip {
                    -a
                } else {
                    a
                },
                h,
            ));
        }
    }
    Err(Error::Accuracy("native root disks could not separate the exact product/phase sheet signs within resource limits".into()))
}
