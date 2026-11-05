//! Exact reduction domains used by physical boundary reuse.
use crate::{
    DifferentialSystem, Error, Precision, Result, family::substitute, kinematics::KinematicPath,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;

type ExactPolynomial = MultivariatePolynomial<IntegerRing, u16>;

/// Opt-in metadata only. `bits=exact` means no numerical precision is in scope.
/// End markers denote scope exit, including an error or early rejection.
pub(crate) struct PhysicalPhaseTiming {
    phase: &'static str,
    index: usize,
    subindex: usize,
    bits: Option<u32>,
    start: Option<std::time::Instant>,
}

impl PhysicalPhaseTiming {
    pub(crate) fn start(
        phase: &'static str,
        index: usize,
        subindex: usize,
        bits: Option<u32>,
    ) -> Self {
        let timer = Self {
            phase,
            index,
            subindex,
            bits,
            start: (std::env::var_os("RUSTFLOW_PHYSICAL_PHASE_TIMINGS").is_some_and(|v| v == "1"))
                .then(std::time::Instant::now),
        };
        if timer.start.is_some() {
            timer.emit("begin", 0.0);
        }
        timer
    }

    fn emit(&self, event: &str, seconds: f64) {
        use std::io::Write;
        let bits = self
            .bits
            .map_or_else(|| "exact".into(), |bits| bits.to_string());
        let _ = writeln!(
            std::io::stderr().lock(),
            "RUSTFLOW_PHYSICAL_PHASE event={event} phase={} index={} subindex={} bits={bits} elapsed_seconds={seconds:.6}",
            self.phase,
            self.index,
            self.subindex
        );
    }
}

impl Drop for PhysicalPhaseTiming {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            self.emit("end", start.elapsed().as_secs_f64());
        }
    }
}

/// Validate the rational domain before canonical cancellation can erase a pole.
fn rational_expression(
    a: AtomView<'_>,
    variables: &BTreeSet<Symbol>,
    poles: &mut Vec<Atom>,
) -> Result<()> {
    match a {
        AtomView::Num(n) if matches!(n.get_coeff_view().to_owned(), Coefficient::Complex(_)) => Ok(()),
        AtomView::Var(v) if variables.contains(&v.get_symbol()) => Ok(()),
        AtomView::Add(v) => v.iter().try_for_each(|a| rational_expression(a, variables, poles)),
        AtomView::Mul(v) => v.iter().try_for_each(|a| rational_expression(a, variables, poles)),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            if let AtomView::Num(n) = exponent
                && let Coefficient::Complex(power) = n.get_coeff_view().to_owned()
                && power.im.is_zero() && power.re.is_integer()
            {
                if power.re < 0 { poles.push(base.to_owned()); }
                return rational_expression(base, variables, poles);
            }
            Err(Error::Unsupported("reduction conditions require integer rational powers".into()))
        }
        _ => Err(Error::Unsupported("reduction conditions require exact rational-complex expressions in declared physical variables and epsilon".into())),
    }
}

fn rational_parts(a: &Atom) -> Result<[Atom; 2]> {
    let rational: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(a)
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::Unsupported(format!("rational reduction condition required: {e}")))?;
    let imaginary = BTreeMap::from([(
        Atom::var(crate::family::imaginary_parameter()),
        Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
    )]);
    Ok([
        rational.numerator.to_expression(),
        rational.denominator.to_expression(),
    ]
    .map(|part| substitute(&part, &imaginary).expand()))
}

/// Split only visible factors, after rational_expression retained all poles.
/// Over rational functions, a product (or a nonzero integer power) is nonzero
/// exactly when every factor (or its base) is nonzero. Do not factor sums.
fn nonzero_factors(a: AtomView<'_>, out: &mut Vec<Atom>) {
    match a {
        AtomView::Mul(v) => v.iter().for_each(|factor| nonzero_factors(factor, out)),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            if let AtomView::Num(n) = exponent
                && let Coefficient::Complex(power) = n.get_coeff_view().to_owned()
                && power.im.is_zero()
                && power.re.is_integer()
                && !power.re.is_zero()
            {
                nonzero_factors(base, out);
            } else {
                out.push(a.to_owned());
            }
        }
        _ => out.push(a.to_owned()),
    }
}

/// All returned conditions are polynomials. Every original denominator is
/// retained independently before rational cancellation. Pure epsilon factors
/// remain valid generic conditions; epsilon is not set to zero here.
pub(crate) fn canonical_conditions(
    conditions: &[Atom],
    variables: &BTreeSet<Symbol>,
) -> Result<Vec<Atom>> {
    let _timing = PhysicalPhaseTiming::start("conditions.canonical", 0, 0, None);
    let mut candidates = conditions.to_vec();
    for condition in conditions {
        rational_expression(condition.as_view(), variables, &mut candidates)?;
    }
    let mut factors = Vec::new();
    for candidate in candidates {
        nonzero_factors(candidate.as_view(), &mut factors);
    }
    let mut result = Vec::new();
    for (index, candidate) in factors.into_iter().enumerate() {
        let _timing = PhysicalPhaseTiming::start("conditions.canonical_parts", index, 0, None);
        for part in rational_parts(&candidate)? {
            if part.is_zero() {
                return Err(Error::InvalidInput(
                    "an identically zero reduction condition cannot define a physical domain"
                        .into(),
                ));
            }
            if matches!(part.as_view(), AtomView::Num(_)) {
                continue;
            }
            result.push(part);
        }
    }
    result.sort_by_cached_key(AtomCore::to_canonical_string);
    result.dedup();
    Ok(result)
}

/// Retain original rational denominator bases before a path expression is
/// simplified. Its value may be zero; only denominator domains are constrained.
pub(crate) fn rational_denominator_conditions(
    expressions: &[Atom],
    variables: &BTreeSet<Symbol>,
) -> Result<Vec<Atom>> {
    let mut poles = Vec::new();
    for expression in expressions {
        rational_expression(expression.as_view(), variables, &mut poles)?;
    }
    canonical_conditions(&poles, variables)
}

pub(crate) fn validate_conditions_at(
    conditions: &[Atom],
    epsilon: Symbol,
    coordinates: &BTreeMap<Symbol, Atom>,
) -> Result<()> {
    let substitutions = coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    for (index, condition) in conditions.iter().enumerate() {
        let _timing = PhysicalPhaseTiming::start("conditions.validate_at", index, 0, None);
        let value = substitute(condition, &substitutions);
        rational_expression(value.as_view(), &BTreeSet::from([epsilon]), &mut Vec::new())?;
        if rational_parts(&value)?.iter().any(|a| a.is_zero()) {
            return Err(Error::InvalidInput(
                "physical point violates a required nonzero reduction condition".into(),
            ));
        }
    }
    Ok(())
}

/// Real path values annihilate an epsilon polynomial iff they annihilate all
/// real and imaginary epsilon coefficients. Native exact polynomial gcd finds
/// that common zero set without approximating epsilon or implementing Q(i).
fn common_real_zero_polynomial(
    polynomial: &Atom,
    epsilon: Symbol,
    parameter: Symbol,
    p: Precision,
    index: usize,
    part: usize,
) -> Result<Atom> {
    let _timing = PhysicalPhaseTiming::start("conditions.common_real", index, part, Some(p.bits));
    let mut gcd: Option<ExactPolynomial> = None;
    let variables = Arc::new(vec![PolyVariable::Symbol(parameter)]);
    let mut component_index = 0;
    for (_, coefficient) in
        crate::coefficient::exact_coefficient_list(&polynomial.expand(), &[Atom::var(epsilon)])?
    {
        let mut real = Atom::new();
        let mut imaginary = Atom::new();
        for (monomial, value) in crate::coefficient::exact_coefficient_list(
            &coefficient.expand(),
            &[Atom::var(parameter)],
        )? {
            let AtomView::Num(n) = value.as_view() else {
                return Err(Error::Unsupported(
                    "guarded physical paths require rational-complex coefficients".into(),
                ));
            };
            let Coefficient::Complex(value) = n.get_coeff_view().to_owned() else {
                return Err(Error::Unsupported(
                    "guarded physical paths require exact coefficients".into(),
                ));
            };
            real += &monomial * Atom::num(value.re);
            imaginary += &monomial * Atom::num(value.im);
        }
        for component in [real, imaginary] {
            let subindex = component_index;
            component_index += 1;
            if component.is_zero() {
                continue;
            }
            let timing = PhysicalPhaseTiming::start(
                "conditions.component_parts",
                index,
                subindex,
                Some(p.bits),
            );
            let rational: RationalPolynomial<IntegerRing, u16> = component
                .try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))
                .map_err(|e| {
                    Error::Unsupported(format!("polynomial physical path required: {e}"))
                })?;
            drop(timing);
            if !rational.denominator.is_constant() {
                return Err(Error::Unsupported(
                    "nonpolynomial condition coefficient after clearing denominators".into(),
                ));
            }
            gcd = Some(if let Some(previous) = gcd {
                let _timing =
                    PhysicalPhaseTiming::start("conditions.gcd", index, subindex, Some(p.bits));
                previous.gcd(&rational.numerator)
            } else {
                rational.numerator
            });
            if gcd.as_ref().is_some_and(|p| p.is_constant()) {
                return Ok(Atom::one());
            }
        }
    }
    Ok(gcd.map_or_else(Atom::new, |polynomial| polynomial.to_expression()))
}

/// Admit the complete real parameter interval [0,1], including both endpoints.
/// Algebraic paths are rejected explicitly until their exact common zero set
/// can be represented by the native rational polynomial machinery.
pub(crate) fn conditions_admit_path(
    conditions: &[Atom],
    epsilon: Symbol,
    path: &KinematicPath,
    p: Precision,
    digits: u32,
) -> Result<bool> {
    let _timing = PhysicalPhaseTiming::start("conditions.path", 0, 0, Some(p.bits));
    path.validate()?;
    if conditions.is_empty() {
        return Ok(true);
    }
    if path.parameter == epsilon {
        return Err(Error::InvalidInput(
            "guard path parameter cannot be epsilon".into(),
        ));
    }
    let substitutions = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    let variables = BTreeSet::from([epsilon, path.parameter]);
    for (index, condition) in conditions.iter().enumerate() {
        let _timing = PhysicalPhaseTiming::start("conditions.restrict", index, 0, Some(p.bits));
        let restricted = substitute(condition, &substitutions);
        let mut candidates = vec![restricted.clone()];
        rational_expression(restricted.as_view(), &variables, &mut candidates)?;
        let mut factors = Vec::new();
        for candidate in candidates {
            nonzero_factors(candidate.as_view(), &mut factors);
        }
        drop(_timing);
        for (candidate_index, candidate) in factors.into_iter().enumerate() {
            let timing = PhysicalPhaseTiming::start(
                "conditions.path_parts",
                index,
                candidate_index,
                Some(p.bits),
            );
            let parts = rational_parts(&candidate)?;
            drop(timing);
            for (part_index, polynomial) in parts.into_iter().enumerate() {
                let common = common_real_zero_polynomial(
                    &polynomial,
                    epsilon,
                    path.parameter,
                    p,
                    index,
                    2 * candidate_index + part_index,
                )?;
                if common.is_zero() {
                    return Ok(false);
                }
                if common.is_one() {
                    continue;
                }
                for endpoint in [Atom::new(), Atom::one()] {
                    if substitute(
                        &common,
                        &BTreeMap::from([(Atom::var(path.parameter), endpoint)]),
                    )
                    .expand()
                    .is_zero()
                    {
                        return Ok(false);
                    }
                }
                let denominator_probe = DifferentialSystem {
                    variable: path.parameter,
                    matrix: vec![vec![Atom::one() / common]],
                };
                let timing = PhysicalPhaseTiming::start(
                    "conditions.probe_compile",
                    index,
                    2 * candidate_index + part_index,
                    Some(p.bits),
                );
                let compiled = denominator_probe.compile(p, &Default::default())?;
                drop(timing);
                let tolerance = p.tolerance(
                    digits
                        .checked_add(8)
                        .ok_or_else(|| Error::Limit("condition precision overflow".into()))?,
                );
                if compiled.poles.iter().any(|root| {
                    root.re >= -tolerance.clone()
                        && root.re <= p.real(1) + &tolerance
                        && root.im >= -tolerance.clone()
                        && root.im <= tolerance
                }) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

pub(crate) fn epsilon_leading_coefficient(polynomial: &Atom, epsilon: Symbol) -> Result<Atom> {
    let variable = Atom::var(epsilon);
    let mut leading: Option<(u32, Atom)> = None;
    for (monomial, coefficient) in crate::coefficient::exact_coefficient_list(
        &polynomial.expand(),
        std::slice::from_ref(&variable),
    )? {
        if coefficient.is_zero() {
            continue;
        }
        let exponent = if monomial.is_one() {
            0
        } else if monomial == variable {
            1
        } else if let AtomView::Pow(power) = monomial.as_view() {
            let (base, exponent) = power.get_base_exp();
            if base != variable.as_view() {
                return Err(Error::Unsupported(
                    "nonpolynomial epsilon denominator".into(),
                ));
            }
            exponent
                .to_string()
                .parse::<u32>()
                .map_err(|_| Error::Unsupported("noninteger epsilon denominator power".into()))?
        } else {
            return Err(Error::Unsupported(
                "nonpolynomial epsilon denominator".into(),
            ));
        };
        if leading.as_ref().is_none_or(|(old, _)| exponent < *old) {
            leading = Some((exponent, coefficient));
        }
    }
    Ok(leading.map_or_else(Atom::new, |(_, coefficient)| coefficient))
}

/// Original matrix poles and algebraic branch bases survive cancellation by a
/// stationary path coordinate. Their leading epsilon coefficient must stay
/// nonzero for one regular Laurent chart along the entire transport path.
/// This is stricter than a generic reduction guard such as s-epsilon: the
/// matrix 1/(s-epsilon) has coefficients 1/s, 1/s^2, ... and excludes s=0.
pub(crate) fn matrix_domain_conditions(
    system: &crate::kinematics::KinematicSystem,
) -> Result<Vec<Atom>> {
    fn bases(a: AtomView<'_>, out: &mut Vec<Atom>) -> Result<()> {
        match a {
            AtomView::Add(v) => v.iter().try_for_each(|a| bases(a, out)),
            AtomView::Mul(v) => v.iter().try_for_each(|a| bases(a, out)),
            AtomView::Pow(v) => {
                let (base, exponent) = v.get_base_exp();
                if let AtomView::Num(n) = exponent
                    && let Coefficient::Complex(value) = n.get_coeff_view().to_owned()
                    && value.im.is_zero()
                {
                    if value.re < 0 || !value.re.is_integer() {
                        out.push(base.to_owned());
                    }
                    return bases(base, out);
                }
                Err(Error::Unsupported(
                    "physical matrix powers require exact rational exponents".into(),
                ))
            }
            _ => Ok(()),
        }
    }
    let mut poles = Vec::new();
    for entry in system.derivatives.values().flatten().flatten() {
        bases(entry.as_view(), &mut poles)?;
    }
    let mut conditions = Vec::new();
    for (index, pole) in poles.into_iter().enumerate() {
        let _timing = PhysicalPhaseTiming::start("conditions.matrix_domain", index, 0, None);
        for part in rational_parts(&pole)? {
            conditions.push(epsilon_leading_coefficient(&part, system.epsilon)?);
        }
    }
    Ok(conditions)
}

#[cfg(test)]
mod factor_guard_tests {
    use super::*;

    fn symbols() -> (Symbol, Symbol, Symbol) {
        symbol!(
            "factor_guard_test::s",
            "factor_guard_test::epsilon",
            "factor_guard_test::t"
        )
    }

    // Pre-change whole-expression conversion, including original denominators.
    // Feeding these expanded polynomials to the path checker exercises its
    // legacy common-zero/GCD route without the visible product optimization.
    fn legacy_conditions(guard: &Atom, variables: &BTreeSet<Symbol>) -> Vec<Atom> {
        let mut candidates = vec![guard.clone()];
        rational_expression(guard.as_view(), variables, &mut candidates).unwrap();
        candidates
            .iter()
            .flat_map(|a| rational_parts(a).unwrap())
            .filter(|a| !matches!(a.as_view(), AtomView::Num(_)))
            .collect()
    }

    #[test]
    fn nested_rational_products_and_integer_powers_keep_separate_guards() {
        let (s, epsilon, _) = symbols();
        let variables = BTreeSet::from([s, epsilon]);
        let s = Atom::var(s);
        let epsilon = Atom::var(epsilon);
        let nested = &s + Atom::one() / (&s + Atom::one());
        let guard = nested.clone().pow(3) * ((&s + Atom::num(2)) / (&s + Atom::num(3))).pow(-2);
        assert_eq!(
            canonical_conditions(&[guard], &variables).unwrap(),
            canonical_conditions(
                &[
                    nested,
                    &s + Atom::one(),
                    &s + Atom::num(2),
                    &s + Atom::num(3)
                ],
                &variables
            )
            .unwrap()
        );
        let guard = epsilon.clone().pow(3) * (&s - &epsilon).pow(2) / (&s + &epsilon).pow(4);
        assert_eq!(
            canonical_conditions(&[guard], &variables).unwrap(),
            canonical_conditions(&[epsilon.clone(), &s - &epsilon, &s + &epsilon], &variables)
                .unwrap()
        );
        assert!(canonical_conditions(&[s.pow((1, 2))], &variables).is_err());
        assert!(canonical_conditions(&[Atom::new()], &variables).is_err());
    }

    #[test]
    fn canceled_denominator_guards_survive_factor_splitting() {
        let (s, epsilon, parameter) = symbols();
        let variables = BTreeSet::from([s, epsilon]);
        let a = Atom::var(s);
        let canceled = ((&a + Atom::one()).pow(2) - Atom::one()) / (&a * (&a + Atom::num(2)));
        assert!(canceled.together().cancel().is_one());
        // Keep the cancellation inside a sum: factor splitting alone cannot
        // retain these poles; the original denominator walk is essential.
        let guard = ((canceled + Atom::one()) * (&a + Atom::num(3)).pow(2)).pow(-2);
        let conditions = canonical_conditions(&[guard.clone()], &variables).unwrap();
        for value in [0, -2, -3] {
            assert!(
                validate_conditions_at(
                    &conditions,
                    epsilon,
                    &BTreeMap::from([(s, Atom::num(value))])
                )
                .is_err()
            );
        }
        validate_conditions_at(&conditions, epsilon, &BTreeMap::from([(s, Atom::num(1))])).unwrap();
        let legacy = legacy_conditions(&guard, &variables);
        for (start, end, expected) in [(-1, 1, false), (-3, -1, false), (1, 2, true)] {
            let path = KinematicPath::straight_line(
                parameter,
                &BTreeMap::from([(s, Atom::num(start))]),
                &BTreeMap::from([(s, Atom::num(end))]),
            )
            .unwrap();
            let p = Precision::decimal(60).unwrap();
            assert_eq!(
                conditions_admit_path(&conditions, epsilon, &path, p, 20).unwrap(),
                expected
            );
            assert_eq!(
                conditions_admit_path(&legacy, epsilon, &path, p, 20).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn factor_guards_match_legacy_exact_domain_path_counterexamples() {
        let (s, epsilon, parameter) = symbols();
        let variables = BTreeSet::from([s, epsilon]);
        let a = Atom::var(s);
        let e = Atom::var(epsilon);
        let i = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
        let cases = [
            (&e * (&a - &e), Atom::num(-1), Atom::num(1), true),
            (&e * &a, Atom::num(-1), Atom::num(1), false),
            ((&a - &e) / (&a + &e), Atom::num(-1), Atom::num(1), true),
            (&a / (&a - &e), Atom::num(-1), Atom::num(1), false),
            (
                (&a - Atom::one()) * (&a + Atom::one()),
                Atom::num(-2),
                Atom::num(2),
                false,
            ),
            (
                (&a - Atom::one()) * (&a + Atom::one()),
                Atom::num(1),
                Atom::num(2),
                false,
            ),
            (
                (&a - Atom::one()) * (&a + Atom::one()),
                Atom::num(2),
                Atom::num(3),
                true,
            ),
            ((&a + &i) * (&a - &i), Atom::num(-2), Atom::num(2), true),
            ((&a + &i) * (&a - &i), -i.clone(), i.clone(), false),
            (
                (&a + &i * &e) * (&a - &i * &e),
                Atom::num(-1),
                Atom::num(1),
                true,
            ),
        ];
        let p = Precision::decimal(60).unwrap();
        for (guard, start, end, expected) in cases {
            let split = canonical_conditions(&[guard.clone()], &variables).unwrap();
            let legacy = legacy_conditions(&guard, &variables);
            let path = KinematicPath::straight_line(
                parameter,
                &BTreeMap::from([(s, start)]),
                &BTreeMap::from([(s, end)]),
            )
            .unwrap();
            assert_eq!(
                conditions_admit_path(&split, epsilon, &path, p, 20).unwrap(),
                expected
            );
            assert_eq!(
                conditions_admit_path(&legacy, epsilon, &path, p, 20).unwrap(),
                expected
            );
        }
        // A nonvanishing formal-epsilon numerator must not hide a pole in a
        // rational coordinate at an interior path parameter.
        let path = KinematicPath {
            parameter,
            coordinates: BTreeMap::from([(
                s,
                Atom::one() / (Atom::var(parameter) - Atom::num((1, 2))),
            )]),
        };
        let guard = &e * (&a - &e);
        for conditions in [
            canonical_conditions(&[guard.clone()], &variables).unwrap(),
            legacy_conditions(&guard, &variables),
        ] {
            assert!(!conditions_admit_path(&conditions, epsilon, &path, p, 20).unwrap());
        }
    }
}
