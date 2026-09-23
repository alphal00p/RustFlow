use crate::{Error, MassMode, Result};
use rustred::algebra::{Coefficient, CoefficientContext};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Integral(pub Vec<i16>);

/// A propagator affine in scalar products: loop-loop upper triangle first,
/// then loop-external products in loop-major order, as in RustRed.
#[derive(Clone, Debug)]
pub struct Propagator {
    pub constant: Atom,
    pub scalar_products: Vec<Atom>,
}

impl Propagator {
    /// Construct (sum a_i l_i + sum b_j p_j)^2 - mass_squared.
    pub fn quadratic(
        loops: &[i64],
        external: &[i64],
        mass_squared: Atom,
        gram: &[Vec<Atom>],
    ) -> Result<Self> {
        if gram.len() != external.len() || gram.iter().any(|r| r.len() != external.len()) {
            return Err(Error::InvalidInput("external Gram dimensions".into()));
        }
        let mut products = Vec::new();
        for (i, &a) in loops.iter().enumerate() {
            for (j, &b) in loops.iter().enumerate().skip(i) {
                products.push(Atom::num(a) * Atom::num(b) * Atom::num(if i == j { 1 } else { 2 }));
            }
        }
        for &a in loops {
            for &b in external {
                products.push(Atom::num(2) * Atom::num(a) * Atom::num(b));
            }
        }
        let mut constant = -mass_squared;
        for (i, &a) in external.iter().enumerate() {
            for (j, &b) in external.iter().enumerate() {
                constant += Atom::num(a) * Atom::num(b) * &gram[i][j];
            }
        }
        Ok(Self {
            constant: constant.together().cancel(),
            scalar_products: products,
        })
    }
}

#[derive(Clone, Debug)]
pub struct IntegralFamily {
    pub name: String,
    pub loops: Vec<String>,
    pub external: Vec<String>,
    pub external_gram: Vec<Vec<Atom>>,
    pub propagators: Vec<Propagator>,
    pub physical_propagators: usize,
    pub epsilon: Symbol,
    pub dimension: i64,
}

#[derive(Clone, Debug, Default)]
pub struct KinematicPoint(pub BTreeMap<Atom, Atom>);
impl KinematicPoint {
    pub fn apply(&self, a: &Atom) -> Atom {
        substitute(a, &self.0)
    }
}

pub fn substitute(a: &Atom, rules: &BTreeMap<Atom, Atom>) -> Atom {
    if rules.is_empty() {
        return a.clone();
    }
    let variables_only = rules
        .keys()
        .all(|key| matches!(key.as_view(), AtomView::Var(_)));
    a.replace_map(|view, _, out| {
        // Native parameter renaming only matches variables. Cloning every
        // polynomial subtree to look it up cannot produce a match in that case.
        if variables_only && !matches!(view, AtomView::Var(_)) {
            return;
        }
        if let Some(r) = rules.get(&view.to_owned()) {
            **out = r.clone();
        }
    })
}

pub(crate) fn scalar_symbols(a: AtomView<'_>, out: &mut BTreeSet<Atom>) -> Result<()> {
    match a {
        AtomView::Var(_) => {
            out.insert(a.to_owned());
        }
        AtomView::Add(v) => {
            for a in v.iter() {
                scalar_symbols(a, out)?;
            }
        }
        AtomView::Mul(v) => {
            for a in v.iter() {
                scalar_symbols(a, out)?;
            }
        }
        AtomView::Pow(v) => {
            for a in v.iter() {
                scalar_symbols(a, out)?;
            }
        }
        AtomView::Num(n) => match n.get_coeff_view().to_owned() {
            symbolica::coefficient::Coefficient::Complex(c) => {
                if !c.im.is_zero() {
                    out.insert(Atom::var(imaginary_parameter()));
                }
            }
            _ => {
                return Err(Error::InvalidInput(
                    "family inputs must be exact rational or rational-complex coefficients".into(),
                ));
            }
        },
        _ => {
            return Err(Error::Unsupported(
                "family coefficients must be rational functions of scalar symbols".into(),
            ));
        }
    }
    Ok(())
}

/// Reserved formal coefficient parameter used to carry exact complex input
/// through RustRed's rational coefficient field. It is specialized to i during
/// numerical evaluation; all returned nonzero conditions still apply.
pub(crate) fn imaginary_parameter() -> Symbol {
    symbol!("symbolica_amflow::imaginary_unit")
}

pub(crate) fn encode_complex(a: &Atom) -> Atom {
    a.replace_map(|view, _, out| {
        if let AtomView::Num(n) = view
            && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
            && !c.im.is_zero()
        {
            **out = Atom::num(c.re) + Atom::var(imaginary_parameter()) * Atom::num(c.im);
        }
    })
}

pub(crate) struct ConvertedFamily {
    pub family: rustred::family::IntegralFamily,
    pub reverse: BTreeMap<Atom, Atom>,
}

impl IntegralFamily {
    pub fn validate_integral(&self, i: &Integral) -> Result<()> {
        if self.physical_propagators > self.propagators.len() {
            return Err(Error::InvalidInput(
                "invalid physical propagator count".into(),
            ));
        }
        if i.0.len() != self.propagators.len() {
            return Err(Error::InvalidInput(
                "integral index count does not match family".into(),
            ));
        }
        if i.0[self.physical_propagators..].iter().any(|&v| v > 0) {
            return Err(Error::InvalidInput(
                "irreducible numerator slots must have nonpositive powers".into(),
            ));
        }
        Ok(())
    }
    /// Content key for the integral value, independent of denominator order,
    /// unused numerator slots and display names. Loop coordinates are retained.
    pub fn integral_key(&self, integral: &Integral) -> Result<String> {
        self.validate_integral(integral)?;
        let substitution = BTreeMap::from([(
            Atom::var(self.epsilon),
            Atom::var(symbol!("symbolica_amflow::canonical_epsilon")),
        )]);
        let canonical = |a: &Atom| substitute(a, &substitution).together().cancel().to_string();
        let mut factors = BTreeMap::<Vec<String>, i64>::new();
        for (propagator, &power) in self.propagators.iter().zip(&integral.0) {
            if power == 0 {
                continue;
            }
            let key = std::iter::once(canonical(&propagator.constant))
                .chain(propagator.scalar_products.iter().map(&canonical))
                .collect();
            *factors.entry(key).or_default() += i64::from(power);
        }
        factors.retain(|_, power| *power != 0);
        let content = format!(
            "{}:{}:{}:{:?}:{factors:?}",
            self.loops.len(),
            self.external.len(),
            self.dimension,
            self.external_gram
                .iter()
                .map(|row| row.iter().map(&canonical).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        );
        Ok(blake3::hash(content.as_bytes()).to_hex().to_string())
    }

    pub fn at(&self, point: &KinematicPoint) -> Self {
        let mut out = self.clone();
        for d in &mut out.propagators {
            d.constant = point.apply(&d.constant);
            for a in &mut d.scalar_products {
                *a = point.apply(a);
            }
        }
        for a in out.external_gram.iter_mut().flatten() {
            *a = point.apply(a);
        }
        out
    }
    /// Intrinsic squared mass of an ordinary rank-one quadratic propagator.
    /// The affine constant alone also contains the external momentum square.
    pub fn mass_squared(&self, index: usize) -> Result<Atom> {
        let propagator = self
            .propagators
            .get(index)
            .ok_or_else(|| Error::InvalidInput("propagator index out of range".into()))?;
        let l = self.loops.len();
        let e = self.external.len();
        let offset = l * (l + 1) / 2;
        if self.external_gram.len() != e || self.external_gram.iter().any(|row| row.len() != e) {
            return Err(Error::InvalidInput("external Gram dimensions".into()));
        }
        if propagator.scalar_products.len() != offset + l * e {
            return Err(Error::InvalidInput("scalar-product dimensions".into()));
        }
        let mut diagonal = 0;
        for i in 0..l {
            let scale = &propagator.scalar_products[diagonal];
            if !scale.is_zero() {
                let mut mass = -propagator.constant.clone();
                for a in 0..e {
                    for b in 0..e {
                        mass += &propagator.scalar_products[offset + i * e + a]
                            * &propagator.scalar_products[offset + i * e + b]
                            * &self.external_gram[a][b]
                            / (Atom::num(4) * scale);
                    }
                }
                return Ok(mass.together().cancel());
            }
            diagonal += l - i;
        }
        Err(Error::Unsupported(
            "propagator has no squared loop momentum".into(),
        ))
    }

    pub fn deform(&self, eta: Symbol, mode: &MassMode) -> Result<(Self, Vec<bool>)> {
        if self.physical_propagators > self.propagators.len() {
            return Err(Error::InvalidInput(
                "invalid physical propagator count".into(),
            ));
        }
        let selected = match mode {
            MassMode::Auto if self.loops.len() == 1 => (0..self.physical_propagators).collect(),
            MassMode::Auto => vec![
                (0..self.physical_propagators)
                    .find(|&i| self.mass_squared(i).is_ok_and(|mass| !mass.is_zero()))
                    .unwrap_or(0),
            ],
            MassMode::All => (0..self.physical_propagators).collect(),
            MassMode::Propagators(v) => v.clone(),
            MassMode::Mass => {
                let mut groups: Vec<(Atom, Vec<usize>)> = Vec::new();
                for i in 0..self.physical_propagators {
                    let mass = self.mass_squared(i)?;
                    if mass.is_zero() {
                        continue;
                    }
                    if let Some((_, indices)) = groups
                        .iter_mut()
                        .find(|(m, _)| (m - &mass).together().cancel().is_zero())
                    {
                        indices.push(i);
                    } else {
                        groups.push((mass, vec![i]));
                    }
                }
                groups
                    .into_iter()
                    .map(|(_, indices)| indices)
                    .min_by_key(|indices| (indices.len(), indices.clone()))
                    .ok_or_else(|| {
                        Error::Unsupported(
                            "mass placement requires a nonzero intrinsic mass".into(),
                        )
                    })?
            }
            MassMode::Propagator | MassMode::Branch | MassMode::Loop => {
                self.topological_placement(mode)?
            }
        };
        if selected.is_empty() || selected.iter().any(|&i| i >= self.physical_propagators) {
            return Err(Error::InvalidInput(
                "mass insertion must select physical propagators".into(),
            ));
        }
        let mut out = self.clone();
        let mut mask = vec![false; out.propagators.len()];
        for i in selected {
            if !mask[i] {
                out.propagators[i].constant = &out.propagators[i].constant - Atom::var(eta);
                mask[i] = true;
            }
        }
        Ok((out, mask))
    }
    fn topological_placement(&self, mode: &MassMode) -> Result<Vec<usize>> {
        self.validate()?;
        let mut branches: Vec<(Vec<Atom>, Vec<usize>)> = Vec::new();
        for (i, propagator) in self.propagators[..self.physical_propagators]
            .iter()
            .enumerate()
        {
            let direction = crate::regions::branch(propagator, self.loops.len())?;
            if let Some((_, indices)) = branches.iter_mut().find(|(b, _)| *b == direction) {
                indices.push(i);
            } else {
                branches.push((direction, vec![i]));
            }
        }
        if branches.is_empty() {
            return Err(Error::InvalidInput("empty physical topology".into()));
        }
        match mode {
            MassMode::Propagator => {
                let (_, indices) = branches
                    .iter()
                    .max_by_key(|(_, indices)| (indices.len(), std::cmp::Reverse(indices[0])))
                    .unwrap();
                Ok(vec![indices[0]])
            }
            MassMode::Branch => Ok(branches
                .iter()
                .map(|(_, indices)| indices)
                .min_by_key(|indices| (indices.len(), *indices))
                .unwrap()
                .clone()),
            MassMode::Loop => {
                // A degree-(L-1) monomial of the first Symanzik polynomial
                // selects an independent branch span. The possible final
                // factors are precisely the branches outside that span.
                let rank = self
                    .loops
                    .len()
                    .checked_sub(1)
                    .ok_or_else(|| Error::InvalidInput("zero-loop family".into()))?;
                let mut choices = vec![Vec::<usize>::new()];
                for _ in 0..rank {
                    let mut next = Vec::new();
                    for choice in choices {
                        let start = choice.last().map_or(0, |i| i + 1);
                        for i in start..branches.len() {
                            let mut candidate = choice.clone();
                            candidate.push(i);
                            let matrix = candidate.iter().map(|&j| branches[j].0.clone()).collect();
                            if crate::algebra::rref(matrix).1.len() == candidate.len() {
                                next.push(candidate);
                                if next.len() > 10000 {
                                    return Err(Error::Limit(
                                        "loop-placement routing budget exhausted".into(),
                                    ));
                                }
                            }
                        }
                    }
                    choices = next;
                }
                choices
                    .into_iter()
                    .filter_map(|choice| {
                        let matrix: Vec<_> =
                            choice.iter().map(|&i| branches[i].0.clone()).collect();
                        let mut selected = Vec::new();
                        for (direction, indices) in &branches {
                            let mut trial = matrix.clone();
                            trial.push(direction.clone());
                            if crate::algebra::rref(trial).1.len() > rank {
                                selected.extend(indices);
                            }
                        }
                        selected.sort();
                        (!selected.is_empty()).then_some(selected)
                    })
                    .min_by_key(|indices| (indices.len(), indices.clone()))
                    .ok_or_else(|| {
                        Error::Unsupported("topology does not span its loop momenta".into())
                    })
            }
            _ => unreachable!(),
        }
    }
    pub fn validate(&self) -> Result<()> {
        self.convert().map(|_| ())
    }
    pub(crate) fn convert(&self) -> Result<ConvertedFamily> {
        self.convert_at_epsilon(None)
    }
    pub(crate) fn convert_at_epsilon(&self, epsilon: Option<&Rational>) -> Result<ConvertedFamily> {
        if self.physical_propagators > self.propagators.len() {
            return Err(Error::InvalidInput("too many physical propagators".into()));
        }
        let dimension = Atom::num(self.dimension) - Atom::num(2) * Atom::var(self.epsilon);
        let mut symbols = BTreeSet::new();
        scalar_symbols(dimension.as_view(), &mut symbols)?;
        for p in &self.propagators {
            scalar_symbols(p.constant.as_view(), &mut symbols)?;
            for a in &p.scalar_products {
                scalar_symbols(a.as_view(), &mut symbols)?;
            }
        }
        for a in self.external_gram.iter().flatten() {
            scalar_symbols(a.as_view(), &mut symbols)?;
        }
        let names = (0..symbols.len())
            .map(|i| format!("amflow_parameter_{i}"))
            .collect::<Vec<_>>();
        let ctx = CoefficientContext::try_new(names.clone())
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        let forward = symbols
            .into_iter()
            .zip(
                names
                    .iter()
                    .map(|n| ctx.parameter(n).unwrap().to_expression()),
            )
            .collect::<BTreeMap<_, _>>();
        let coefficient = |a: &Atom| -> Result<Coefficient> {
            let a = if let Some(value) = epsilon {
                substitute(
                    a,
                    &BTreeMap::from([(Atom::var(self.epsilon), Atom::num(value.clone()))]),
                )
            } else {
                a.clone()
            };
            let a = substitute(&encode_complex(&a), &forward);
            let v = a
                .try_to_rational_polynomial(&Q, &Z, Some(ctx.one().get_variables().clone()))
                .map_err(|e| Error::InvalidInput(e.to_string()))?;
            if !ctx.contains(&v) {
                return Err(Error::InvalidInput(
                    "undeclared coefficient variables".into(),
                ));
            }
            Ok(v)
        };
        let denominators = self
            .propagators
            .iter()
            .map(|p| {
                Ok(rustred::family::AffineDenominator::new(
                    coefficient(&p.constant)?,
                    p.scalar_products
                        .iter()
                        .map(&coefficient)
                        .collect::<Result<_>>()?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let gram = self
            .external_gram
            .iter()
            .map(|r| r.iter().map(&coefficient).collect::<Result<Vec<_>>>())
            .collect::<Result<Vec<_>>>()?;
        let dim = coefficient(&dimension)?;
        let shifts = vec![ctx.zero(); denominators.len()];
        let family = rustred::family::IntegralFamily::new(
            self.name.clone(),
            self.loops.clone(),
            self.external.clone(),
            ctx,
            dim,
            denominators,
            gram,
            shifts,
        )
        .map_err(|e| Error::InvalidInput(e.to_string()))?;
        Ok(ConvertedFamily {
            family,
            reverse: forward.into_iter().map(|(a, b)| (b, a)).collect(),
        })
    }
}

pub fn eta_derivative(integral: &Integral, shifted: &[bool]) -> Result<Vec<(Integral, Atom)>> {
    if integral.0.len() != shifted.len() {
        return Err(Error::InvalidInput("deformation mask length".into()));
    }
    let mut out = Vec::new();
    for (i, (&n, &q)) in integral.0.iter().zip(shifted).enumerate() {
        if q && n != 0 {
            let mut target = integral.clone();
            target.0[i] = n
                .checked_add(1)
                .ok_or_else(|| Error::Limit("integral index overflow".into()))?;
            out.push((target, Atom::num(n as i64)));
        }
    }
    Ok(out)
}
