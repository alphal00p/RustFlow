use crate::{Error,Result,IntegralFamily};
use rustred::algebra::{Coefficient,CoefficientContext};
use std::collections::{BTreeMap,BTreeSet};
use symbolica::prelude::*;
pub use symbolica_amflow::family::substitute;
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
pub(crate) struct ConvertedFamily {pub family:rustred::family::IntegralFamily,pub reverse:BTreeMap<Atom,Atom>}
pub(crate) trait DraftConvert {fn convert(&self)->Result<ConvertedFamily>;fn convert_at_epsilon(&self,epsilon:Option<&Rational>)->Result<ConvertedFamily>;}
impl DraftConvert for IntegralFamily {
    fn convert(&self) -> Result<ConvertedFamily> {
        self.convert_at_epsilon(None)
    }
    fn convert_at_epsilon(&self, epsilon: Option<&Rational>) -> Result<ConvertedFamily> {
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
