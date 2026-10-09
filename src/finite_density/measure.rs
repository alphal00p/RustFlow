//! Distribution-weighted identities, derived before energy-shell integration.
//!
//! `C_n(s)=(-1)^(n-1) delta^(n-1)(s)/(n-1)!` for required cuts.
//! `H_0(h)=theta(h)` and `H_n(h)=C_n(h)` for n>=1. Thus dH_0/dh=H_1,
//! whereas dH_n/dh=-n H_(n+1). Upper occupation and lower positive-energy
//! endpoints are separate H factors. Never also differentiate an integrated
//! on-shell energy when using these identities.
use super::InversePropagatorBasis;
use super::guarded::{GuardedIdentity, GuardedIdentityTerm};
use crate::{Error, Result};
use rustred::solver::guarded::{IndexBounds, IndexDomain, IndexRole};
use std::collections::BTreeMap;
use symbolica::prelude::*;

/// Polynomial factors in independent Minkowski scalar/medium coordinates.
/// This owner generates sources; the caller must establish their common
/// contour, support and vanishing total-derivative boundary assumptions.
/// The caller must also establish that products of shell, occupation and
/// virtual factors are well-defined in that prescription. In particular,
/// coincident shell deltas and virtual poles at occupied endpoints are not
/// admitted merely by constructing these sources. The generated domains split
/// occupation indices but leave ordinary indices unrestricted; intersect them
/// with any additional physically admitted index domains before reduction.
pub struct WeightedMeasure<const N: usize> {
    coordinates: Vec<Atom>,
    factors: [Atom; N],
    roles: [IndexRole; N],
    numerator_basis: InversePropagatorBasis,
    numerator_slots: Vec<usize>,
}

impl<const N: usize> WeightedMeasure<N> {
    pub fn new(coordinates: Vec<Atom>, factors: [Atom; N], roles: [IndexRole; N]) -> Result<Self> {
        let numerator_slots = (0..N)
            .filter(|&i| roles[i] != IndexRole::Occupation)
            .collect::<Vec<_>>();
        let physical = numerator_slots
            .iter()
            .map(|&i| factors[i].clone())
            .collect::<Vec<_>>();
        let labels = (0..physical.len() + coordinates.len())
            .map(|i| {
                Atom::parse(
                    &format!("rustflow_weighted::rho_{i}"),
                    "rustflow_weighted",
                    Default::default(),
                )
                .map_err(|e| Error::InvalidInput(e.to_string()))
            })
            .collect::<Result<Vec<_>>>()?;
        let numerator_basis = InversePropagatorBasis::new(coordinates.clone(), physical, labels)?;
        if numerator_basis.slots().len() != numerator_slots.len() {
            return Err(Error::InvalidInput(
                "weighted factors need explicit ordinary numerator-completion slots".into(),
            ));
        }
        for factor in &factors {
            // Includes occupation endpoints: no rational/nonpolynomial weight
            // may silently become an opaque coefficient of a weighted source.
            numerator_basis.convert(factor, &vec![0; numerator_slots.len()])?;
        }
        Ok(Self {
            coordinates,
            factors,
            roles,
            numerator_basis,
            numerator_slots,
        })
    }

    pub fn roles(&self) -> &[IndexRole; N] {
        &self.roles
    }

    fn multiply(
        &self,
        expression: &Atom,
        shift: [i16; N],
        out: &mut BTreeMap<[i16; N], Atom>,
    ) -> Result<()> {
        for (integral, coefficient) in self
            .numerator_basis
            .convert(expression, &vec![0; self.numerator_slots.len()])?
        {
            let mut shifted = shift;
            for (&slot, degree) in self.numerator_slots.iter().zip(integral.0) {
                shifted[slot] = shifted[slot]
                    .checked_add(degree)
                    .ok_or_else(|| Error::Limit("weighted source shift overflow".into()))?;
            }
            *out.entry(shifted).or_default() += coefficient;
        }
        Ok(())
    }

    fn source(
        id: String,
        terms: BTreeMap<[i16; N], Atom>,
        bounds: [IndexBounds; N],
    ) -> Result<Option<GuardedIdentity<N>>> {
        let terms = terms
            .into_iter()
            .filter_map(|(shift, coefficient)| {
                let coefficient = coefficient.together().cancel();
                (!coefficient.is_zero()).then_some(GuardedIdentityTerm { shift, coefficient })
            })
            .collect::<Vec<_>>();
        if terms.is_empty() {
            return Ok(None);
        }
        Ok(Some(GuardedIdentity {
            id,
            terms,
            domain: IndexDomain::new(bounds).map_err(|e| Error::InvalidInput(e.to_string()))?,
            nonzero_conditions: vec![],
        }))
    }

    /// Total derivative with a supplied polynomial vector field. `direction`
    /// contains its action on each scalar coordinate, and `divergence` its exact
    /// divergence in the integration measure. Sources split occupation zero
    /// from positive indices; negative occupation indices are never admitted.
    pub fn ibp(
        &self,
        id: &str,
        indices: &[Symbol; N],
        direction: &[Atom],
        divergence: &Atom,
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        if indices
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != N
            || indices.iter().any(|&index| {
                self.coordinates
                    .iter()
                    .chain(self.factors.iter())
                    .chain(self.numerator_basis.labels())
                    .any(|a| !a.derivative(index).is_zero())
            })
        {
            return Err(Error::InvalidInput("weighted index symbols must be distinct from each other and the measure coordinates/parameters".into()));
        }
        if direction.len() != self.coordinates.len() {
            return Err(Error::InvalidInput(
                "weighted vector-field coordinate count".into(),
            ));
        }
        let occupations = (0..N)
            .filter(|&i| self.roles[i] == IndexRole::Occupation)
            .collect::<Vec<_>>();
        let cases = 1_usize
            .checked_shl(
                u32::try_from(occupations.len())
                    .map_err(|_| Error::Limit("occupation domain count".into()))?,
            )
            .ok_or_else(|| Error::Limit("occupation domain count".into()))?;
        if cases > domain_budget {
            return Err(Error::Limit(
                "weighted occupation domain budget exhausted".into(),
            ));
        }
        let mut sources = Vec::new();
        for case in 0..cases {
            let mut bounds = *IndexDomain::for_roles(&self.roles).bounds();
            for (j, &slot) in occupations.iter().enumerate() {
                bounds[slot] = if (case >> j) & 1 == 0 {
                    IndexBounds::fixed(0)
                } else {
                    IndexBounds::new(Some(1), None)
                        .map_err(|e| Error::InvalidInput(e.to_string()))?
                };
            }
            let mut terms = BTreeMap::new();
            self.multiply(divergence, [0; N], &mut terms)?;
            for slot in 0..N {
                let mut derivative = Atom::new();
                for (coordinate, action) in self.coordinates.iter().zip(direction) {
                    let AtomView::Var(v) = coordinate.as_view() else {
                        unreachable!()
                    };
                    derivative += self.factors[slot].derivative(v.get_symbol()) * action;
                }
                let coefficient = if self.roles[slot] == IndexRole::Occupation
                    && bounds[slot] == IndexBounds::fixed(0)
                {
                    Atom::num(1)
                } else {
                    -Atom::var(indices[slot])
                };
                let mut shift = [0; N];
                shift[slot] = 1;
                self.multiply(&(coefficient * derivative), shift, &mut terms)?;
            }
            if let Some(source) =
                Self::source(format!("{id}/occupation-case-{case}"), terms, bounds)?
            {
                sources.push(source);
            }
        }
        Ok(sources)
    }

    /// Exact polynomial/distribution multiplication identities. These include
    /// h*delta(h)=0 and h*C_n(h)=C_(n-1)(h), n>=2; h*theta(h) is not zero.
    pub fn multiplication_sources(&self) -> Result<Vec<GuardedIdentity<N>>> {
        let mut sources = Vec::new();
        for slot in 0..N {
            let mut cases = Vec::new();
            let mut default_bounds = *IndexDomain::for_roles(&self.roles).bounds();
            if self.roles[slot] == IndexRole::Occupation {
                default_bounds[slot] = IndexBounds::fixed(1);
                cases.push((default_bounds, false));
                default_bounds[slot] = IndexBounds::new(Some(2), None)
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                cases.push((default_bounds, true));
            } else {
                cases.push((default_bounds, true));
            }
            for (case, (bounds, lower)) in cases.into_iter().enumerate() {
                let mut terms = BTreeMap::new();
                self.multiply(&self.factors[slot], [0; N], &mut terms)?;
                if lower {
                    let mut shift = [0; N];
                    shift[slot] = -1;
                    *terms.entry(shift).or_default() -= Atom::num(1);
                }
                if let Some(source) =
                    Self::source(format!("multiplication/{slot}/{case}"), terms, bounds)?
                {
                    sources.push(source);
                }
            }
        }
        Ok(sources)
    }

    /// Lorentz IBPs in Gram coordinates followed by `u.q_i`, with `u.u=1`.
    /// The energies remain independent shell variables throughout this step.
    pub fn lorentz_ibps(
        &self,
        loops: usize,
        dimension: &Atom,
        indices: &[Symbol; N],
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        if self.coordinates.len() != loops * (loops + 1) / 2 + loops {
            return Err(Error::InvalidInput(
                "Lorentz IBP scalar/medium coordinate count".into(),
            ));
        }
        let mut pairs = Vec::new();
        for i in 0..loops {
            for j in i..loops {
                pairs.push((i, j));
            }
        }
        let scalar = |a: usize, b: usize| {
            self.coordinates[pairs
                .iter()
                .position(|&(i, j)| (i, j) == (a.min(b), a.max(b)))
                .unwrap()]
            .clone()
        };
        let energy = |a: usize| self.coordinates[pairs.len() + a].clone();
        let mut sources = self.multiplication_sources()?;
        for differentiated in 0..loops {
            for vector in 0..=loops {
                let contraction = |a| {
                    if vector == loops {
                        energy(a)
                    } else {
                        scalar(a, vector)
                    }
                };
                let mut direction = pairs
                    .iter()
                    .map(|&(a, b)| {
                        Atom::num(i32::from(a == differentiated)) * contraction(b)
                            + Atom::num(i32::from(b == differentiated)) * contraction(a)
                    })
                    .collect::<Vec<_>>();
                for a in 0..loops {
                    direction.push(
                        Atom::num(i32::from(a == differentiated))
                            * if vector == loops {
                                Atom::num(1)
                            } else {
                                energy(vector)
                            },
                    );
                }
                let divergence = if vector == differentiated {
                    dimension.clone()
                } else {
                    Atom::new()
                };
                sources.extend(self.ibp(
                    &format!("lorentz/{differentiated}/{vector}"),
                    indices,
                    &direction,
                    &divergence,
                    domain_budget,
                )?);
            }
        }
        Ok(sources)
    }
}

/// The occupied one-line energy residue before the spatial integral, retaining
/// all upper-surface derivatives. No theta(0) convention is substituted.
#[derive(Clone, Debug)]
pub struct EnergyResidue {
    /// Coefficient multiplying theta(mu-E).
    pub bulk: Atom,
    /// Coefficients multiplying delta^(k)(mu-E), derivatives of their argument.
    pub upper_surface: Vec<Atom>,
}

impl EnergyResidue {
    pub fn simple(energy: Symbol) -> Self {
        Self {
            bulk: -Atom::num(1) / (Atom::num(2) * Atom::var(energy)),
            upper_surface: vec![],
        }
    }
    /// Generate one additional Euclidean propagator power by -1/n d/dm^2
    /// at fixed spatial momentum and chemical potential, where dE/dm^2=1/(2E).
    /// This integrated-shell representation must not be combined with a second
    /// differentiation of the independent-energy shell distribution.
    pub fn raise(&self, energy: Symbol, power: u16) -> Result<Self> {
        if power == 0 {
            return Err(Error::InvalidInput(
                "cannot raise a zero-power residue".into(),
            ));
        }
        let scale = Atom::num(1) / (Atom::num(2 * u32::from(power)) * Atom::var(energy));
        let mut surface = vec![Atom::new(); self.upper_surface.len() + 1];
        surface[0] += &scale * &self.bulk;
        for (k, coefficient) in self.upper_surface.iter().enumerate() {
            surface[k] -= &scale * coefficient.derivative(energy);
            surface[k + 1] += &scale * coefficient;
        }
        Ok(Self {
            bulk: (-&scale * self.bulk.derivative(energy)).together().cancel(),
            upper_surface: surface.into_iter().map(|a| a.together().cancel()).collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn raised_residue_has_nonzero_surface_with_correct_sign() {
        let e = symbol!("fd_residue_e");
        let x = Atom::var(e);
        let raised = EnergyResidue::simple(e).raise(e, 1).unwrap();
        assert!(
            (raised.bulk + Atom::num(1) / (Atom::num(4) * x.clone().pow(3)))
                .together()
                .cancel()
                .is_zero()
        );
        assert!(
            (&raised.upper_surface[0] + Atom::num(1) / (Atom::num(4) * x.pow(2)))
                .together()
                .cancel()
                .is_zero()
        );
    }
    #[test]
    fn theta_zero_is_not_a_cut_zero_and_lower_endpoint_is_retained() {
        let g = parse!("fd_weighted_g");
        let e = parse!("fd_weighted_e");
        let measure = WeightedMeasure::new(
            vec![g.clone(), e.clone()],
            [
                g - parse!("fd_weighted_m"),
                e.clone(),
                parse!("fd_weighted_mu") - &e,
                e,
            ],
            [
                IndexRole::RequiredCut,
                IndexRole::Ordinary,
                IndexRole::Occupation,
                IndexRole::Occupation,
            ],
        )
        .unwrap();
        let indices = [
            symbol!("fd_n0"),
            symbol!("fd_n1"),
            symbol!("fd_n2"),
            symbol!("fd_n3"),
        ];
        let sources = measure
            .lorentz_ibps(1, &parse!("4-2*fd_eps"), &indices, 4)
            .unwrap();
        let bulk = sources
            .iter()
            .find(|s| s.id == "lorentz/0/1/occupation-case-0")
            .unwrap();
        assert_eq!(bulk.domain.bounds()[2], IndexBounds::fixed(0));
        assert!(
            bulk.terms
                .iter()
                .any(|term| term.shift[2] == 1 && term.coefficient == Atom::num(-1))
        );
        assert!(
            bulk.terms
                .iter()
                .any(|term| term.shift[3] == 1 && term.coefficient == Atom::num(1))
        );
        assert!(measure.lorentz_ibps(1, &parse!("4"), &indices, 3).is_err());
        let mut colliding = indices;
        colliding[0] = symbol!("fd_weighted_g");
        assert!(
            measure
                .lorentz_ibps(1, &parse!("4"), &colliding, 4)
                .is_err()
        );
        assert!(
            sources
                .iter()
                .filter(|s| s.id.starts_with("multiplication/2/"))
                .all(|s| !s.domain.bounds()[2].contains(0))
        );
    }
}
