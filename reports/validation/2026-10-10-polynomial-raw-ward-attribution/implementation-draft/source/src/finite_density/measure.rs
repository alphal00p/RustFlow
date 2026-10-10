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
    // Only actual physical factors are stored. N is backend storage capacity.
    factors: Vec<Atom>,
    roles: [IndexRole; N],
    numerator_basis: InversePropagatorBasis,
    numerator_slots: Vec<usize>,
}

impl<const N: usize> WeightedMeasure<N> {
    /// Exact-physical-arity compatibility constructor.
    pub fn new(coordinates: Vec<Atom>, factors: [Atom; N], roles: [IndexRole; N]) -> Result<Self> {
        Self::from_physical(coordinates, factors.into(), roles.into())
    }

    /// Physical factors are never padded; only emitted index/guard arrays use N.
    pub fn from_physical(
        coordinates: Vec<Atom>,
        factors: Vec<Atom>,
        physical_roles: Vec<IndexRole>,
    ) -> Result<Self> {
        let physical_arity = factors.len();
        if physical_arity == 0 || physical_arity > N || physical_roles.len() != physical_arity {
            return Err(Error::InvalidInput(
                "weighted physical arity exceeds storage or role count".into(),
            ));
        }
        let roles = std::array::from_fn(|i| {
            physical_roles
                .get(i)
                .copied()
                .unwrap_or(IndexRole::Ordinary)
        });
        let numerator_slots = (0..physical_arity)
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

    /// The complete physical factor list, without backend padding factors.
    pub fn factors(&self) -> &[Atom] {
        &self.factors
    }

    /// Number of physical inverse-propagator and distribution factors.
    pub fn physical_arity(&self) -> usize {
        self.factors.len()
    }

    /// Native guard coordinates beyond the physical measure are exactly zero.
    fn default_bounds(&self) -> [IndexBounds; N] {
        let mut bounds = *IndexDomain::for_roles(&self.roles).bounds();
        bounds[self.physical_arity()..].fill(IndexBounds::fixed(0));
        bounds
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
        &self,
        id: String,
        terms: BTreeMap<[i16; N], Atom>,
        bounds: [IndexBounds; N],
    ) -> Result<Option<GuardedIdentity<N>>> {
        if bounds[self.physical_arity()..]
            .iter()
            .any(|b| *b != IndexBounds::fixed(0))
            || terms
                .keys()
                .any(|shift| shift[self.physical_arity()..].iter().any(|&n| n != 0))
        {
            return Err(Error::InvalidInput(
                "weighted source escapes its physical coordinate image".into(),
            ));
        }
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

    /// Generate only the guarded polynomial raw Ward pair. There is no rational
    /// vector-field seed and no extension of the admitted ordinary index box.
    /// The typed evidence checks the complete exact measure before conversion.
    pub(crate) fn raw_singleton_ward_sources(
        &self,
        evidence: &super::massless_endpoint::SingletonRawWardEvidence,
        dimension: &Atom,
        indices: &[Symbol; N],
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        evidence.validate_measure(&self.coordinates, &self.factors, &self.roles)?;
        if domain_budget < 2 {
            return Err(Error::Limit(
                "raw Ward source domain budget exhausted".into(),
            ));
        }
        if indices
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != N
            || indices.iter().any(|&s| {
                !dimension.derivative(s).is_zero()
                    || self
                        .coordinates
                        .iter()
                        .chain(self.factors.iter())
                        .chain(self.numerator_basis.labels())
                        .any(|a| !a.derivative(s).is_zero())
            })
            || self.coordinates.iter().any(|c| {
                let AtomView::Var(v) = c.as_view() else {
                    return true;
                };
                !dimension.derivative(v.get_symbol()).is_zero()
            })
        {
            return Err(Error::InvalidInput(
                "raw Ward index/dimension symbol collision".into(),
            ));
        }
        // E is affine in the complete factor basis. Convert the actual Atom,
        // never infer one special factor slot or assume a unit coefficient.
        let mut energy_terms = BTreeMap::new();
        let mut upper_shift = [0_i16; N];
        upper_shift[evidence.upper()] = 1;
        self.multiply(evidence.energy(), upper_shift, &mut energy_terms)?;
        // An affine energy has at most (#coordinates + 1) factor monomials.
        // Two rows add one diagonal term each. Fail atomically on any excess.
        let per_row_terms = self
            .coordinates
            .len()
            .checked_add(2)
            .ok_or_else(|| Error::Limit("raw Ward term budget overflow".into()))?;
        let total_bound = per_row_terms
            .checked_mul(2)
            .ok_or_else(|| Error::Limit("raw Ward total term budget overflow".into()))?;
        if energy_terms
            .len()
            .checked_add(1)
            .is_none_or(|n| n > per_row_terms)
        {
            return Err(Error::Limit(
                "raw Ward affine energy conversion exceeded term bound".into(),
            ));
        }
        let mut common = self.default_bounds();
        for slot in evidence.physical_slots()..evidence.input_slots() {
            common[slot] =
                IndexBounds::new(None, Some(0)).map_err(|e| Error::InvalidInput(e.to_string()))?;
        }
        for &slot in evidence.energy_dependent_slots() {
            common[slot] = IndexBounds::fixed(0);
        }
        common[evidence.cut()] = IndexBounds::fixed(1);
        common[evidence.lower()] = IndexBounds::fixed(0);
        let mut result = Vec::with_capacity(2);
        let mut total_terms = 0_usize;
        for positive in [false, true] {
            let mut bounds = common;
            bounds[evidence.upper()] = if positive {
                IndexBounds::new(Some(1), None).map_err(|e| Error::InvalidInput(e.to_string()))?
            } else {
                IndexBounds::fixed(0)
            };
            let mut terms = BTreeMap::from([([0_i16; N], dimension - Atom::num(2))]);
            let coefficient = if positive {
                Atom::var(indices[evidence.upper()])
            } else {
                Atom::num(-1)
            };
            for (shift, value) in &energy_terms {
                *terms.entry(*shift).or_default() += &coefficient * value;
            }
            let source = self
                .source(
                    format!(
                        "raw-polynomial-singleton-Ward-v1/upper-{}",
                        if positive { "positive" } else { "zero" }
                    ),
                    terms,
                    bounds,
                )?
                .ok_or_else(|| {
                    Error::InvalidInput("raw Ward identity unexpectedly empty".into())
                })?;
            total_terms = total_terms
                .checked_add(source.terms.len())
                .ok_or_else(|| Error::Limit("raw Ward emitted term count overflow".into()))?;
            if total_terms > total_bound {
                return Err(Error::Limit("raw Ward total term budget exhausted".into()));
            }
            result.push(source);
        }
        Ok(result)
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
        self.ibp_partitioned(id, indices, direction, divergence, domain_budget, false)
    }

    fn ibp_partitioned(
        &self,
        id: &str,
        indices: &[Symbol; N],
        direction: &[Atom],
        divergence: &Atom,
        domain_budget: usize,
        split_inactive: bool,
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
        if direction
            .iter()
            .chain(std::iter::once(divergence))
            .any(|value| {
                indices[self.physical_arity()..]
                    .iter()
                    .any(|&symbol| !value.derivative(symbol).is_zero())
            })
        {
            return Err(Error::InvalidInput(
                "physical vector field depends on a storage-only index".into(),
            ));
        }
        // A factor annihilated by this vector field contributes no derivative
        // in any of its distribution branches. Keep that occupation axis at
        // its full admitted range rather than splitting an inactive theta face.
        let derivatives = self
            .factors
            .iter()
            .map(|factor| {
                self.coordinates
                    .iter()
                    .zip(direction)
                    .fold(Atom::zero(), |sum, (coordinate, action)| {
                        let AtomView::Var(variable) = coordinate.as_view() else {
                            unreachable!()
                        };
                        sum + factor.derivative(variable.get_symbol()) * action
                    })
                    .expand()
                    .together()
                    .cancel()
            })
            .collect::<Vec<_>>();
        let occupations = (0..self.physical_arity())
            .filter(|&i| {
                self.roles[i] == IndexRole::Occupation
                    && (split_inactive || !derivatives[i].is_zero())
            })
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
            let mut bounds = self.default_bounds();
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
            for slot in 0..self.physical_arity() {
                let derivative = &derivatives[slot];
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
                self.source(format!("{id}/occupation-case-{case}"), terms, bounds)?
            {
                sources.push(source);
            }
        }
        Ok(sources)
    }

    /// Exact polynomial IBPs tangent to occupied shells. The first
    /// `compact_loops` Gram momenta are real compact shell variables. Spatial
    /// angular fields also annihilate every occupation factor and are emitted
    /// first; the radial shell-tangent fields retain their explicit Fermi flux.
    /// These are direct total derivatives supplied to the same native guarded
    /// source owner, without changing its integral ordering or elimination.
    pub fn compact_tangent_ibps(
        &self,
        loops: usize,
        compact_loops: usize,
        dimension: &Atom,
        indices: &[Symbol; N],
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        if compact_loops > loops {
            return Err(Error::InvalidInput("compact tangent loop count".into()));
        }
        let mut sources = Vec::new();
        for angular in [true, false] {
            for differentiated in 0..compact_loops {
                for vector in 0..loops {
                    if angular && differentiated == vector {
                        continue;
                    }
                    let (direction, divergence) =
                        self.tangent_action(loops, dimension, differentiated, vector, angular)?;
                    let kind = if angular {
                        "angular-tangent"
                    } else {
                        "shell-tangent"
                    };
                    sources.extend(self.ibp(
                        &format!("{kind}/{differentiated}/{vector}"),
                        indices,
                        &direction,
                        &divergence,
                        domain_budget,
                    )?);
                }
            }
        }
        Ok(sources)
    }

    fn tangent_action(
        &self,
        loops: usize,
        dimension: &Atom,
        differentiated: usize,
        vector: usize,
        angular: bool,
    ) -> Result<(Vec<Atom>, Atom)> {
        if self.coordinates.len() != loops * (loops + 1) / 2 + loops
            || differentiated >= loops
            || vector >= loops
            || (angular && differentiated == vector)
        {
            return Err(Error::InvalidInput("compact tangent Gram geometry".into()));
        }
        let pairs = (0..loops)
            .flat_map(|a| (a..loops).map(move |b| (a, b)))
            .collect::<Vec<_>>();
        let scalar = |a: usize, b: usize| {
            self.coordinates[pairs
                .iter()
                .position(|&(i, j)| (i, j) == (a.min(b), a.max(b)))
                .unwrap()]
            .clone()
        };
        let energy = |a: usize| self.coordinates[pairs.len() + a].clone();
        let i = differentiated;
        let j = vector;
        let ei = energy(i);
        let ej = energy(j);
        let gij = scalar(i, j);
        let flux = &ei * &ej - &gij;
        // u²=1. In the angular case, v=(Ei²-gii)(qj-Ej u)
        // -(Ei Ej-gij)(qi-Ei u); otherwise v=Ei qj-gij u.
        let contraction = |a| {
            if angular {
                (ei.pow(2) - scalar(i, i)) * (scalar(a, j) - &ej * energy(a))
                    - &flux * (scalar(a, i) - &ei * energy(a))
            } else {
                &ei * scalar(a, j) - &gij * energy(a)
            }
            .expand()
        };
        let mut direction = pairs
            .iter()
            .map(|&(a, b)| {
                (Atom::num(i32::from(a == i)) * contraction(b)
                    + Atom::num(i32::from(b == i)) * contraction(a))
                .expand()
            })
            .collect::<Vec<_>>();
        for a in 0..loops {
            direction.push(if !angular && a == i {
                flux.clone().expand()
            } else {
                Atom::zero()
            });
        }
        let divergence = if angular {
            (Atom::num(2) - dimension) * flux
        } else if i == j {
            (dimension - Atom::one()) * ei
        } else {
            Atom::zero()
        };
        Ok((direction, divergence.expand()))
    }

    /// Shell-normal total derivatives v_i=u/(2E_i), lowered through a
    /// certified completion f_slot=c E_i. The caller must prove E_i nonzero on
    /// the complete distribution support and admit all integer powers there.
    /// This is the ordinary u identity at a_slot+1, multiplied by c/2; it does
    /// not introduce a second elimination or differentiate an integrated shell.
    pub fn compact_normal_ibps(
        &self,
        loops: usize,
        completions: &[(usize, usize, Rational)],
        indices: &[Symbol; N],
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        let pairs = (0..loops)
            .flat_map(|a| (a..loops).map(move |b| (a, b)))
            .collect::<Vec<_>>();
        if self.coordinates.len() != pairs.len() + loops {
            return Err(Error::InvalidInput("shell-normal Gram geometry".into()));
        }
        let mut sources = Vec::new();
        let mut used = std::collections::BTreeSet::new();
        for &(i, slot, ref coefficient) in completions {
            if i >= loops
                || slot >= self.physical_arity()
                || self.roles[slot] != IndexRole::Ordinary
                || coefficient.is_zero()
                || !used.insert((i, slot))
            {
                return Err(Error::InvalidInput(
                    "shell-normal completion certificate".into(),
                ));
            }
            let energy = &self.coordinates[pairs.len() + i];
            if !(&self.factors[slot] - Atom::num(coefficient.clone()) * energy)
                .expand()
                .together()
                .cancel()
                .is_zero()
            {
                return Err(Error::InvalidInput(
                    "shell-normal factor differs from certified energy".into(),
                ));
            }
            let mut direction = pairs
                .iter()
                .map(|&(a, b)| {
                    Atom::num(i32::from(a == i)) * &self.coordinates[pairs.len() + b]
                        + Atom::num(i32::from(b == i)) * &self.coordinates[pairs.len() + a]
                })
                .collect::<Vec<_>>();
            direction.extend((0..loops).map(|a| Atom::num(i32::from(a == i))));
            let mut shifted_sources = self.ibp(
                &format!("shell-normal/{i}/{slot}"),
                indices,
                &direction,
                &Atom::zero(),
                domain_budget,
            )?;
            let index = Atom::var(indices[slot]);
            let factor = Atom::num(coefficient.clone()) / Atom::num(2);
            for source in &mut shifted_sources {
                for term in &mut source.terms {
                    term.shift[slot] = term.shift[slot]
                        .checked_add(1)
                        .ok_or_else(|| Error::Limit("shell-normal shift overflow".into()))?;
                    term.coefficient = (&factor
                        * term
                            .coefficient
                            .replace(index.clone())
                            .with(&index + Atom::one()))
                    .together()
                    .cancel();
                }
                let mut shift = [0; N];
                shift[slot] = 1;
                source.domain = source
                    .domain
                    .pullback(&shift)
                    .map_err(|error| Error::InvalidInput(error.to_string()))?;
            }
            sources.extend(shifted_sources);
        }
        Ok(sources)
    }

    /// Exact polynomial/distribution multiplication identities. These include
    /// h*delta(h)=0 and h*C_n(h)=C_(n-1)(h), n>=2; h*theta(h) is not zero.
    pub fn multiplication_sources(&self) -> Result<Vec<GuardedIdentity<N>>> {
        let mut sources = Vec::new();
        for slot in 0..self.physical_arity() {
            let mut cases = Vec::new();
            let mut default_bounds = self.default_bounds();
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
                    self.source(format!("multiplication/{slot}/{case}"), terms, bounds)?
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
        self.lorentz_ibps_partitioned(loops, dimension, indices, domain_budget, false)
    }

    /// Equivalent Lorentz identities with every occupation face partitioned.
    /// Native discovery may cover different domains for equivalent source
    /// presentations, so the validated legacy presentation remains available.
    pub(crate) fn lorentz_ibps_legacy(
        &self,
        loops: usize,
        dimension: &Atom,
        indices: &[Symbol; N],
        domain_budget: usize,
    ) -> Result<Vec<GuardedIdentity<N>>> {
        self.lorentz_ibps_partitioned(loops, dimension, indices, domain_budget, true)
    }

    fn lorentz_ibps_partitioned(
        &self,
        loops: usize,
        dimension: &Atom,
        indices: &[Symbol; N],
        domain_budget: usize,
        split_inactive: bool,
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
                sources.extend(self.ibp_partitioned(
                    &format!("lorentz/{differentiated}/{vector}"),
                    indices,
                    &direction,
                    &divergence,
                    domain_budget,
                    split_inactive,
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

    fn tangent_measure() -> WeightedMeasure<9> {
        let coordinates = (0..5)
            .map(|i| Atom::var(symbol!(format!("fd_tangent::c_{i}"))))
            .collect::<Vec<_>>();
        let c = &coordinates;
        WeightedMeasure::new(
            coordinates.clone(),
            [
                &c[0] - Atom::num(2),
                &c[2] - Atom::num(3),
                &c[0] + &c[2] - Atom::num(2) * &c[1] - Atom::num(5),
                c[3].clone(),
                c[4].clone(),
                Atom::num(7) - &c[3],
                c[3].clone(),
                Atom::num(11) - &c[4],
                c[4].clone(),
            ],
            [
                IndexRole::RequiredCut,
                IndexRole::RequiredCut,
                IndexRole::Ordinary,
                IndexRole::Ordinary,
                IndexRole::Ordinary,
                IndexRole::Occupation,
                IndexRole::Occupation,
                IndexRole::Occupation,
                IndexRole::Occupation,
            ],
        )
        .unwrap()
    }

    #[test]
    fn tangent_gram_actions_and_divergences_match_independent_cartesian_fields() {
        let measure = tangent_measure();
        for dimension in 2..=5 {
            let variables = (0..2)
                .map(|i| {
                    (0..dimension)
                        .map(|mu| symbol!(format!("fd_tangent_cartesian::q_{dimension}_{i}_{mu}")))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let q = variables
                .iter()
                .map(|row| row.iter().map(|&s| Atom::var(s)).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let dot = |i: usize, j: usize| {
                (1..dimension).fold(&q[i][0] * &q[j][0], |s, mu| s - &q[i][mu] * &q[j][mu])
            };
            let gram = [
                dot(0, 0),
                dot(0, 1),
                dot(1, 1),
                q[0][0].clone(),
                q[1][0].clone(),
            ];
            let substitution = measure
                .coordinates
                .iter()
                .cloned()
                .zip(gram.iter().cloned())
                .collect::<BTreeMap<_, _>>();
            for angular in [false, true] {
                for i in 0..2 {
                    for j in 0..2 {
                        if angular && i == j {
                            continue;
                        }
                        let (action, divergence) = measure
                            .tangent_action(2, &Atom::num(dimension as i64), i, j, angular)
                            .unwrap();
                        let ri2 = q[i][0].pow(2) - dot(i, i);
                        let flux = &q[i][0] * &q[j][0] - dot(i, j);
                        let cartesian = (0..dimension)
                            .map(|mu| {
                                let time = Atom::num(i32::from(mu == 0));
                                if angular {
                                    &ri2 * (&q[j][mu] - &q[j][0] * &time)
                                        - &flux * (&q[i][mu] - &q[i][0] * &time)
                                } else {
                                    &q[i][0] * &q[j][mu] - dot(i, j) * time
                                }
                            })
                            .collect::<Vec<_>>();
                        let shell_action = (1..dimension)
                            .fold(&q[i][0] * &cartesian[0], |s, mu| {
                                s - &q[i][mu] * &cartesian[mu]
                            });
                        assert!(shell_action.expand().is_zero());
                        if angular {
                            assert!(cartesian[0].expand().is_zero());
                        }
                        let cartesian_divergence = (0..dimension).fold(Atom::zero(), |s, mu| {
                            s + cartesian[mu].derivative(variables[i][mu])
                        });
                        assert!(
                            (cartesian_divergence
                                - crate::family::substitute(&divergence, &substitution))
                            .expand()
                            .is_zero()
                        );
                        for (coordinate, abstract_action) in gram.iter().zip(action) {
                            let cartesian_action = (0..dimension).fold(Atom::zero(), |s, mu| {
                                s + coordinate.derivative(variables[i][mu]) * &cartesian[mu]
                            });
                            assert!(
                                (cartesian_action
                                    - crate::family::substitute(&abstract_action, &substitution))
                                .expand()
                                .is_zero()
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn inactive_occupations_are_unsplit_and_tangent_sources_do_not_raise_shells() {
        let measure = tangent_measure();
        let indices = std::array::from_fn(|i| symbol!(format!("fd_tangent::a_{i}")));
        let sources = measure
            .compact_tangent_ibps(2, 2, &Atom::num(4), &indices, 4)
            .unwrap();
        let angular = sources
            .iter()
            .filter(|s| s.id.starts_with("angular-tangent/"))
            .collect::<Vec<_>>();
        assert_eq!(angular.len(), 2);
        for source in &sources {
            assert!(
                source
                    .terms
                    .iter()
                    .all(|term| term.shift[0] <= 0 && term.shift[1] <= 0)
            );
        }
        for source in angular {
            assert!(
                source.domain.bounds()[5..]
                    .iter()
                    .all(|b| *b == IndexBounds::new(Some(0), None).unwrap())
            );
            assert!(
                source
                    .terms
                    .iter()
                    .all(|term| term.shift[5..].iter().all(|&s| s == 0))
            );
        }
        // Each ordinary Lorentz field acts on only one loop's two endpoints,
        // so four cases suffice even though there are four occupations total.
        let lorentz = measure.lorentz_ibps(2, &Atom::num(4), &indices, 4).unwrap();
        assert!(
            lorentz
                .iter()
                .filter(|s| s.id.starts_with("lorentz/0/"))
                .all(|s| s.domain.bounds()[7..]
                    .iter()
                    .all(|b| *b == IndexBounds::new(Some(0), None).unwrap()))
        );
    }

    #[test]
    fn normalized_shell_sources_match_direct_distribution_product_rule() {
        let g = parse!("fd_normal::g");
        let e = parse!("fd_normal::e");
        let indices = std::array::from_fn(|i| symbol!(format!("fd_normal::a_{i}")));
        let measure = WeightedMeasure::new(
            vec![g.clone(), e.clone()],
            [
                g - Atom::num(2),
                -Atom::num(2) * &e,
                Atom::num(3) - &e,
                e.clone(),
            ],
            [
                IndexRole::RequiredCut,
                IndexRole::Ordinary,
                IndexRole::Occupation,
                IndexRole::Occupation,
            ],
        )
        .unwrap();
        let sources = measure
            .compact_normal_ibps(1, &[(0, 1, Rational::from(-2))], &indices, 4)
            .unwrap();
        let bulk = sources
            .iter()
            .find(|s| {
                s.domain.bounds()[2] == IndexBounds::fixed(0)
                    && s.domain.bounds()[3] == IndexBounds::fixed(0)
            })
            .unwrap();
        let terms = bulk
            .terms
            .iter()
            .map(|t| (t.shift, t.coefficient.clone()))
            .collect::<BTreeMap<_, _>>();
        let expected = BTreeMap::from([
            ([1, 0, 0, 0], -Atom::var(indices[0])),
            (
                [0, 2, 0, 0],
                -Atom::num(2) * (Atom::var(indices[1]) + Atom::one()),
            ),
            ([0, 1, 1, 0], Atom::one()),
            ([0, 1, 0, 1], -Atom::one()),
        ]);
        assert_eq!(terms.len(), expected.len());
        for (shift, coefficient) in expected {
            assert!((&terms[&shift] - coefficient).expand().is_zero());
        }
        // Independent Cartesian divergence, shell action and temporal flux of
        // v=u/(2E), in several integer spacetime dimensions.
        for dimension in 2..=5 {
            let q = (0..dimension)
                .map(|mu| symbol!(format!("fd_normal_cartesian::q_{dimension}_{mu}")))
                .collect::<Vec<_>>();
            let energy = Atom::var(q[0]);
            let shell = (1..dimension).fold(energy.pow(2) - Atom::num(2), |s, mu| {
                s - Atom::var(q[mu]).pow(2)
            });
            let field = Atom::one() / (Atom::num(2) * &energy);
            assert!(
                (field.derivative(q[0]) + Atom::one() / (Atom::num(2) * energy.pow(2)))
                    .together()
                    .cancel()
                    .is_zero()
            );
            assert!(
                (shell.derivative(q[0]) * &field - Atom::one())
                    .together()
                    .cancel()
                    .is_zero()
            );
            assert!(
                ((Atom::num(3) - &energy).derivative(q[0]) * &field + &field)
                    .expand()
                    .is_zero()
            );
        }
        assert!(
            measure
                .compact_normal_ibps(1, &[(0, 1, Rational::from(1))], &indices, 4)
                .is_err()
        );
    }

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

    #[test]
    fn storage_capacity_preserves_every_physical_source_presentation() {
        let g = parse!("fd_storage::g");
        let e = parse!("fd_storage::e");
        let coordinates = vec![g.clone(), e.clone()];
        let factors = [g - Atom::num(2), e.clone(), Atom::num(3) - &e, e];
        let roles = [
            IndexRole::RequiredCut,
            IndexRole::Ordinary,
            IndexRole::Occupation,
            IndexRole::Occupation,
        ];
        let exact = WeightedMeasure::new(coordinates.clone(), factors.clone(), roles).unwrap();
        let padded =
            WeightedMeasure::<6>::from_physical(coordinates, factors.to_vec(), roles.to_vec())
                .unwrap();
        assert_eq!(exact.factors(), padded.factors());
        assert_eq!(padded.physical_arity(), 4);
        assert_eq!(padded.roles()[4..], [IndexRole::Ordinary; 2]);
        let indices = std::array::from_fn(|i| symbol!(format!("fd_storage::a_{i}")));
        let exact_indices = std::array::from_fn(|i| indices[i]);
        let dimension = parse!("4-2*fd_storage::epsilon");
        let exact_sources = [
            exact.multiplication_sources().unwrap(),
            exact
                .lorentz_ibps(1, &dimension, &exact_indices, 4)
                .unwrap(),
            exact
                .lorentz_ibps_legacy(1, &dimension, &exact_indices, 4)
                .unwrap(),
            exact
                .compact_tangent_ibps(1, 1, &dimension, &exact_indices, 4)
                .unwrap(),
            exact
                .compact_normal_ibps(1, &[(0, 1, Rational::one())], &exact_indices, 4)
                .unwrap(),
        ];
        let padded_sources = [
            padded.multiplication_sources().unwrap(),
            padded.lorentz_ibps(1, &dimension, &indices, 4).unwrap(),
            padded
                .lorentz_ibps_legacy(1, &dimension, &indices, 4)
                .unwrap(),
            padded
                .compact_tangent_ibps(1, 1, &dimension, &indices, 4)
                .unwrap(),
            padded
                .compact_normal_ibps(1, &[(0, 1, Rational::one())], &indices, 4)
                .unwrap(),
        ];
        for (physical_sources, stored_sources) in exact_sources.into_iter().zip(padded_sources) {
            assert!(!physical_sources.is_empty());
            assert_eq!(physical_sources.len(), stored_sources.len());
            for (physical, stored) in physical_sources.into_iter().zip(stored_sources) {
                assert_eq!(physical.id, stored.id);
                assert_eq!(physical.domain.bounds(), &stored.domain.bounds()[..4]);
                assert!(
                    stored.domain.bounds()[4..]
                        .iter()
                        .all(|b| *b == IndexBounds::fixed(0))
                );
                assert_eq!(physical.nonzero_conditions, stored.nonzero_conditions);
                assert_eq!(physical.terms.len(), stored.terms.len());
                for (physical_term, stored_term) in physical.terms.into_iter().zip(stored.terms) {
                    assert_eq!(physical_term.shift, stored_term.shift[..4]);
                    assert_eq!(stored_term.shift[4..], [0, 0]);
                    assert!(
                        (physical_term.coefficient - &stored_term.coefficient)
                            .together()
                            .cancel()
                            .is_zero()
                    );
                    assert!(
                        indices[4..]
                            .iter()
                            .all(|&index| stored_term.coefficient.derivative(index).is_zero())
                    );
                }
            }
        }
        // A source author cannot smuggle a storage index into a physical vector field.
        assert!(
            padded
                .ibp(
                    "invalid-tail-field",
                    &indices,
                    &[Atom::var(indices[4]), Atom::zero()],
                    &Atom::zero(),
                    4
                )
                .is_err()
        );
        assert!(
            padded
                .ibp(
                    "invalid-tail-divergence",
                    &indices,
                    &[Atom::zero(), Atom::zero()],
                    &Atom::var(indices[5]),
                    4
                )
                .is_err()
        );
        assert!(
            padded
                .compact_normal_ibps(1, &[(0, 4, Rational::one())], &indices, 4)
                .is_err()
        );
    }

    #[test]
    fn physical_factor_constructor_rejects_invalid_capacity_and_roles() {
        let x = parse!("fd_storage_bad::x");
        assert!(
            WeightedMeasure::<1>::from_physical(
                vec![x.clone()],
                vec![x.clone(), x.clone()],
                vec![IndexRole::Ordinary; 2]
            )
            .is_err()
        );
        assert!(
            WeightedMeasure::<2>::from_physical(vec![x.clone()], vec![x.clone()], vec![]).is_err()
        );
        assert!(WeightedMeasure::<2>::from_physical(vec![x], vec![], vec![]).is_err());
    }
}
