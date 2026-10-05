//! Exact common-mass AMF connections from supplied canonical physical systems.
//!
//! The supplied connection is a scientific input, just like a reduction table.
//! Its association with the declared physical integrals is the caller's contract.
//! This module certifies homogeneity, the gauge identity, root cancellation, and
//! the exact deformation; it does not supply numerical boundary constants.
use crate::algebraic::{CanonicalAlgebraicSystem, normalized_root_entry};
use crate::engine::SuppliedAuxiliarySystem;
use crate::family::{imaginary_parameter, substitute};
use crate::kinematics::KinematicPath;
use crate::reduction::ReducedSystem;
use crate::{
    DifferentialSystem, Error, Integral, IntegralFamily, KinematicPoint, Progress, Result,
    RunContext,
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

/// An ordinary integral basis related to a supplied canonical system by I = M F,
/// up to a common normalization independent of every physical coordinate.
///
/// The family has unit nonzero masses. Every canonical coordinate has mass
/// dimension two; propagator coefficients are kinematics independent and external Gram
/// products have degree one in these coordinates. All physical massive lines
/// receive the same auxiliary mass; irreducible numerator slots remain massless.
pub struct HomogeneousCanonicalBasis<'a> {
    pub family: &'a IntegralFamily,
    pub integrals: &'a [Integral],
    pub canonical: &'a CanonicalAlgebraicSystem,
    pub ordinary_from_canonical: &'a [Vec<Atom>],
    pub canonical_from_ordinary: &'a [Vec<Atom>],
    /// Provenance of the supplied physical connection and exact basis map.
    pub provenance: &'a str,
}

impl HomogeneousCanonicalBasis<'_> {
    fn validate(&self, variable: Symbol) -> Result<Vec<bool>> {
        self.family.validate()?;
        let n = self.canonical.dimension();
        if self.provenance.trim().is_empty()
            || self.family.epsilon != self.canonical.epsilon()
            || variable == self.family.epsilon
            || variable == imaginary_parameter()
            || self.canonical.variables().contains(&variable)
            || self.canonical.roots().iter().any(|r| r.symbol == variable)
            || self.integrals.len() != n
            || self.integrals.iter().collect::<BTreeSet<_>>().len() != n
            || [self.ordinary_from_canonical, self.canonical_from_ordinary]
                .iter()
                .any(|m| m.len() != n || m.iter().any(|row| row.len() != n))
        {
            return Err(Error::InvalidInput(
                "common-mass basis has incompatible dimensions, symbols, or repeated integrals"
                    .into(),
            ));
        }
        for integral in self.integrals {
            self.family.validate_integral(integral)?;
        }
        let allowed = self
            .canonical
            .variables()
            .iter()
            .copied()
            .chain([self.family.epsilon, imaginary_parameter()])
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        for expression in self.family.external_gram.iter().flatten().chain(
            self.family
                .propagators
                .iter()
                .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
        ) {
            let mut found = BTreeSet::new();
            crate::family::scalar_symbols(expression.as_view(), &mut found)?;
            if !found.is_subset(&allowed) {
                return Err(Error::InvalidInput(
                    "common-mass family contains undeclared physical parameters".into(),
                ));
            }
        }
        let homogeneous = |expression: &Atom, degree: i64| {
            let derivative = self
                .canonical
                .variables()
                .iter()
                .fold(Atom::Zero, |sum, &v| {
                    sum + Atom::var(v) * expression.derivative(v)
                });
            (derivative - expression * degree)
                .together()
                .cancel()
                .is_zero()
        };
        if self
            .family
            .external_gram
            .iter()
            .flatten()
            .any(|a| !homogeneous(a, 1))
        {
            return Err(Error::Unsupported(
                "common-mass pullback requires degree-one external scalar products".into(),
            ));
        }
        let mut mask = Vec::with_capacity(self.family.propagators.len());
        for (index, propagator) in self.family.propagators.iter().enumerate() {
            // The existing region owner rejects linear and non-rank-one loop
            // quadratics; do not silently apply quadratic homogeneity to them.
            let branch = crate::regions::branch(propagator, self.family.loops.len())?;
            let pivot = branch.iter().position(|a| !a.is_zero()).ok_or_else(|| {
                Error::Unsupported("common-mass pullback requires a quadratic loop branch".into())
            })?;
            let offset = self.family.loops.len() * (self.family.loops.len() + 1) / 2;
            for (i, coefficient) in branch.iter().enumerate() {
                for a in 0..self.family.external.len() {
                    if !(&propagator.scalar_products[offset + i * self.family.external.len() + a]
                        - coefficient
                            * &propagator.scalar_products
                                [offset + pivot * self.family.external.len() + a]
                            / &branch[pivot])
                        .together()
                        .cancel()
                        .is_zero()
                    {
                        return Err(Error::Unsupported(
                            "common-mass external routing must share the quadratic loop branch"
                                .into(),
                        ));
                    }
                }
            }
            let mass = self.family.mass_squared(index)?;
            if (!mass.is_zero() && !mass.is_one())
                || (index >= self.family.physical_propagators && !mass.is_zero())
                || !homogeneous(&(&propagator.constant + &mass), 1)
                || propagator.scalar_products.iter().any(|a| {
                    self.canonical
                        .variables()
                        .iter()
                        .any(|&v| !a.derivative(v).is_zero())
                })
            {
                return Err(Error::Unsupported(
                    "common-mass pullback requires homogeneous quadratic propagators with unit physical masses and massless numerator slots".into(),
                ));
            }
            mask.push(index < self.family.physical_propagators && !mass.is_zero());
        }
        if !mask.iter().any(|&shifted| shifted) {
            return Err(Error::Unsupported(
                "common-mass pullback needs a massive physical line".into(),
            ));
        }
        Ok(mask)
    }

    /// Derive dI/deta at an exact point with z=1+eta and x(eta)=x(0)/z.
    ///
    /// If h_i=L*(D0/2-eps)-sum(a_i), then I_i(eta)=z^h_i M_ij(x(eta),eps) F_j.
    /// The common fractional epsilon power cancels from the gauge conjugation;
    /// h_i-h_j is always integral. The canonical connection is exactly linear in
    /// epsilon, so its order-one pullback retains its full regulator dependence.
    pub fn pullback(
        &self,
        point: &KinematicPoint,
        variable: Symbol,
        context: &RunContext,
    ) -> Result<SuppliedAuxiliarySystem> {
        context.cancellation.check()?;
        let deformation_mask = self.validate(variable)?;
        let epsilon = self.family.epsilon;
        let z = Atom::one() + Atom::var(variable);
        let variables = self
            .canonical
            .variables()
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        if point.0.keys().cloned().collect::<BTreeSet<_>>()
            != variables
                .iter()
                .copied()
                .map(Atom::var)
                .collect::<BTreeSet<_>>()
        {
            return Err(Error::InvalidInput(
                "common-mass point must specify every canonical coordinate exactly once".into(),
            ));
        }
        let path = KinematicPath {
            parameter: variable,
            coordinates: variables
                .iter()
                .map(|&v| (v, &point.0[&Atom::var(v)] / &z))
                .collect(),
        };
        context.emit(Progress::Stage {
            name: "pulling the canonical connection along the common auxiliary mass".into(),
        })?;
        let pulled = self.canonical.pullback(&path, 1)?;
        let roots = &pulled.roots;
        let substitutions = path
            .coordinates
            .iter()
            .map(|(&v, a)| (Atom::var(v), a.clone()))
            .collect();
        let mut physical = variables;
        physical.insert(epsilon);
        physical.insert(imaginary_parameter());
        let expressions = self
            .ordinary_from_canonical
            .iter()
            .flatten()
            .chain(self.canonical_from_ordinary.iter().flatten())
            .cloned()
            .collect::<Vec<_>>();
        let mut conditions = pulled.nonzero_conditions;
        let mut source_conditions = crate::algebraic::registered_domain_conditions(
            &expressions,
            self.canonical.roots(),
            &physical,
        )?;
        source_conditions.extend(crate::physical_conditions::rational_denominator_conditions(
            &self
                .family
                .external_gram
                .iter()
                .flatten()
                .chain(
                    self.family
                        .propagators
                        .iter()
                        .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products)),
                )
                .cloned()
                .collect::<Vec<_>>(),
            &physical,
        )?);
        for condition in source_conditions {
            let leading =
                crate::physical_conditions::epsilon_leading_coefficient(&condition, epsilon)?;
            if substitute(&leading, &substitutions)
                .together()
                .cancel()
                .is_zero()
            {
                return Err(Error::InvalidInput(
                    "common-mass path changes a gauge condition's generic epsilon valuation".into(),
                ));
            }
            conditions.push(substitute(&condition, &substitutions));
        }
        conditions.push(z.clone());
        let allowed = [variable, epsilon, imaginary_parameter()]
            .into_iter()
            .chain(roots.iter().map(|r| r.symbol))
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        let normalize = |entry: &Atom| -> Result<Atom> {
            let normalized = normalized_root_entry(entry, roots, &allowed)?;
            Ok(normalized.expression(roots))
        };
        let restrict_matrix = |matrix: &[Vec<Atom>]| -> Result<Vec<Vec<Atom>>> {
            matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|a| normalize(&substitute(a, &substitutions)))
                        .collect()
                })
                .collect()
        };
        let forward = restrict_matrix(self.ordinary_from_canonical)?;
        let inverse = restrict_matrix(self.canonical_from_ordinary)?;
        context.emit(Progress::Stage {
            name: "certifying the exact physical/canonical inverse".into(),
        })?;
        for (i, row) in crate::algebra::matmul(&forward, &inverse)
            .iter()
            .enumerate()
        {
            context.cancellation.check()?;
            for (j, entry) in row.iter().enumerate() {
                if !normalize(&(entry - Atom::num(i64::from(i == j))))?.is_zero() {
                    return Err(Error::InvalidInput(
                        "common-mass physical/canonical inverse identity failed".into(),
                    ));
                }
            }
        }
        let n = self.integrals.len();
        let canonical = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        &pulled.system.matrices[0][i][j]
                            + Atom::var(epsilon) * &pulled.system.matrices[1][i][j]
                    })
                    .collect()
            })
            .collect::<Vec<Vec<Atom>>>();
        let mut differentiated = forward
            .iter()
            .map(|row| {
                row.iter()
                    .map(|entry| {
                        let mut derivative = entry.derivative(variable);
                        for root in roots {
                            derivative += entry.derivative(root.symbol)
                                * root.radicand.derivative(variable)
                                / (Atom::var(root.symbol) * 2);
                        }
                        normalize(&derivative)
                    })
                    .collect()
            })
            .collect::<Result<Vec<Vec<Atom>>>>()?;
        context.emit(Progress::Stage {
            name: "transforming the common-mass differential connection to ordinary integrals"
                .into(),
        })?;
        let coupled = crate::algebra::matmul(&forward, &canonical);
        for (row, source) in differentiated.iter_mut().zip(coupled) {
            for (entry, coupling) in row.iter_mut().zip(source) {
                *entry = normalize(&(&*entry + coupling))?;
            }
        }
        context.cancellation.check()?;
        let mut unscaled = crate::algebra::matmul(&differentiated, &inverse);
        for (i, row) in unscaled.iter_mut().enumerate() {
            context.cancellation.check()?;
            for (j, entry) in row.iter_mut().enumerate() {
                let normalized = normalized_root_entry(entry, roots, &allowed)?;
                if normalized.terms.keys().any(|indices| !indices.is_empty()) {
                    return Err(Error::Unsupported(format!(
                        "ordinary common-mass connection retains algebraic roots at ({i},{j}) after exact quotient reduction: {}",
                        normalized.expression(roots),
                    )));
                }
                for coefficient in normalized.decoded_inverse_coefficients() {
                    conditions.extend(crate::physical_conditions::rational_denominator_conditions(
                        &[coefficient],
                        &BTreeSet::from([variable, epsilon, imaginary_parameter()]),
                    )?);
                }
                *entry = normalized.expression(roots);
            }
        }
        // Recheck the final rational gauge, without its proposed inverse.
        context.emit(Progress::Stage {
            name: "certifying root cancellation and the ordinary differential gauge".into(),
        })?;
        for (row, expected) in crate::algebra::matmul(&unscaled, &forward)
            .iter()
            .zip(&differentiated)
        {
            context.cancellation.check()?;
            for (entry, expected) in row.iter().zip(expected) {
                if !normalize(&(entry - expected))?.is_zero() {
                    return Err(Error::InvalidInput(
                        "common-mass differential gauge identity failed".into(),
                    ));
                }
            }
        }
        let sums = self
            .integrals
            .iter()
            .map(|integral| integral.0.iter().map(|&a| i64::from(a)).sum::<i64>())
            .collect::<Vec<_>>();
        let loops = i64::try_from(self.family.loops.len())
            .map_err(|_| Error::Limit("loop count overflows homogeneity".into()))?;
        let dimension = loops
            .checked_mul(self.family.dimension)
            .ok_or_else(|| Error::Limit("dimension overflows homogeneity".into()))?;
        for (i, row) in unscaled.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry *= z.pow(sums[j] - sums[i]);
                if i == j {
                    *entry += (Atom::num(Rational::from((dimension, 2)))
                        - Atom::num(sums[i])
                        - Atom::var(epsilon) * loops)
                        / &z;
                }
                *entry = entry.together().cancel();
            }
        }
        conditions.extend(crate::physical_conditions::rational_denominator_conditions(
            &unscaled.iter().flatten().cloned().collect::<Vec<_>>(),
            &BTreeSet::from([variable, epsilon, imaginary_parameter()]),
        )?);
        let conditions = crate::physical_conditions::canonical_conditions(
            &conditions,
            &BTreeSet::from([variable, epsilon, imaginary_parameter()]),
        )?;
        crate::physical_conditions::validate_conditions_at(
            &conditions,
            epsilon,
            &BTreeMap::from([(variable, Atom::Zero)]),
        )?;
        let system = DifferentialSystem {
            variable,
            matrix: unscaled,
        };
        system.validate()?;
        let targets = self
            .integrals
            .iter()
            .map(|integral| BTreeMap::from([(integral.clone(), Atom::one())]))
            .collect::<Vec<_>>();
        let reduced = ReducedSystem {
            basis: self.integrals.to_vec(),
            matrix: system.matrix.clone(),
            candidates: self
                .integrals
                .iter()
                .cloned()
                .zip(targets.iter().cloned())
                .collect(),
            targets,
            nonzero_conditions: conditions,
            transformations: Vec::new(),
        };
        Ok(SuppliedAuxiliarySystem {
            variable,
            reduced,
            deformation_mask,
            provenance: format!("common-mass-v1; {}", self.provenance),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gg_hg::{HiggsJetIntegralSystem, PluginFamilyKind};

    fn complete_published_connection(kind: PluginFamilyKind) -> Result<()> {
        use std::{sync::Arc, time::Instant};
        let started = Instant::now();
        let input = HiggsJetIntegralSystem::load(kind, "common_mass_complete_certificate")?;
        let family = input.map.family()?;
        let n = kind.dimension();
        let dense = |rows: &[Vec<crate::gg_hg::CanonicalTerm>]| {
            let mut matrix = vec![vec![Atom::Zero; n]; n];
            for (i, row) in rows.iter().enumerate() {
                for term in row {
                    matrix[i][term.integral] = term.coefficient.clone();
                }
            }
            matrix
        };
        let forward = dense(&input.map.canonical_to_publication);
        let inverse = dense(&input.map.publication_to_canonical);
        let values = if kind == PluginFamilyKind::Planar {
            [(1, 10), (1, 25), (1, 50)]
        } else {
            [(1, 10), (1, 5), (1, 1)]
        };
        let point = KinematicPoint(
            input
                .map
                .coordinates
                .into_iter()
                .map(Atom::var)
                .zip(values.map(|(a, b)| Atom::num(Rational::from((-a, b)))))
                .collect(),
        );
        let context = RunContext {
            progress: Some(Arc::new(move |progress| {
                eprintln!("{:?} {progress:?}", started.elapsed())
            })),
            ..Default::default()
        };
        let basis = HomogeneousCanonicalBasis {
            family: &family,
            integrals: &input.map.integrals,
            canonical: &input.canonical,
            ordinary_from_canonical: &forward,
            canonical_from_ordinary: &inverse,
            provenance: &input.map.evidence.source_sha256,
        };
        let eta = symbol!("common_mass_complete_certificate::eta");
        let connection = basis.pullback(&point, eta, &context)?;
        let nonzero = connection
            .reduced
            .matrix
            .iter()
            .flatten()
            .filter(|a| !a.is_zero())
            .count();
        eprintln!(
            "CERTIFIED {kind:?}: dimension={n}, nonzero={nonzero}, conditions={}, elapsed={:?}",
            connection.reduced.nonzero_conditions.len(),
            started.elapsed()
        );
        assert_eq!(connection.reduced.basis, input.map.integrals);
        assert_eq!(
            connection.deformation_mask,
            vec![false, false, false, false, false, true, true, false, false]
        );
        if kind == PluginFamilyKind::Planar {
            assert!(
                (&connection.reduced.matrix[0][0]
                    + Atom::var(input.map.epsilon) / (Atom::one() + Atom::var(eta)))
                .together()
                .cancel()
                .is_zero()
            );
            assert!(
                connection.reduced.matrix[0][1..]
                    .iter()
                    .all(|a| a.is_zero())
            );
        }
        Ok(())
    }

    #[test]
    #[ignore = "complete 48-master symbolic connection certificate"]
    fn complete_planar_connection_is_rational_and_gauge_certified() -> Result<()> {
        complete_published_connection(PluginFamilyKind::Planar)
    }

    #[test]
    #[ignore = "complete 61-master symbolic connection certificate"]
    fn complete_nonplanar_connection_is_rational_and_gauge_certified() -> Result<()> {
        complete_published_connection(PluginFamilyKind::Nonplanar)
    }

    #[test]
    fn first_planar_master_has_the_independent_tadpole_mass_exponent() -> Result<()> {
        let input = HiggsJetIntegralSystem::load(PluginFamilyKind::Planar, "common_mass_first_pl")?;
        let family = input.map.family()?;
        // The first canonical component is a closed factorized sector. Its
        // massless bubble is fixed along AMF; the squared massive tadpole scales
        // as (1+eta)^(-epsilon), independently of the canonical gauge formula.
        assert!(
            input
                .canonical
                .constant_matrices()
                .iter()
                .all(|m| m[0][1..].iter().all(|a| a.is_zero()))
        );
        assert_eq!(input.map.canonical_to_publication[0].len(), 1);
        assert_eq!(input.map.publication_to_canonical[0].len(), 1);
        let matrices = input
            .canonical
            .constant_matrices()
            .iter()
            .map(|m| vec![vec![m[0][0].clone()]])
            .collect::<Vec<_>>();
        let canonical = CanonicalAlgebraicSystem::new(
            input.map.epsilon,
            input.canonical.variables(),
            input.canonical.letters(),
            &matrices,
            input.canonical.roots().to_vec(),
        )?;
        let forward = vec![vec![
            input.map.canonical_to_publication[0][0].coefficient.clone(),
        ]];
        let inverse = vec![vec![
            input.map.publication_to_canonical[0][0].coefficient.clone(),
        ]];
        let point = KinematicPoint(
            input
                .map
                .coordinates
                .into_iter()
                .map(Atom::var)
                .zip([(1, 10), (1, 25), (1, 50)].map(|(a, b)| Atom::num(Rational::from((-a, b)))))
                .collect(),
        );
        let eta = symbol!("common_mass_first_pl::eta");
        let basis = HomogeneousCanonicalBasis {
            family: &family,
            integrals: &input.map.integrals[..1],
            canonical: &canonical,
            ordinary_from_canonical: &forward,
            canonical_from_ordinary: &inverse,
            provenance: &input.map.evidence.source_sha256,
        };
        let connection = basis.pullback(&point, eta, &RunContext::default())?;
        assert_eq!(
            connection.deformation_mask,
            vec![false, false, false, false, false, true, true, false, false]
        );
        assert!(
            (&connection.reduced.matrix[0][0]
                + Atom::var(input.map.epsilon) / (Atom::one() + Atom::var(eta)))
            .together()
            .cancel()
            .is_zero()
        );
        assert!(
            connection
                .reduced
                .nonzero_conditions
                .iter()
                .any(|c| c == &(Atom::one() + Atom::var(eta)))
        );
        Ok(())
    }

    #[test]
    fn common_mass_rejects_an_incorrect_exact_inverse_and_nonunit_masses() -> Result<()> {
        let input = HiggsJetIntegralSystem::load(PluginFamilyKind::Planar, "common_mass_invalid")?;
        let family = input.map.family()?;
        let matrices = input
            .canonical
            .constant_matrices()
            .iter()
            .map(|m| vec![vec![m[0][0].clone()]])
            .collect::<Vec<_>>();
        let canonical = CanonicalAlgebraicSystem::new(
            input.map.epsilon,
            input.canonical.variables(),
            input.canonical.letters(),
            &matrices,
            input.canonical.roots().to_vec(),
        )?;
        let forward = vec![vec![Atom::one()]];
        let inverse = vec![vec![Atom::num(2)]];
        let point = KinematicPoint(
            input
                .map
                .coordinates
                .into_iter()
                .map(Atom::var)
                .zip([(1, 10), (1, 25), (1, 50)].map(|(a, b)| Atom::num(Rational::from((-a, b)))))
                .collect(),
        );
        let eta = symbol!("common_mass_invalid::eta");
        let basis = HomogeneousCanonicalBasis {
            family: &family,
            integrals: &input.map.integrals[..1],
            canonical: &canonical,
            ordinary_from_canonical: &forward,
            canonical_from_ordinary: &inverse,
            provenance: &input.map.evidence.source_sha256,
        };
        assert!(
            matches!(basis.pullback(&point, eta, &RunContext::default()), Err(Error::InvalidInput(message)) if message.contains("inverse identity"))
        );
        let mut unequal_masses = family.clone();
        unequal_masses.propagators[5].constant -= Atom::one();
        let basis = HomogeneousCanonicalBasis {
            family: &unequal_masses,
            ..basis
        };
        assert!(
            matches!(basis.validate(eta), Err(Error::Unsupported(message)) if message.contains("unit physical masses"))
        );
        Ok(())
    }
}
