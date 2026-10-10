//! Native, guarded reduction of identities derived by the weighted-measure owner.
//!
//! This adapter does not derive identities or certify their physical validity.
//! RustRed certifies exact consequences of the supplied corpus, including its
//! integer domains and nonzero conditions. Uncovered domains remain explicit.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rustred::algebra::{Coefficient, CoefficientPolynomial};
use rustred::persistence::BinaryIoLimits;
pub use rustred::solver::guarded::{
    GuardedApplicationFailure, GuardedReductionLimits, GuardedUnresolved, IndexBounds, IndexDomain,
    IndexRole,
};
use rustred::solver::guarded::{GuardedProgram, GuardedSource, GuardedSourceSystem};
use rustred::solver::{Integral as NativeIntegral, SearchOptions, Term};
use serde::{Deserialize, Serialize};
use symbolica::prelude::*;

use crate::{Error, Result};

/// Caller-owned descriptions of the actual distribution and deformation.
/// Every field participates in persisted-program identity. Descriptions must
/// include physical parameter assignments, rather than only benchmark labels.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GuardedMeasureIdentity {
    pub measure: String,
    pub support: String,
    pub orientation: String,
    pub normalization: String,
    pub branch: String,
    pub deformation: String,
}

impl GuardedMeasureIdentity {
    fn native_identity(&self) -> Result<String> {
        if [
            &self.measure,
            &self.support,
            &self.orientation,
            &self.normalization,
            &self.branch,
            &self.deformation,
        ]
        .iter()
        .any(|part| part.trim().is_empty())
        {
            return Err(Error::InvalidInput(
                "guarded reduction requires measure, support, orientation, normalization, branch and deformation identities".into(),
            ));
        }
        let metadata = serde_json::to_string(self).map_err(|e| Error::Cache(e.to_string()))?;
        Ok(format!(
            "rustflow-weighted-sources-v1:{}:{}:{metadata}",
            env!("RUSTRED_SOURCE_DIGEST"),
            env!("DEPENDENCY_SOURCE_DIGEST"),
        ))
    }
}

/// A coefficient multiplying `I(n + shift)` in a homogeneous identity.
#[derive(Clone, Debug)]
pub struct GuardedIdentityTerm<const N: usize> {
    pub shift: [i16; N],
    pub coefficient: Atom,
}

#[derive(Clone, Debug)]
pub struct GuardedIdentity<const N: usize> {
    pub id: String,
    pub terms: Vec<GuardedIdentityTerm<N>>,
    pub domain: IndexDomain<N>,
    pub nonzero_conditions: Vec<Atom>,
}

/// Bounded native discovery, independent of the physical loop count.
#[derive(Clone, Copy, Debug)]
pub struct GuardedDiscoveryOptions {
    pub max_depth: u32,
    pub max_domains: usize,
    pub sample_seed: u64,
}

impl Default for GuardedDiscoveryOptions {
    fn default() -> Self {
        Self {
            max_depth: 2,
            max_domains: 256,
            sample_seed: 0,
        }
    }
}

#[derive(Debug)]
pub struct GuardedContext<const N: usize> {
    sources: Arc<GuardedSourceSystem<N>>,
    physical_arity: usize,
    dummy_symbols: Vec<Symbol>,
}

#[derive(Debug)]
pub struct GuardedDiscovery<const N: usize> {
    pub program: GuardedReductionProgram<N>,
    /// These domains are not certified covered and are never called masters.
    pub unresolved: Vec<GuardedUnresolved<N>>,
}

#[derive(Debug)]
pub struct GuardedReductionProgram<const N: usize> {
    native: GuardedProgram<N>,
    physical_arity: usize,
    dummy_symbols: Vec<Symbol>,
}

#[derive(Clone, Debug)]
pub struct GuardedAtomUnresolved<const N: usize> {
    pub integral: [i64; N],
    pub coefficient: Atom,
    pub reason: GuardedApplicationFailure,
}

#[derive(Clone, Debug)]
pub struct GuardedAtomReduction<const N: usize> {
    pub terms: BTreeMap<[i64; N], Atom>,
    pub unresolved: Vec<GuardedAtomUnresolved<N>>,
    pub nonzero_conditions: Vec<Atom>,
    pub rule_applications: usize,
}

impl<const N: usize> GuardedContext<N> {
    /// Attach concrete measure-owner zero evidence before discovery. The native
    /// owner binds these exact domains into application and persisted replay;
    /// no ordinary sector census or algebraic rank inference is involved.
    pub(crate) fn with_measure_zero_domains(self, domains: Vec<IndexDomain<N>>) -> Result<Self> {
        for domain in &domains {
            validate_storage_domain(domain, self.physical_arity)?;
        }
        let sources = Arc::try_unwrap(self.sources).map_err(|_| {
            Error::InvalidInput(
                "measure zero evidence must be attached before guarded discovery".into(),
            )
        })?;
        let sources = sources
            .with_zero_domains(domains)
            .map_err(|error| Error::Reduction(error.to_string()))?;
        Ok(Self {
            sources: Arc::new(sources),
            physical_arity: self.physical_arity,
            dummy_symbols: self.dummy_symbols,
        })
    }

    /// Lower native Symbolica expressions onto one exact coefficient map.
    /// Index symbols occur first; parameters follow in the supplied order.
    /// Rational source coefficients are cleared exactly, retaining all poles.
    pub fn new(
        measure: GuardedMeasureIdentity,
        roles: [IndexRole; N],
        index_symbols: [Symbol; N],
        parameter_symbols: Vec<Symbol>,
        identities: Vec<GuardedIdentity<N>>,
    ) -> Result<Self> {
        Self::new_with_physical_arity(
            measure,
            roles,
            index_symbols,
            parameter_symbols,
            identities,
            N,
        )
    }

    /// Embed actual physical coordinates in fixed inline storage. No factor is
    /// introduced for a storage axis: its index is identically zero, all source
    /// shifts vanish there, and its symbol cannot occur in physical algebra.
    pub fn new_with_physical_arity(
        measure: GuardedMeasureIdentity,
        roles: [IndexRole; N],
        index_symbols: [Symbol; N],
        parameter_symbols: Vec<Symbol>,
        identities: Vec<GuardedIdentity<N>>,
        physical_arity: usize,
    ) -> Result<Self> {
        validate_storage_arity::<N>(physical_arity)?;
        if roles[physical_arity..]
            .iter()
            .any(|role| *role != IndexRole::Ordinary)
        {
            return Err(Error::InvalidInput(
                "storage tail roles must be ordinary".into(),
            ));
        }
        let dummy_symbols = index_symbols[physical_arity..].to_vec();
        for identity in &identities {
            validate_storage_domain(&identity.domain, physical_arity)?;
            for term in &identity.terms {
                if term.shift[physical_arity..].iter().any(|&shift| shift != 0) {
                    return Err(Error::InvalidInput(
                        "guarded source shifts a storage tail coordinate".into(),
                    ));
                }
                validate_storage_expression(&term.coefficient, &dummy_symbols)?;
            }
            for condition in &identity.nonzero_conditions {
                validate_storage_expression(condition, &dummy_symbols)?;
            }
        }
        let symbols: Vec<_> = index_symbols.into_iter().chain(parameter_symbols).collect();
        let distinct: BTreeSet<_> = symbols.iter().copied().collect();
        if distinct.len() != symbols.len() {
            return Err(Error::InvalidInput(
                "guarded coefficient/index symbols must be distinct".into(),
            ));
        }
        let variables = Arc::new(symbols.into_iter().map(PolyVariable::Symbol).collect());
        let mut sources = Vec::with_capacity(identities.len());
        for identity in identities {
            // Preserve the domains of the supplied Atom syntax before the
            // rational polynomial owner cancels removable factors. A caller
            // that simplified earlier must provide those conditions explicitly.
            let raw_poles = crate::physical_conditions::rational_denominator_conditions(
                &identity
                    .terms
                    .iter()
                    .map(|term| term.coefficient.clone())
                    .chain(identity.nonzero_conditions.iter().cloned())
                    .collect::<Vec<_>>(),
                &distinct,
            )?;
            let coefficients = identity
                .terms
                .iter()
                .map(|term| native_coefficient(&term.coefficient, &variables))
                .collect::<Result<Vec<_>>>()?;
            let mut conditions = Vec::new();
            for condition in identity.nonzero_conditions.iter().chain(&raw_poles) {
                let converted = native_coefficient(condition, &variables)?;
                if converted.numerator.is_zero() {
                    return Err(Error::InvalidInput(
                        "guarded source declares an identically zero nonzero condition".into(),
                    ));
                }
                conditions.push(converted.numerator);
                if !converted.denominator.is_one() {
                    conditions.push(converted.denominator);
                }
            }
            for coefficient in &coefficients {
                if !coefficient.denominator.is_one() {
                    conditions.push(coefficient.denominator.clone());
                }
            }
            let mut row = Vec::with_capacity(identity.terms.len());
            for (position, term) in identity.terms.into_iter().enumerate() {
                let mut coefficient = coefficients[position].numerator.clone();
                for (other, rational) in coefficients.iter().enumerate() {
                    if other != position {
                        coefficient = &coefficient * &rational.denominator;
                    }
                }
                row.push(Term {
                    integral: NativeIntegral::symbolic(term.shift)
                        .map_err(|e| Error::Unsupported(e.to_string()))?,
                    coefficient,
                });
            }
            sources.push(
                GuardedSource::new(identity.id, row, identity.domain)
                    .with_nonzero_conditions(conditions),
            );
        }
        let mut native_identity = measure.native_identity()?;
        if physical_arity != N {
            native_identity.push_str(&format!(
                ":zero-tail-storage-v1:physical={physical_arity}:capacity={N}"
            ));
        }
        let sources =
            GuardedSourceSystem::new(native_identity, roles, std::array::from_fn(|i| i), sources)
                .map_err(|e| Error::Reduction(e.to_string()))?;
        Ok(Self {
            sources: Arc::new(sources),
            physical_arity,
            dummy_symbols,
        })
    }

    pub fn physical_arity(&self) -> usize {
        self.physical_arity
    }

    pub fn sources(&self) -> &Arc<GuardedSourceSystem<N>> {
        &self.sources
    }

    /// Terminal indices are explicit stopping points for this lower-level
    /// service. They do not alter native ordering or establish independence.
    pub fn discover(
        &self,
        domains: Vec<IndexDomain<N>>,
        terminals: impl IntoIterator<Item = [i64; N]>,
        options: GuardedDiscoveryOptions,
    ) -> Result<GuardedDiscovery<N>> {
        self.discover_with_priority_points(domains, terminals, &[], options)
    }

    /// Spend the same bounded native search budget on domains containing these
    /// indices first. Hints do not become terminals or establish coverage:
    /// every unvisited box is still returned as an unresolved discovery gap.
    pub fn discover_with_priority_points(
        &self,
        domains: Vec<IndexDomain<N>>,
        terminals: impl IntoIterator<Item = [i64; N]>,
        priority_points: &[[i64; N]],
        options: GuardedDiscoveryOptions,
    ) -> Result<GuardedDiscovery<N>> {
        self.discover_scheduled(domains, terminals, priority_points, options, None)
    }

    /// Allocate a separate bounded native search to each supplied domain.
    /// The total is a conservative allocation cap, not a measured visit count.
    /// All rules are assembled and source-replayed once; unvisited domains
    /// remain native discovery gaps, including when the total budget is zero.
    pub fn discover_partitioned(
        &self,
        domains: Vec<IndexDomain<N>>,
        terminals: impl IntoIterator<Item = [i64; N]>,
        priority_points: &[[i64; N]],
        options: GuardedDiscoveryOptions,
        max_domains_per_domain: usize,
    ) -> Result<GuardedDiscovery<N>> {
        self.discover_scheduled(
            domains,
            terminals,
            priority_points,
            options,
            Some(max_domains_per_domain),
        )
    }

    fn discover_scheduled(
        &self,
        domains: Vec<IndexDomain<N>>,
        terminals: impl IntoIterator<Item = [i64; N]>,
        priority_points: &[[i64; N]],
        options: GuardedDiscoveryOptions,
        per_domain: Option<usize>,
    ) -> Result<GuardedDiscovery<N>> {
        for domain in &domains {
            validate_storage_domain(domain, self.physical_arity)?;
        }
        let terminals = terminals.into_iter().collect::<Vec<_>>();
        for terminal in &terminals {
            validate_storage_label(terminal, self.physical_arity)?;
        }
        for point in priority_points {
            validate_storage_label(point, self.physical_arity)?;
        }
        let search = SearchOptions {
            max_depth: Some(options.max_depth),
            sample_seed: options.sample_seed,
            ..Default::default()
        };
        let found = if let Some(per_domain) = per_domain {
            self.sources.solve_domains_partitioned(
                domains,
                priority_points,
                search,
                per_domain,
                options.max_domains,
            )
        } else {
            self.sources.solve_domains_with_priority_points(
                domains,
                priority_points,
                search,
                options.max_domains,
            )
        }
        .map_err(|e| Error::Reduction(e.to_string()))?;
        let native = GuardedProgram::new(self.sources.clone(), found.rules, terminals)
            .map_err(|e| Error::Reduction(e.to_string()))?;
        validate_storage_program(&native, self.physical_arity)?;
        Ok(GuardedDiscovery {
            program: GuardedReductionProgram {
                native,
                physical_arity: self.physical_arity,
                dummy_symbols: self.dummy_symbols.clone(),
            },
            unresolved: found.unresolved,
        })
    }

    /// Loading replays the stored rules against this exact source corpus and
    /// complete measure identity. No serialized rule is accepted on trust.
    pub fn decode(
        &self,
        bytes: &[u8],
        limits: BinaryIoLimits,
    ) -> Result<GuardedReductionProgram<N>> {
        let native = GuardedProgram::decode_generated(bytes, self.sources.clone(), limits)
            .map_err(|e| Error::Cache(e.to_string()))?;
        validate_storage_program(&native, self.physical_arity)?;
        Ok(GuardedReductionProgram {
            native,
            physical_arity: self.physical_arity,
            dummy_symbols: self.dummy_symbols.clone(),
        })
    }
}

impl<const N: usize> GuardedReductionProgram<N> {
    pub fn physical_arity(&self) -> usize {
        self.physical_arity
    }

    pub fn native(&self) -> &GuardedProgram<N> {
        &self.native
    }

    pub fn encode(&self, limits: BinaryIoLimits) -> Result<Vec<u8>> {
        self.native
            .encode_native(limits)
            .map_err(|e| Error::Cache(e.to_string()))
    }

    /// Replace the stopping set and replay the unchanged native proof corpus.
    /// Prior provisional terminals must not suppress residual search requests.
    pub fn with_terminals_replayed(
        self,
        terminals: impl IntoIterator<Item = [i64; N]>,
        max_rules: usize,
    ) -> Result<Self> {
        let terminals = terminals.into_iter().collect::<Vec<_>>();
        for terminal in &terminals {
            validate_storage_label(terminal, self.physical_arity)?;
        }
        let native = self
            .native
            .with_terminals_replayed(terminals, max_rules)
            .map_err(|e| Error::Reduction(e.to_string()))?;
        validate_storage_program(&native, self.physical_arity)?;
        Ok(Self {
            native,
            physical_arity: self.physical_arity,
            dummy_symbols: self.dummy_symbols,
        })
    }

    /// Retain previously proved rules as fallbacks, replaying the complete
    /// union natively against one exact source corpus. Only the supplied
    /// stopping set survives; prior terminals never become implicit masters.
    pub fn union_replayed(
        self,
        fallback: Self,
        terminals: impl IntoIterator<Item = [i64; N]>,
        max_rules: usize,
    ) -> Result<Self> {
        if self.physical_arity != fallback.physical_arity
            || self.dummy_symbols != fallback.dummy_symbols
        {
            return Err(Error::InvalidInput(
                "cannot reuse guarded rules with different physical storage mappings".into(),
            ));
        }
        let terminals = terminals.into_iter().collect::<Vec<_>>();
        for terminal in &terminals {
            validate_storage_label(terminal, self.physical_arity)?;
        }
        let native = self
            .native
            .union_replayed(fallback.native, terminals, max_rules)
            .map_err(|e| Error::Reduction(e.to_string()))?;
        validate_storage_program(&native, self.physical_arity)?;
        Ok(Self {
            native,
            physical_arity: self.physical_arity,
            dummy_symbols: self.dummy_symbols,
        })
    }

    pub fn reduce(
        &self,
        target: [i64; N],
        limits: GuardedReductionLimits,
    ) -> Result<GuardedAtomReduction<N>> {
        validate_storage_label(&target, self.physical_arity)?;
        let reduced = self
            .native
            .reduce(target, limits)
            .map_err(|e| Error::Reduction(e.to_string()))?;
        let result = GuardedAtomReduction {
            terms: reduced
                .terms
                .into_iter()
                .map(|(indices, coefficient)| (indices, coefficient.to_expression()))
                .collect(),
            unresolved: reduced
                .unresolved
                .into_iter()
                .map(|term| GuardedAtomUnresolved {
                    integral: term.integral,
                    coefficient: term.coefficient.to_expression(),
                    reason: term.reason,
                })
                .collect(),
            nonzero_conditions: reduced
                .nonzero_conditions
                .iter()
                .map(CoefficientPolynomial::to_expression)
                .collect(),
            rule_applications: reduced.rule_applications,
        };
        for (indices, coefficient) in &result.terms {
            validate_storage_label(indices, self.physical_arity)?;
            validate_storage_expression(coefficient, &self.dummy_symbols)?;
        }
        for unresolved in &result.unresolved {
            validate_storage_label(&unresolved.integral, self.physical_arity)?;
            validate_storage_expression(&unresolved.coefficient, &self.dummy_symbols)?;
        }
        for condition in &result.nonzero_conditions {
            validate_storage_expression(condition, &self.dummy_symbols)?;
        }
        Ok(result)
    }
}

pub(crate) fn validate_storage_arity<const N: usize>(physical_arity: usize) -> Result<()> {
    if physical_arity == 0 || physical_arity > N {
        return Err(Error::InvalidInput(format!(
            "physical arity {physical_arity} does not fit guarded storage capacity {N}"
        )));
    }
    Ok(())
}

pub(crate) fn validate_storage_label<const N: usize>(
    indices: &[i64; N],
    physical_arity: usize,
) -> Result<()> {
    validate_storage_arity::<N>(physical_arity)?;
    if indices[physical_arity..].iter().any(|&index| index != 0) {
        return Err(Error::InvalidInput(
            "nonzero storage tail is not a physical integral".into(),
        ));
    }
    Ok(())
}

pub(crate) fn validate_storage_domain<const N: usize>(
    domain: &IndexDomain<N>,
    physical_arity: usize,
) -> Result<()> {
    validate_storage_arity::<N>(physical_arity)?;
    if domain.bounds()[physical_arity..]
        .iter()
        .any(|bound| *bound != IndexBounds::fixed(0))
    {
        return Err(Error::InvalidInput(
            "guarded storage tail domains must be fixed at zero".into(),
        ));
    }
    Ok(())
}

fn validate_storage_expression(expression: &Atom, dummy_symbols: &[Symbol]) -> Result<()> {
    if dummy_symbols
        .iter()
        .any(|&symbol| expression.contains_symbol(symbol))
    {
        return Err(Error::InvalidInput(
            "physical algebra depends on a storage tail symbol".into(),
        ));
    }
    Ok(())
}

fn validate_storage_program<const N: usize>(
    program: &GuardedProgram<N>,
    physical_arity: usize,
) -> Result<()> {
    validate_storage_arity::<N>(physical_arity)?;
    for terminal in program.terminals() {
        validate_storage_label(terminal, physical_arity)?;
    }
    for rule in program.rules() {
        validate_storage_domain(rule.domain(), physical_arity)?;
        for integral in std::iter::once(&rule.candidate().target)
            .chain(rule.candidate().rhs.iter().map(|term| &term.integral))
        {
            if integral.powers()[physical_arity..]
                .iter()
                .any(|power| power.is_symbolic() || power.value() != 0)
            {
                return Err(Error::Reduction(
                    "native rule escaped fixed-zero storage tail".into(),
                ));
            }
        }
    }
    Ok(())
}

fn native_coefficient(atom: &Atom, variables: &Arc<Vec<PolyVariable>>) -> Result<Coefficient> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(atom.as_view(), &mut symbols)?;
    for symbol in symbols {
        let AtomView::Var(variable) = symbol.as_view() else {
            unreachable!()
        };
        if !variables.contains(&PolyVariable::Symbol(variable.get_symbol())) {
            return Err(Error::InvalidInput(format!(
                "undeclared guarded coefficient variable {symbol}"
            )));
        }
    }
    atom.try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))
        .map_err(|e| Error::InvalidInput(format!("invalid guarded rational coefficient: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn measure() -> GuardedMeasureIdentity {
        GuardedMeasureIdentity {
            measure: "synthetic-adapter-regression-v1".into(),
            support: "occupation-zero-is-bulk".into(),
            orientation: "positive".into(),
            normalization: "unit".into(),
            branch: "rational-real".into(),
            deformation: "none".into(),
        }
    }

    fn context(identity: GuardedMeasureIdentity) -> GuardedContext<2> {
        let a = symbol!("guarded_test_a");
        let b = symbol!("guarded_test_b");
        let x = symbol!("guarded_test_x");
        let x_atom = Atom::var(x);
        GuardedContext::new(
            identity,
            [IndexRole::RequiredCut, IndexRole::Occupation],
            [a, b],
            vec![x],
            vec![GuardedIdentity {
                id: "bulk-surface-test-relation".into(),
                terms: vec![
                    GuardedIdentityTerm {
                        shift: [0, 0],
                        coefficient: Atom::num(1) / &x_atom,
                    },
                    GuardedIdentityTerm {
                        shift: [-1, 1],
                        coefficient: Atom::num(-1) / &x_atom,
                    },
                ],
                domain: IndexDomain::new([
                    IndexBounds::new(Some(2), None).unwrap(),
                    IndexBounds::fixed(0),
                ])
                .unwrap(),
                nonzero_conditions: Vec::new(),
            }],
        )
        .unwrap()
    }

    #[test]
    fn roles_exact_replay_conditions_and_measure_bound_persistence() {
        let context = context(measure());
        assert!(!context.sources().is_zero(&[1, 0]));
        assert!(context.sources().is_zero(&[0, 0]));
        assert!(!context.sources().valid_indices(&[1, -1]));
        let found = context
            .discover(
                vec![
                    IndexDomain::new([
                        IndexBounds::new(Some(2), None).unwrap(),
                        IndexBounds::fixed(0),
                    ])
                    .unwrap(),
                ],
                [[1, 1]],
                Default::default(),
            )
            .unwrap();
        assert!(found.unresolved.is_empty());
        let reduced = found.program.reduce([2, 0], Default::default()).unwrap();
        assert!(reduced.unresolved.is_empty());
        assert_eq!(reduced.terms[&[1, 1]], Atom::num(1));
        assert!(
            reduced
                .nonzero_conditions
                .iter()
                .any(|condition| { condition == &Atom::var(symbol!("guarded_test_x")) })
        );
        let bytes = found.program.encode(Default::default()).unwrap();
        let restored = context.decode(&bytes, Default::default()).unwrap();
        assert_eq!(
            restored.reduce([2, 0], Default::default()).unwrap().terms,
            reduced.terms
        );
        let mut other_measure = measure();
        other_measure.deformation = "a-different-fixed-shell-path".into();
        assert!(
            self::context(other_measure)
                .decode(&bytes, Default::default())
                .is_err()
        );
        let uncovered = restored.reduce([1, 0], Default::default()).unwrap();
        assert!(uncovered.terms.is_empty());
        assert_eq!(uncovered.unresolved.len(), 1);
        let invalid = restored.reduce([1, -1], Default::default()).unwrap();
        assert_eq!(invalid.unresolved.len(), 1);
    }

    #[test]
    fn priority_hints_preserve_storage_and_occupation_domains() {
        let domain = IndexDomain::new([
            IndexBounds::new(Some(2), None).unwrap(),
            IndexBounds::fixed(0),
            IndexBounds::fixed(0),
        ])
        .unwrap();
        let context = GuardedContext::new_with_physical_arity(
            measure(),
            [
                IndexRole::RequiredCut,
                IndexRole::Occupation,
                IndexRole::Ordinary,
            ],
            [
                symbol!("priority_cut"),
                symbol!("priority_occupation"),
                symbol!("priority_tail"),
            ],
            vec![],
            vec![GuardedIdentity {
                id: "padded-priority-source".into(),
                terms: vec![
                    GuardedIdentityTerm {
                        shift: [0, 0, 0],
                        coefficient: Atom::one(),
                    },
                    GuardedIdentityTerm {
                        shift: [-1, 1, 0],
                        coefficient: -Atom::one(),
                    },
                ],
                domain: domain.clone(),
                nonzero_conditions: vec![],
            }],
            2,
        )
        .unwrap();
        let discover = |points: &[[i64; 3]]| {
            context.discover_with_priority_points(
                vec![domain.clone()],
                [[1, 1, 0]],
                points,
                Default::default(),
            )
        };
        assert!(discover(&[[2, 0, 1]]).is_err());
        assert!(discover(&[[2, -1, 0]]).is_err());
        let found = discover(&[[2, 0, 0]]).unwrap();
        let restored = context
            .decode(
                &found.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let reduced = restored.reduce([2, 0, 0], Default::default()).unwrap();
        assert!(reduced.unresolved.is_empty());
        // The hinted index was reduced; it was not added as a terminal.
        assert_eq!(reduced.terms, BTreeMap::from([([1, 1, 0], Atom::one())]));
        let uncovered = restored.reduce([1, 0, 0], Default::default()).unwrap();
        assert!(uncovered.terms.is_empty());
        assert!(!uncovered.unresolved.is_empty());
    }

    #[test]
    fn uncancelled_source_and_condition_denominators_remain_required() {
        let a = symbol!("guarded_raw_poles_a");
        let b = symbol!("guarded_raw_poles_b");
        let x = symbol!("guarded_raw_poles_x");
        let variable = Atom::var(x);
        let coefficient = (variable.pow(2) - Atom::one()) / (&variable - Atom::one());
        let condition = (variable.pow(2) - Atom::num(4)) / (&variable - Atom::num(2));
        assert!(
            (coefficient.together().cancel() - (&variable + Atom::one()))
                .expand()
                .is_zero()
        );
        let context = GuardedContext::new(
            measure(),
            [IndexRole::RequiredCut, IndexRole::Occupation],
            [a, b],
            vec![x],
            vec![GuardedIdentity {
                id: "removable-source-and-condition-poles".into(),
                terms: vec![
                    GuardedIdentityTerm {
                        shift: [0, 0],
                        coefficient: coefficient.clone(),
                    },
                    GuardedIdentityTerm {
                        shift: [-1, 1],
                        coefficient: -coefficient,
                    },
                ],
                domain: IndexDomain::new([
                    IndexBounds::new(Some(2), None).unwrap(),
                    IndexBounds::fixed(0),
                ])
                .unwrap(),
                nonzero_conditions: vec![condition],
            }],
        )
        .unwrap();
        for forbidden in [1, 2] {
            let replacements = BTreeMap::from([(variable.clone(), Atom::num(forbidden))]);
            assert!(
                context.sources().sources()[0]
                    .nonzero_conditions
                    .iter()
                    .any(|condition| {
                        crate::family::substitute(&condition.to_expression(), &replacements)
                            .is_zero()
                    })
            );
        }
        let found = context
            .discover(
                vec![
                    IndexDomain::new([
                        IndexBounds::new(Some(2), None).unwrap(),
                        IndexBounds::fixed(0),
                    ])
                    .unwrap(),
                ],
                [[1, 1]],
                Default::default(),
            )
            .unwrap();
        let reduced = found.program.reduce([2, 0], Default::default()).unwrap();
        assert!(reduced.unresolved.is_empty());
        for forbidden in [1, 2] {
            let replacements = BTreeMap::from([(variable.clone(), Atom::num(forbidden))]);
            assert!(reduced.nonzero_conditions.iter().any(|condition| {
                crate::family::substitute(condition, &replacements).is_zero()
            }));
        }
    }

    #[test]
    fn generated_ball_sources_reduce_surface_moments_to_volume() {
        use crate::finite_density::measure::WeightedMeasure;

        let x = symbol!("guarded_ball_x");
        let radius = symbol!("guarded_ball_R");
        let a = symbol!("guarded_ball_a");
        let b = symbol!("guarded_ball_b");
        let x_atom = Atom::var(x);
        let r = Atom::var(radius);
        let roles = [IndexRole::Ordinary, IndexRole::Occupation];
        let weighted = WeightedMeasure::new(
            vec![x_atom.clone()],
            [x_atom.clone(), &r - x_atom.clone().pow(2)],
            roles,
        )
        .unwrap();
        let mut identities = weighted
            .ibp("ball-dilation", &[a, b], &[x_atom], &Atom::num(1), 2)
            .unwrap();
        identities.extend(weighted.multiplication_sources().unwrap());
        // This validation integral contains polynomial x insertions only.
        // Positive a would introduce an unprescribed pole inside the ball;
        // forbid those source instances in the actual replay corpus.
        let polynomial_domain = IndexDomain::new([
            IndexBounds::new(None, Some(0)).unwrap(),
            IndexBounds::new(Some(0), None).unwrap(),
        ])
        .unwrap();
        for identity in &mut identities {
            identity.domain = identity.domain.intersection(&polynomial_domain).unwrap();
        }

        // Independent integration checks the physical sources: for integer
        // k>=0, I(-2k,0)=2 R^(k+1/2)/(2k+1) and
        // I(-2k,1)=R^(k-1/2). Odd moments vanish by x -> -x symmetry.
        // Each generated bulk IBP and first-surface multiplication identity
        // vanishes after these exact integrated values are substituted.
        for k in 0..=4_i64 {
            let index = -2 * k;
            for identity in &identities {
                for occupation in [0, 1] {
                    if !identity.domain.contains(&[index, occupation]) {
                        continue;
                    }
                    if identity
                        .terms
                        .iter()
                        .any(|term| occupation + i64::from(term.shift[1]) > 1)
                    {
                        continue;
                    }
                    let substitution = BTreeMap::from([
                        (Atom::var(a), Atom::num(index)),
                        (Atom::var(b), Atom::num(occupation)),
                    ]);
                    let mut integrated = Atom::new();
                    for term in &identity.terms {
                        let power = index + i64::from(term.shift[0]);
                        let surface = occupation + i64::from(term.shift[1]);
                        assert!(power <= 0 && power % 2 == 0);
                        let k = -power / 2;
                        let moment = if surface == 0 {
                            Atom::num(2) * r.clone().pow(Atom::num(Rational::from((2 * k + 1, 2))))
                                / Atom::num(2 * k + 1)
                        } else {
                            r.clone().pow(Atom::num(Rational::from((2 * k - 1, 2))))
                        };
                        integrated +=
                            crate::family::substitute(&term.coefficient, &substitution) * moment;
                    }
                    assert!(integrated.expand().together().cancel().is_zero());
                }
            }
        }
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "dx on the real line; x^(-a) H_b(R-x^2)".into(),
                support: "R>0; -sqrt(R)<=x<=sqrt(R)".into(),
                orientation: "increasing x; positive delta(R-x^2)".into(),
                normalization: "unit Lebesgue measure".into(),
                branch: "positive sqrt(R); a nonpositive even integer".into(),
                deformation: "none; exact compact seed source validation".into(),
            },
            roles,
            [a, b],
            vec![radius],
            identities,
        )
        .unwrap();
        let targets = [[0, 1], [-2, 1]];
        let found = context
            .discover(
                targets
                    .iter()
                    .map(|target| IndexDomain::new(target.map(IndexBounds::fixed)).unwrap())
                    .collect(),
                [[0, 0], [-2, 0]],
                GuardedDiscoveryOptions {
                    max_depth: 3,
                    ..Default::default()
                },
            )
            .unwrap();
        for (target, expected) in targets.into_iter().zip([
            Atom::num(1) / (Atom::num(2) * &r),
            Atom::num(Rational::from((1, 2))),
        ]) {
            let reduced = found.program.reduce(target, Default::default()).unwrap();
            assert!(reduced.unresolved.is_empty(), "{:?}", reduced.unresolved);
            // Native ordering closes on two polynomial bulk moments here.
            // Their independently integrated values imply the exact basis
            // change I(-2,0)=R/3 I(0,0). Closure is not a minimal-master claim.
            let coefficient = reduced
                .terms
                .iter()
                .map(|(indices, coefficient)| match indices {
                    [0, 0] => coefficient.clone(),
                    [-2, 0] => coefficient * &r / Atom::num(3),
                    unexpected => panic!("unexpected native bulk terminal {unexpected:?}"),
                })
                .fold(Atom::new(), |sum, coefficient| sum + coefficient);
            assert!((coefficient - expected).together().cancel().is_zero());
        }
    }
}
