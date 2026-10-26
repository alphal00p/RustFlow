//! Optional attribution of accepted equations; no search or dispatch logic.
use crate::{Integral, NativeReductionOutcome, NativeRuleKind, NativeRuleProvenance, Progress};
use rustred::algebra::CoefficientPolynomial;
use rustred::solver::{IntegralOrder, NumericResult, RuleCandidate, SectorRule, SourceSystem};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};
use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
use symbolica::prelude::*;

fn hash_field(hash: &mut blake3::Hasher, value: &str) {
    hash.update(&(value.len() as u64).to_le_bytes());
    hash.update(value.as_bytes());
}

// Native coefficients and guards are rational in scalar symbols. Rebuild the
// exact rational pair in a variable order chosen by caller names, not IDs or
// native parameter slots. Atom canonical printing alone does not normalize
// equivalent fractions such as 1/(z-a) and -1/(a-z).
fn caller_rational(value: &Atom) -> RationalPolynomial<IntegerRing> {
    let mut symbols = BTreeSet::new();
    crate::family::scalar_symbols(value.as_view(), &mut symbols)
        .expect("native provenance contains only rational scalar coefficients");
    let mut variables = symbols
        .into_iter()
        .map(|atom| {
            let AtomView::Var(var) = atom.as_view() else {
                unreachable!()
            };
            (
                atom.to_canonical_string(),
                symbolica::poly::PolyVariable::Symbol(var.get_symbol()),
            )
        })
        .collect::<Vec<_>>();
    variables.sort_by(|a, b| a.0.cmp(&b.0));
    let variables = Arc::new(variables.into_iter().map(|(_, var)| var).collect());
    let rational = value
        .try_to_rational_polynomial(&Q, &Z, Some(variables))
        .expect("native provenance coefficients admit exact rational conversion");
    RationalPolynomial::from_num_den(rational.numerator, rational.denominator, &Z, true)
}

fn canonical_rational(value: &Atom) -> String {
    caller_rational(value).to_expression().to_canonical_string()
}

fn canonical_guard(value: &Atom) -> String {
    let rational = caller_rational(value);
    assert!(
        rational.denominator.is_constant(),
        "native guard must be polynomial"
    );
    let mut polynomial = rational.numerator;
    if !polynomial.is_zero() {
        polynomial = polynomial.make_primitive();
        if polynomial.lcoeff().is_negative() {
            polynomial = -polynomial;
        }
    }
    polynomial.to_expression().to_canonical_string()
}

/// Separate from cache identities: internal parameter assignment and Atom
/// registration order must not enter the caller-chart rule identity.
pub(super) fn family_scope(
    family: &crate::family::ConvertedFamily,
    physical_propagators: usize,
    cuts: Option<&crate::cuts::CutMetadata>,
) -> String {
    let mut hash = blake3::Hasher::new();
    hash_field(&mut hash, "rustflow-native-family-scope-v1");
    for names in [
        family.family.loop_momenta(),
        family.family.external_momenta(),
    ] {
        hash_field(&mut hash, &names.len().to_string());
        for name in names {
            hash_field(&mut hash, name);
        }
    }
    hash_field(&mut hash, &physical_propagators.to_string());
    // CutMetadata contains ordered integer slots, rational routings and enums,
    // with no Atom ordering or internal coefficient names.
    hash_field(&mut hash, &format!("{cuts:?}"));
    hash_field(&mut hash, &family.family.denominator_count().to_string());
    hash_field(&mut hash, &family.family.external_gram().len().to_string());
    for group in std::iter::once(std::slice::from_ref(family.family.dimension()))
        .chain(family.family.denominators().iter().flat_map(|denominator| {
            [
                std::slice::from_ref(denominator.constant()),
                denominator.coefficients(),
            ]
        }))
        .chain(family.family.external_gram().iter().map(Vec::as_slice))
        .chain(std::iter::once(family.family.power_shifts()))
    {
        hash_field(&mut hash, &group.len().to_string());
        for coefficient in group {
            let caller = crate::family::substitute(&coefficient.to_expression(), &family.reverse);
            hash_field(&mut hash, &canonical_rational(&caller));
        }
    }
    hash.finalize().to_hex().to_string()
}

pub(super) struct Attribution<'a, const N: usize> {
    sources: &'a SourceSystem<N>,
    order: &'a IntegralOrder<N>,
    family_scope: &'a str,
    reverse: &'a BTreeMap<Atom, Atom>,
    index_map: BTreeMap<Atom, Atom>,
    index_symbols: Vec<Atom>,
}

impl<'a, const N: usize> Attribution<'a, N> {
    pub fn new(
        sources: &'a SourceSystem<N>,
        order: &'a IntegralOrder<N>,
        family_scope: &'a str,
        reverse: &'a BTreeMap<Atom, Atom>,
    ) -> Self {
        let mut forbidden = BTreeSet::new();
        let mut occupied = BTreeSet::new();
        // Include every caller symbol, including unused parameters, and the
        // native names so the two substitution stages cannot capture either.
        for atom in reverse.keys().chain(reverse.values()).cloned().chain(
            sources
                .coefficient_variables()
                .iter()
                .map(|var| var.to_atom()),
        ) {
            let symbols = atom.get_all_symbols(true);
            forbidden.extend(symbols.iter().map(|s| s.get_name().to_owned()));
            occupied.extend(symbols);
        }
        let mut index_map = BTreeMap::new();
        let mut index_symbols = Vec::with_capacity(N);
        for (axis, &variable) in sources.index_variables().iter().enumerate() {
            let (name, symbol) = (0u64..)
                .find_map(|suffix| {
                    let name = format!("symbolica_amflow_provenance_v1::n{axis}_{suffix}");
                    if forbidden.contains(&name) {
                        return None;
                    }
                    let symbol = symbol!(&name);
                    // Also reject aliases of occupied symbols.
                    (!occupied.contains(&symbol) && symbol.get_name() == name)
                        .then_some((name, symbol))
                })
                .unwrap();
            let index = Atom::var(symbol);
            forbidden.insert(name);
            occupied.insert(symbol);
            index_map.insert(
                sources.coefficient_variables()[variable].to_atom(),
                index.clone(),
            );
            index_symbols.push(index);
        }
        Self {
            sources,
            order,
            family_scope,
            reverse,
            index_map,
            index_symbols,
        }
    }

    pub fn polynomial(&self, value: &CoefficientPolynomial) -> Atom {
        // Alpha-rename native indices before caller parameters are restored:
        // callers may legally use RustRed's own native index names.
        let indexed = crate::family::substitute(&value.to_expression(), &self.index_map);
        crate::family::substitute(&indexed, self.reverse)
    }

    pub fn bank_rule(&self, rule: &SectorRule<N>) -> Arc<NativeRuleProvenance> {
        self.rule(
            &rule.candidate,
            &rule.exceptions.branches,
            if rule.candidate.case.is_numerical() {
                NativeRuleKind::BankFixedCase
            } else {
                NativeRuleKind::BankSymbolicIndices
            },
        )
    }

    fn rule(
        &self,
        candidate: &RuleCandidate<N>,
        exceptions: &[Vec<CoefficientPolynomial>],
        kind: NativeRuleKind,
    ) -> Arc<NativeRuleProvenance> {
        let affine_equalities = candidate
            .case
            .affine()
            .map(|case| {
                case.equations()
                    .iter()
                    .map(|p| self.polynomial(p))
                    .collect()
            })
            .unwrap_or_default();
        let mut provenance = NativeRuleProvenance {
            identity: String::new(),
            family_scope: self.family_scope.to_owned(),
            kind,
            fixed_indices: candidate.case.fixed().to_vec(),
            index_symbols: self.index_symbols.clone(),
            affine_equalities,
            excluded_zero_branches: exceptions
                .iter()
                .map(|branch| branch.iter().map(|p| self.polynomial(p)).collect())
                .collect(),
            source_conditions: self
                .sources
                .conditions()
                .iter()
                .map(|p| self.polynomial(p))
                .collect(),
            denominators: candidate
                .rhs
                .iter()
                .filter(|term| !term.coefficient.denominator.is_constant())
                .map(|term| self.polynomial(&term.coefficient.denominator))
                .collect(),
            direct_hit: candidate.stats.direct_hit,
            exact_trace_rows: candidate.stats.exact_trace_rows,
            source_rows: candidate.sources.len(),
        };
        // Rational pairs use caller-name variable order; polynomial guards
        // are primitive with positive leading coefficient in that order.
        // Length framing separates every identity field.
        // RHS coefficients enter only this digest, not the exported event.
        let mut hash = blake3::Hasher::new();
        let mut field = |value: &str| hash_field(&mut hash, value);
        field("rustflow-native-rule-v1");
        field(self.family_scope);
        field(&format!("{:?}", self.order));
        field(&format!("{:?}", candidate.case.fixed()));
        field(&format!(
            "{:?}",
            candidate
                .target
                .powers()
                .map(|p| (p.is_symbolic(), p.value()))
        ));
        field(&provenance.index_symbols.len().to_string());
        for atom in &provenance.index_symbols {
            field(&atom.to_canonical_string());
        }
        for group in [&provenance.affine_equalities, &provenance.source_conditions] {
            let mut guards = group.iter().map(canonical_guard).collect::<Vec<_>>();
            guards.sort();
            field(&group.len().to_string());
            for guard in guards {
                field(&guard);
            }
        }
        let mut branches = provenance
            .excluded_zero_branches
            .iter()
            .map(|branch| {
                let mut guards = branch.iter().map(canonical_guard).collect::<Vec<_>>();
                guards.sort();
                guards
            })
            .collect::<Vec<_>>();
        branches.sort();
        field(&branches.len().to_string());
        for branch in branches {
            field(&branch.len().to_string());
            for guard in branch {
                field(&guard);
            }
        }
        field(&candidate.rhs.len().to_string());
        for term in &candidate.rhs {
            let powers = term.integral.powers().map(|p| (p.is_symbolic(), p.value()));
            field(&format!("{powers:?}"));
            let numerator = self.polynomial(&term.coefficient.numerator);
            let denominator = self.polynomial(&term.coefficient.denominator);
            field(&canonical_rational(&(&numerator / &denominator)));
            // Preserve the original pole guard even if the fraction cancels.
            field(&canonical_guard(&denominator));
        }
        provenance.identity = hash.finalize().to_hex().to_string();
        Arc::new(provenance)
    }

    pub fn numeric_events(&self, solved: &NumericResult<N>, max_depth: u32) -> Vec<Progress> {
        let mut events = Vec::new();
        for rule in &solved.rules {
            let provenance = self.rule(
                rule,
                &[],
                if rule.stats.direct_hit {
                    NativeRuleKind::ConcreteDirect
                } else {
                    NativeRuleKind::ConcreteElimination
                },
            );
            let nonzero_conditions = provenance
                .source_conditions
                .iter()
                .chain(&provenance.denominators)
                .cloned()
                .collect();
            events.push(Progress::NativeReduction {
                target: Integral(rule.target.powers().iter().map(|p| p.value()).collect()),
                sector: self.order.sector().to_vec(),
                outcome: NativeReductionOutcome::Rule {
                    provenance,
                    rhs: rule
                        .rhs
                        .iter()
                        .map(|term| {
                            Integral(term.integral.powers().iter().map(|p| p.value()).collect())
                        })
                        .collect(),
                    nonzero_conditions,
                },
            });
        }
        for case in &solved.residuals {
            events.push(Progress::NativeReduction {
                target: Integral(
                    case.fixed()
                        .iter()
                        .map(|p| p.expect("numeric solver residuals are fully fixed"))
                        .collect(),
                ),
                sector: self.order.sector().to_vec(),
                outcome: NativeReductionOutcome::Residual { max_depth },
            });
        }
        events
    }
}
