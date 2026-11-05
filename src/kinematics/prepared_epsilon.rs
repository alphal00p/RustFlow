use super::*;
use crate::diffexp::EpsilonSystem;
use crate::physical_conditions::{PhysicalPhaseTiming, rational_denominator_conditions};

/// Immutable epsilon expansion of a rational physical connection, independent
/// of any requested point or path. Each partial connection must be regular at
/// epsilon=0. Path-dependent cancellation of epsilon poles remains available
/// through the legacy `KinematicSystem::pullback` route.
#[derive(Clone, Debug)]
pub struct PreparedKinematicEpsilon {
    coefficients: Vec<KinematicSystem>,
    original_poles: Vec<Atom>,
    expansion_poles: Vec<Atom>,
}

fn decode_complex(value: &Atom) -> Atom {
    value
        .replace(Atom::var(crate::family::imaginary_parameter()))
        .with(Atom::i())
}

impl KinematicSystem {
    /// Expand all rational partial matrices once through `max_order` inclusive.
    /// For a solution requested on [leading,last], prepare at least last-leading.
    /// No numerical precision, boundary credit or cache identity is changed.
    pub fn prepare_epsilon(&self, max_order: usize) -> Result<PreparedKinematicEpsilon> {
        let _timing = PhysicalPhaseTiming::start("transport.epsilon_prepare", 0, 0, None);
        self.validate()?;
        if max_order >= 1024 {
            return Err(Error::Limit("at most 1024 epsilon orders supported".into()));
        }
        let mut variables: BTreeSet<_> = self.derivatives.keys().copied().collect();
        variables.insert(self.epsilon);
        let expressions = self
            .derivatives
            .values()
            .flatten()
            .flatten()
            .cloned()
            .collect::<Vec<_>>();
        // Keep original denominators even above the retained order or in a
        // stationary coordinate; expansion/cancellation must not erase them.
        let original_poles = rational_denominator_conditions(&expressions, &variables)?;
        let mut coefficients = vec![
            KinematicSystem {
                epsilon: self.epsilon,
                derivatives: BTreeMap::new(),
            };
            max_order + 1
        ];
        let mut expansion_poles = Vec::new();
        for (&variable, matrix) in &self.derivatives {
            let expanded = EpsilonSystem::from_differential_system(
                &DifferentialSystem {
                    variable,
                    matrix: matrix.clone(),
                },
                self.epsilon,
                max_order,
            )?;
            for (coefficient, mut matrix) in coefficients.iter_mut().zip(expanded.matrices) {
                for entry in matrix.iter_mut().flatten() {
                    *entry = decode_complex(entry);
                    exact_expression(entry.as_view(), &variables)?;
                    if entry.contains_symbol(self.epsilon) {
                        return Err(Error::Unsupported(
                            "epsilon remained in a prepared coefficient".into(),
                        ));
                    }
                }
                coefficient.derivatives.insert(variable, matrix);
            }
            for entry in matrix.iter().flatten() {
                let rational: RationalPolynomial<IntegerRing, u16> =
                    crate::family::encode_complex(entry)
                        .try_to_rational_polynomial(&Q, &Z, None)
                        .map_err(|e| {
                            Error::Unsupported(format!("rational epsilon system required: {e}"))
                        })?;
                // A zero constant denominator after specialization can change
                // epsilon valuation, even when the specialized legacy function
                // is regular (e.g. epsilon/(s+epsilon) at s=0). Reject it rather
                // than silently using the generic Taylor coefficients.
                let constant = decode_complex(&rational.denominator.to_expression())
                    .replace(Atom::var(self.epsilon))
                    .with(Atom::zero())
                    .expand();
                if constant.is_zero() {
                    return Err(Error::Unsupported(
                        "partial matrix has an epsilon pole; supply an epsilon-regular basis"
                            .into(),
                    ));
                }
                if !matches!(constant.as_view(), AtomView::Num(_)) {
                    expansion_poles.push(constant);
                }
            }
        }
        expansion_poles.sort_by_cached_key(AtomCore::to_canonical_string);
        expansion_poles.dedup();
        Ok(PreparedKinematicEpsilon {
            coefficients,
            original_poles,
            expansion_poles,
        })
    }
}

impl PreparedKinematicEpsilon {
    /// Maximum retained connection order (inclusive), independent of the
    /// solution's common Laurent leading power.
    pub fn max_order(&self) -> usize {
        self.coefficients.len() - 1
    }

    pub(crate) fn check_order(&self, order: usize) -> Result<()> {
        if order > self.max_order() {
            return Err(Error::InvalidInput(format!(
                "prepared epsilon depth {} does not cover requested order {order}",
                self.max_order()
            )));
        }
        Ok(())
    }

    /// Pull back only epsilon-independent coefficient matrices. Original poles
    /// and exceptional epsilon specializations are checked before Jacobians or
    /// truncation can hide them. Exceptional lines are rejected; callers may
    /// explicitly try the legacy route to obtain their specialized expansion.
    pub fn pullback(&self, path: &KinematicPath, order: usize) -> Result<EpsilonSystem> {
        self.check_order(order)?;
        self.coefficients[0].validate_path(path)?;
        let substitutions = path
            .coordinates
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect();
        for (poles, expansion) in [(&self.original_poles, false), (&self.expansion_poles, true)] {
            for pole in poles {
                // These are retained polynomials: expansion suffices, with no
                // rational (epsilon,parameter) GCD or epsilon series per path.
                if substitute(pole, &substitutions).expand().is_zero() {
                    return Err(if expansion {
                        Error::Unsupported("path has an exceptional epsilon specialization; use the legacy pullback".into())
                    } else {
                        Error::InvalidInput("path lies in an original kinematic matrix pole".into())
                    });
                }
            }
        }
        let mut matrices = Vec::with_capacity(order + 1);
        for coefficient in &self.coefficients[..=order] {
            let pulled = coefficient.pullback_validated(path)?;
            pulled.validate()?;
            let mut matrix = pulled.matrix;
            // Match EpsilonSystem's formal complex-coefficient convention.
            for entry in matrix.iter_mut().flatten() {
                *entry = crate::family::encode_complex(entry).together().cancel();
            }
            matrices.push(matrix);
        }
        Ok(EpsilonSystem {
            variable: path.parameter,
            matrices,
        })
    }
}
