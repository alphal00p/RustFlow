use crate::{Error, Result};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecursionMode {
    Amf,
    Ft,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Prescription {
    PlusI0,
    MinusI0,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MassMode {
    /// All lines at one loop; one massive line (or the first line) at higher loops.
    Auto,
    All,
    /// Shift the smallest group of lines with an equal nonzero intrinsic mass.
    Mass,
    /// Shift a single line in a branch containing the most propagators.
    Propagator,
    /// Shift the branch containing the fewest propagators.
    Branch,
    /// Shift the smallest set of lines outside an independent (L-1)-branch span.
    Loop,
    Propagators(Vec<usize>),
}

#[derive(Debug, Clone)]
pub struct FlowOptions {
    pub digits: u32,
    pub guard_digits: u32,
    pub series_order: usize,
    pub max_steps: usize,
    pub max_precision_attempts: usize,
    pub workers: usize,
    pub dimension: i64,
    pub recursion: RecursionMode,
    pub mass_mode: MassMode,
    pub prescription: Prescription,
    /// Factor dimension dependence before numerical specialization. This takes
    /// precedence over sampled_reduction and retains symbolic epsilon in IBPs.
    pub refine_basis: bool,
    pub skip_reduction: bool,
    pub sampled_reduction: bool,
    pub cache_directory: Option<std::path::PathBuf>,
}
impl Default for FlowOptions {
    fn default() -> Self {
        Self {
            digits: 20,
            guard_digits: 40,
            series_order: 80,
            max_steps: 1000,
            max_precision_attempts: 3,
            workers: 1,
            dimension: 4,
            recursion: RecursionMode::Amf,
            mass_mode: MassMode::Auto,
            prescription: Prescription::PlusI0,
            refine_basis: false,
            skip_reduction: false,
            sampled_reduction: true,
            cache_directory: None,
        }
    }
}
impl FlowOptions {
    pub fn validate(&self) -> Result<()> {
        if self.digits == 0
            || self.workers == 0
            || self.series_order < 8
            || self.max_steps == 0
            || self.max_precision_attempts == 0
        {
            return Err(Error::InvalidInput(
                "positive digits, workers, limits and series_order >= 8 required".into(),
            ));
        }
        if self.digits.checked_add(self.guard_digits).is_none() {
            return Err(Error::InvalidInput("working precision overflow".into()));
        }
        if self.series_order > 10000 || self.max_precision_attempts > 16 || self.workers > 256 {
            return Err(Error::Limit(
                "limits are 10000 series terms, 16 refinement attempts, and 256 workers".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default)]
pub struct CancellationToken(Arc<AtomicBool>);
impl CancellationToken {
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }
    pub fn check(&self) -> Result<()> {
        if self.0.load(Ordering::Relaxed) {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
}

#[derive(Clone, Debug)]
pub enum Progress {
    Reduction {
        integrals: usize,
    },
    Substitution {
        rules: usize,
    },
    ReductionFrontier {
        searched: usize,
        remaining: usize,
    },
    DifferentialClosure {
        round: usize,
        requested: usize,
        basis_size: usize,
        new_derivatives: usize,
    },
    SectorReduction {
        active_lines: usize,
        integrals: usize,
        visited: usize,
    },
    ParametricReduction {
        rays: usize,
        applied: usize,
        uncovered: usize,
        domains: usize,
        domain_applied: usize,
        elapsed_ms: u128,
    },
    SymmetryReduction {
        candidates: usize,
        automorphisms: usize,
        applied: usize,
        uncovered: usize,
        transport_failures: usize,
        search_limited: bool,
        elapsed_ms: u128,
    },
    SectorReduced {
        integrals: usize,
        seeds: usize,
        rows: usize,
        exact_trace_rows: usize,
        elapsed_ms: u128,
        exact_ms: u128,
    },
    BoundaryPlan {
        basis_size: usize,
        region_series: usize,
        coefficients: usize,
        max_half_order: usize,
    },
    Prepared {
        basis_size: usize,
        blocks: Vec<usize>,
    },
    Step {
        index: usize,
    },
    Sample {
        index: usize,
        total: usize,
    },
}

#[derive(Clone, Default)]
pub struct RunContext {
    pub cancellation: CancellationToken,
    pub progress: Option<Arc<dyn Fn(Progress) + Send + Sync>>,
}
impl RunContext {
    pub fn emit(&self, event: Progress) -> Result<()> {
        self.cancellation.check()?;
        if let Some(f) = &self.progress {
            f(event);
        }
        Ok(())
    }
}
