//! Bounded independent whole-unit work. No internal reduction queue is changed.
use super::*;
use rayon::prelude::*;

/// Bounds the input window and the number of private-pool workers. This bounds
/// concurrent reductions, not CAS scratch or decoded coefficient heap size.
#[derive(Clone, Copy, Debug)]
pub struct GuardedUnitBatchLimits {
    pub max_workers: usize,
    pub max_in_flight: usize,
}
impl Default for GuardedUnitBatchLimits {
    fn default() -> Self {
        Self {
            max_workers: 1,
            max_in_flight: 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct GuardedUnitBatchStats {
    pub input_calls: usize,
    pub prefetched_native_calls: usize,
    pub fallback_native_calls: usize,
    pub native_calls_returning_error: usize,
    pub cache_hits: usize,
    /// Includes misses whose reduction or insertion returned an error.
    pub cache_misses: usize,
    pub private_pool_workers: usize,
    pub simultaneous_native_calls_bound: usize,
    /// Known application counts from every native Ok result, whether or not
    /// that result was eventually returned to the caller. Native Err currently
    /// carries no partial application count; the boolean below exposes this.
    pub actual_known_rule_applications: u128,
    pub actual_application_count_complete: bool,
    pub returned_logical_rule_applications: u128,
    pub returned_memoized_rule_applications: u128,
    pub returned_uncached_rule_applications: u128,
    pub unused_prefetch_known_rule_applications: u128,
}

/// This explicit batch is not fail-fast. Every item has the same logical unit
/// result/limits as an ordered reduce_memoized call, including native failures
/// and codec errors. The caller must inspect every item. Do not confuse the
/// executed batch with a sequential caller's shorter prefix before an error.
#[derive(Debug)]
pub struct GuardedUnitBatch<const N: usize> {
    pub items: Vec<Result<GuardedMemoizedReduction<N>, SolverError>>,
    /// Ordered snapshots after each item's cache action, including errors.
    /// This allows a caller to retain its consumed-prefix peak semantics.
    pub memo_usage_after: Vec<Option<GuardedUnitMemoUsage>>,
    pub stats: GuardedUnitBatchStats,
}
struct Job<const N: usize> {
    target: [i64; N],
    result: Option<Result<GuardedReduction<N>, SolverError>>,
}

impl<const N: usize> GuardedProgram<N> {
    /// Borrow one immutable verified program. The coordinator owns the memo
    /// lock; workers use only plain reduce and never acquire that lock. Cache
    /// hit/recency/eviction/insertion occur in the original input order.
    ///
    /// No entry or rule is cloned to construct worker programs. Plain reduce's
    /// ordering, conditions, application/pending limits and failures are intact.
    /// Pool construction errors happen before native work or cache mutations.
    /// There is no cancellation poll inside a native unit call, matching reduce.
    pub fn reduce_batch_memoized(
        &self,
        targets: &[[i64; N]],
        limits: GuardedReductionLimits,
        batch: GuardedUnitBatchLimits,
    ) -> Result<GuardedUnitBatch<N>, SolverError> {
        if batch.max_workers == 0
            || batch.max_in_flight == 0
            || batch.max_workers > batch.max_in_flight
            || targets.len() > batch.max_in_flight
        {
            return Err(invalid(
                "invalid or exhausted independent unit batch bounds",
            ));
        }
        let mut memo = self
            .unit_memo
            .as_ref()
            .map(|m| {
                m.lock()
                    .map_err(|_| invalid("native unit memo lock poisoned"))
            })
            .transpose()?;
        let mut jobs = Vec::<Job<N>>::new();
        let mut indices = Vec::<Option<usize>>::with_capacity(targets.len());
        let mut first_missing = BTreeMap::new();
        for target in targets {
            let present = memo.as_ref().is_some_and(|m| {
                m.slots.iter().any(|entry| {
                    entry
                        .as_ref()
                        .is_some_and(|e| e.target == *target && same_limits(e.limits, limits))
                })
            });
            if present {
                indices.push(None);
                continue;
            }
            // Only collapse planned duplicate misses. The actual ordered cache
            // lookup still decides every call. If the first result cannot be
            // retained or is evicted, a later miss executes again, as serial does.
            if memo.is_some() {
                if let Some(&index) = first_missing.get(target) {
                    indices.push(Some(index));
                    continue;
                }
            }
            let index = jobs.len();
            first_missing.insert(*target, index);
            indices.push(Some(index));
            jobs.push(Job {
                target: *target,
                result: None,
            });
        }
        let workers = batch.max_workers.min(jobs.len());
        let pool = if workers > 1 {
            Some(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .map_err(|error| {
                        invalid(&format!("cannot build bounded native unit pool: {error}"))
                    })?,
            )
        } else {
            None
        };
        if let Some(pool) = &pool {
            pool.install(|| {
                jobs.par_iter_mut().for_each(|job| {
                    job.result = Some(self.reduce(job.target, limits));
                })
            });
        } else {
            for job in &mut jobs {
                job.result = Some(self.reduce(job.target, limits));
            }
        }
        let mut stats = GuardedUnitBatchStats {
            input_calls: targets.len(),
            prefetched_native_calls: jobs.len(),
            private_pool_workers: if pool.is_some() { workers } else { 0 },
            simultaneous_native_calls_bound: workers,
            actual_application_count_complete: true,
            ..Default::default()
        };
        for job in &jobs {
            match job.result.as_ref().expect("each bounded job executed") {
                Ok(result) => {
                    stats.actual_known_rule_applications += result.rule_applications as u128
                }
                Err(_) => {
                    stats.native_calls_returning_error += 1;
                    stats.actual_application_count_complete = false;
                }
            }
        }
        let mut items = Vec::with_capacity(targets.len());
        let mut memo_usage_after = Vec::with_capacity(targets.len());
        for (target, index) in targets.iter().zip(indices) {
            let hit = match memo.as_mut() {
                Some(m) => m.hit(*target, limits),
                None => Ok(None),
            };
            let result = match hit {
                Err(error) => Err(error), // Owned codec errors never become misses.
                Ok(Some(reduction)) => {
                    stats.cache_hits += 1;
                    Ok(GuardedMemoizedReduction {
                        reduction,
                        cache_hit: true,
                    })
                }
                Ok(None) => {
                    if memo.is_some() {
                        stats.cache_misses += 1;
                    }
                    let staged = index.and_then(|index| jobs[index].result.take());
                    let reduced = match staged {
                        Some(reduced) => reduced,
                        None => {
                            // Snapshot hit evicted by an earlier ordered insert,
                            // or duplicate result was not retained. Preserve the
                            // serial native call instead of an unbounded side cache.
                            stats.fallback_native_calls += 1;
                            stats.simultaneous_native_calls_bound =
                                stats.simultaneous_native_calls_bound.max(1);
                            let r = self.reduce(*target, limits);
                            match &r {
                                Ok(r) => {
                                    stats.actual_known_rule_applications +=
                                        r.rule_applications as u128
                                }
                                Err(_) => {
                                    stats.native_calls_returning_error += 1;
                                    stats.actual_application_count_complete = false;
                                }
                            }
                            r
                        }
                    };
                    reduced.and_then(|reduction| {
                        if let Some(m) = memo.as_mut() {
                            m.insert(*target, limits, &reduction)?;
                        }
                        Ok(GuardedMemoizedReduction {
                            reduction,
                            cache_hit: false,
                        })
                    })
                }
            };
            if let Ok(result) = &result {
                let n = result.reduction.rule_applications as u128;
                stats.returned_logical_rule_applications += n;
                if result.cache_hit {
                    stats.returned_memoized_rule_applications += n;
                } else {
                    stats.returned_uncached_rule_applications += n;
                }
            }
            memo_usage_after.push(memo.as_ref().map(|m| m.usage()));
            items.push(result);
        }
        for job in jobs {
            if let Some(Ok(result)) = job.result {
                stats.unused_prefetch_known_rule_applications += result.rule_applications as u128;
            }
        }
        Ok(GuardedUnitBatch {
            items,
            memo_usage_after,
            stats,
        })
    }
}
