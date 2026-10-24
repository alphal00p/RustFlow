//! Grow native recursion stacks without replacing the caller's resource scope.

use crate::{Error, ResourceLimitError, Result};
use symbolica::utils::resources;

// Stack allocation granularity and headroom, not limits on recursive work.
const RED_ZONE: usize = 256 * 1024;
const STACK_BYTES: usize = 2 * 1024 * 1024;

pub(crate) fn run<T>(operation: impl FnOnce() -> Result<T>) -> Result<T> {
    resources::checkpoint()?;
    let result = if stacker::remaining_stack().is_some_and(|remaining| remaining >= RED_ZONE) {
        operation()
    } else {
        resources::check_allocation(STACK_BYTES as u64)?;
        stacker::maybe_grow(RED_ZONE, STACK_BYTES, || {
            // Stacker's unsupported-target backend can call through without
            // switching stacks. Verify headroom instead of claiming it grew.
            if stacker::remaining_stack().is_none_or(|remaining| remaining < RED_ZONE) {
                return Err(Error::Resource(ResourceLimitError::MonitorUnavailable {
                    resource: "call stack",
                    message:
                        "native stack growth cannot establish recursive headroom on this target"
                            .into(),
                }));
            }
            resources::checkpoint()?;
            operation()
        })
    };
    let checkpoint = resources::checkpoint();
    match result {
        Err(error) => Err(error),
        Ok(value) => {
            checkpoint?;
            Ok(value)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ResourceBudget, ResourceLimits};

    #[test]
    fn an_exhausted_outer_budget_refuses_work_before_entry() {
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1),
                ..Default::default()
            },
            || Ok(2),
        );
        assert!(matches!(
            budget.run(|| run(|| -> Result<()> { panic!("work must not start") })),
            Err(Error::Resource(ResourceLimitError::Memory { .. }))
        ));
    }

    #[test]
    fn stack_growth_keeps_the_thread_and_resource_context() -> Result<()> {
        fn descend(depth: usize, thread: std::thread::ThreadId) -> Result<()> {
            run(|| {
                let live_frame = [0_u8; 8192];
                assert_eq!(std::thread::current().id(), thread);
                resources::checkpoint()?;
                if depth != 0 {
                    descend(depth - 1, thread)?;
                }
                std::hint::black_box(&live_frame);
                Ok(())
            })
        }
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(u64::MAX),
                ..Default::default()
            },
            || Ok(0),
        );
        budget.run(|| descend(384, std::thread::current().id()))
    }
}
