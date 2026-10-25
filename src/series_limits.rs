//! Conservative allocation estimates for long numerical series.
//! Shared resource policies are optional; representation checks always apply.
use crate::{Error, Precision, Result};

pub(crate) fn complex_bytes(p: Precision) -> u128 {
    2 * (u128::from(p.bits).div_ceil(8) + 64)
}

pub(crate) fn admit_bytes(order: usize, bytes: u128, operation: &str) -> Result<()> {
    order
        .checked_add(1)
        .ok_or_else(|| Error::Limit(format!("{operation} series length overflow")))?;
    i64::try_from(order).map_err(|_| {
        Error::Limit(format!(
            "{operation} derivative index exceeds i64 representation"
        ))
    })?;
    let bytes = u64::try_from(bytes).map_err(|_| {
        Error::Limit(format!(
            "{operation} allocation estimate exceeds u64 representation"
        ))
    })?;
    symbolica::utils::resources::check_allocation(bytes)?;
    Ok(())
}

pub(crate) fn admit_complex(
    p: Precision,
    order: usize,
    values: u128,
    operation: &str,
) -> Result<()> {
    admit_bytes(order, values.saturating_mul(complex_bytes(p)), operation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn admission_uses_optional_shared_resources_and_checked_representation() {
        let p = Precision::decimal(24_000).unwrap();
        assert!(admit_complex(p, 27_000, 12 * 4 * 27_001, "Taylor chart").is_ok());
        assert!(admit_complex(p, 200_001, 1, "Taylor chart").is_ok());
        assert!(
            admit_complex(Precision { bits: u32::MAX }, 10_001, 10_002, "Taylor chart").is_ok()
        );
        assert!(matches!(
            admit_complex(p, 27_000, u128::MAX, "Taylor chart"),
            Err(Error::Limit(_))
        ));
        assert!(matches!(
            admit_complex(p, usize::MAX, 1, "Taylor chart"),
            Err(Error::Limit(_))
        ));
        let budget = crate::ResourceBudget::with_memory_monitor(
            crate::ResourceLimits {
                max_process_rss_bytes: Some(1024),
                ..Default::default()
            },
            || Ok(0),
        );
        assert!(matches!(
            budget.run(|| admit_complex(p, 27_000, 12 * 4 * 27_001, "Taylor chart")),
            Err(Error::Resource(crate::ResourceLimitError::Memory { .. }))
        ));
    }
}
