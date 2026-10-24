use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    InvalidInput(String),
    #[error("unsupported input: {0}")]
    Unsupported(String),
    #[error("IBP reduction failed: {0}")]
    Reduction(String),
    #[error("incomplete reduction or boundary data: {0}")]
    IncompleteReduction(String),
    #[error("numerical failure: {0}")]
    Numerical(String),
    #[error("resource limit: {0}")]
    Limit(String),
    #[error("{0}")]
    Resource(#[from] symbolica::utils::resources::ResourceLimitError),
    #[error("requested accuracy was not reached: {0}")]
    Accuracy(String),
    /// Working precision cannot resolve endpoint arithmetic conditioning.
    /// Re-evaluation must use fresh source/boundary data at the suggested bits.
    #[error("insufficient working precision (at least {minimum_bits} bits): {context}")]
    InsufficientPrecision { minimum_bits: u32, context: String },
    #[error("calculation cancelled")]
    Cancelled,
    #[error("cache error: {0}")]
    Cache(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}
pub type Result<T> = std::result::Result<T, Error>;

/// Preserve shared resource causes before adding a consumer's contextual message.
pub(crate) fn preserve_resource(error: &(dyn std::error::Error + 'static)) -> Option<Error> {
    let mut cause = Some(error);
    while let Some(error) = cause {
        if let Some(resource) =
            error.downcast_ref::<symbolica::utils::resources::ResourceLimitError>()
        {
            return Some(Error::Resource(resource.clone()));
        }
        if let Some(Error::Resource(resource)) = error.downcast_ref::<Error>() {
            return Some(Error::Resource(resource.clone()));
        }
        cause = error.source();
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as StdError;

    #[test]
    fn resource_cause_remains_visible_to_standard_error_walkers() {
        let original = symbolica::utils::resources::ResourceLimitError::Memory {
            rss_bytes: 11,
            additional_bytes: 2,
            limit_bytes: 10,
        };
        let error = Error::Resource(original.clone());
        let cause = error
            .source()
            .expect("resource wrapper must expose its cause");
        assert_eq!(
            cause.downcast_ref::<symbolica::utils::resources::ResourceLimitError>(),
            Some(&original)
        );
    }
}
