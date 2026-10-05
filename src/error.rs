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
