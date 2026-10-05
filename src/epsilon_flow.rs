//! Epsilon-rescaled physical flows with authenticated boundary conversion.
use crate::diffexp::EpsilonBoundary;
use crate::transport_cache::{BoundaryIdentity, CachedBoundary, EpsilonRange};
use crate::{EpsilonShearing, Error, Precision, Result, RustFlow};

/// A rational physical connection regularized by exact diagonal epsilon powers.
/// The original and transformed cache identities remain distinct. Numerical
/// transport and cache selection use the ordinary shared [`RustFlow`] engine.
pub struct EpsilonShearedFlow {
    pub(crate) original_identity: BoundaryIdentity,
    pub(crate) flow: RustFlow,
    pub(crate) shearing: EpsilonShearing,
}

impl EpsilonShearedFlow {
    pub fn flow(&self) -> &RustFlow {
        &self.flow
    }
    pub fn shearing(&self) -> &EpsilonShearing {
        &self.shearing
    }
    pub fn original_identity(&self) -> &BoundaryIdentity {
        &self.original_identity
    }

    /// Rectangular transformed range sufficient to restore every original
    /// component through the requested last order. The common leading bound
    /// must be justified for the original solution, not inferred from a
    /// boundary coefficient that happens to vanish at one point.
    pub fn required_range(&self, original: EpsilonRange) -> Result<EpsilonRange> {
        let (leading, last) = self
            .shearing
            .required_sheared_range(original.leading, original.last)?;
        EpsilonRange::new(leading, last)
    }

    fn source_coefficients(
        &self,
        boundary: &CachedBoundary,
        identity: &BoundaryIdentity,
    ) -> Result<EpsilonBoundary> {
        boundary.validate()?;
        if boundary.identity.key() != identity.key() {
            return Err(Error::InvalidInput(
                "epsilon boundary conversion requires the matching declared basis identity".into(),
            ));
        }
        // Epsilon reindexing does not use a path coordinate. Physical point,
        // branch, kind and original error provenance are copied separately.
        Ok(EpsilonBoundary {
            point: Precision {
                bits: boundary.accuracy.working_bits(),
            }
            .zero(),
            leading: boundary.range.leading,
            coefficients: boundary.coefficients.clone(),
        })
    }

    /// Convert a boundary belonging to the original system into the exact
    /// rescaled basis. No cache is mutated. Missing positive source orders are
    /// rejected, and all coefficient errors and evidence caps are retained.
    pub fn to_sheared_boundary(
        &self,
        source: &CachedBoundary,
        last: i32,
    ) -> Result<CachedBoundary> {
        let original = self.source_coefficients(source, &self.original_identity)?;
        let (shifted, errors) =
            self.shearing
                .to_sheared(&original, source.accuracy.comparison_errors(), last)?;
        let result = CachedBoundary {
            identity: self.flow.identity().clone(),
            point: source.point.clone(),
            kind: source.kind,
            range: EpsilonRange::new(shifted.leading, last)?,
            coefficients: shifted.coefficients,
            accuracy: source.accuracy.reindexed(
                errors,
                "Exact diagonal epsilon rescaling of the original boundary",
            )?,
        };
        result.validate()?;
        Ok(result)
    }

    /// Restore an explicitly requested original range. This neither increases
    /// the source's accuracy cap nor supplies omitted epsilon coefficients.
    pub fn to_original_boundary(
        &self,
        source: &CachedBoundary,
        range: EpsilonRange,
    ) -> Result<CachedBoundary> {
        let sheared = self.source_coefficients(source, self.flow.identity())?;
        let (restored, errors) = self.shearing.to_original(
            &sheared,
            source.accuracy.comparison_errors(),
            range.leading,
            range.last,
        )?;
        let result = CachedBoundary {
            identity: self.original_identity.clone(),
            point: source.point.clone(),
            kind: source.kind,
            range,
            coefficients: restored.coefficients,
            accuracy: source
                .accuracy
                .reindexed(errors, "Exact restoration of the original epsilon basis")?,
        };
        result.validate()?;
        Ok(result)
    }
}
