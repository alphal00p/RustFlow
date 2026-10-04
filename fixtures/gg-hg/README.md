# Supplied mixed QCD–EW Higgs-plus-jet systems

These scientific data independently serialize the complete 48-dimensional
planar and 61-dimensional nonplanar canonical systems from
`mg5_higgs_ew_plugin` commit `89f64b93d0bdbd8ee90eb85737021189229b034a`.
The original computational implementation is loaded externally for optional
reference regeneration; no original solver code is included here.

`full-systems.json` contains exact logarithmic letters, rational constant
matrices, square-root radicands, and sixteen supplied boundary/endpoint cases:
four invariant permutations for each of the W and Z mass-ratio contexts in
each family. Every case includes epsilon orders zero through four. All
coordinates are exact rational strings. `projections.json` contains exact
canonical-master weights for the four form factors and the supplied
`preCan_functions` diagnostic coefficient rows. These rows retain their
explicit epsilon indices; they do not establish a physical propagator/basis
map that is absent from the source archive.

The upstream source and endpoint files, equation exports, reference solver,
and independent serialization scripts are hash-pinned in the fixture metadata
and validation report. The serialized canonical connection is
`dI = epsilon * sum(C_k dlog(letter_k)) I`; the native
`CanonicalAlgebraicSystem` supplies the epsilon factor.

Root sheets reproduce the original threshold replacement, with an explicit
germ in every case. In particular, `sqrt(-b)` above `b=0` is opposite to the
principal root. Negative nonthreshold multivariate radicands retain the
original principal convention. Explicit infinitesimals in form-factor roots
are taken to their exact analytic side limits; no finite numerical regulator
appears in the native coefficients.

The source accuracy is capped at 24 digits and additionally carries the
recorded inherited integration error, rounded upward to a decimal envelope.
Decimal serialization keeps 100 significant digits and adds an absolute
`1e-80` reserve to each inexact input error. Long original mantissas are not
additional evidence of integration accuracy. Native returned error estimates
use `error <= 10^-digits * max(1, |coefficient|)`; for coefficients below one
this is an absolute criterion. Form-factor reports therefore also give
relative errors after projection and normalization.

The initial W and Z cases are separate regular physical contexts. Their
normalization is stated in exact `MH^2=1` units, so `MV^2=1/b` and the final
factor is `-b^2/(4*pi)^4`, applied once. These cases do not constitute a coherent
W–Z interference or HEPKit tensor/helicity contraction at one shared physical
point, nor do they exercise automatic boundary generation or IBP reduction.

The ordinary offline test runs one complete planar case. The full test runs
all 4,360 coefficients twice with higher precision/order and verifies binary
cache restart; Mathematica is not required:

```sh
cargo test --release --test gg_hg complete_planar_gg_hg_system_matches_original_and_restarts
cargo test --release --test gg_hg complete_crossed_gg_hg_systems_refine_and_match_original -- --ignored
```

Reference regeneration uses the unchanged upstream application, a compatible
Mathematica installation, the pinned `atilde`, boundary and form-factor files,
and the independently authored export/check scripts recorded in the report.
The original code remains outside the MIT implementation. Mathematical input
data, computed outputs and provenance are retained as scientific fixtures.
