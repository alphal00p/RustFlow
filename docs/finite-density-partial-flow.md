# Partial fixed-shell occupied flows

The direct `PreparedOccupiedFlow` Rust API accepts an explicit
`MassMode::Propagators` mask for a separately certified class with one virtual
loop and two occupied massless shells. Its proof and boundary wiring passed 50
tests in one combined crate. This is not yet a closed native-flow or numerical
amplitude validation. The complete `PreparedDensityFlow` interface retains its
All/Auto policy.

The deformation is the existing fixed-shell one: native uncut denominators
`D_i-eta` on the selected physical slots, or Euclidean `rho_i+eta`. Cut shells,
positive-energy support, occupied spatial regions, physical-point numerator
indices and coefficients stay fixed. The selected slots are canonicalized and
must be distinct uncut physical denominators. A cut or completion slot cannot
be shifted. The evaluated mask must equal the prepared mask.

## Structural admission

All physical masses must be zero, occupied chemical magnitudes positive, and
completion factors polynomial. The constructor reconstructs the exact cut
geometry from the original graph and routing. One common rational translation
of virtual momenta must remove compact offsets from all unshifted virtual
denominators. Every pure compact physical denominator must be shifted and be a
square of a rational multiple of a compact momentum difference. The complete
physical momentum rows, including factors absent from a particular target, are
part of the certificate.

For every positive-denominator support, the proof uses native off-shell Gaussian
data and exact parameter charts. Independent positive regulators precede null
shell specialization. Rank-deficient supports require a genuine unrestricted
virtual null direction. Full-rank supports retain the actual regulated quadratic
form and nonnegative compact invariants. Unsupported signs, residual invariants,
translations or exhausted proof budgets cause an explicit error. There is no
graph-name or benchmark-ID dispatch. The underlying algebraic proof owners have
a wider rank-at-most-two scope, but that does not grant wider numerical admission.

## Identities, endpoints and domains

Weighted sources retain the cut and occupation distributions, including upper
Fermi-surface derivatives. Nonpositive cut indices and zero occupation indices
keep their different meanings. Lower-origin zero jets are continued from a
sufficiently large real dimension for each finite label; no single dimension is
claimed for an unbounded symbolic index ray. Optional virtual-zero boxes come
from the bound origin proof. Partial placements grant no singleton raw-Ward
identity.

Before discovery, an independent class checker inspects every actual emitted
source row. Its domain has nonpositive completion indices, nonnegative
occupation indices and zero storage tails. A source image that crosses an
integer boundary is permitted only when its exact combined coefficient vanishes
there. All fixed guard coordinates are specialized, and raw conditions survive
row cancellation. The certificate binds source rows, maps, guards, conditions and
zero boxes. Native seed, recentering and replay preserve this class; the checker
does not independently establish the physical validity of the identities.

Original targets, the final basis, target and candidate rows, and raw
eta-derivative labels receive finite degree and distribution-jet audits. Raw
numerator and coefficient poles are captured before mass substitution or cut
zeros. Matrix, target and candidate denominator conditions are also attached to
the actual `ReducedSystem` before its audit signature is sealed. Transport checks
these conditions even for an empty basis.

The endpoint argument concerns the complete physical target. Native weighted
reconstruction and the actual convergent regular-singular Frobenius expansion
remain required. Rational weights are multiplied into the power/log series
before the shared endpoint projector is applied. Terms whose exponents differ
by integers are combined before continuity is used. On an open high-dimensional
domain, continuity removes any surviving divergent or logarithmic contribution;
the result is then continued meromorphically in dimension. An integer shift from
a rational target weight cannot turn a dimension-dependent exponent into an
identically constant zero exponent. Insufficient series order, an uncancelled
divergence or a nongeneric sampled indicial system remains an error. Transport
is restricted to the certified right half of the eta plane.

## Large-mass boundary

The ordinary region owner enumerates compact-preserving regions. All-hard
virtual factors use ordinary recursive AMF with bubble shortcuts disabled and
the tadpole-only terminal policy. Explicit parent slot numbers are not forwarded
as slot numbers of ordinary children; those children use `MassMode::All`.

A region retaining the virtual soft loop must obtain an exact native ordinary
zero certificate after its common translation, or an unrestricted-polynomial
null-direction proof. The certificate binds the actual expression, coordinate
map, distributions and parent continuation identity. This does not turn an
arbitrary nonpolynomial compact factor into a moment. Raw parameter poles are
captured in the complete parent coordinates before projection or cancellation,
checked at the sample, and retained in provenance. A removable compact-only pole
still has no ordinary-zero authority. Other soft factors require weighted
recursion and are rejected.

The retained infinity exponents remain symbolic in epsilon. Ordinary hard seeds
and certified virtual zeros supply the needed coefficients to the existing
Frobenius matching owner; no reference value supplies a boundary.

## Validation and use

For the example `massless_three_loop_chain.json`, the direct cut `[0,3]` admits
masks `[1,2]` and `[2,4]` structurally. These are examples of the tests, not special
cases in the implementation. The runtime occupied-flow harness accepts
`RUSTFLOW_WEIGHTED_SHIFTED_SLOTS='[1,2]'`; its usual input, cut, work-budget and
report settings remain required. This option does not change the complete
assembly interface.

The frozen [combined-crate report](../reports/validation/2026-10-10-partial-flow-wiring-draft/README.md)
records type-checking and all 50 selected tests, including actual hard/soft
region integration. The [continuation bundle report](../reports/validation/2026-10-10-partial-continuation-integration-draft/README.md)
records the underlying proof checks and their limitations. Work limits cover
retained expressions, supports, charts, labels and source images; they do not
bound scratch allocation or interrupt one internal CAS operation. Numerical
closure, refinement and independent-reference gates remain separate.
