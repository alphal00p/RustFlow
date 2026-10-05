# Published gg → Hg inputs and the automatic seed gate

The exact map between ordinary physical integrals and the plugin's 48-planar
and 61-nonplanar canonical masters is now available. The correct primary source
is [arXiv:2112.07578v1](https://arxiv.org/abs/2112.07578v1), whose ancillary files
provide both physical-power-vector replacement rules and the canonical dlog
matrices used by the plugin. This resolves the earlier ambiguity with the
separate 64-component basis in the 2020 paper.

Exact basis equivalence is a prerequisite, not a numerical AMF seed validation.
The full native seed acceptance still requires independent Euclidean checks
and production seed evaluation; no upstream numerical boundary is substituted.

## Plugin 48/61 physical-to-canonical map

`PluginBasisMap::load(kind, namespace)` parses all coefficients as exact
Symbolica atoms in a caller-selected namespace. It exposes `s`, `t`, `b`,
`eps`, and the same `root1` through `root8` names used by the native canonical
systems. `family()` constructs an ordinary RustFlow family from the published
routing, while `canonical_rows()` supplies one exact projection per master.

The scalar variables are `s=(p1+p2)^2`, `t=(p1-p3)^2`, and
`b=(p1+p2-p3)^2`, with three massless external momenta. This is the scattering
convention, distinct from the decay convention of `PublishedFamily` below.
The source variables are specialized exactly by `s12=s`, `s13=-t`,
`s23=s+t-b`, `mZ=1`. The ancillary variable `mZ` denotes the squared boson mass.
All source signed powers are retained, including the two numerator slots.
The paper uses `q²-m²` denominators, so no denominator-sign conversion occurs.

| Map | Physical requests | Forward nonzeros | Inverse nonzeros | Roots |
| --- | ---: | ---: | ---: | ---: |
| Planar | 48 | 157 | 201 | 2 |
| Nonplanar | 61 | 360 | 453 | 8 |

Writing the source rule as `J = M F`, the native canonical result is

```text
F = M^(-1) exp(2 epsilon EulerGamma) I,
```

where `I` uses the ordinary RustFlow measure. The mass scale and renormalization
scale are both one here, as in the plugin's dimensionless systems.
`factors(root_germ)` supplies the exact common normalization to
`solve_integral_projections_normalized`. Named roots and the EulerGamma factor
are evaluated separately at every finite epsilon and working precision before
Laurent fitting. Root sheets are explicit input; the map does not infer a
continuation path or turn an approximate root into an exact coefficient.

The extraction verified exact reconstruction of every primary replacement
rule, exact `M M^(-1)=1`, and equality of every published dlog-matrix entry to
the corresponding plugin entry after `mmH=b`. These are symbolic comparisons,
not numerical fits or matches of index counts. `verify_inverse()` independently
checks the inverse using the native Symbolica matrix implementation and exact
polynomial cancellation. Admitting this certified mathematical input does not
claim successful numerical seeds.

## Exact physical families

`gg_hg::PublishedFamily::new` constructs ordinary RustFlow `IntegralFamily`
objects using the existing `Propagator::quadratic` implementation and native
RustRed validation. Its arguments are exact Symbolica atoms for
`s=(p1+p2)^2`, `t=(p1+p3)^2`, `u=(p2+p3)^2`, and the squared vector-boson mass.
The three external momenta are massless and the Higgs virtuality is `s+t+u`.
There are two loops, seven physical denominators, and two numerator slots.
No new graph, momentum-basis, or tensor algorithm is introduced.

| Kind | Definition | Stored integral requests |
| --- | --- | --- |
| `Planar2018` | arXiv:1810.05138v1, section 2 | Appendix A, T1 through T48, in that order |
| `Planar2020` | arXiv:2007.09813v2, table 1 | 45 distinct planar requests from appendix A |
| `Nonplanar2020` | arXiv:2007.09813v2, table 1 | 18 nonplanar requests and 3 additional DE requests from appendix A |

The last count is an input list, not a claim of 21 linearly independent masters.
Reduction and differential closure still belong to the native reducer.

All constructor outputs use `q²-m²` denominators and the ordinary measure
`d^Dk/(i π^(D/2))`. Publication conventions remain explicit metadata:

- The 2018 paper uses the negative of each stored denominator. Its exact
  `denominator_factor` is `(-1)^sum(powers)`, including signed numerator powers.
  Its additional per-loop measure is
  `exp(epsilon EulerGamma) (mB²/mu²)^epsilon`.
- The 2020 paper uses the stored denominator sign and an additional
  `1/Gamma(1+epsilon)` per loop.

Callers must apply those factors when comparing published values. The family
constructor does not silently change RustFlow's integral convention.

## Independently serialized published canonical data

`PublishedCanonicalBasis::load` reads the 64-component vector in the 2020
ancillary file. It retains 64 distinct integral references, 191 nonzero matrix
entries, and seven separately named square roots. The extraction checks exact
linear reconstruction of every source expression before serialization. The
loader uses Symbolica's exact rational polynomial conversion to reject weights
outside the declared scalar variables and root generators.

The scalar variables are `sh=mh²`, `y=-t/mh²`, `z=-u/mh²`, and
`rho=-mV²/mh²`. The common factor relative to ordinary two-loop integrals is

```text
(-mh²)^(2 epsilon) / Gamma(1+epsilon)^2.
```

`PublishedNormalization::evaluate` evaluates that factor at each exact epsilon
sample using existing Symbolica/MPFR operations. For positive real `mh²`, the
physical `+i0` convention uses `log(-mh²-i0)`. Its opposite prescription is the
complex conjugate. Zero Higgs virtuality is rejected. This helper deliberately
accepts real rational virtualities; it does not infer a complex homotopy.

The source includes `PLx12` and `PLx123` integral heads. They are retained as
distinct typed labels. Their names alone do not prove a permutation of physical
momenta or of invariant variables. `require_resolved_crossings()` rejects these
unresolved references. In particular, a naive swap of the first two external
momenta must not be inferred from the text `x12`.

## Generic physical form factors

`HiggsJetFormFactors::load(namespace)` loads the four finite scalar form factors
as functions of the exact dimensionless invariants `s`, `t`, `b`. The data
contain 2,470 grouped canonical-master coefficients and fourteen separately
registered coefficient roots. They are not restricted to a stored grid point.
All four expressions were checked symbolically against the primary 2112
ancillary expressions after the variable rename and removal of explicit
infinitesimals. Every grouped term was also reconstructed exactly, including
its root factors and the original infinitesimal before taking a side limit.

`evaluate` accepts exact physical `s`, `t`, `mH²`, `mV²` and four cached planar
plus four cached nonplanar boundaries. Their order is `(s,u,b)`, `(s,t,b)`,
`(u,t,b)`, `(t,s,b)`, with `u=b-s-t` in boson-mass units. It checks each boundary's
canonical identity, exact endpoint, epsilon range and physical root germ.
Wrong crossings, basis identities or root sheets are rejected. The final
physical normalization is `-1/(mV²)^2/(4*pi)^4`, applied once after projection.

The coefficient roots retain their exact radicand and the coefficient of the
positive imaginary infinitesimal. A negative radicand approached from below
therefore has a negative imaginary square root; a root without a prescribed
infinitesimal keeps the principal convention. No finite numerical regulator is
introduced. This physical form-factor interface currently accepts exact real
rational inputs and positive squared Higgs/vector masses; it does not infer
complex continuation paths for the observable.

Absolute master-component errors are propagated with the modulus of each
weight and the common normalization, with an additional arithmetic reserve.
Relative digits are assessed for each resulting form factor and capped by the
input evidence. An unresolved zero receives no relative-digit claim. These
remain empirical bounds conditional on the supplied master errors, rather
than outward-rounded interval enclosures. Tests compare all generic weights
at the independent W and Z point specializations and check root side limits.

The production fixture `fixtures/gg-hg/generic-form-factors.json` contains only
exact expressions and source hashes. Existing numerical projection files are
read by tests only, never as production boundaries or coefficient inputs.

## Remaining automatic-seed acceptance

The plugin's `preCan_functions_*` records alone do not identify physical
integrals; the 2112 ancillary replacement rules supply that missing scientific
information. The separate 2020 64-component vector remains a useful independent
input and must not be substituted for the plugin's 61-vector.

Acceptance requires independent AMF checks at the analytic Euclidean anchors.
Production seeds then come from native AMF at each of the 16 existing physical
source cases; the checked local physical paths transport those values. The
original numerical grid remains comparison data only. Reduction-search,
recursive-boundary and precision limits must return typed failures, never a
cached success. An exact basis certificate alone supplies no accuracy estimate
or performance claim for this numerical workflow.

Module tests exercise family signs/masses/numerators, source-normalization
branches, exact data loading and native exact inverse checks. Automatic full
48/61 seed evaluation is a separate acceptance gate.

## Python notebook and long acceptance

The community notebook is `examples/hep/gg_hg.py`. Its controller owns one
native HEPKit model, a seed-only boundary bank, and a progressively growing
transport bank. Loading seeds merges them into the transport bank; it preserves
existing intermediate points. Forced recomputation bypasses completed numerical
samples and verified boundaries while retaining reusable exact reductions.
Only completed raw sample vectors are checkpointed, separately from boundaries
that passed Laurent and precision checks.

Opening the Marimo notebook does not run the two-loop boundary calculation.
With the shared community extension built in release mode, run the separate
long acceptance from the community checkout:

```sh
python examples/hep/gg_hg_acceptance.py --directory /path/to/empty-cache --interrupt-after-samples 1
```

Use `--resume` with the same directory after interruption. The runner compares
4,360 coherent reference coefficients, checks the observable reference allowances,
tests binary reload and exact repeated hits, forces fresh boundaries at higher
precision, and evaluates nearby kinematics using accumulated points. The optional
interruption probe cancels after complete sample files exist, verifies their
payloads survive, and creates a fresh controller before resuming. It measures
cancellation latency; reducer calls can delay cancellation until a batch ends.
Cold, resumed, warm, forced, and nearby stages are timed separately in
`acceptance.json`. References are loaded only after native results exist.

The amplitude component test uses archived transported form factors and checks
three observables against independent original-code calculations to 30 relative
digits. Its input allowances support 19 physical digits for the EW square and
20 for interference. This component test does not satisfy the empty-cache native
boundary gate. The notebook refines native seed precision if propagated
observable uncertainty falls short of its requested target.

## Provenance

Mathematical definitions are independently serialized in
`fixtures/gg-hg/published-physical-families.json` and
`fixtures/gg-hg/published-canonical-64.json`, and
`fixtures/gg-hg/plugin-physical-map.json`. No upstream integration or reduction
implementation is included. The primary papers and their authors retain
attribution:

- Matteo Becchetti, Francesco Moriello and Armin Schweitzer,
  [Two-loop amplitude for mixed QCD-EW corrections to gg → Hg](https://arxiv.org/abs/2112.07578v1).
- Matteo Becchetti, Roberto Bonciani, Valerio Casconi, Vittorio Del Duca and
  Francesco Moriello,
  [Planar master integrals for the two-loop light-fermion electroweak corrections
  to Higgs plus jet production](https://arxiv.org/abs/1810.05138v1).
- Marco Bonetti, Erik Panzer, Vladimir A. Smirnov and Lorenzo Tancredi,
  [Two-loop mixed QCD-EW corrections to gg → Hg](https://arxiv.org/abs/2007.09813v2).

Downloaded primary source archive SHA-256 values:

```text
1810.05138v1: 9749d594218e7eb1052057021b527a5be0b4655d8138bd862aec361ea6e441a8
2007.09813v2: 82d7c206ba734b21ab61de1eef8763b085fb41e5f684f3196df7060beb47d922
2112.07578v1: 7026b05b0d8950c052a4642c3985cf03f71b838d57e25913a4f15a81fa823eb6
```

The 2020 `canonical_basis_NP_system.m` SHA-256 is
`cf673b527f1fcd25fc2df9f6c04998d5027b0909cd9d64b32e758d86eabbf39c`.
The 2018 `CanonicalBasisHggg.mx` SHA-256 is
`0cf655d17ed5696d633ff7bd6710fced7c5f7df7c5d31969f05fcc9e4c47c847`;
the accompanying `CanonicalBasisHgggNote.txt` SHA-256 is
`603a3e71506beca563339b2f7160b5beff3cf4a2f1960eaf99b84db9ac9d33e7`.
Raw downloads and extraction diagnostics are retained locally under
`target/gg-hg-basis-work`.

The 2112 source was retrieved from `https://export.arxiv.org/src/2112.07578v1`.
Replacement rules have SHA-256
`a6bd80f1edeb056f89ee1613561be3066bb3b3ea1fe7f4ab65b928891fa74380`
(planar) and
`4a65e20f6244405152a8d87290daca63bb3a80aca7f63e024efa54cb9d8dcccd`
(nonplanar). The independently serialized map also records source routing,
published dlog and plugin dlog hashes, and all exact extraction assertions.
