# Three-loop boundary and transport audit

This is a read-only source and saved-proof audit of the optional
`PolynomialClosure` integration. It contains no numerical predictions, no new
production rules, and no imported boundary values. The coherent release and
full-amplitude runtime gates remain separate.

## Source-equivalence gate

`tests/finite_density_polynomial_sources.rs` was independently reviewed before
source freeze. The test reads the historical v2 metadata only to reconstruct
the source rows, roles, guards, coefficient variables and zero domains. It
discards all 19 saved rules. The public factory uses a newly bound origin
permit and the new policy identity; it never imports the historical program
as reduction evidence. All 52 legacy rows must agree exactly after explicit
index renaming. The ten additions may differ by one nonzero rational constant
per row, with every retained condition checked in both directions up to
nonzero rational associates and constant conditions. The empty-program
roundtrip tests the new source context, not a claim of closed reduction.
No API or mathematical-contract blocker was found in this test review.

## Actual closed cut-[0] system

The evidence is the completed unchanged-native baseline result at
`../2026-10-10-native-source-portfolio-experiment/active-resume/polynomial62-baseline/result.json`.
Its seven basis labels, in physical slot order followed by zero completion,
occupation and storage-tail indices, are:

| Basis | Positive physical slots | Diagonal derivative coefficient |
| --- | --- | --- |
| 0 | 0,3,4 | (D-2)/eta |
| 1 | 0,2,4 | (D-2)/eta |
| 2 | 0,1,4 | (D-2)/eta |
| 3 | 0,1,3 | (D-2)/eta |
| 4 | 0,1,3,4 | (D-3)/eta |
| 5 | 0,1,2 | (D-2)/eta |
| 6 | 0,1,2,4 | (D-3)/eta |

Every positive power is one. The occupied slot 0 is C1; both occupation
indices and all completion indices are zero. The saved two target rows and
all seven actual weighted derivative rows have been replayed. Target
coefficients contain poles through eta^-3 and dimension-dependent conditions.
The diagonal Euler system does not justify dropping those target weights or
substituting a numerical dimension into the indicial classes.

This is one standalone closure, not proof that the public preparation will
necessarily select the same basis under other search budgets. Its full
symbolic target conditions must also survive the public adapter.

## Exact hard families expected from the production expansion

Use the actual input routing
`[q, k, l, q-l, k-l]`, where q is the occupied momentum. In the all-virtual-hard
region set k=sqrt(eta) K and l=sqrt(eta) L, keeping q fixed and soft. The
leading native ordinary factors are

```
slot 1: A = K^2 - 1
slot 2: B = L^2 - 1
slot 3: B = L^2 - 1
slot 4: C = (K-L)^2 - 1.
```

The leading coefficient in the seven rows is consequently
`1/(BC), 1/(BC), 1/(AC), 1/(AB), 1/(ABC), 1/(AB), 1/(ABC)`, respectively,
multiplied by the same raw native compact C1/H0/H0 moment. The first, second,
third, fourth and sixth hard integrals become products of two unit-mass
tadpoles by unimodular real linear loop changes. The fifth and seventh have
the same connected equal-mass two-loop vacuum hard family. This describes
integrands only; no tadpole or vacuum value is supplied here.

The eta powers, including the two hard loop measures, are D-2 and D-3,
exactly as in the source-derived connection. Their infinity powers are
2-D and 3-D. For this diagonal basis, leading all-hard coefficients suffice
to match the constants. All remaining regions leave an unrestricted virtual
soft direction; after every uncut physical denominator is expanded in eta,
its coefficient is polynomial in that direction and integrates to a
dimensionally scaleless zero. The occupied polynomial factor is retained.

This routing calculation follows the factors consumed by
`OccupiedCutFamily::region_family`, `regions::expand_region`, and
`integrand::projected_factor_region`. A standalone call through those public
owners can record their exact chosen coordinates after the coherent binary
is linked. Until that replay is saved, the table is an independent exact
routing derivation, not a claim that the production mapping has been run.

## Coverage and remaining gates by occupied sector

| Cut | Existing sealed admission | Boundary requirements | Outstanding numerical gate |
| --- | --- | --- | --- |
| [0] | SingletonGerm, two virtual directions, no external-only uncut physical factor | The seven-label closure above needs C1 polynomial moments, tadpole products and the connected equal-mass two-loop vacuum hard coefficient | Public-factory closure, integrated boundary, shared transport and physical projection have not yet run together for this three-loop sector |
| [3] | SingletonGerm after the exact occupied routing, again two virtual directions | The scalar has the same virtual geometry; the raised original slot 0 is now virtual, so its closed basis cannot be inferred from [0] | No completed new-policy source closure is currently saved |
| [0,3] | RankOneBlocks: one virtual loop k, denominators k^2 and (k-p)^2, plus the pure spacelike transfer p=qA-qB | One ordinary hard loop, two compact shells and any finite upper/shell jets actually present in its eventual basis | No completed new-policy source closure is currently saved; this sector supplies the nonzero physical result |

For [0,3], the all-hard coefficient is polynomial in compact momenta after
the deformed pure-transfer factor is expanded. The all-virtual-soft region
has an unrestricted polynomial virtual integral and vanishes. Thus the
current polynomial boundary owner is structurally sufficient; it needs no
retained nonpolynomial soft denominator or new weighted recursion for this
all-uncut-shifted fixture. The actual requested orders and moments must still
be derived from the final native basis.

The new raw Ward pair is deliberately restricted to singleton C1,
energy-free source bases. It is absent for the double cut. Both valid shifted
polynomial source groups remain available there. The cut-[0] closure therefore
does not establish that the new policy closes the other two sectors.

The 677 retained conditions in the completed cut-[0] baseline were also
specialized by exact sparse rational polynomial arithmetic at epsilon
`-5/4`, `-19/6`, `1/8`, `1/1000`, and `1/2000`. None vanishes identically;
each resulting eta polynomial is a nonzero monomial. The input hashes and
counts are retained in `cut0-saved-condition-samples.json`. This finite
sample check neither covers a full Laurent grid nor replaces the public
connection's checks of its own eventual conditions.

## Shared-owner checks and explicit limits

`PreparedOccupiedFlow::prepare_with_source_options` audits the union of the
final basis, original targets, reduced targets, retained candidate keys and
their right-hand-side labels against the bound endpoint/origin permit. The
boundary revalidates the actual basis and family. No reciprocal massless
energy completion is admitted by this policy.

`OccupiedFlowBoundary` retains symbolic dimensional eta powers and uses
`region_orders` to derive the needed depth. It declares zero explicit log
degree at fixed generic symbolic D: all hard denominators are single-scale
massive factors, and Taylor expansion leaves polynomial soft weights.
Laurent-expansion logarithms from eta-dependent dimensional powers must not
be confused with an explicit fixed-D boundary logarithm. Generic-indicial
and exact-condition checks reject exceptional samples rather than silently
changing these classes.

The integrated owner uses `RecursiveTerminalPolicy::TadpolesOnly` and
`bubble_subloops=false`. A two-hard child has fewer loops than the weighted
three-loop parent. Its connected equal-mass vacuum coefficient must be
obtained by the existing ordinary recursive AMF/FT machinery. Current
occupied multi-hard coverage explicitly exercises tadpole products, while
the terminal-policy regression exercises a single-mass sunset. Neither is
an already completed numerical test of this particular connected equal-mass
hard coefficient inside a weighted three-loop boundary.

Limits remain explicit: half-order at most min(series_order,100), 10000
coefficient/product terms, distribution order 32, and tensor rank 32. They
do not obstruct the saved seven-label basis, which has no high distribution
orders or numerator tensors. The unknown [3] and [0,3] bases still require an
actual audit; their mere mathematical finite-label admission does not
guarantee these numerical budgets suffice.

`ConnectionTransport` checks all retained nonzero conditions before even
the empty-basis case, keeps paths in Re(eta)>=0, and multiplies each target's
rational weights into the endpoint series before projecting the physical
constant. In particular the eta^-3 target weights above must remain joint
with their basis modes. The singleton Euler modes carry nonzero symbolic
epsilon slope, so their dimensionally regulated endpoint is zero under the
sealed singleton theorem; the implementation still computes the native
boundary and transport instead of installing this zero as an answer.

That last fact limits what a successful full-amplitude comparison proves
about the singleton boundary constants: every mode of this saved singleton
connection has nonzero symbolic epsilon slope, so the projector returns
zero independently of their numerical constants. Completing the actual
singleton calculation exercises the recursive boundary and transport path,
but its zero physical value is not an independent numerical accuracy check
of the connected hard coefficient. The nonzero double-cut contribution is
the decisive physical-value comparison. An additional finite-eta or hard
coefficient comparison would be a separate validation gate if needed.

Failures can be localized without changing source: the assembly emits
`evaluating occupied cut [...]` before each sector, and the shared connection
then emits compilation, infinity, native boundary matching, transport,
endpoint and physical projection stages. The resource log therefore records
the active sector and numerical stage even though the outer `failure.json`
uses `complete evaluation`. Each successful prediction retains the full
`OccupiedBoundaryProvenance` through the runtime report's boundary field;
this includes region transformations, exponents and requested orders. The
report does not currently serialize each ordinary hard integrand or child
value, and an interrupted full evaluation does not save a partial amplitude.

Full acceptance needs the vacuum and all three occupied cuts in one complete
report, including the raised original numerator and its upper-surface terms,
followed by independent fixed-dimension and Laurent comparisons with working
precision, series depth and starting-point refinement. A closed cut-[0]
program alone establishes none of those numerical outcomes.
