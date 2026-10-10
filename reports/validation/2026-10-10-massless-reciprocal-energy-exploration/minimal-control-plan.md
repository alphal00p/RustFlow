# Minimal controls for a reciprocal-energy source extension

These are validation designs, not implemented production admission. The physical inputs contain only polynomial numerators. Native proofs must be generated afresh from the newly identified reciprocal-energy source family; previous polynomial-domain programs cannot be reused as certificates for the enlarged family.

## One-loop graph and exact radial owner

Use one vertex and one charged massless self-loop, routing q, chemical potential mu>0, and original shifted-Euclidean targets `(n,N)=(1,1),(2,1),(2,u1^2)`. The corresponding occupied coordinates are `g=q^2` and `E=u.q`; use factors `[g,E,mu-E,E]`, with roles `[RequiredCut,Ordinary,Occupation,Occupation]`. The second factor is the individual energy completion, not an additional physical propagator. Its domain is all integer powers only under the proposed joint radial prescription. Occupations remain nonnegative, and positive lower contacts use that joint continuation rather than empty support.

There is no genuine fixed-shell auxiliary-mass flow for this minimal physical graph: after cutting its only line, there is no uncut physical denominator to shift. The correct one-loop gate is native source/replay plus signed compact boundary/terminal evaluation. A zero or fictitious eta connection must not be reported as an AMF validation. The existing `MasslessFlowEvidence::new` also requires at least two loops; it should not be bypassed to manufacture a one-loop flow permit. A future additive compact-origin certificate can express the elementary radial proof independently of virtual endpoint evidence.

At spatial dimension d and for the original real-energy monomial E^r with no radial insertion, the Euclidean occupied value is

```text
I_n(r) = -mu^(d+r+1-2n)
         / [(4*pi)^(d/2) Gamma(n) Gamma(d/2-n+1) (d+r+1-2n)].
```

This is the complete binomial expression after its angular Gamma cancellation, initially derived on its high-D convergence domain and then continued. A literal original P0^r supplies its fixed Wick factor i^r. Keep the binomial form when combining possible removable factors before assigning a dimension; do not separately evaluate singular Gamma/ratio factors. This formula is a validation reference only, not a production shortcut for virtual integrals.

For raw spatial-normalized shell moments, write `M_n(r)=int C_n E^r H_upper,0 H_lower,0` and `S_n(r)=int C_n E^r H_upper,1 H_lower,0`. Then `I_n(r)=(-1)^n M_n(r)`. The native shell-normal identity at n=1 is

```text
M_2(r) = (r-1)/2 M_1(r-2) - 1/2 S_1(r-1),
```

after the jointly defined lower contact vanishes. For r=0 it forces both reciprocal-energy terms:

```text
M_2(0) = -1/2 M_1(-2) - 1/2 S_1(-1).
```

For r=2 it gives `M_2(2)=M_1(0)/2-S_1(1)/2`; the original Euclidean numerator u1^2 contributes the additional minus Wick phase. Each identity must be verified by native source replay and by independently integrated bulk and upper moments. The upper term is nonzero and cannot be omitted. The existing order can prefer eliminating a positive energy completion over a raised-cut label, so the gate should reconstruct the complete identity without imposing a particular pivot orientation.

At D=5, d=4, mu=1, `A_d=1/(8*pi^2)`. Independent references are:

| Object | Euclidean spatial value |
|---|---|
| M_1(-2) | 1/(16*pi^2) |
| S_1(-1) | 1/(16*pi^2) |
| M_2(0)=I_2(0) | -1/(16*pi^2) |
| physical full scalar n=1 | -1/(48*pi^2) |
| physical full scalar n=2 | -1/(16*pi^2) |
| physical full n=2, original u1^2 | +1/(48*pi^2) |

The vacuum terms are massless scaleless integrals. These physical values are already checked by the existing polynomial terminal regression; the new test is reconstruction through explicitly admitted reciprocal-energy intermediates with the same values. Convert raw spatial moments to the native occupied normalization exactly once by `(2*pi)^d/pi^(D/2)`, then apply the whole-graph measure and Wick factors as already prescribed.

Negative checks must include the true pole of `M_1(-2)` at D=4, a positive-energy removable binomial/ratio example from the existing seed contract, mu=0, a mixed energy sum, an inverse virtual energy, negative occupations, and source-context mismatch. A numerical pole in an auxiliary moment cannot be silently treated as zero even if another representation of a complete target could cancel it.

## Genuine flow control

The minimal existing physical flow is the massless two-loop sunset with cuts `[0]`, `[1]`, and `[0,1]`, already admitted by the singleton-germ/pure-transfer proof. Preserve all uncut physical eta shifts and the original scalar/raised-medium targets. Regenerate each source context with the reciprocal completion proof only where the exact family contains c*E_i; audit every basis, candidate and target label with the strengthened radial bounds. Require complete native derivative closure, integrated occupied boundary matching, shared transport and physical endpoint reconstruction. Compare the unchanged polynomial physical amplitudes against the saved independent two-loop references with precision/order/start/epsilon refinement. A smaller frontier or a direct source application alone is not this gate.

## Three-loop diagnostic before implementation

For the current three-loop singleton cut0, the proposed native control must derive its completion mapping from `PreparedDensityInput::occupied_cut([0],...)` and verify exactly that the claimed slot is the ordinary input completion `c*E_0` with c nonzero rational. Save the actual factor, medium coordinate, shell loop and equality residual. The candidate slot6 must not be accepted merely because a previous report used that index. Keep the original full source rows and measure factors, bind a new origin/source identity, enlarge only the certified completion and its matching lower/free-virtual zero domains, and discover/replay every rule fresh. Compare identical point and native-budget controls before any discovery-performance claim.
