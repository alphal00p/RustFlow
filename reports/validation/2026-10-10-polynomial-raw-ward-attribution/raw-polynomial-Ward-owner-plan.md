# Raw polynomial singleton-shell Ward source: proof and owner plan

This is a read-only design, with separate report-only native attribution probes.
It does not enable a production source or reciprocal-energy domain. The source
is the previously derived C1 normalized shell boost before applying occupation
multiplication identities.

For one future massless occupied momentum q=E(1,n), a simultaneous Lorentz boost
of every loop in the spatial direction n has occupied-shell action E*d/dE.
All Gram products, scalar virtual factors and scalar deformations stay fixed.
Virtual Lorentz measures have zero divergence; the shell measure
E^(D−3)dE dOmega supplies divergence D−2. Scalar Gram numerator degree therefore
does not enter this coefficient. The argument uses the existing joint
high-dimensional origin prescription at fixed positive eta; an origin cutoff
vanishes before meromorphic continuation. It does not multiply a reciprocal
energy distribution by an apex contact.

## Raw rows and guards

With C1, all energy-dependent ordinary-factor indices fixed zero, and lower
occupation index zero:

    upper s=0:   (D−2) I_s − E I_(s+1) = 0,
    upper s>=1:  (D−2) I_s + s E I_(s+1) = 0.

`E I` means a polynomial numerator multiplying that distribution-weighted
integral. If an ordinary factor is exactly E, it is index shift −1. It is not
an inverse factor. Upper differentiation has minus sign for theta and plus s
for the normalized positive-index Hs. The lower flux E H1(E) vanishes in the
same origin prescription.

An upper-index base pullback s=H−1 gives:

    H=1:   (D−2) I_(H−1) − E I_H = 0,
    H>=2:  (D−2) I_(H−1) + (H−1) E I_H = 0.

The upper guards must become exactly H=1 and H>=2; all other source-base guards
remain unchanged. In particular the base energy index stays zero even though
the row contains an energy-numerator term. A native recurrence targeted at
energy index −1 can consequently be valid; that does not enlarge the source
base domain to arbitrary energy powers.

Multiplication gives E H1(mu−E)=mu H1(mu−E) and, for s>=1,
E H_(s+1)(mu−E)=mu H_(s+1)(mu−E)−H_s(mu−E). Thus the raw rows are exactly
equivalent to the earlier mu rows. Their different first-hit behavior can
reflect native source presentation alone. It cannot establish that reciprocal
energy admission was needed.

The narrow source does not cover raised cuts: q² Cn=C(n−1) survives for n>1.
It also does not cover positive cut mass, multiple occupied shells, or an
energy-dependent ordinary factor with nonzero source index. Those require
additional derivations. Do not import a shifted old boost source from a
forbidden inverse-energy seed to claim this new source's provenance.

## Generic production placement, if validated

1. The geometric/physical admission should be in `preparation.rs`, bound to the
   existing sealed `MasslessFlowEvidence`. Require a typed singleton-germ
   variant and exact input/family/deformation/source-option binding; one shell,
   assigned mass zero, mu>0 and the fixed-T/eta joint origin prescription.
   A typed private evidence accessor can expose the occupied loop and shell
   slots. Do not test a serialized classification string or graph name.
2. Put polynomial row construction in `measure.rs`, using its exact
   `multiply` and `source` owners. The numerator `E` must be converted through
   the actual inverse-propagator basis; never assume slot6 or unit coefficient.
   This handles a rational rescaling or an affine completion basis while
   preserving all generated physical shifts and exact coefficients.
3. Detect every non-Lorentz-scalar ordinary factor structurally and fix its
   source index zero. A conservative test is exact independence of all medium
   energy coordinates; verify all physical virtual factors are invariant.
   Scalar Gram completions retain their original nonpositive bounds. Intersect
   the whole row domain with the existing admitted box and freeze storage tails.
   Merely recognizing completion names beginning with `u` is insufficient.
4. Build the two upper cases with cut index exactly1 and lower index exactly0.
   Lower-positive sectors keep only existing origin-zero evidence. Preserve
   theta at index0 and all original occupation conventions. Budget the emitted
   rows/converted terms with checked arithmetic; failure provides no partial
   source claim.
5. Keep default source behavior unchanged initially. An explicit presentation
   option may append/prepend this pair, with a new theorem/presentation identity
   bound into source/cache provenance. It must not set
   `positive_compact_energy_powers` or widen any original zero domain.
6. A whole polynomial boost shifted by +1 on a completion equal to cE is a
   separate presentation of an already admitted identity only on the complete
   pulled-back source domain. From original index<=0 this gives new index<=−1.
   Shift coefficient index symbols, every row shift and source conditions too;
   leave the original support-zero boxes unchanged. Do not confuse that source
   transformation with the separately proved C1 raw Ward pair.

No numerical boundary adaptation is needed for these polynomial source rows.
Any new derivative-closed native basis still passes the existing label,
boundary and endpoint audits. This is a finite-eta source identity, not an
endpoint-zero value or an analytic subdiagram evaluation.

## Required tests before production use

* Exact raw and multiplied radial distribution checks for upper H0 and several
  positive Hs, including nonzero upper surfaces; equality of the raw/mu rows
  after the existing multiplication source.
* Typed admission rejects mass>0, raised Cn, multiple occupied loops, mu=0,
  missing/mismatched origin evidence, any nonzero energy-factor index, negative
  occupation and nonzero storage tails. Native rule-domain tests distinguish
  source base E-index0 from the legitimate E-numerator image index−1.
* Exact coordinate/routing shear and a nonunit rational energy completion cE
  give the same polynomial row under the mapped basis. Energy-dependent mixed
  completions receive zero-index guards; arbitrary Gram polynomials do not.
* Original domains, nonzero conditions and support-zero premises survive native
  sealing, encoding/decoding and full-context replay. Recentered and
  unrecentered variants have exact pulled-back upper guards.
* Matched native controls compare raw Ward only, shifted polynomial boost only,
  both orders together, and the widened reciprocal experiment. Use identical
  bounded search and full proof accounting before attributing any gain.
* Only a subsequent complete fixed-D/precision/Laurent and target-reconstruction
  gate can support numerical acceptance; shorter local RHS is insufficient.
