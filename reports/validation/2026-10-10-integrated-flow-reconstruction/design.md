# Integrated large-eta reconstruction: fixed exploratory protocol

This is a report-only, isolated exploration. It uses integrated large-eta coefficients supplied by the separate moment-expansion owner, ordinary exact scalar algebra, and the existing ordinary ODE transport. It does not invoke occupation-aware IBP, weighted closure, saved weighted reductions, numerical reference values, or a topology-specific evaluated formula.

## Coefficient contract

The producer exposes the complete known prefix of

    I(eta,D) = sum_(sector,ell,k) M_k(D) eta^alpha_sector(D)
               (log eta)^ell sum_(n>=0) c_(sector,ell,k,n)(D) t^n,
    t = eta^(-1/q).

Start with q=2. A change to q=1 requires a producer proof of the parity/grading, not a fit that ignores odd terms. Preserve raw target eta offsets, normalization/phases, regulator order, exact physical input/cut/deformation identity, all retained parameter conditions, and the source of every integrated coefficient. An absent coefficient is unknown unless the producer explicitly proved it zero. Each M_k is independent of n and has a stable exact ordinary hard-master/compact normalization identity. Scalar coefficients belong to Q(epsilon), or Q after exact specialization. No decimal rationalization is permitted. Rational D alone does not imply rational integrated values: Gamma factors and hard periods must remain explicit.

Fit all relevant exact channels jointly. Treating distinct period symbols as formal independent channels is a sufficient coefficientwise annihilation test, not a claim that the periods are independent or that the resulting order is minimal. Merge overlapping exponent/normalization contributions exactly when their equality is proved. Keep exponent and log sectors explicit, including resonances at special D. D=13/2 (epsilon=-5/4 for D=4-2epsilon) is the first fixed value; symbolic epsilon is a separate result and cannot follow by inference from a few fixed dimensions.

## Budgets fixed before fitting

- Producer census/pilot: at most 24 coefficients and 180 seconds. Extension, only if feasible: at most 64 coefficients and 900 seconds. Existing half-order<=100 and tensor-rank limits remain in force; no silent owner extension.
- The intended main fit reserves indices 0..47 as training and 48..63 as untouched heldout data. If only a shorter prefix exists, report that fact; do not silently relabel training data as heldout or claim the main protocol passed. A smaller diagnostic fit must disclose its exact split and satisfy the same excess-equation requirement.
- Recurrence shift order r<=4 and polynomial degree d<=4. Direct theta-ODE derivative order r<=4 and polynomial t-degree s<=4. At most 25 ansatz coefficients, at least 8 more training equations than the projective unknown count U-1. Candidate enumeration is deterministic by increasing U, then order and degree. All identically zero/startup rows are recorded; matrix rank and nullity are reported.
- At most 32 exact coefficient channels, log power<=2, 2048 matrix rows, 25 columns and 24 shape attempts. These are conservative admission bounds, not mathematical limits.
- Fit stage: cumulative 300 seconds. Certificate stage: separate cumulative 300 seconds. No budget increase after a miss.
- Retained exact algebra payload: at most 200,000 monomials and 128 MiB of serialized canonical algebra data. Reject oversized inputs before fitting and check retained outputs. This is not a bound on internal CAS scratch or peak process memory. Timeout and process RSS are recorded separately; a killed CAS operation is an incomplete result.

Every coefficient prefix, split, proposed shape, attempted rank, resource rejection and candidate digest is saved. Heldout data cannot choose the preferred nullspace vector or shape: choose a candidate using training alone, freeze its JSON and SHA256 before reading the heldout file, then evaluate that one frozen candidate. A heldout failure ends this dataset attempt; do not tune or search a higher order using the same heldout coefficients. Multiple surviving nullspace directions remain nonuniqueness, not a unique reconstructed equation. A zero sequence needs a producer zero certificate; fitting a zero prefix is insufficient.

## Exact fitting and conversion

Use pinned Symbolica Matrix::from_nested_vec and Matrix::row_reduce, over Q at fixed D and RationalPolynomialField over Q(epsilon) if affordable. The small adapter constructs the linear equations and extracts free-column nullspace vectors; it does not implement a new elimination engine. Canonicalize only by an explicit nonzero rational/rational-function scale and retain its denominator conditions.

For a recurrence candidate

    sum_(j=0)^r p_j(n) c_(n+j) = 0,

training rows use only training coefficients; heldout rows must contain at least one heldout coefficient and no unknown coefficient. Preserve singular leading-coefficient indices. The corresponding generating-function equation generally has startup forcing:

    L F = B(t),   L = sum_j t^(r-j) p_j(theta-j),   theta=t*d/dt,

where B is computed exactly from the initial coefficients and has degree<r. Do not silently discard B. Its differential order is <=d, which differs from the recurrence shift order. If homogenization is requested, left multiplication by product_(k=0)^(r-1)(theta-k) has order<=d+r<=8; otherwise preserve the inhomogeneous polynomial forcing and its bounded augmentation explicitly.

Direct theta-ODE fitting avoids this conversion ambiguity. For L=sum_(s,j) a_(s,j) t^s theta^j and a channel t^beta sum c_n t^n, each coefficient row is sum_(s,j) a_(s,j)(n-s+beta)^j c_(n-s)=0. Log channels use the exact action theta[t^lambda(log t)^ell]=lambda*t^lambda(log t)^ell+ell*t^lambda(log t)^(ell-1), including log(eta)=-q log(t). Initial absent negative-index terms are zero only for a declared one-sided series.

A certified scalar operator can be converted to a rational companion system in t using state (f,theta f,...). Divide by the leading theta coefficient only on declared nonzero patches; retain apparent singularities, local branches and initial data. DifferentialSystem::compile and CompiledSystem::transport are the existing transport owners. No artificial ordinary-integral basis is installed through SuppliedAuxiliarySystem: its constructor is a declared exact-identity trust boundary, not a proof of a finite fit.

## What would certify an equation

A finite fit, unique nullspace, heldout success, or stable numerical continuation is not an identity certificate. Two acceptable routes are:

1. An exact all-n recurrence proof from the producer's generic moment representation. It must cover the actual sum of all regions/channels and all indices, with explicit startup/singular-index handling and justified conversion from the integral to the represented germ.
2. A generic parametric telescoping certificate for the original scalar parameter integrand, after the producer supplies an admissible fixed-domain representation. For hyperexponential terms T, verify exactly that sum_j p_j(n) T(n+j)/T(n) equals sum_i [d_i R_i + R_i d_i(log T)], or the analogous differential-parameter identity. Clear only declared nonzero denominators and prove the resulting polynomial identity with Symbolica. Every boundary term must be evaluated, vanish in a common convergence region, or appear as an explicit inhomogeneous component; analytic continuation then uses the same regulator prescription.

The initial certificate search is deliberately small: at most 4 integration variables; candidate denominators use only at most 8 existing irreducible integrand factors, each exponent<=2 and total denominator degree<=12; numerator total integration degree<=3 and n-degree<=2; at most 420 certificate unknowns and 2048 exact coefficient equations. Existing factors and their allowed powers are fixed before search. Retained algebra caps and the 300-second limit also apply. Nonrational shift/log-derivative ratios, unsupported finite sums, moving-domain boundary data, or a required larger ansatz return explicit unsupported/incomplete results. A generic telescoping discovery API was not found in the inspected pinned Symbolica owners, so this is a bounded exact certificate adapter proposal, not a claim that a complete telescoper already exists.

A proof about formal coefficients alone still needs an actual analytic-germ or asymptotic-uniqueness argument before it becomes an equation for the physical integral. Beyond-all-orders terms can share every inverse-power coefficient. The producer must state whether its expansion converges in a declared domain or is only asymptotic; expansion by regions alone must not be relabeled a convergent Taylor germ.

## Direct continuation fallback

Without a certified ODE, the existing ODE residual checks cannot validate continuation of a guessed equation. Padé approximants of the coefficient prefix may be reported as diagnostics, with order/sampling/precision comparisons, but are not accepted as the physical solution solely from those comparisons.

A controlled alternative needs an independent analytic remainder bound. For example, an exact local analytic germ with a proved disk bound M on radius R gives a geometric Taylor tail bound M(r/R)^N/(1-r/R); overlapping chart continuation requires such bounds at each step and the prescribed branch path. Another option is direct evaluation of the exact scalar parameter representation with certified quadrature, UV/origin subtractions and contour control. Neither bound nor a generic quadrature owner is supplied by a finite fit. If those owners are unavailable, the honest output is a bounded candidate with explicit missing proof, not a numerical acceptance equation.

## Executed fitter scope

The first isolated executable accepts exact rational, log-free channels at D=13/2 only. It rejects logs and any incomplete prefix rather than silently projecting them. The broader symbolic-epsilon/log certificate design above remains unimplemented. Its recurrence output applies to each coefficient core c_n; channels with different beta require different conjugations theta -> theta-beta when converted, so no common full-integral scalar ODE is claimed from that recurrence alone. Direct theta fitting includes beta in every exact matrix row and avoids that ambiguity.

The CLI receives separate training and heldout files. Fit does not receive the heldout path. Validate checks the caller-supplied frozen BLAKE3 of the candidate and the candidate's training digest before opening heldout data. Outputs use create_new and cannot overwrite an earlier candidate. The protocol requires recording both SHA256 and BLAKE3 of a selected candidate before validation; the mathematical discipline also requires not using a rejected holdout to choose another operator family or shape. A caller must choose the family or a deterministic training-only cross-family selection order before fitting. No physical family selection has occurred yet.

Eight synthetic tests passed, including a false extension with the same 48 training terms, common formal channels, exact rational exponents and the candidate-hash-before-holdout check. One first compile failed on test-only owned String moves; its source/log/build binding are retained. These tests establish adapter behavior, not an identity or numerical result for a physical integral.

## Pre-fit grading amendment

The parent approved q=1 after the producer's source-bound hard-reversal/parity and raw-target-offset proof. The count remains 64 meaningful integer coefficients per fractional branch, with 48 training and 16 heldout; it is not padded by 32 structural half-grade zeros. The entire fractional eta prefactor survives as beta=-q*alpha in each t^beta channel. `grading_proof_identity` is mandatory for q=1 and must match between training and heldout. This amendment was recorded before any physical fit. Direct theta is the predetermined primary physical operator family; a heldout rejection will not select a different family against the same holdout. All original order/degree/time/algebra caps remain unchanged.
