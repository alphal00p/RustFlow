# Frozen bounded physical coefficient fit

The scalar target has a heldout-validated recurrence candidate. The raised mixed-medium target has no candidate within the prospective caps. No equation has been certified, so this attempt does not establish a three-loop solution or a Laurent flow.

The frozen protocol selected a separate scalar operator for each original target. Within each target, all formal period channels were fitted jointly. The actual data have one period channel and the source-proved integer grading q=1, with physical branch beta=-1/4. Each target uses 48 training coefficients and 16 heldout coefficients; theta order/degree and recurrence shift/degree are each at most 4, with at most 25 unknowns and at least eight excess equations. Theta training was attempted first, then recurrence training only on a theta miss. All candidates and misses were frozen before any holdout was opened. Neither a numerical reference nor an analytic answer was used.

| Target | Training outcome | Heldout outcome |
|---|---|---|
| Scalar | Theta miss in 20 shapes; recurrence shift 2 / degree 4, rank 14 of 15 unknowns from 46 rows | All 16 genuine heldout equations pass |
| Raised mixed-medium | Theta miss in 20 shapes; recurrence miss in 20 shapes | Unopened |

The aggregate fitting, digest and scalar holdout validation took 0.845 seconds. The scalar candidate SHA256 is `252a5aa37c821fc0dad292675e7e088cfa4074b9216908cfd00281ed3d00c35a`. The exact selected operator, every attempted shape, frozen candidate identities, commands and resource counters are saved alongside this file. The original bounded fitter and prospective per-target protocol are in the sibling reconstruction and target-protocol reports; they were not modified.

For `sum_j p_j(n)c_(n+j)=0`, the exact generating-function conversion is

    sum_j t^(2-j) p_j(theta-beta-j) g(t) = t^beta B(t),
    g(t)=t^beta sum_n c_n t^n, beta=-1/4, B(t)=-14/117.

It is a fourth-order *inhomogeneous* ODE despite having recurrence shift order two. The startup forcing and fractional branch cannot be discarded. Conversion, ordinary transport and the conditional eta-zero diagnostic are in `../2026-10-10-integrated-flow-conditional-transport/`; those numerical results do not promote the candidate to an identity.

## Certificate scope

The earlier certificate preflight remains accurate: the producer supplies exact integer-index coefficients, but a generic machine-ready all-index telescoping adapter has not been implemented. No such search is claimed here.

A separate necessary-condition check now rules out the prospectively bounded direct-eta rational certificate ansatz for this scalar candidate. In the actual source kernel, `K=eta+V/U^2` (projective gauge U=1) occurs with exponent 5/4. The candidate's leading eta-derivative coefficient is `-32/1053*eta^2*(eta+1)*(eta+4)`. Therefore `L(T)/T` has a generically nonzero fourth-order pole on K=0. A certificate vector component with K-denominator exponent at most two produces a pole of at most order three after divergence and logarithmic differentiation. Direct exact specialization of the actual U/F/V polynomials gives a nonzero fourth-order residue; the other denominator remains nonzero. `certificate-pole-check.py` and its JSON retain that calculation.

This is an obstruction to that bounded direct-eta ansatz, not proof that no differential equation, recurrence certificate, larger-pole certificate or boundary-coupled representation exists. No budget was enlarged after fitting, and the raised target's heldout data were not inspected. Finite-prefix agreement, numerical refinement and a convergent source germ do not by themselves establish an all-coefficient identity.
