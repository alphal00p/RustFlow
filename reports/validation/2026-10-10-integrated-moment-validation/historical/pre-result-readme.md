# Integrated-moment independent validation

The finite-eta scalar and raised reference evaluator is prepared, but no reference numbers have been generated. Its driver refuses execution until the first integrated producer coefficients and native conversion metadata are frozen and bound by hashes. The producer receives neither the reference formula nor its values.

`native-conversion-audit.json` checks the actual census: all 39 factor/mass/deformation equalities, three signed target terms and 108 exact original-numerator/raised-jet identities pass. For this fixture only, the numerator uses unchanged completion factors and therefore `N_eta=N`. Physical raising acts on the entire original polynomial, kernel, shell measure and Fermi support; the upper surface is retained separately.

`finite_eta_reference.py` uses independent compact Gauss-Jacobi and virtual-parameter Gauss-Legendre quadrature with mpmath. It emits scalar, raised total, and the measure/numerator/kernel/surface pieces. `run-reference.py` enforces an individual 180-second and cumulative 600-second numerical budget, preserves failed attempts, and passes no producer values into the evaluator. Start at eta 8 with node/precision profiles 16/50, 24/50, 32/50 and 32/80. The prospective eta 2 and 1/2 amendment is already bound but is not yet released for execution.

`assess-reference.py` checks saved independent node and arithmetic refinements and the virtual parameter cross-method check. Its error estimate is empirical. There is no ODE, physical finite-eta comparison, production change, weighted closure, or three-/four-loop acceptance result here. See `plan.md` for branch, normalization, heldout and physical-mass checks; `implementation-readiness.json` records the prepared scripts.
