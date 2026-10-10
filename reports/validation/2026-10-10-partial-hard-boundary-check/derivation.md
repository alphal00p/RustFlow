The proposed fixed-shell partial placements leave a one-mass two-loop vacuum period among their all-hard coefficients. This validation uses the ordinary native family `(K²-1, L², (K-L)²)` with powers `(1,1,1)`, at `D=12/5` (`epsilon=4/5`). It is distinct from the earlier equal-mass boundary check. The ordinary provider uses `RecursiveTerminalPolicy::TadpolesOnly` and `bubble_subloops=false`; its preparation counter must be nonzero. No analytic sunset or bubble value is supplied to the native calculation.

For an independent reference, Wick rotate to positive Euclidean denominators with measures `d^Dk/pi^(D/2)`. Gaussian integration gives

`J_E = integral_0^infinity da db dc exp(-a) (ab+ac+bc)^(-D/2)`.

Set `b=a*x`, `c=a*y`. The Jacobian is `a²`; integrating `a` gives `Gamma(3-D)`. The remaining integral is

`integral dx dy (x+y+xy)^(-D/2)`.

The `y` integral equals `x^(1-D/2)/[(D/2-1)*(1+x)]`. The final beta integral therefore gives

`J_E = Gamma(D/2-1)^2 Gamma(2-D/2) Gamma(3-D) / Gamma(D/2)`.

All these integrals converge absolutely for `2<Re(D)<3`, including `D=12/5`. Native Minkowski normalization has the factor `(-1)^3`, so the expected native value is `-J_E`. This sign follows from the actual three inverse propagators, not a blanket finite-density normalization.

The reference implementation evaluates only positive-argument gamma functions through shifted Stirling series with exact rational Bernoulli numbers. The first omitted term bounds each logarithmic gamma truncation on the positive axis. Decimal precision is independently varied through 35, 60 and 85 digits. Native predictions use three separate `(digits,series_order)` profiles `(18,60)`, `(28,60)` and `(28,80)`, with 40 guard digits. Reference and native prediction programs neither import one another's outputs nor read oracle records. Comparison occurs only after prediction files are saved.

This is a nonzero ordinary hard-period check. It neither admits partial occupied placements nor validates their complete region factorization, source closure, physical endpoint, or assembled three-loop amplitude.
