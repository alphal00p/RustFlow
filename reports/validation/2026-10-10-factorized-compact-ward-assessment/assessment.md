The factorized compact-leg extension is justified under the existing joint dimensional prescription, with the following exact source-base guards. This is a mathematical assessment, not production admission or a native closure result.

For an occupied loop i require its physical shell to be exactly `q_i²`, its cut index to be one, its lower factor to be `E_i` with index zero, and its upper factor to be `mu_i−E_i` with `mu_i>0`. Every ordinary factor that depends on `E_i` or any incident Gram coordinate `g_ij` must have index exactly zero. This includes physical denominators and numerator completions, including negative numerator powers. Every other cut and occupation factor must be exactly independent of those coordinates. The existing polynomial-completion and frozen-storage domains remain in force.

These are sufficient, deliberately conservative conditions. Test dependence on the actual transformed factors, not source IDs, graph names, slots, positive supports alone, or a restriction to the shell. A factor such as `g_ii+g_jj` still counts as dependent. Nonzero powers of a dependent factor cannot be admitted merely because its on-shell value looks simpler. Source guards constrain bases; the resulting row is allowed to contain the explicit polynomial energy numerator on its right-hand side.

Under these guards the actual Cartesian measure factorizes, independently of whether its scalar-coordinate Gram Jacobian looks coupled. Apart from a common angular/normalization constant and the untouched remainder, the selected leg is

`M_s = integral_0^infinity E^(D−3) H_s(mu−E) dE`.

The total radial derivative of `E^(D−2) H_s(mu−E)` gives

* `(D−2) M_0 − M_1[E] = 0` for the upper theta;
* `(D−2) M_s + s M_(s+1)[E] = 0` for s positive.

Here `H_s=(-1)^(s−1) delta^(s−1)/(s−1)!`; consequently its derivative with respect to its argument is `−s H_(s+1)`. The extra minus sign from `mu−E` explains the positive coefficient in the second row. No division by mu or E remains in either identity.

For the selected C1 leg the lower radial flux vanishes on the open domain `Re D>2`. Positive upper derivatives have support at the separated endpoint mu, so they introduce no new lower-origin singularity. The other occupied legs may have arbitrary finite raised cuts and upper derivatives already admitted by their own sealed origin prescription. Choose a common sufficiently high nonresonant dimension for their finite jets, perform the fixed-positive-eta and fixed-T regulator removal in the established order, and continue the complete tensor-product identity meromorphically. Virtual UV continuation is applied to the independent remainder. This needs no new pointwise product at the origin, no virtual endpoint-zero theorem, and no ordinary-vacuum zero classification of an occupied integral.

The same rows must not be used for a raised selected cut: the `q_i² C_n` contact no longer vanishes for n>1. They must also reject a selected energy numerator or incident Gram numerator; such insertions alter radial degree or couple the angular factor. The assessment supplies explicit nonzero residual counterexamples for the former two exclusions. Other legs' raised indices remain unchanged, so this is not a differentiation of the C1 Ward row into an unsupported Cn identity.

For the actual three-loop chain with cuts `[0,3]`, rational completion gives `q0=P1`, `q1=P1−P3`, `q2=P2`. The uncut factors are `q2²−eta`, `(q0−q1)²−eta`, and `(q2−q0+q1)²−eta`. Thus either compact-leg row requires physical slots 2 and 4 to have index zero; slot 1 is independent. Additional completion guards must still be derived from the actual factor map. This confines the proposed benefit to factorized lower sectors and says nothing about reducing the original coupled positive-power sector.

`check.py` derives these physical routing guards from the input JSON and performs exact rational product-moment checks with a second cut of orders one through three, nonzero upper contacts, and an independent second-leg energy numerator. It reads no native predictions or supplied answers. The test checks signs and cross-cut factorization, not native search performance.

A minimal future implementation would use a separate sealed factorized-leg certificate, bound to the exact family/deformation and existing origin evidence. It should reuse the raw-energy row conversion with its stronger per-leg base domain, leaving the current singleton theorem intact. Necessary native tests include own C2 rejection, dependent positive and negative ordinary powers rejected, other-leg C2/upper contacts allowed, mixed energy routing, unchanged zero/tail domains, both upper branches, full original-source replay, and a measured comparison on the relevant lower-sector requests before making any closure claim.
