# Independent validation of integrated moment series

This is an isolated validation design. Weighted closure remains paused. No native closure, production edit, finite-eta reference evaluation, ODE inference or reference comparison has been run here. The first executable check is exact rational routing/basis algebra only. Numerical reference generation must wait for a producer artifact containing its saved predictions and a frozen equation/training specification.

## Match the actual deformation first

Use the unchanged `massless_three_loop_chain` input and double cut `[0,3]`, with **all uncut physical slots `[1,2,4]` shifted**. This differs from the interrupted partial-placement control. In Euclidean notation let `qA=P1`, `qB=P1-P3`, `p=P3=qA-qB`, virtual momentum `k=P2`, and `h=p^2`. The denominators become `(k^2+eta)*((k-p)^2+eta)*(h+eta)`. Shells and Fermi supports remain fixed; eta is not a physical cut mass.

The producer must save the actual prepared independent-mass basis, converted indexed targets, Wick factors, physical mass assignments and shifted factor substitution. The independent Fraction calculation in `basis-crosscheck.json` gives physical rank five and completions `[g1_2,u1,u2,u3]`. For this particular input the original numerator is `rho5+rho6*rho7`, so it uses only unshifted completions and its conversion coefficients do not depend on any physical mass. Thus its deformation is expected to satisfy `N_eta=N`. Acceptance of that statement still requires the producer's actual prepared metadata to match. No such equality is assumed for a generic numerator. In general reconstruct the eta-dependent polynomial from the frozen indexed combination before comparing anything.

Physical mass raising must precede mass assignment. If an indexed unraised target is `sum c_I(a) J_I(a,eta)`, its raised version includes both `-c_I'(a) J_I` and the denominator/distribution derivatives. A helper that discards `c_I'(a)` is unsound in general. In this fixture the coefficient derivative is zero, but the shell-substituted original polynomial still depends on `a`; that dependence is essential. Direct C2H0 integration is valid only if its distribution owner acts on the entire remaining polynomial and support. C2H0 does not imply that the upper surface is absent.

## Reference-only finite-eta integral

Take `d=D-1`, `alpha=(D-2)/2`, `nu=D/2-2` and

```
A = area(S^(D-2))/(2*pi)^(D-1),
w_alpha(t) = t^(alpha-1)*(1-t)^(alpha-1)/B(alpha,alpha),
C = (4*pi)^(-D/2)*Gamma(2-D/2),
I_eta(h) = integral_0^1 [eta+x*(1-x)*h]^nu dx,
K_eta(h) = C*I_eta(h)/(eta+h).
```

The Feynman-parameter Gaussian formula is initially justified with convergent virtual powers and continued meromorphically in D. At real positive eta and massless future shells, `h=4*r1*r2*t>=0`, so its compact integrand is smooth away from the ordinary radial measure endpoints. D=13/2 is the first proposed witness, chosen by the existing high-D mass-jet proof before any producer outcome; D=31/3 is an independent generic-D check if cost permits. Neither is a virtual Gamma pole. This is a UV-continued reference, not a claim that the unregulated high-D virtual momentum integral converges.

The scalar contribution is

```
S(eta) = A^2/4 * integral_[0,mu]^2 dr1 dr2
          r1^(D-3)*r2^(D-3) * average_t K_eta(4*r1*r2*t).
```

For the raised target introduce `a=m0^2` at fixed original spatial vectors **before** differentiation:

```
E1=sqrt(r1^2+a), E2=r2,
h=-a+2*(E1*E2-r1*r2*z), z=1-2*t,
h_a=-1+E2/E1,
Nbar=(h-a)/4-E1^2/2+E1*E2/2,
Nbar_a=-1+E2/(2*E1).
```

The first virtual moment is `k=p/2` by the exact change of variable `k -> p-k`; equal virtual masses eta preserve it. This observation is confined to reference generation and must not be passed into the generic producer as a subloop shortcut or equation ansatz.

With `Wi=A*r_i^(D-2)/(2*Ei) dr_i`, the raised bulk is

```
integral W1 W2 average_t [
 K_eta(h)*(Nbar/(2*E1^2)-Nbar_a) - K_eta'(h)*Nbar*h_a
] at a=0.
```

The **positive moving upper-surface contribution**, saved separately, is

```
A*mu^(D-4)/4 * integral W2 average_t [K_eta(h)*Nbar] at E1=mu,a=0.
```

At zero physical mass the bulk can be evaluated without subtracting singular symbolic pieces:

```
K_eta(h)*[3/4+(r2/r1)*(t/2-1/4)]
- K_eta'(h)*[r1^2/2-r1*r2*(t+1)+r2^2*(t+1/2)].
```

Keep the measure, original-numerator derivative, kernel derivative and upper surface as separately reported components even when their stable sum uses this expression. Differentiate the reference kernel directly:

```
K_eta'(h)=C*(nu*integral x*(1-x)*[eta+x*(1-x)*h]^(nu-1) dx/(eta+h)
                 -I_eta(h)/(eta+h)^2).
```

For fixed positive eta this compact raised integrand is integrable at D=13/2. Conservative high-D origin bounds apply before continuation, and no pointwise multiplication of a singular lower-origin jet is introduced. The pre-existing eta=0 reference can be checked only after the finite-eta producer and equation are frozen; it does not supply initial conditions or train the equation.

## Branch and cross-method checks

For Re(eta)>0, both `eta+h` and `eta+x*(1-x)*h` remain in the right half-plane for the massless reference, so use the branch continued from positive real eta. A finite physical-mass difference check uses `0<a<min(mu^2,Re(eta))/16`; its possible negative h cannot cross either kernel cut. The radial upper limit is `sqrt(mu^2-a)`. Differentiate with the spatial radius fixed, then include that moving boundary; differentiating at fixed E without its Jacobian/support is a different operation.

An optional reference accelerator is

```
I_eta(h)=eta^nu * 2F1(-nu,1;3/2;-h/(4*eta)).
```

It follows by transforming the same one-dimensional parameter integral; its principal branch is justified by the same right-half-plane condition. Cross-check this accelerator against direct x-quadrature before using it to remove one integration dimension. NIST's Euler integral states the parameter/branch constraints, and its notation uses the regularized function, which must not be confused with ordinary 2F1. [DLMF 15.6.1](https://dlmf.nist.gov/15.6.E1). The angular and parameter moments use Euler's beta integral. [DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1).

## Affordable ordered validation gates

1. **Artifact and physical matching.** Hash the unchanged input and actual emitted factors/targets, then verify the all-uncut placement, original coefficients and physical raising algebra. No ODE/reference values participate in this gate.
2. **Freeze producer outputs before references.** Save the moment normalization, exact dimension, training interval, equation coefficients/order/degree, singularity list, initial data and predicted values. Freeze every producer source/command/dependency hash. The reference process reads only the saved output after this point and cannot write to producer inputs.
3. **Unused moment coefficients.** Reserve at least 16 consecutive coefficients beyond the training range and another 8 at a separated higher range; increase that count if it is less than twice the guessed free-coefficient count. Generate them independently from the parameter/radial/angular integral, rather than asking the fitted ODE for them. Exact rational-D normalized moments are preferable. Preserve unsuccessful orders/ansatz choices; do not select the held-out range after seeing errors.
4. **Finite-eta values.** Initial fixed real points are `eta/mu^2 = 8, 2, 1/2`, with both targets and all raised components. These cover a large-eta series point and points beyond its immediate geometric convergence disk. Add `2+i` and `2-i` as a paired branch/path check once the real controls pass. A smaller `1/8` point is a separate bounded refinement, not needed to choose the first equation.
5. **Independent quadrature refinement.** Start with tensor Gauss-Jacobi radial/angular and Gauss-Legendre x orders 16/24/32 at 50 digits, then repeat the largest order at 80 digits. The raised r2/r1 term may use its own radial Jacobi weight. Surface terms have one fewer radial dimension. Use two successive order changes plus an independent precision change; do not infer precision from the final order alone. If direct four-dimensional quadrature is too expensive, retain it at a bounded subset and use the independently checked parameter-integral accelerator for denser three-dimensional grids.
6. **Physical-mass/surface cross-check.** At eta=2, compare the separate analytic jet plus surface with one-sided mass differences of the original unraised polynomial integral at decreasing positive a. Reparameterize the changing radius only for numerical quadrature, retaining its exact a dependence. Establish empirical convergence; do not assume an integer Richardson order when fractional origin remainders may occur. A finite real-temperature Fermi profile is an optional additional check of the positive upper delta sequence after the same fixed-eta regulator order.
7. **Producer truncation/precision/path controls.** Independently increase retained moment order, arithmetic precision, ODE Taylor order and start eta, changing one at a time. Transport from one common large positive basepoint to the same endpoint along the positive ray and a path entirely in Re(eta)>0. Require agreement with direct finite-eta reference values and between paths; conjugation alone is not an independent reference.

The initial numerical comparison criterion remains relative `1e-12`, switching to absolute `1e-25` below magnitude `1e-20`, with `1e-25` imaginary-zero tolerance only for real positive eta. References must first show at least two extra stable digits under independent order/precision refinement. If that accuracy is not reached, report an unresolved reference accuracy limit rather than relaxing the acceptance threshold after inspecting predictions. Refinement estimates are empirical, not rigorous error bounds.

Finite held-out agreement is evidence for a guessed ODE, not a proof of the differential equation. Production promotion needs an exact all-order recurrence or annihilator/integration certificate, including endpoint/contact terms, or a clearly bounded numerical-only contract. A finite-dimensional fitted equation at one D does not establish symbolic-D continuation, Laurent coefficients, or the physical eta endpoint.

## Generic four-loop requirement and current limits

The producer should receive graph routing, propagator powers, original polynomial, shell/occupation distributions and deformation. Its moment objects must retain every compact Gram invariant and virtual tensor contraction. The three-loop reference happens to reduce to one h; the producer must not special-case that reduction or call this bubble formula. A future four-loop structural gate needs at least two independent compact invariants and nontrivial virtual tensors, with the same generic moment and distribution owners. An independent multi-parameter Gaussian/Schwinger and compact-angle reference can then be prepared under its own budget. No four-loop finite-eta reference or ODE claim follows from the present single-bubble control.

Immediate remaining dependencies: actual producer basis/target metadata, saved reference-free predictions, and its proposed generic equation/recurrence contract. No expensive reference calculation is needed before those artifacts exist.
