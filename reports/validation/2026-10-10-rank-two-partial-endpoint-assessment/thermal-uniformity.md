# Uniform real-temperature pairing after regulator removal

This is an additive clarification of the reviewed proposal; it does not modify the frozen origin source or claim a new production consumer. Independent line regulators are removed first at each fixed T>0 inside the Fermi pole strip. Only the resulting **real** Fermi weights are considered uniformly for 0<T<=T0.

Let the angularly dominated finite-jet kernel and its endpoint difference have the bounds proved in the proposal. On finite energy intervals, each needed derivative is bounded by an integrable radial majorant; the difference has the same majorant times rho^sigma. At large energies the majorant grows by a fixed finite polynomial. At each compact origin, the strengthened strict lower-contact inequality gives a positive residual power after every required origin jet. All constants may depend on the finite label, the fixed generic D neighborhood, positive chemical magnitudes and T0, but not on eta, sufficiently small line regulators after their prescribed removal, or real T in (0,T0].

For H0 the real occupation f_T(E-mu) lies between 0 and 1 and, for E>=mu, is at most exp(-(E-mu)/T0). It therefore supplies a common tail majorant. Real lower-support smoothing is also bounded by 1; the separate positive origin power controls its limiting endpoint.

For upper H_s, s>=1, transfer exactly **s-1** derivatives by integration by parts onto the kernel, leaving the positive logistic delta approximation

    delta_T(E-mu) = 1/[4T cosh^2((E-mu)/(2T))].

Its L1 norm on the half-line is at most 1. This uses the existing s-1 jet count; transferring all s derivatives to a primitive theta weight would unnecessarily ask for one more jet. The original factorials and signs remain those of the existing H_s owner and do not affect the estimate.

On a fixed neighborhood of mu>0 the finite-jet kernel difference is bounded uniformly, so pairing with delta_T costs only its bounded total mass. Outside that neighborhood the delta_T tails admit a T-independent exponential bound: for x=|E-mu|>=a>0,

    T^-1 exp(-x/T) <= C(a,T0) exp(-x/(2T0)).

Consequently finite polynomial energy growth is uniformly integrable. This proves a uniform O(rho^sigma) endpoint difference for the upper surface pairing without a hidden T^-s loss. The same argument applies successively to several compact loops after the angular Holder majorants are taken.

For a lower H_ell, ell>=1, transfer ell-1 derivatives and retain delta_T centered at zero. The remaining kernel jet is O(E^a) for some a>0 on a common high-D neighborhood. Its pairing tends to zero uniformly in eta: the local moment is O(T^a), while the remote tail is exponentially small. The endpoint difference obeys the corresponding rho^sigma majorant. This is a regulated limiting statement; it does not define a pointwise delta derivative times a singular origin factor.

If upper and lower derivatives occur on the same loop, use a smooth partition near E=0, near E=mu and on the complement. The two centers are separated by fixed mu>0. In each localized piece the other derivative weight has exponentially small overlap, which dominates every fixed inverse power of T. The total number of transferred derivatives never exceeds (s-1)_++(ell-1)_+, already included in J. The high-D radial estimates control the remaining kernel derivatives. This also explains why no coincident-mu0 endpoint is admitted.

These bounds justify commuting the real-T limit with the certified eta endpoint on the stated high-D neighborhood, after the independent complex line regulators have been removed in their required fixed-T order. They do **not** justify an unrestricted joint complex-regulator/T limit. Meromorphic continuation in D follows only after this common high-D pairing construction.
