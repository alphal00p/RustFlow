This follow-up is a source-presentation assessment only. It neither changes production admission nor reports native closure.

For a massless label with L loop momenta, cut indices n_c≥1, ordinary physical powers a_j, and homogeneous polynomial completion factors F_b of momentum degree d_b, sum the existing local dilation IBPs over all L loops. Assume all uncut physical factors are D_j=q_j²−eta, the cut factors are q_c², the upper factors are mu_c−E_c, and lower occupations have index zero. The lower terms vanish in the existing jointly continued origin prescription. Ordinary completion indices remain nonpositive.

Write s_c for an upper-occupation index and k(s)=−1 for s=0, k(s)=s for s>0. Before applying upper multiplication identities the exact row is

`[L D − 2 sum_c n_c − 2 sum_j a_j − sum_b d_b a_b] I`

`− 2 eta sum_j a_j I_(a+e_j) + sum_c k(s_c) I_(s+e_c)[E_c] = 0`.

Using E_c=mu_c−h_c, h_c H_1=0 and h_c H_(s+1)=H_s for s≥1 gives

`[L D − 2 sum_physical a_j − sum_b d_b a_b − sum_c s_c] I`

`− 2 eta partial_eta I + sum_c k(s_c) mu_c I_(s+e_c) = 0`.

Here partial_eta I is exactly the existing weighted derivative sum over shifted uncut slots; it is not a new source variable or external formula. Equivalently the last sum is −sum_c mu_c partial_(mu_c) I, with the other parameters and original factor powers held fixed. The cut term is valid for every finite n_c: it uses q_c² C_(n_c+1)=C_(n_c), retaining the raised-cut distribution contact. No C1-only simplification was made. Negative ordinary physical indices are allowed; their finite polynomial numerator degree is included by the same a_j coefficients.

This identity holds for coupled compact and virtual momenta. Compact-only dilation does not by itself have this simple form: a mixed factor (A_compact+A_virtual)²−eta acquires cross Gram insertions. Adding the virtual dilations cancels that complication and yields the global formula. If some physical masses are nonzero, or some uncut factors are not shifted, the corresponding extra mass/shift terms must be retained. A completion must have an exact Euler degree; an inhomogeneous completion cannot silently be assigned degree one or two.

As a potential source presentation, the summed row is finite, polynomial, and a literal linear combination of existing sources. Occupation branches must be intersected consistently, and every condition/zero-domain premise retained. It gives one relation among all upper surfaces. With two cuts it does not determine each surface separately. Its eta-derivative terms may still be the native leading terms, yielding a recurrence into surfaces rather than reducing them. Thus no claim of taming occupation growth follows from the identity alone; any additional presentation would need matched native measurements. There is no new mathematical source content or new endpoint theorem here.

Recompleting ordinary numerator coordinates in the occupied local Gram/energy coordinates is mathematically legitimate. `InversePropagatorBasis::new` already performs exact rational affine completion while preserving every physical slot, including dependent slots. Passing the unchanged routed physical factors and the occupied coordinates would select local independent Gram coordinates and each local energy, rather than inherit mixed original-loop energies. The routing determinant, shells, chemical endpoints, physical-edge identities and cut indices do not change merely because the completion is changed.

The necessary conversion is a complete linear map on represented targets. For each old term, retain its physical denominator powers, reconstruct the polynomial product of old completion factors, and convert that product once in the new physical-endpoint basis. Keep the old Wick coefficient and routing normalization exactly once. Terms whose required-cut index becomes nonpositive vanish by the existing cut-role identity. Do not interpret them as ordinary integrals. Preserve raw coefficient nonzero conditions before cancellations, and keep all new completion indices nonpositive. Dependent physical slots still need their original slot identity and power convention.

Conversion must be done before eta deformation and with independent physical masses still available. If a change of Gram completion produces physical inverse factors, the new fixed-index extension can differ from the old one at finite eta even though their physical endpoints agree. A pure linear change among energy completions has no such eta issue. In either case it is incorrect to reuse an old differential system or finite-eta rule corpus solely from endpoint equality. Original numerator mass differentiation must include all converted coefficient derivatives, or equivalently be performed in the old represented family and then mapped exactly. Replacing physical masses early can conceal those derivatives.

Implementation would need an explicit, versioned completion policy and a bound exact rebase map. The current sealed evidence compares the supplied family against the canonical occupied family, so a rebased family will correctly fail until that proof constructor validates the policy/map as well as the physical geometry. Source/cache identities, basis/target audits, storage arity, region conversion and boundary reconstruction must bind the new factors. Native closure must be rediscovered. Minimal tests would include rational shears/scalings, mixed occupied energies, a raised original medium numerator with nonzero mass-derivative coefficient terms, required-cut cancellation, both-way polynomial equivalence, and fresh old/new physical predictions. No topology name or special period value is needed.

For the current chain, the parent's separate rerouted-input control applies both an exact unimodular loop change and a physical-edge permutation. It also changes which numerator coordinates the original input completion selects. It is useful evidence for a representation effect, but cannot isolate completion from routing or edge-order changes.
