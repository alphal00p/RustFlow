Partial fixed-shell placements are algebraically supported by the weighted source/deformation owners, but are not currently admitted by the occupied flow/boundary proof. Two narrow candidates in the actual three-loop family are worth a bounded source-only comparison before implementation: singleton cut [3] with shifted physical slot [0], and double cut [0,3] with shifted slots [1,2] or [2,4]. These are diagnostic candidates selected from explicit factor geometry, not a topology-name dispatch or a claim of numerical validity already implemented.

The current restrictions occur in three distinct places. `flow.rs::validate_options` rejects MassMode other than All/Auto, and preparation constructs the mask from every uncut physical slot. `flow_boundary.rs::constants` independently requires exactly that mask. `massless_endpoint.rs::MasslessFlowEvidence` binds the all-uncut channel certificate and uses the shifted list as the complete physical row list in structural and finite-label bounds. The last assumption is significant: simply allowing a smaller mask would incorrectly omit the undeformed denominator powers from P/jet/support audits. `OccupiedCutFamily::deformed_measure`, the native derivative owner, ordinary region expansion and compact-compatible region enumeration already accept partial masks at the algebraic level.

The practical boundary restriction is equally real. All-uncut deformation leaves polynomial soft coefficients; the existing integrated occupied owner can integrate those or prove an unrestricted polynomial virtual loop scaleless. Partial placement can retain soft denominators. At present even a translated massless ordinary soft tadpole is rejected as nonpolynomial. More general partial placements retain coupled compact/virtual soft amplitudes, requiring decreasing weighted recursion or new generic compact integration owners. They cannot be treated as zero merely because a loop is soft.

For original cut [3], take q=P1−P3, t=P3 and l=P2. The uncut factors are

* slot0: (t+q)²;
* slot1: l²;
* slot2: t²;
* slot4: (l−t)².

Shifting only slot0 makes every undeformed denominator independent of q in these coordinates. At large eta, the all-hard coefficient is an ordinary two-loop vacuum with one massive and three massless factors, multiplied by compact polynomial moments. Every region with an actual virtual soft denominator has only a homogeneous massless virtual vacuum factor after a rational translation; the compact variables enter its numerator polynomial only. Such a factor is scaleless in the ordinary virtual measure. This is a candidate extension of the existing boundary zero proof, not permission to classify the full occupied integral with ordinary vacuum heuristics. The current production virtual coordinates need not equal (t,l), so the removal of compact dependence must be proved by an exact virtual affine shift; it must not be inferred from a name or from q²=0 alone.

There is a useful exact finite-eta germ check. Give slots0,2,1,4 Schwinger parameters a,b,c,d. Direct Gaussian elimination gives

`U=(a+b)(c+d)+cd`,

`V=a[b(c+d)+cd]`,

`F=eta*a*U + q_E²*V`.

The identity `a*U−V=a²(c+d)` proves `0≤V/(aU)≤1` on the positive interior. Thus the normalized q_E²/eta germ has a uniform nonzero unit for a small complex disk, despite the mass polynomial vanishing on parameter faces. A production permit must obtain and replay the actual native full-kinematic U/F and account for these massless faces through the same UV/IR meromorphic parameter continuation used by the existing singleton germ, with independent line regulators and fixed-T removal order. This Gaussian calculation is an independent assessment, not a replacement for that owner.

For finite labels with positive slot0 power and full virtual rank, the same polynomial inequality restricts to every active support. If slot0 has no positive power, the remaining virtual denominators can be made q-independent and have no mass scale; polynomial numerator factors from slot0 do not restore a virtual scale. Those sectors need an explicit ordinary-virtual scaleless certificate, not accidental master promotion. Rank-deficient active supports still require the existing free-virtual-direction proof. All positive undeformed powers must enter the finite-label degree and parameter-face bounds.

Shifting only slot2 is less immediately suitable. Its mass polynomial is bU, whereas V contains acd. The ratio V/(bU) is unbounded on b→0 with a,c,d positive. This does not prove the placement impossible, but it loses the simple uniform external-cone germ and leaves a genuinely q-dependent massless soft subgraph. It should not inherit the slot0 proof.

For double cut [0,3], use q=P1, r=P1−P3, p=q−r=P3 and k=P2. The only uncut physical factors are slot1 k², slot2 p², and slot4 (k−p)². Write h=−p²≥0. Shifting the transfer slot2 and either one virtual line gives a compact factor −(h+eta), a bubble with one mass eta and one zero mass, and the unchanged occupied shells.

For shifted slots [1,2], with bubble parameters x,y on slots1,4,

`U=x+y`, `V=xy`, `F=eta*x*U+h*V`,

and `x*U−V=x²` proves `0≤V/(xU)≤1`. The [2,4] case exchanges x and y. The virtual simplex massless endpoint remains and requires a finite-label high-D/meromorphic treatment; it must not be called a uniform massive gap on the closed simplex. Nonetheless the bounded ratio provides a direct candidate germ around h=0. The shifted pure transfer preserves a strict eta gap at compact origins. If the shifted virtual line has nonpositive power, the remaining positive virtual support is a translated single massless propagator, hence an ordinary scaleless virtual sector.

At infinity these two placements have only a one-loop massive/massless hard vacuum coefficient times compact polynomials. Their virtual-soft coefficient is a translated massless tadpole times compact polynomials and is zero by the same narrowly scoped virtual-vacuum proof. This is the smallest apparent boundary extension: recognize certified massless rational virtual-soft vacuum zeros, while continuing to reject retained compact denominators and all nonzero weighted soft amplitudes.

By contrast, shifting only the two virtual lines [1,4] leaves h to a negative power in every relevant compact soft coefficient. Shifting only the transfer [2] leaves a massless virtual bubble with compact external p. Neither fits the present polynomial compact boundary owner. They would require a genuine compact-transfer moment/weighted-recursion extension; evaluating a known bubble formula in production is not an acceptable substitute.

A generic partial-placement policy could first require that every undeformed positive denominator, after an exact rational virtual translation, is a homogeneous massless quadratic independent of compact momenta. It would then audit every generated region/support: either a surviving hard factor goes through ordinary RecursiveBoundary, a compact soft coefficient is polynomial, or a separated ordinary virtual soft factor has a replayable scaleless proof. All remaining cases fail explicitly. The class-specific contour/origin evidence must also certify the actual weighted mass polynomial `U*sum_shifted(alpha)` and full kinematic coefficients for every consumed positive support; keep shifted slots separate from the complete physical row inventory. The old pure-power/log assertion should be retained only after that region proof, with formal epsilon-dependent exponents unchanged.

A fair next reduction-only experiment can regenerate original polynomial sources for each candidate mask, bind a fresh context/deformation identity, retain all original polynomial/completion/storage guards, and compare original-target plus actual weighted eta-derivative closure under identical budgets. Derivatives raise only the selected shifted slots. Old rules cannot be imported. Lower-contact zero domains may be retained only after the partial-placement origin proof has been explicitly reviewed/bound; alternatively a preliminary algebra-only test must retain those contacts and disclose the weaker source corpus. Native closure alone would not authorize numerical flow. Before a full extension, required gates include exact native U/F replay, finite-label audits covering absent shifted factors, the rational virtual-soft zero proof, systematic region rank/depth matching, and independent complete lower-loop amplitudes with raised original numerators and surfaces.

For already certified strictly massive heavy-edge or positive-sunset domains, adding eta to only a subset of uncut masses preserves the underlying nonnegative real contour contribution. That observation does not solve their retained soft boundary factors. The current three-loop massless candidates need the stronger weighted-mass/origin analysis above; no blanket statement that every subset is admissible is justified.
