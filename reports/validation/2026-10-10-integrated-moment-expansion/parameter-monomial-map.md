# Exact parameter-monomial map to native vacuum integrals

Let h be the virtual rank, a_e the positive virtual indices, P=sum a_e, and lambda=hD/2-P. The actual HEPKit convention is F=eta U S+V, S=sum alpha_e for equally shifted masses. A numerator group with one Gaussian mean is written N_p(alpha,q)/U^p, with p=0 or1 in this prototype and homogeneous alpha degree hp.

At Taylor order n, expand N_p V^n into alpha monomials alpha^m times compact polynomials. Each such monomial has sum m=(h+1)n+hp. The coefficient is mapped to the unit-mass vacuum family at

    D' = D+2(n+p),     a'_e = a_e+m_e,

multiplied by

    (-1)^(n-sum m) / n! * product_e Gamma(a_e+m_e)/Gamma(a_e).

Indeed P'-hD'/2=P-hD/2+n, and the projective U power is -D/2-n-p, so both Schwinger kernels agree exactly. The gamma ratios are finite rising products of integers, not evaluated integral periods. The remaining eta exponent is lambda-n. Applying the generic Gram determinant dimension shift yields the fixed-D ordinary integral submitted to RustRed. The pure-compact shifted factors expand separately by their exact integer binomial series and are convolved; each contributes its explicit sign(-1)^a and eta^-a.

No coefficient in this derivation depends on there being only one compact invariant, one virtual denominator pair, or a known solution for the integrated graph. Virtual numerator degree two or above needs the additional Gaussian covariance terms; the prototype rejects it explicitly.
