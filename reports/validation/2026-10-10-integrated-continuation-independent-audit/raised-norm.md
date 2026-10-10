# A separate absolute C2 component norm

This derives a conservative bound for the actual raised mixed-medium target at D=13/2, mu0=mu1=1. It uses the native source polynomial, not the independent finite-eta reference. All quantities below are stripped of the common native hard-master and compact angular normalization, as in the producer.

Let x be the second Schwinger parameter on the unit simplex. The emitted mean is k=x(q0-q1). The complete original numerator q0.k+E0 k0 therefore becomes

    N=x[a0+E0^2-2E0 E1+r0 r1 c],

where a0=q0^2, E_i=sqrt(r_i^2+a_i), and c is the relative spatial cosine. At zero masses,

    N0=x[r0^2+r0 r1(c-2)],
    N_a=x[2-r1/r0],
    h_a=-1+r1/r0,  h=-(q0-q1)^2.

The first raw mass derivative of the shell contribution has bulk

    [-N0/(2r0^2)+N_a] Phi + N0 h_a Phi_h

and upper surface -N0 Phi/(2r0) at r0=1. The lower term uses the existing joint high-D origin continuation. These expressions retain the complete original numerator and the residue/upper derivatives; an overall native target sign is immaterial to an absolute bound.

For |w|=r<1, define

    B0=(1+r)^(5/4)/(1-r),
    B1=(r/4)[(5/16)(1+r)^(1/4)/(1-r)+(1+r)^(5/4)/(1-r)^2].

The source-transformed kernel satisfies |Phi|<=B0 and |Phi_h|<=B1: its virtual h coefficient is x(1-x)<=1/4, and the pure-transfer coefficient is1. The compact polynomial bounds are

    |N0|<=x(r0^2+3r0r1),
    |N_a|<=x(2+r1/r0),
    |h_a|<=1+r1/r0.

The positive radial measure is r0^(7/2) r1^(7/2)dr0dr1, and the angular measure is normalized. Put m_e=1/(9/2+e). Integration of the absolute bulk/surface terms gives

    C_bulk=(5/2)(m_0^2+m_-1*m_1)=1580/6237,
    C_surface=(m_0+3m_1)/2=38/99,
    C_derivative=4(m_2*m_0+m_1^2)=3808/14157.

Multiplying by integral_0^1 x dx=1/2 and the native ordinary I2/I1 ratio9/4 gives

    M_raised(r) <= (1987/2772) B0 + (476/1573) B1.

All inverse-radial factors are integrable at this dimension. This bound does not infer positivity from the signed leading coefficient. As an internal sign check, setting w=0, using the normalized angular average c=0 and keeping the signed bulk/surface terms yields

    (9/8)[(3/2)m_0^2 - m_0/2 + m_1] = 43/264,

exactly the frozen leading raised coefficient. This finite check is supplementary; the absolute bound follows from the displayed component inequalities.

For N known transformed coefficients and 0<=w<r<1, an absolute tail bound is M_raised(r)(w/r)^N/(1-w/r). If the exact rational partial sum S_N satisfies |S_N|>T_N for a rigorous upper tail enclosure T_N, then relative error is at most T_N/(|S_N|-T_N). The common physical normalization and R^alpha w^-alpha cancel from this relative ratio. At w=1 there is still no geometric bound. No evaluation or endpoint acceptance is asserted by this note.
