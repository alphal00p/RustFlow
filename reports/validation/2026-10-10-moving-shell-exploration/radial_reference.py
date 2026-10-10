"""Independent quadrature checks for the exploratory fixed-spatial-cut flow.

Validation only: no production/oracle values, no fitted AMF boundary data.
Angular volume and (2*pi)^(-d) are omitted consistently in every radial value.
"""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp


def encode(x):
    if isinstance(x, (mp.mpf, mp.mpc)):
        return {"re": mp.nstr(mp.re(x), 95), "im": mp.nstr(mp.im(x), 95)}
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def integral(d, a, b, eta, moment=0):
    # r=sqrt(b)*t fixes the integration domain before mass differentiation.
    if b <= 0:
        return mp.mpf(0)
    return b ** (d / 2) / 2 * mp.quad(
        lambda t: t ** (d - 1) * (a + eta + b * t*t) ** (moment - mp.mpf("0.5")), [0, 1]
    )


def surface(d, a, b, eta):
    # H_1(b-r^2), not delta(mu-E): the radial Jacobian is included here.
    return b ** ((d - 2) / 2) / (4 * mp.sqrt(a + eta + b))


def asymptotic(d, a, b, eta, order):
    aa = a + eta
    prefactor = b ** (d / 2) / (2 * mp.sqrt(aa))
    total = mp.mpf(0)
    coefficient = mp.mpf(1)
    for n in range(order):
        total += coefficient * (b / aa) ** n / (d + 2*n)
        coefficient *= -(mp.mpf(n) + mp.mpf("0.5")) / (n + 1)
    return prefactor * total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    checks = []
    profiles = []
    for digits in (60, 90):
        mp.mp.dps = digits
        a, b = mp.mpf(1)/4, mp.mpf(3)/4
        local = []
        for d_string in ("7/5", "3", "9/2"):
            numerator, _, denominator = d_string.partition("/")
            d = mp.mpf(numerator) / mp.mpf(denominator or 1)
            for eta in (mp.mpf(0), mp.mpf(1)/3, mp.mpf(16), mp.mpc(1, mp.mpf(1)/2)):
                aa = a + eta
                value = integral(d, a, b, eta)
                upper = surface(d, a, b, eta)
                derivative = mp.diff(lambda e: integral(d, a, b, e), eta)
                ode = (d-1)/(2*aa)*value - b/aa*upper
                second = mp.diff(lambda e: surface(d, a, b, e), eta)
                k = integral(d, a, b, eta, 1)
                hyper = b**(d/2)/(2*d*mp.sqrt(aa))*mp.hyp2f1(mp.mpf("0.5"), d/2, 1+d/2, -b/aa)
                identity_errors = {
                    "radial_ibp_with_upper_surface": abs(derivative-ode),
                    "upper_surface_ode": abs(second+upper/(2*(aa+b))),
                    "medium_moment_reduction": abs(k-(2*b*(aa+b)*upper+aa*value)/(d+1)),
                    "independent_hypergeometric_check": abs(value-hyper),
                }
                assert max(identity_errors.values()) < mp.mpf(10)**(-digits+8)
                local.append(value)
                records.append(encode({"digits": digits, "spatial_dimension": d_string, "mass_squared": a,
                    "radius_squared": b, "eta": eta, "I": value, "S": upper, "I_prime": derivative,
                    "errors": identity_errors}))
            # At eta=0 vary the physical mass at fixed mu and original N=g+u^2.
            eta = mp.mpf(0)
            value = integral(d, a, b, eta)
            upper = surface(d, a, b, eta)
            eta_derivative = mp.diff(lambda e: integral(d, a, b, e), eta)
            physical = mp.diff(lambda mass: integral(d, mass, 1-mass, 0), a)
            medium = mp.diff(lambda mass: mass*integral(d, mass, 1-mass, 0)
                             + integral(d, mass, 1-mass, 0, 1), a)
            medium_bulk = mp.mpf("1.5")*value + a*eta_derivative
            medium_surface = -(a+1)*upper
            assert abs(physical-(eta_derivative-upper)) < mp.mpf(10)**(-digits+8)
            assert abs(medium-(medium_bulk+medium_surface)) < mp.mpf(10)**(-digits+8)
            # An eta derivative alone is demonstrably not physical raising.
            assert abs(physical-eta_derivative) > mp.mpf("0.01")
            checks.append(encode({"digits": digits, "spatial_dimension": d_string,
                "physical_scalar_mass_derivative": physical, "eta_scalar_derivative": eta_derivative,
                "scalar_upper_surface": -upper, "physical_medium_mass_derivative": medium,
                "medium_bulk": medium_bulk, "medium_upper_surface": medium_surface,
                "original_numerator": "g+u^2; g=q^2 Minkowski, u=q^0",
                "fixed_input_coefficient_eta_extension": "a*I + K, K=int E_eta/2",
                "eta_only_raising_is_incorrect": True}))
        profiles.append(local)
    mp.mp.dps = 90
    refinements = [abs(a-b)/max(abs(b),mp.mpf(1)) for a,b in zip(*profiles)]
    assert max(refinements) < mp.mpf("1e-55")
    a,b,d=mp.mpf(1)/4,mp.mpf(3)/4,mp.mpf(3)
    boundaries=[]
    for start in (8, 16):
        exact=integral(d,a,b,start)
        for order in (16, 32, 64):
            estimate=asymptotic(d,a,b,start,order)
            # Positive-real alternating tail decreases in magnitude.
            coefficient=mp.binomial(2*order,order)/4**order
            bound=b**(d/2)/(2*mp.sqrt(a+start))*coefficient*(b/(a+start))**order/(d+2*order)
            error=abs(estimate-exact)
            assert error <= bound
            boundaries.append(encode({"start_eta":start,"order":order,"estimate":estimate,
                "exact_quadrature":exact,"absolute_error":error,"alternating_tail_bound":bound}))
    # A massless physical derivative requires a stronger initial convergence strip.
    d=mp.mpf(9)/2
    value=integral(d,mp.mpf(0),mp.mpf(1),mp.mpf(0))
    assert abs(value-1/(2*(d-1))) < mp.mpf("1e-80")
    derivative_bulk=-mp.quad(lambda r:r**(d-4),[0,1])/4
    massless={"spatial_dimension":"9/2","I_at_eta_zero":value,
        "physical_mass_derivative":derivative_bulk-mp.mpf(1)/4,
        "known_from_direct_convergent_radial_integration":-(d-2)/(4*(d-3)),
        "order_of_limits":"first eta->0 at Re(d)>3 for one physical mass derivative, then meromorphic dimension continuation"}
    assert abs(massless["physical_mass_derivative"]-massless["known_from_direct_convergent_radial_integration"]) < mp.mpf("1e-80")
    # The mixed sunset's leading hard subgraph retains external invariant p^2=1.
    sunset=[]
    D=mp.mpf(12)/5; d=D-1; a=b=mp.mpf(1)/4; c=mp.mpf(1); radius=1-a
    exponent=D/2-2
    prefactor=mp.gamma(2-D/2)/(4*mp.pi)**(D/2)
    leading=prefactor*mp.quad(lambda x:(1-x*(1-x))**exponent,[0,1])
    vacuum_leading=prefactor
    assert abs(leading-vacuum_leading) > mp.mpf("0.01")*abs(leading)
    for eta in (16,64,256,1024):
        bubble=prefactor*mp.quad(lambda x:(x*(b+eta)+(1-x)*(c+eta)-x*(1-x)*(a+eta))**exponent,[0,1])
        weighted=integral(d,a,radius,eta)*bubble
        expected=radius**(d/2)/(2*d)*leading
        scaled=weighted*mp.mpf(eta)**(-(D-5)/2)
        sunset.append(encode({"eta":eta,"mixed_scalar_radial_value":weighted,
            "scaled_value":scaled,"leading_coefficient":expected,"relative_remainder":abs(scaled/expected-1)}))
    errors=[mp.mpf(row["relative_remainder"]["re"]) for row in sunset]
    assert all(a>b for a,b in zip(errors,errors[1:]))
    output={"schema":1,"status":"passed","scope":"validation-only exploratory one-loop flow and scalar mixed sunset asymptotic; no production acceptance",
        "angular_normalization":"radial integrals omit Omega_(d-1)/(2*pi)^d; the virtual bubble retains (4*pi)^(-D/2)",
        "mpmath_version":mp.__version__,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "one_loop":records,"physical_mass_derivatives":checks,"max_precision_refinement_error":encode(max(refinements)),
        "convergent_large_eta_boundaries":boundaries,"massless_endpoint":encode(massless),
        "mixed_sunset":sunset,"leading_hard_coefficient":encode(leading),"incorrect_zero_external_momentum_coefficient":encode(vacuum_leading),
        "production_reference_values_loaded":False,"four_loop_acceptance":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"status":"passed","one_loop_points":len(records),"physical_mass_checks":len(checks),
        "boundary_checks":len(boundaries),"mixed_sunset_scalings":len(sunset)}))


if __name__ == "__main__":
    main()
