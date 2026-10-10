"""Independent all-uncut finite-eta reference; never imported by the producer.

The only target-specific reduction here is the equal-mass virtual first moment.
It is validation only. The driver must first bind saved producer coefficients.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dec(value):
    return mp.nstr(value, mp.mp.dps)


def radial_rule(n, exponent, mu):
    nodes, weights = mp.gauss_quadrature(n, "jacobi", 0, exponent)
    return [(mu * (nodes[i] + 1) / 2,
             weights[i] * (mu / 2) ** (exponent + 1)) for i in range(n)]


def angular_rule(n, alpha):
    nodes, weights = mp.gauss_quadrature(n, "jacobi", alpha - 1, alpha - 1)
    norm = 2 ** (2 * alpha - 1) * mp.beta(alpha, alpha)
    return [((nodes[i] + 1) / 2, weights[i] / norm) for i in range(n)]


def parameter_rule(n):
    nodes, weights = mp.gauss_quadrature(n, "legendre")
    return [((nodes[i] + 1) / 2, weights[i] / 2) for i in range(n)]


def quadrature_moments(radial, angular, parameter, exponent, alpha, mu):
    errors = []
    # Independent analytic moments check the transformed, normalized rules.
    for p in range(5):
        exact_r = mu ** (exponent + p + 1) / (exponent + p + 1)
        exact_t = mp.beta(alpha + p, alpha) / mp.beta(alpha, alpha)
        exact_x = mp.mpf(1) / (p + 1)
        errors += [abs(mp.fsum(w*r**p for r,w in radial) / exact_r - 1),
                   abs(mp.fsum(w*t**p for t,w in angular) / exact_t - 1),
                   abs(mp.fsum(w*x**p for x,w in parameter) / exact_x - 1)]
    worst = max(errors)
    assert worst < mp.mpf(10) ** (8 - mp.mp.dps), worst
    return {"checks": len(errors), "max_relative_error": dec(worst)}


def kernel(h, eta, nu, gaussian, parameter):
    values, derivatives = [], []
    for x,w in parameter:
        v = x*(1-x)
        base = eta + v*h
        assert base > 0  # Initial real-positive reference contract only.
        if nu == mp.mpf(5)/4:
            power = base * mp.sqrt(mp.sqrt(base))
        else:
            power = base ** nu
        values.append(w*power)
        derivatives.append(w*v*power/base)
    integral = mp.fsum(values)
    derivative = nu*mp.fsum(derivatives)
    denominator = eta+h
    return (gaussian*integral/denominator,
            gaussian*(derivative/denominator-integral/denominator**2))


def evaluate(dimension, eta, mu, n):
    D, eta, mu = map(mp.mpf, (dimension, eta, mu))
    assert D == mp.mpf(13)/2 and mu == 1 and eta > 0
    alpha, nu = (D-2)/2, D/2-2
    # The r^(D-4) weight makes all raised integrands polynomial times a smooth
    # kernel; multiplying r1*r2 restores the scalar r^(D-3) measure.
    radial = radial_rule(n, D-4, mu)
    angular = angular_rule(n, alpha)
    parameter = parameter_rule(n)
    checks = quadrature_moments(radial, angular, parameter, D-4, alpha, mu)
    area = 2*mp.pi**((D-1)/2)/mp.gamma((D-1)/2)
    A = area/(2*mp.pi)**(D-1)
    gaussian = (4*mp.pi)**(-D/2)*mp.gamma(2-D/2)
    scalar, measure, numerator, transfer = [], [], [], []
    # Exchange symmetry halves kernel evaluations. The two original target
    # contributions are averaged algebraically, not assumed separately equal.
    for i,(r1,w1) in enumerate(radial):
        for j in range(i,n):
            r2,w2 = radial[j]
            pair_weight = w1*w2*(1 if i == j else 2)
            product, squares = r1*r2, r1*r1+r2*r2
            for t,wt in angular:
                K,Kh = kernel(4*product*t, eta, nu, gaussian, parameter)
                w = pair_weight*wt
                scalar.append(w*K*product)
                measure.append(w*K*(squares*(t+mp.mpf('0.5'))/4-product/4))
                numerator.append(w*K*(product-squares/4))
                transfer.append(-w*Kh*product*(t+1)*(r1-r2)**2/2)
    normalization = A*A/4
    pieces = {"measure": normalization*mp.fsum(measure),
              "original_numerator_derivative": normalization*mp.fsum(numerator),
              "kernel_derivative": normalization*mp.fsum(transfer)}
    surface = []
    for r,w in radial:
        for t,wt in angular:
            K,_ = kernel(4*mu*r*t, eta, nu, gaussian, parameter)
            original_numerator = mu*r*(t+mp.mpf('0.5'))-mu*mu/2
            surface.append(w*wt*r*K*original_numerator)
    pieces["upper_surface"] = A*A*mu**(D-4)/8*mp.fsum(surface)
    scalar_value = normalization*mp.fsum(scalar)
    raised_value = mp.fsum(pieces.values())
    # This independently checks the direct parameter integral and its derivative.
    # Errors are reported, not silently absorbed into a changed tolerance.
    cross = []
    for h in (mp.mpf(0),mp.mpf(1),mp.mpf(4)):
        direct,direct_h = kernel(h,eta,nu,gaussian,parameter)
        analytic = lambda z: gaussian*eta**nu*mp.hyp2f1(-nu,1,mp.mpf(3)/2,-z/(4*eta))/(eta+z)
        value,derivative = analytic(h),mp.diff(analytic,h)
        cross.append({"h":dec(h),"value_relative_error":dec(abs(direct/value-1)),
                      "derivative_relative_error":dec(abs(direct_h/derivative-1))})
    return {"scalar":dec(scalar_value),"raised":dec(raised_value),
            "raised_components":{k:dec(v) for k,v in pieces.items()},
            "normalization":{"spatial_angular_measure":dec(A),"virtual_gaussian":dec(gaussian),
                             "raw_euclidean":True,"cut_residue_jacobian":"1/(2*E) for each cut",
                             "routing_determinant":"1"},
            "quadrature_moment_checks":checks,"parameter_cross_checks":cross}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    receipt=json.loads(args.receipt.read_text())
    assert receipt['prediction_artifacts_frozen_before_reference'] is True
    profile=receipt['profile']
    assert profile['nodes'] in (16,24,32) and profile['digits'] in (50,80)
    assert profile['dimension']=='6.5' and profile['mu']=='1'
    assert profile['eta'] in ('8','2','0.5')
    assert receipt['reference_evaluator_sha256']==sha(__file__)
    assert not args.output.exists()
    mp.mp.dps=profile['digits']
    started=time.monotonic()
    result=evaluate(profile['dimension'],profile['eta'],profile['mu'],profile['nodes'])
    result.update({"schema":1,"status":"reference quadrature completed; refinement not yet assessed",
                   "scope":"validation only; no producer equation or coefficient used",
                   "profile":profile,"receipt_sha256":sha(args.receipt),
                   "evaluator_sha256":sha(__file__),"mpmath_version":mp.__version__,
                   "elapsed_seconds":time.monotonic()-started,
                   "empirical_accuracy_only":True,"values":'raw Euclidean two-cut contribution'})
    temporary=args.output.with_suffix('.tmp')
    temporary.write_text(json.dumps(result,indent=2)+'\n')
    temporary.replace(args.output)


if __name__=='__main__':
    main()
