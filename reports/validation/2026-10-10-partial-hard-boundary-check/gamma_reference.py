"""Gamma via shifted Stirling with exact Bernoulli numbers; reused validated implementation."""
from decimal import Decimal as D
from fractions import Fraction as Q

def bernoulli(n):
    a = []
    for m in range(n + 1):
        a.append(Q(1, m + 1))
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
    return a[0]

def atan_inverse(q, tolerance):
    x = D(1) / q
    power = x
    total = power
    for k in range(1, 10000):
        power *= -x * x
        term = power / (2 * k + 1)
        total += term
        if abs(term) < tolerance:
            return total
    raise ValueError('pi series exhausted')

def gamma_positive(z, digits):
    tolerance = D(10) ** (-digits - 10)
    pi = 16 * atan_inverse(D(5), tolerance) - 4 * atan_inverse(D(239), tolerance)
    w = z + 128
    value = (w - D('0.5')) * w.ln() - w + (2 * pi).ln() / 2
    previous = D('Infinity')
    for k in range(1, 100):
        b = bernoulli(2 * k)
        term = (D(b.numerator) / b.denominator) / (2 * k * (2 * k - 1) * w ** (2 * k - 1))
        if abs(term) >= previous:
            raise ValueError('Stirling remainder stopped decreasing')
        if abs(term) < tolerance:
            remainder = abs(term)
            break
        value += term
        previous = abs(term)
    else:
        raise ValueError('Gamma work limit exhausted')
    for j in range(128):
        value -= (z + j).ln()
    return value.exp(), remainder, k
