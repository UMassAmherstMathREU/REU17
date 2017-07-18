from ptdt_package import *

prec = 20
R.<a,b> = ZZ[]
c = - a - b
coeffs = (0, 0, 1)
FF = FractionField(R)
L = LaurentSeriesRing(FF, 'q', default_prec=prec)
RatFunc = FractionField(R['q'])

def PT(shape, powers=()):
    return L(weighted_sum(coeffs, powers, shape, 'pt', prec))

def continued_frac(series):
    neg_tail = series.truncate(1)
    pos_tail = series - neg_tail
    q = RatFunc.gen()
    rat_tail = sum(c * q^e for c, e in zip(neg_tail.coefficients(),
                                           neg_tail.exponents()))
    if pos_tail != 0:
        rest = continued_frac(1/pos_tail)
        return rat_tail + 1 / rest
    else:
        return rat_tail

def euler_polynomial(n):
    if n not in NN:
        raise ValueError("n must be a natural number")
    if n == 0:
        return L(1)
    else:
        prev = euler_polynomial(n - 1)
        return (q * (1 - q) * prev.derivative()
                + prev * (1 + (n - 1) * q))

### Demonstrate formula for PT(box, n)
q = L.gen()
for n in range(1,10):
    pt1 = PT((1,), n)
    pt2 = c^n * q * euler_polynomial(n) / (q - 1)^(n+2)
    assert pt1 == pt2, "Formula fails for n=%d" % n

### Formula for PT(box, 1^n) = https://oeis.org/A154283
def strange_poly(n):
    s = -(q - 1)^(2*n + 1) * sum((k * (k + 1) / 2)^n * q^(k-1)
                                for k in range(2*n))
    return s.truncate(2*n-1)

for n in range(1,10):
    pt1 = PT((1,), (1,) * n)
    pt2 = -c^n * q^2 * strange_poly(n) / (q - 1)^(2*n+1)
    assert pt1 == pt2

### PT((1, 1), m)
PS = L.power_series_ring()
for n in range(1, 10):
    pt11 = PS(PT((1, 1), n))
    pt1 = PS(PT((1,), n))
    rem = a^n * q * (1 + 2*q) / (1 - q) / (1 - q^2)^2
    assert pt11 == pt1.V(2) + pt1 / (1 - q) + rem

### PT((1, 1, 1), m)
for n in range(1, 10):
    pt111 = PS(PT((1, 1, 1), n))
    pt1 = PS(PT((1,), n))
    A = (pt1 / (1 - q^2) / (1 - q)
         + pt1.V(2) / (1 - q)
         + pt1.V(3) * (-q^4 + 2*q^3 - 2*q + 1) / (1 - q^2) / (1 - q)^2)
    assert A == pt111

### experimenting
p = sum(j*q^(i+j) for j in range(20) for i in range(j+1)) + O(q^20)
rat_func = continued_frac(p)
if rat_func != 0:
    print rat_func.numerator().factor()
    print rat_func.denominator().factor()
else:
    print 0
