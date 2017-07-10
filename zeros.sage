from ptdt_package import *

def check_formula(n, shape, prec=6):
    R.<a,b,c> = ZZ[]
    DT1 = weighted_sum((a, b, c), [0] * n, shape, 'dt', prec)
    DT2 = sum(k *
              weighted_sum((a, b, c), [0] * i, [], 'dt', prec) *
              weighted_sum((a, b, c), [0] * j, shape, 'pt', prec)
              for (i, j), k in binomial_coefficients(n).items())
    assert DT1 == DT2
