from ptdt_package import *
"""
This file contains all the formulas from the paper that we saw.
It implements the weights and functions in a different way from what
we did before, so be careful if loading this at the same time as another
file.
"""

def chern_character(part, m, coeffs=None, invert=False):
    """ Compute a component of Chern character for a given partition

    Check that this actually satisfies the definition given:
    sage: from ptdt_package import *
    sage: R.<a,b> = QQ[]; c = -a-b
    sage: PS.<q> = R[[]]
    sage: part = Tableau([[4, 2, 1], [2, 1]])
    sage: f = (sum(chern_character(part, k)*q^k for k in range(20))
    ....:      + O(q^20))
    sage: p = prod(1 - exp(t * q) for t in [a, b, c])
    sage: T = lambda i, j, k: a * i + b * j + c * k
    sage: s = sum(exp(T(i,j,k) * q)
    ....:         for i, j in part.cells()
    ....:         for k in range(part[i][j]))
    sage: f == 1 - p * s
    True
    """
    if m == 0:
        return 1
    if coeffs is None:
        R = PolynomialRing(QQ, 'a,b')
        a, b = R.gens()
        coeffs = (a, b, -a-b)
    a, b, c = coeffs
    if invert:
        T = lambda i, j, k: a * i + b * j + c * (-1-k)
    else:
        T = lambda i, j, k: a * i + b * j + c * k
    return -sum((-1) ** len(s) * (T(i, j, k) + sum(s)) ** m
                for i, row in enumerate(part)
                for j, col in enumerate(row)
                if col is not None
                for k in range(col)
                for s in powerset(coeffs)) / factorial(m)

def chern_product(part, ks, coeffs=None, invert=False):
    if ks in ZZ:
        ks = (ks,)
    return prod(chern_character(part, k, coeffs, invert) for k in ks)

def new_weighted_sum(P, ks, invert=False, coeffs=None, prec=8):
    if coeffs is None:
        R = PolynomialRing(QQ, 'a,b')
        a, b = R.gens()
        coeffs = (a, b, -a-b)
    
    R = parent(sum(coeffs))
    q = R[['q']].gen()
    return sum(chern_product(pi, ks, coeffs) * (-q) ** size
               for size in range(prec)
               for pi in P.graded_component(size)) + O(q ** prec)

def PT(shape, ks=(), coeffs=None, prec=8):
    return new_weighted_sum(ReversePlanePartitions(shape),
                            ks, True, coeffs, prec)

def DT(shape, ks=(), coeffs=None, prec=8):
    return new_weighted_sum(SkewPlanePartitions(shape),
                            ks, False, coeffs, prec)

def DTn(shape, ks=(), coeffs=None, prec=8):
    """ Compute the normalized DT series, or DT' """
    return DT(shape, ks, coeffs, prec) / DT([], (), coeffs, prec)

# This runs infinitely, use C-c to terminate
# Just noticed that PT is always 0 when defined this way.
# Same goes for DT of anything non-empty
def check_all_0():
    for ks, part in cantor_product(Partitions(), Partitions()):
        if not ks:
            continue
        print "Checking k={}, lambda={}".format(ks, part)
        assert PT(part, ks, (1, -1, 0)) == 0

def odd_eisenstein(n, prec=8):
    if n < 0 or n % 2 == 0:
        raise ValueError("{} is not positive and odd".format(n))
    q = QQ[['q']].gen()
    return sum(d ** (n-1) * q ** i
               for i in range(1, prec)
               for d in divisors(i)) + O(q**prec)

def check_formulas():
    # make sure we work in the right ring
    R = QQ['a,b']
    PS = R[['q']]
    a, b = R.gens()
    q = PS.gen()
    c = -a-b
    D = (a + b) * (a + c) * (b + c)
    
    E = lambda n: PS(odd_eisenstein(n))
    monomial = SymmetricFunctions(QQ).monomial()
    m = lambda p: monomial(p).expand(3)(a,b,c)
    
    assert DTn([], 3) == -D * E(3)(-q), "DT3"
    assert DTn([], 4) == 0, "DT4"
    assert DTn([], 5) == -D/12 * (m([2]) + m([1,1])) * E(5)(-q), "DT5"
    assert DTn([], 6) == D**2 /2 * E(3)(-q) ** 2, "DT6"
    assert DTn([], 7) == (
        -D/360 * (m([4]) + 2*m([3,1]) + 3*m([2,2]) + 5*m([2,1,1]))
        * E(7)(-q)), "DT7"

inds = [0, 1, 2, 3]
varnames = ['x%s_%s' % (ind1, ind2)
            for ind1, _ in enumerate(inds)
            for ind2, _ in enumerate(inds)]
R1.<a,b> = QQ[]
R2 = PolynomialRing(R1, varnames)
gens = (R2(a),R2(b),R2(-a-b))
shape = Partition([2, 1, 1])
empty = Partition([])

def get_var(i, j):
    return R2('x%s_%s' % (i, j))

def get_eq(shape):
    shape = Partition(shape)
    DT1 = DT(shape, 3, gens)
    DT2 = sum(get_var(i, j) * DT(empty, ti, gens, 6) * PT(shape, tj, gens, 6)
              for i, ti in enumerate(inds)
              for j, tj in enumerate(inds))
    return DT1 - DT2

terms = [term
         for shape in Partitions(3)
         for term in get_eq(shape)]

unknowns = [get_var(i, j)
            for i, ti in enumerate(inds)
            for j, tj in enumerate(inds)]

mat = matrix([[t[g] for g in unknowns]
              for t in terms])
vec = vector(-t.constant_coefficient() for t in terms)

print "Equations: %s" % mat.nrows()
print "Rank: %s" % mat.rank()
print "Unknown: %s" % len(unknowns)
print "Solving..."
print mat.solve_right(vec)
