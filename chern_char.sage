from ptdt_package import *
"""
This file contains all the formulas from the paper that we saw.
It implements the weights and functions in a different way from what
we did before, so be careful if loading this at the same time as another
file.
"""

def column_sum(value, i, j, m, coeffs, invert):
    a, b, c = coeffs
    if value is None:
        # Infinite column of DT
        T = a * i + b * j
        return sum((-1)**len(s) * (T + sum(s))**m
                   for s in powerset([a, b]))
    elif not invert:
        # Finite part of DT
        T = lambda k: a * i + b * j + c * k
        return sum((-1)**len(s) * (T(k) + sum(s))**m
                   for k in range(value)
                   for s in powerset([a, b, c]))
    else:
        # Infinite part of PT
        T = a * i + b * j - (2 + value) * c
        return -sum((-1)**len(s) * (T + sum(s)) ** m
                    for s in powerset([a, b]))

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
    return -sum(column_sum(val, i, j, m, coeffs, invert)
                for i, row in enumerate(part)
                for j, val in enumerate(row)) / factorial(m)
    
def chern_product(part, ks, coeffs=None, invert=False):
    if ks in ZZ:
        ks = (ks,)
    return prod(chern_character(part, k, coeffs, invert) for k in ks)

@cached_method
def all_pps(P, size):
    return tuple(P.graded_component(size))

@cached_method
def new_weighted_sum(P, ks, invert=False, coeffs=None, prec=8):
    if coeffs is None:
        R = PolynomialRing(QQ, 'a,b')
        a, b = R.gens()
        coeffs = (a, b, -a-b)
    
    R = parent(sum(coeffs))
    q = R[['q']].gen()
    return sum(chern_product(pi, ks, coeffs, invert) * (-q) ** size
               for size in range(prec)
               for pi in all_pps(P, size)) + O(q ** prec)

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

### Solve matrix stuff
@cached_method
def ptdt_prod(shape, k1, k2, n, gens, prec):
    a, b, c = gens
    l = n - sum(k1) - sum(k2)
    if l < 0:
        return 0
    coeff = (-1) ** len(k2) * (2 * c) ** l / factorial(l)
    dt = DT([], k1, prec=prec).change_ring(c.parent())
    pt = PT(shape, k2, prec=prec).change_ring(c.parent())
    return coeff * dt * pt

def solve_matrix(solve_for, prec=8):
    inds = [p for s in range(sum(solve_for)+1)
            for p in Partitions(s, min_part=2,
                                max_length=len(solve_for))]
    indstr = ["".join(str(x) for x in ind) if ind else "0"
              for ind in inds]
    varnames = ['x%s_%s' % (ind1, ind2)
                for ind1 in indstr
                for ind2 in indstr]

    known = {(i, j): 0
             for i, ti in enumerate(inds)
             for j, tj in enumerate(inds)
             if sum(ti) + sum(tj) > sum(solve_for)}
    
    for i, ti in enumerate(inds):
        for j, tj in enumerate(inds):
            if DT([], ti) == 0 or PT([1,1], tj) == 0 or 2 in ti:
                known[i, j] = 0

    R1 = PolynomialRing(QQ, varnames)
    R2.<a,b> = R1[]
    gens = (R2(a),R2(b),R2(-a-b))
    empty = Partition([])

    def get_var(i, j):
        if (i, j) in known:
            return known[i, j]
        return R1('x%s_%s' % (indstr[i], indstr[j]))

    def get_eq(shape):
        print "Getting equations for {}".format(shape)
        shape = Partition(shape)
        DT1 = DT(shape, solve_for, gens, prec)
        DT2 = sum(get_var(i, j) *
                  ptdt_prod(shape, ti, tj, sum(solve_for), gens, prec)
                  for i, ti in enumerate(inds)
                  for j, tj in enumerate(inds)
                  if (i, j) not in known)
        return DT1 - DT2

    terms = [R1(t)
             for shape in [(), (1,), (1, 1), (2, 1)]
             for term in get_eq(shape)
             for t, _ in term]
    unknowns = [get_var(i, j)
                for i, ti in enumerate(inds)
                for j, tj in enumerate(inds)
                if (i, j) not in known]
    mat = matrix([[t[g] for g in unknowns]
                  for t in terms])
    vec = vector(-t.constant_coefficient() for t in terms)

    print "Equations: %s" % mat.nrows()
    print "Rank: %s" % mat.rank()
    print "Unknown: %s" % len(unknowns)
    print "Solving..."
    solution = mat.solve_right(vec)
    if mat.rank() < len(unknowns):
        raise ValueError("Rank is too low")

    for i, a in enumerate(solution):
        if a != 0:
            print "{}: {}".format(unknowns[i], a)

    return {str(unk): a
            for unk, a in zip(unknowns, solution)
            if a != 0}

def find_all():
    # Find all formulas where the sum is under 10
    parts_to_solve = [p for i in range(1,10)
                      for p in Partitions(i, min_part=2)]
    f = open("dt_formulas", "w")
    for p in parts_to_solve:
        print "Working on {}".format(p)
        solution = solve_matrix(p)
        f.write("DT({}) = {}\n".format(p, solution))
    f.close()

### Check formulas
def known_formulas(shape):
    R1.<a,b> = QQ[]
    R2.<q> = R1[[]]
    assert DTn(shape) == PT(shape), "No insertions"
    assert DTn(shape, 1) == 0, "DT(1)"
    assert PT(shape, 1) == 0, "PT(1)"
    assert DTn(shape, 2) == -PT(shape, 2), "DTn(2)"
    assert DTn(shape, 3) == (-2 * c * PT(shape, 2)
                             - PT(shape, 3)
                             + DTn([], 3) * PT(shape)), "DTn(3)"
    assert DTn(shape, 4) == (-2 * c^2 * PT(shape, 2)
                             -2 * c * PT(shape, 3)
                             -PT(shape, 4)), "DTn(4)"
    assert DTn(shape, 5) == (-4/3 * c^3 * PT(shape, 2)
                             -2 * c^2 * PT(shape, 3)
                             -2 * c * PT(shape, 4)
                             -PT(shape, 5)
                             -DTn([], 3) * PT(shape, 2)
                             +DTn([], 5) * PT(shape)), "DTn(5)"
    assert DTn(shape, 6) == (-2/3 * c^4 * PT(shape, 2)
                             -4/3 * c^3 * PT(shape, 3)
                             -2 * c^2 * PT(shape, 4)
                             -2 * c * PT(shape, 5)
                             -PT(shape, 6)
                             -2 * c * DTn([], 3) * PT(shape, 2)
                             -DTn([], 3) * PT(shape, 3)
                             +DTn([], 6) * PT(shape)), "DTn(6)"
    assert DTn(shape, 7) == (
        sum(-(2*c)^n/factorial(n) * PT(shape, 7-n) for n in range(7))
        + DTn([], 7) * PT(shape)
        #- DTn([], 6) * PT(shape, 1) # == 0
        - DTn([], 5) * PT(shape, 2)
        #- DTn([], 4) * PT(shape, 3) # == 0
        - DTn([], 3) * PT(shape, 4)
        #- DTn([], 2) * PT(shape, 5) # == 0
        #- DTn([], 1) * PT(shape, 6) # == 0
        #- 2 * c * DTn([], 5) * PT(shape, 1) # == 0
        #- 2 * c * DTn([], 4) * PT(shape, 2) # == 0
        - 2 * c * DTn([], 3) * PT(shape, 3)
        #- 2 * c * DTn([], 2) * PT(shape, 4) # == 0
        #- 2 * c * DTn([], 1) * PT(shape, 5) # == 0
        #- 2 * c^2 * DTn([], 4) * PT(shape, 1) # == 0
        - 2 * c^2 * DTn([], 3) * PT(shape, 2)
        #- 2 * c^2 * DTn([], 2) * PT(shape, 3) # == 0
        #- 2 * c^2 * DTn([], 1) * PT(shape, 4) # == 0
        )
    # Coefficients for PT are 2^n/n!
    # Same in fact for DT
    # I think we have a conjecture!

# General case -- this is really cool!
def check_pt_dt_formula(shape, n, prec=10):
    assert DT(shape, n) == (
        DT([], n, prec=prec) * PT(shape, prec=prec)
        + sum(-(2*c)^(n-i-j)/factorial(n-i-j)
              * DT([], i, prec=prec) * PT(shape, j, prec=prec)
              for i in range(0, n+1)
              for j in range(1, n-i+1)))
