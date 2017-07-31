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
        T = a * i + b * j - (value) * c
        return sum((-1)**len(s) * (T + sum(s)) ** m
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

def odd_eisenstein(n, k=0, prec=8):
    if n < 0 or n % 2 == 0:
        raise ValueError("{} is not positive and odd".format(n))
    q = QQ[['q']].gen()
    E =  sum(d ** (n-1) * q ** i
             for i in range(1, prec)
             for d in divisors(i)) + O(q**prec)
    for i in range(k):
        E = q * E.derivative()
    return E

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
    k2 = tuple(k for k in k2 if k != 0)
    if l < 0:
        return 0
    dt = DT([], k1, prec=prec).change_ring(c.parent())
    pt = PT(shape, k2, prec=prec).change_ring(c.parent())
    return dt * pt

def solve_matrix(solve_for, use_all_inds=False, remainder=False, prec=8):
    inds = [p for s in range(sum(solve_for) + 1)
            for p in Partitions(s, min_part=2)]
    indstr = ["".join(str(x) for x in ind) if ind else "0"
              for ind in inds]
    varnames = ['x%s_%s' % (ind1, ind2)
                for ind1 in indstr
                for ind2 in indstr]

    known = {(i, j): 0
             for i, ti in enumerate(inds)
             for j, tj in enumerate(inds)
             if sum(ti) + sum(tj) != sum(solve_for)}

    print solve_for
    print len(known)
    for i, ti in enumerate(inds):
        for j, tj in enumerate(inds):
            if DT([], ti, prec=prec) == 0 or (PT([2,1], tj, prec=prec) == 0 and
                                              PT([1,1], tj, prec=prec) == 0 and
                                              PT([1], tj, prec=prec) == 0 and
                                              PT([], tj, prec=prec) == 0):
                known[i, j] = 0
    print len(known)

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
        if remainder:
            n1, n2 = solve_for
            DT1 -= sum(ptdt_prod(shape, (k1, k2), (n1-k1, n2-k2), n1+n2, gens, prec)
                       for k1 in range(n1+1)
                       for k2 in range(n2+1))
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
    print len(known)
    print len(inds)
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
    parts_to_solve = [p for i in range(1,20)
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
    c = -a-b
    R2.<q> = R1[[]]
    assert DTn(shape) == PT(shape), "No insertions"
    assert DTn(shape, 1) == 0, "DT(1)"
    assert PT(shape, 1) == 0, "PT(1)"
    assert DTn(shape, 2) == -PT(shape, 2), "DTn(2)"
    assert DTn(shape, 3) == (-PT(shape, 3)
                             + DTn([], 3) * PT(shape)), "DTn(3)"
    assert DTn(shape, 4) == -PT(shape, 4), "DTn(4)"
    assert DTn(shape, 5) == (-PT(shape, 5)
                             -DTn([], 3) * PT(shape, 2)
                             +DTn([], 5) * PT(shape)), "DTn(5)"
    assert DTn(shape, 6) == (-PT(shape, 6)
                             -DTn([], 3) * PT(shape, 3)
                             +DTn([], 6) * PT(shape)), "DTn(6)"
    # Coefficients for PT are 2^n/n!
    # Same in fact for DT
    # I think we have a conjecture!

# General case -- this is really cool!
def check_pt_dt_formula(shape, n, prec=10):
    assert DT(shape, n) == (
        sum(DT([], k, prec=prec) * PT(shape, n-k, prec=prec)
            for k in range(n+1)))

def check_n_3_formula(shape, n, prec=10):
    R.<a,b> = QQ[]
    gens = (a,b,-a-b)
    DTPT = lambda k1, k2: ptdt_prod(shape, k1, k2, n+3, gens, prec)
    assert DT(shape, (n, 3)) == sum(DTPT((k,3), (n-k,)) + DTPT((k,), (n-k,3))
                                    for k in range(n+1))

def check_z_series(shape):
    qprec = 12
    zprec = 12
    R1.<a,b> = QQ[]
    R2.<q, z1, z2> = R1[[]]
    DT0 = sum(DT([], (k1, k2), prec=qprec) * z1^k1 * z2^k2
              for k1 in range(zprec) for k2 in range(zprec)) + O(z1^zprec + z2^zprec)
    DTl = sum(DT(shape, (k1, k2), prec=qprec) * z1^k1 * z2^k2
              for k1 in range(zprec) for k2 in range(zprec)) + O(z1^zprec + z2^zprec)
    PTl = sum(PT(shape, (k1, k2), prec=qprec) * z1^k1 * z2^k2
              for k1 in range(zprec) for k2 in range(zprec)) + O(z1^zprec + z2^zprec)
    print PTl
    print DTl - PTl * DT0

@cached_method
def get_remainder(shape, n1, n2, gens, prec=10):
    expected = sum(DT([], (k1, k2), prec=prec) * PT(shape, (n1-k1, n2-k2), prec=prec)
                   for k1 in range(n1 + 1) for k2 in range(n2 + 1))
    actual = DT(shape, (n1, n2), prec=prec)
    remainder = actual - expected
    a,b,c = gens
    D = (a+b)*(a+c)*(b+c)
    return remainder

def solve_remainder(n1, n2, prec=10, exclude=None):
    if exclude is None:
        exclude = []
    shapes = [(), (1,), (1,1), (2,1), (2,1,1), (2,2)]
    dt_params = [p for s in range(n1 + n2 + 1)
                 for p in Partitions(s, min_part=3)
                 if 4 not in p or len(p) > 1 and p[1] != 3
                 if p not in exclude]
    pt_params = [p for s in range(n1 + n2 + 1)
                 for p in Partitions(s, min_part=2)]
    params = [(p, dt, pt) for dt in dt_params for pt in pt_params
              if sum(dt) + sum(pt) <= n1 + n2
              for p in Partitions(n1+n2 - sum(dt) - sum(pt),
                                  max_length=1)
              if p != [1] and p != [4] and (p != [2] or dt != [4,4])]
    tup2str = lambda t: "_".join(str(i) for i in t)
    var_names = ["x%s__%s__%s" % (tup2str(p), tup2str(dt), tup2str(pt))
                 for p, dt, pt in params]
    V = PolynomialRing(QQ, var_names)
    R.<a,b> = V[]
    P.<q> = R[[]]
    c = -a-b
    gens = (a, b, c)
    
    print "Computing remainders"
    remainders = [get_remainder(shape, n1, n2, gens, prec) for shape in shapes]
    
    print remainders
    var_syms = [V(name) for name in var_names]
    print "Computing sum"
    m = SymmetricFunctions(QQ).monomial()
    sums = [sum(P(v) * m(p).expand(3)(a, b, c)
                * DT([], dt, prec=prec)
                * PT(shape, pt, prec=prec)
                for v, (p, dt, pt) in zip(var_syms, params))
            for shape in shapes]
    print "Computing difference"
    differences = [r - s for r, s in zip(remainders, sums)]
    terms = [t for d in differences for s in d for t, _ in R(s)]
    print "Building matrix"
    mat = matrix([[t[g] for g in var_syms]
                  for t in terms])
    vec = vector([-t.constant_coefficient() for t in terms])
    print 
    print "rank", mat.rank()
    print "unknowns", len(var_syms)
    print "Solving"
    if mat.rank() != len(var_syms):
        for i, v in zip(mat.right_kernel().matrix()[0], var_names):
            print v, i
    
    solution = mat.solve_right(vec)
    return {name: value for name, value in zip(var_names, solution)
            if value != 0}
    
def find_dt_relations(n, prec):
    params = [([n-s], dt)
              for s in range(n+1)
              if n - s != 1
              for dt in Partitions(s, min_part=3)
              if 4 not in dt or len(dt) > 1 and dt[1] != 3]
    m = SymmetricFunctions(QQ).monomial()
    tup2str = lambda t: "_".join(str(i) for i in t)
    var_names = ["x%s__%s" % (tup2str(p), tup2str(dt))
                 for p, dt in params]
    V = PolynomialRing(QQ, var_names)
    R.<a,b> = V[]
    P.<q> = R[[]]
    c = -a-b
    var_syms = [V(name) for name in var_names]
    m = SymmetricFunctions(QQ).monomial()
    print "Building equations"
    s = sum(P(v) * m(p).expand(3)(a, b, c) * DT([], dt, prec=prec)
            for v, (p, dt) in zip(var_syms, params))
    terms = [t for qterm in s for t, _ in qterm]
    mat = matrix([[t[g] for g in var_syms]
                  for t in terms])
    print "Finding kernel"
    K = mat.right_kernel()
    for row in K.matrix():
        print
        print "Found relation:"
        for (p, dt), value in zip(params, row):
            if value != 0:
                print "{} * m_{} * DT_0({})".format(value, p[0], dt)
        print "sums to 0"

def find_4_relation(n, prec):
    V.<x,y,z,w> = QQ[]
    R.<a,b> = V[]
    c = -a-b
    gens = (a, b, c)
    m2 = a^2 + b^2 + c^2
    s1 = m2 * sum(DT([], (k, n-k), gens, prec) * (-1)^k
                  for k in range(n+1))
    s2 = sum(DT([], (k, n-k+2), gens, prec) * (-1)^k
             for k in range(n+3))
    d1 = m2 * DT([], (4, n-4), gens, prec)
    d2 = DT([], (4, n-2), gens, prec)
    d3 = (a^4 + b^4 + c^4) * DT([], (4, n-6), gens, prec)
    eq = s1 + x * s2 + y * d1 + z * d2 + w * d3
    terms = [t for qterm in eq for t, _ in qterm]
    mat = matrix([[t[x], t[y], t[z], t[w]] for t in terms])
    vec = vector([-t.constant_coefficient() for t in terms])
    print mat.rank()
    print mat.solve_right(vec)

def check_4_formula(shape, n, prec, gens=None):
    actual = DTn(shape, (4, n), gens, prec=prec)
    nice_part = sum(DTn([], (k1, k2), gens, prec=prec)
                    * PT(shape, (4-k1, n-k2), gens, prec=prec)
                    for k1 in range(5) for k2 in range(n+1))
    if gens is None:
        R.<a,b> = QQ[]; c = -a-b
    else:
        a, b, c = gens
        R = QQ
    S.<q> = R[[]]
    m = lambda p: SymmetricFunctions(QQ).monomial()(p).expand(3)(a,b,c)
    OE = lambda i, j: odd_eisenstein(i, j, prec)(-q)
    F = [OE(3,1),
         0,
         m([2]) * OE(5,1) / 24,
         m([3]) * OE(3,1) * OE(3,0) / 3,
         m([4]) * OE(7,1) / 720,
         m([5]) * (OE(5,1) * OE(3,0) + OE(5,0) * OE(3, 1)) / 60,
         m([6]) * OE(9,1) / 40320
         + m([2,2,2]) * (OE(3,1) * OE(3,0)^2 - OE(9,1) / 8640) / 2,
         0,0,0]
    remainder = sum(k * PT(shape, k, gens, prec=prec) * F[n-2-k]
                    for k in range(n-1))
    print (actual - nice_part) / m([2,2,2])
    print
    print remainder
    print
    return (actual - nice_part) / m([2,2,2]) - remainder

rem = check_4_formula([1], 11, 13, (1, 1, -2)) / (1 + 1 + (-2)^7) / PT([1], 2, (1, 1, -2), 13)
oes = [q * prod(odd_eisenstein(i, prec=13) for i in p).derivative()
 for p in Partitions(10)
 if all(i % 2 == 1 for i in p)]

V = PolynomialRing(QQ, ['x%d' % d for d, _ in enumerate(oes)])
oes = [oe.change_ring(V) for oe in oes]
rem = rem.change_ring(V)
eq = rem - sum(V('x%d' % d) * oe for d, oe in enumerate(oes))
mat = matrix([[t[g] for g in V.gens()]
              for t in eq])
vec = vector([-t.constant_coefficient() for t in eq])
sol = mat.solve_right(vec)
