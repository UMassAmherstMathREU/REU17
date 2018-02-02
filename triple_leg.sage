P.<s1, s2, s3> = QQ[]
# Power series ring, with coefficients in Q
PI.<z> = P[[]]
PS.<q> = PI[[]]

# Get the smallest monomial ideal with the given localizations
def minimal_ideal(px, py, pz):
    # Define Ix,Iy,Iz, corresponding to monomial ideals from legs
    # TODO: Think about conventions for which way partitions are oriented
    Ix = P.ideal([s2^a * s3^b for (a,b) in Partition(px).outside_corners()])
    Iy = P.ideal([s1^a * s3^b for (a,b) in Partition(py).outside_corners()])
    Iz = P.ideal([s1^a * s2^b for (a,b) in Partition(pz).outside_corners()])
    # Define I to be the intersection of these ideals:
    # This could be computed directly, using lcms, but sage has it built-in
    I = Ix.intersection(Iy, Iz)
    return I

def normalized_length(I):
    # Restrict to a box sufficiently large
    # "sufficiently large" means includes all generators of I
    xd = max(g.degree(s1) for g in I.gens())
    yd = max(g.degree(s2) for g in I.gens())
    zd = max(g.degree(s3) for g in I.gens())
    # Count the number of boxes (monomials not in I)
    count = sum(1
                for a in range(xd)
                for b in range(yd)
                for c in range(zd)
                if s1^a * s2^b * s3^c not in I)
    # Count the sizes of each infinite leg
    # Determine this by checking a monomial just outside the box
    xcyl_count = sum(xd for b in range(yd) for c in range(zd)
                     if s1^xd * s2^b * s3^c not in I)
    ycyl_count = sum(yd for a in range(xd) for c in range(zd)
                     if s1^a * s2^yd * s3^c not in I)
    zcyl_count = sum(zd for a in range(xd) for b in range(yd)
                     if s1^a * s2^b * s3^zd not in I)
    return count - xcyl_count - ycyl_count - zcyl_count

def multiplicity_of_box(i, j, k, ileg, jleg, kleg):
    m = 0
    if (j, k) in ileg.cells():
        m += 1
    if (i, k) in jleg.cells():
        m += 1
    if (i, j) in kleg.cells():
        m += 1
    return 1 - m if m >= 1 else 0

# When added to the infinite legs, this should give the generating function
# for the plane partition defined by I
def finite_base(I, px, py, pz):
    px = Partition(px)
    py = Partition(py)
    pz = Partition(pz)

    xd = max(g.degree(s1) for g in I.gens())
    yd = max(g.degree(s2) for g in I.gens())
    zd = max(g.degree(s3) for g in I.gens())

    return sum(multiplicity_of_box(i, j, k, px, py, pz) * x^i * y^j * z^k
               for i in range(xd)
               for j in range(yd)
               for k in range(zd))

def check_length_agrees(px, py, pz):
    I = minimal_ideal(px, py, pz)
    f = finite_base(I, px, py, pz)
    return f(1,1,1) == normalized_length(I)

def add_at_corner(G, g):
    G = [f for f in G if f != g]
    Irem = P.ideal(G)
    if g * s1 not in Irem:
        G.append(g * s1)
    if g * s2 not in Irem:
        G.append(g * s2)
    if g * s3 not in Irem:
        G.append(g * s3)
    return frozenset(G)

def chern_char(F, coeff, params):
    print "bar"
    g = (1 - (1 - s1) * (1 - s2) * (1 - s3) * F)
    # Composing power series is SLOW
    # We can avoid doing this if we don't expand e^(si * z)
    # But I think we have to if s1 + s2 + s3 != 0
    prec = coeff + 1
    R.<z> = parent(sum(params))[[]]
    print "foo"
    expanded = g([(t * z).exp(prec) for t in params])
    r = expanded[coeff]
    print parent(r)
    return r

def multiple_chern_char(F, insertions, params):
    return prod(chern_char(F, k, params) for k in insertions)

# f - finite part of Qa
# legs = [l1,l2,l3] - lists of generators for the 3 legs
def equiv_vertex_measure(f, legs, params=(s1,s2,s3)):
    print "Start evm"
    LP.<t1, t2, t3> = LaurentPolynomialRing(ZZ)
    t = [t1, t2, t3]
    F = P(f)(t1, t2, t3)
    L = [sum(t2^a * t3^b for (a, b) in legs[0].cells()),
         sum(t1^a * t3^b for (a, b) in legs[1].cells()),
         sum(t1^a * t2^b for (a, b) in legs[2].cells())]
    bar = lambda a: LP(a)(1/t1, 1/t2, 1/t3)
    # After doing the cancellation with the infinite parts,
    # this is what's left
    print "F = {}".format(F)
    V = (F - (bar(F) - F * bar(F) * (1-t1) * (1-t2) * (1-t3)) / (t1 * t2 * t3)
         + sum((L[i] * bar(F) - t[i] * F * bar(L[i]))
               * (1 - t[j]) * (1 - t[k]) / (t1 * t2 * t3)
               for i in range(3) for j in range(3) for k in range(3)
               if i != j and i != k and j < k)
         - sum(L[i] * bar(L[j]) * (1 - t[k]) / (t[i] * t[k])
               for i in range(3) for j in range(3) for k in range(3)
               if i != j and i != k and j != k))
    print "Almost done with evm"
    return prod((i * params[0] + j * params[1] + k * params[2])^(-V[i,j,k])
                for (i, j, k) in V.exponents())

# Given the shapes for the legs at infinity, calculate the
# non-normalized DT vertex
def vertex_series(px,py,pz, num_terms=5, insertions=(), params=(s1,s2,s3)):
    # Work directly with the generators, since sages's implementation
    # of ideals doesn't work well in a set/hashmap
    px = Partition(px)
    py = Partition(py)
    pz = Partition(pz)

    R.<q> = parent(sum(params))[[]]

    legs = [px, py, pz]
    I = minimal_ideal(px,py,pz)
    l = normalized_length(I)
    S = { (frozenset(I.interreduced_basis()), finite_base(I, px, py, pz)) }
    W = q^l
    for k in range(l+1, l+num_terms):
        S = { (add_at_corner(G, g), f + g) for (G, f) in S for g in G }
        W += sum(equiv_vertex_measure(f, legs, params) * multiple_chern_char(f, insertions, params)
                 for (G, f) in S) * q^k
    return W + O(q^(l+num_terms))

def check_macmahon_equiv(prec=6, params=(s1,s2,s3)):
    """ Check that DT = M(-q)^D """
    a, b, c = params
    D = -(a + b) * (a + c) * (b + c) / (a * b * c)
    # Compute M(-q)^D = prod 1/(1-(-q)^k)^Dk using binomial series
    print "Computing predicted result"
    predicted = prod(sum(binomial(-k * D, l) * (-(-q)^k)^l
                         for l in range(prec)
                         if l*k < prec)
                     for k in range(prec)) + O(q^prec)
    print "Computing actual result"
    actual = vertex_series([], [], [], prec, (), params)
    #print "predicted = {}".format(predicted)
    #print "actual = {}".format(actual)
    return predicted == actual

def check_1_box(prec=6, params=(s1,s2,s3)):
    a,b,c = params
    D = -(a + b) * (a + c) * (b + c) / (a * b * c)
    C = (b + c) / a
    M = prod(sum(binomial(-k * D, l) * (-(-q)^k)^l
                 for l in range(prec)
                 if l*k < prec)
             for k in range(prec)) + O(q^prec)
    N = sum(binomial(C, k) * q^k
            for k in range(prec)) + O(q^prec)
    predicted = N * M
    actual = vertex_series([1], [], [], prec, (), params)
    return predicted == actual

def check_minimal(px,py,pz):
    px = Partition(px)
    py = Partition(py)
    pz = Partition(pz)
    I = minimal_ideal(px,py,pz)
    G = I.interreduced_basis()
    return equiv_vertex_measure(G, [px,py,pz])
