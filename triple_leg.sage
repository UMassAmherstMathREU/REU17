# 3-variable polynomial ring
P.<x,y,z> = QQ[]
# Quotient, used for chern character generating series
Q = P.quotient(x * y * z - 1)
# Power series ring, with coefficients in Q
PS.<q> = Q[[]]

# Get the smallest monomial ideal with the given localizations
def minimal_ideal(px, py, pz):
    # Define Ix,Iy,Iz, corresponding to monomial ideals from legs
    # TODO: Think about conventions for which way partitions are oriented
    Ix = P.ideal([y^a*z^b for (a,b) in Partition(px).outside_corners()])
    Iy = P.ideal([x^a*z^b for (a,b) in Partition(py).outside_corners()])
    Iz = P.ideal([x^a*y^b for (a,b) in Partition(pz).outside_corners()])
    # Define I to be the intersection of these ideals:
    # This could be computed directly, using lcms, but sage has it built-in
    I = Ix.intersection(Iy, Iz)
    return I

def normalized_length(I):
    # Restrict to a box sufficiently large
    # "sufficiently large" means includes all generators of I
    xd = max(g.degree(x) for g in I.gens())
    yd = max(g.degree(y) for g in I.gens())
    zd = max(g.degree(z) for g in I.gens())
    # Count the number of boxes (monomials not in I)
    count = sum(1
                for a in range(xd)
                for b in range(yd)
                for c in range(zd)
                if x^a * y^b * z^c not in I)
    # Count the sizes of each infinite leg
    # Determine this by checking a monomial just outside the box
    xcyl_count = sum(xd for b in range(yd) for c in range(zd)
                     if x^xd * y^b * z^c not in I)
    ycyl_count = sum(yd for a in range(xd) for c in range(zd)
                     if x^a * y^yd * z^c not in I)
    zcyl_count = sum(zd for a in range(xd) for b in range(yd)
                     if x^a * y^b * z^zd not in I)
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

    xd = max(g.degree(x) for g in I.gens())
    yd = max(g.degree(y) for g in I.gens())
    zd = max(g.degree(z) for g in I.gens())

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
    if g * x not in Irem:
        G.append(g * x)
    if g * y not in Irem:
        G.append(g * y)
    if g * z not in Irem:
        G.append(g * z)
    return frozenset(G)

# Just a matter of inclusion/exclusion
def chern_char(G):
    # Pass to quotient ring Q, since we're working without EVM
    return P(sum(G)
             - sum(lcm(gi, gj)
                   for i, gi in enumerate(G)
                   for j, gj in enumerate(G)
                   if i < j)
             + sum(lcm(gi, lcm(gj, gk))
                   for i, gi in enumerate(G)
                   for j, gj in enumerate(G)
                   for k, gk in enumerate(G)
                   if i < j and j < k))

def poincare_poly2(partition, var1, var2):
    return (sum(var1^a * var2^b for (a,b) in partition.outside_corners())
            - sum(var1^a * var2^b for (a,b) in partition.inside_corners())
            + 0 * var1 * var2) # hack to make sure it's the right ring

# f - finite part of Qa
# legs = [l1,l2,l3] - lists of generators for the 3 legs
def equiv_vertex_measure(f, legs):
    

# Given the shapes for the legs at infinity, calculate the
# non-normalized DT vertex
def vertex_series(px,py,pz, num_terms=5):
    # Work directly with the generators, since sages's implementation
    # of ideals doesn't work well in a set/hashmap
    px = Partition(px)
    py = Partition(py)
    pz = Partition(pz)

    legs = [px, py, pz]
    I = minimal_ideal(px,py,pz)
    l = normalized_length(I)
    S = { (frozenset(I.interreduced_basis()), finite_base(I, px, py, pz)) }
    W = (-q)^l
    for k in range(l+1, l+1+num_terms):
        S = { (add_at_corner(G, g), f + g) for (G, f) in S for g in G }
        W += sum(equiv_vertex_measure(f, legs) * chern_char(G) for (G, f) in S) * (-q)^k
    return W + O(q^prec)

def check_minimal(px,py,pz):
    px = Partition(px)
    py = Partition(py)
    pz = Partition(pz)
    I = minimal_ideal(px,py,pz)
    G = I.interreduced_basis()
    return equiv_vertex_measure(G, [px,py,pz])
