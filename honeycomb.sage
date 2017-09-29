from ptdt_package import *

# Makes nice-ish ascii tables of PT(z) and DT(z) generating functions
# The number in each slot is the coefficient of e^Tz, where T are coordinates
# on the grid.

R.<x,y> = LaurentPolynomialRing(ZZ)

def symbolic_chern(partition):
    s = sum(x^(i-k)*y^(j-k)
            for i, row in enumerate(partition)
            for j, entry in enumerate(row)
            if entry is not None
            for k in range(entry))
    t = sum(x^i*y^j
            for i, row in enumerate(partition)
            for j, entry in enumerate(row)
            if entry is None)
    return R(1) - (1 - x) * (1 - y) * (t + (1 - x^(-1)*y^(-1)) * s)

def symbolic_chern_pt(partition):
    s = sum(x^(i+k)*y^(j+k)
            for i, row in enumerate(partition)
            for j, k in enumerate(row))
    return R(1) - (1 - x) * (1 - y) * s

def chern_total(shape, n, pt=False):
    P = ReversePlanePartitions(shape, n) if pt else SkewPlanePartitions(shape, n)
    chern = symbolic_chern_pt if pt else symbolic_chern
    return sum(chern(part) for part in P)

def ascii_table(shape, n, pt=False):
    """
    Print an ascii table representing [q^n]PT(q,z) or [q^n]DT(q,z).

    Each entry in the table corresponds to a coefficient of e^Tz, where
    T is the coordinates as a combination of t1, t2, t3.

    shape is an array/tuple/partition representing the shape of the infinite leg.
    Pass shape=[] to get DT_empty.

    n is the coefficient of q^n to pick, i.e. the size of the partitions to consider.

    If pt=True, this generates the table for PT, if not, it does for DT.
    """
    p = chern_total(shape, n, pt)
    top = min(i + j for i, j in p.exponents())
    bottom = max(i + j for i, j in p.exponents())
    nrows = bottom - top + 1
    left = min(j - i for i, j in p.exponents())
    right = max(j - i for i, j in p.exponents())
    ncols = right - left + 1
    table = [[str(p[(r-c+top-left)/2,(r+c+top+left)/2]) if (r + c + top + left) % 2 == 0 else ""
              for c in range(ncols)]
             for r in range(nrows)]
    maxlen = max(len(s) for row in table for s in row)
    table = [[s + " " * (maxlen - len(s))
              for s in row]
             for row in table]
    for row in table:
        print "  ".join(row)

print "EXAMPLES:"
print "q^4 coefficient of DT_empty:"
ascii_table([], 4, False)

print "q^3 coefficient of PT_1:"
ascii_table([1], 3, True)

print "q^5 coefficient of DT_2,1:"
ascii_table([2,1], 5, False)
