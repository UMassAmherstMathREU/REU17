from ptdt_package import *

inds = [(), 0, 1, 2]
varnames = ['x%s_%s' % (ind1, ind2)
            for ind1, _ in enumerate(inds)
            for ind2, _ in enumerate(inds)]

known = {}

R1 = QQ
R2 = PolynomialRing(R1, varnames)
gens = (R2(1), R2(-1), R2(0))
shape = Partition([2, 1, 1])
empty = Partition([])
prec = 8

@cached_method
def PT(shape, k=()):
    print "PT %s, %s" % (shape, k)
    return weighted_sum(gens, k, shape, 'pt', prec)

@cached_method
def DT(shape, k=()):
    print "DT %s, %s" % (shape, k)
    return weighted_sum(gens, k, shape, 'dt', prec)

def get_var(i, j):
    tup = inds[i], inds[j]
    if tup in known:
        return known[tup]
    return R2('x%s_%s' % (i, j))

def get_eq(shape):
    shape = Partition(shape)
    DT1 = DT(shape, 2)
    DT2 = sum(get_var(i, j) * DT(empty, ti) * PT(shape, tj)
              for i, ti in enumerate(inds)
              for j, tj in enumerate(inds))
    return DT1 - DT2

terms = [term
         for shape in Partitions(3)
         for term in get_eq(shape)]

unknowns = [get_var(i, j)
            for i, ti in enumerate(inds)
            for j, tj in enumerate(inds)
            if (ti, tj) not in known]

mat = matrix([[t[g] for g in unknowns]
              for t in terms])
vec = vector(-t.constant_coefficient() for t in terms)

print "Equations: %s" % mat.nrows()
print "Rank: %s" % mat.rank()
print "Unknown: %s" % len(unknowns)
print "Solving..."
print mat.solve_right(vec)
