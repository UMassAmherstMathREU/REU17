"""
This file contains methods to compute the weights of reverse plane
partitions, and skew plane partitions.  Given fixed integers a, b, c
such that a + b + c = 0, and some integer m, the weight of a partition
pi is the sum of (ai)^m + (bj)^m + (ck)^m over all (i, j, k) in pi.
The index of (i, j, k) starts at (0, 0, 0).
"""

from sage.combinat.partition import Partition
from sage.combinat.skew_tableau import SkewTableau
from hillman_grassl import Inverse_HG

def partition_weight(coefficients, powers, partition):
    """Calculate the weight of the partition."""
    a, b, c = coefficients
    try:
        m = Integer(powers)
    except TypeError:
        # powers is not an integer, try a tuple
        try:
            powers_tuple = tuple(powers)
        except TypeError:
            raise TypeError(
                "Powers must be an integer or iterable")
        # product over all elements of powers
        return prod(partition_weight(coefficients, m, partition)
                    for m in powers_tuple)

    partition_tab = SkewTableau(partition)
    # m is now the only power
    return sum((a * i) ** m + (b * j) ** m + (c * k) ** m
               for i, j in partition_tab.cells()
               for k in range(partition_tab[i][j]))

def ReversePlanePartitions(shape, n):
    """Generate a list of all reverse plane partitions"""
    shape_part = Partition(shape)
    cells = shape_part.cells()
    weights = [shape_part.hook_length(i, j) for i, j in cells]
    for vec in WeightedIntegerVectors(n, weights):
        hg = []
        s = 0
        for p in shape_part:
            hg += [vec[s:(s+p)]]
            s += p
        yield Inverse_HG(hg)

def pt_weighted_sum(coefficients, powers, shape, prec, base = None):
    """Calculate the weighted sum over weight*q^size"""
    if base is None:
        base = ZZ
    R = PowerSeriesRing(base, 'q')
    q = R.gen()
    return sum(
        partition_weight(coefficients, powers, part) * q ** size
        for size in range(prec)
        for part in ReversePlanePartitions(shape, size)) + R([], prec)

def pt_weighted_sum_abc(powers, shape, prec):
    """Calculate the sum, but leave in symbols a, b, c"""
    base = PolynomialRing(ZZ, 'a,b,c')
    return pt_weighted_sum(base.gens(), powers, shape, prec, base)
