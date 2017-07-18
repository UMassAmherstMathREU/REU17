# -*- mode: sage -*-
"""
This file contains methods to compute the weights of reverse plane
partitions, and skew plane partitions.  Given fixed integers a, b, c
such that a + b + c = 0, and some integer m, the weight of a partition
pi is the sum of (ai)^m + (bj)^m + (ck)^m over all (i, j, k) in pi.
The index of (i, j, k) starts at (0, 0, 0).
"""

from sage.combinat.partition import Partition
from sage.combinat.skew_tableau import SkewTableau
from sage.structure.element import parent
from sage.rings.all import ZZ, Integer, O
from sage.misc.all import prod
from .reverse_plane_partition import *
from .skew_plane_partition import *

def partition_weight(coefficients, powers, partition, invert=False):
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

    if not isinstance(partition, Tableau):
        partition = SkewTableau(partition)

    # m is now the only power
    if invert:
        return sum((a * i + b * j + c * (-1-k)) ** m
                   for i, j in partition.cells()
                   for k in range(partition[i][j]))
    else:
        return sum((a * i + b * j + c * k) ** m
                   for i, j in partition.cells()
                   for k in range(partition[i][j]))

def weighted_sum(coefficients, powers, shape, domain='pt', prec=6):
    """Calculate the weighted sum over weight*q^size
    EXAMPLES::
        sage: from ptdt_package import *
        sage: R.<a,b,c> = ZZ[]
        sage: Z1 = weighted_sum((a,b,c), [], [], domain='dt', prec=6)
        sage: R.<q> = ZZ[[]]
        sage: Z2 = prod(1 / (1 - q^k)^k for k in range(1, 6)) + O(q^6)
        sage: Z1 == Z2
        True
        sage: weighted_sum((a,b,c), 1, [2, 1], prec=3)
        (a + b - 2*c)*q + (3*a + 3*b - 8*c)*q^2 + O(q^3)
    """
    if domain == 'pt':
        P = ReversePlanePartitions(shape)
        invert = True
    elif domain == 'dt':
        P = SkewPlanePartitions(shape)
        invert = False
    else:
        raise ValueError("Unknown domain (use pt or dt): %s" % domain)
    base_ring = parent(sum(coefficients))
    R = base_ring[['q']]
    q = R.gen()
    return sum(partition_weight(coefficients, powers, part, invert)
               * q ** size
               for size in range(prec)
               for part in P.graded_component(size)) + O(q ** prec)
