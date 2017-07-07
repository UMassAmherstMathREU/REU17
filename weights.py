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
from reverse_plane_partition import ReversePlanePartitions

def partition_weight(coefficients, powers, partition, side='dt'):
    """Calculate the weight of the partition. Argument side defaults to 'dt' and should be
       set to 'dt' for SkewPlanePartitions and 'pt' for ReversePlanePartitions. 

       Args:
           coefficients: A three element tuple (a,b,c) for changing the weight formula
           powers: A collection k of integers k_1, ..., k_n for computing the weight
           partition: A partition to compute the weight of
           side: A string specifying whether to compute the weight for a ReversePlanePartition('pt')
                 or a SkewPlanePartition('dt')

       Returns:
           The weight of the partition (an integer) given the powers and coefficients.
    """
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
        return prod(partition_weight(coefficients, m, partition, side=side)
                    for m in powers_tuple)

    partition_tab = SkewTableau(partition)
    # m is now the only power
    if side == 'dt':
        return dt_weight_formula(coefficients, m, partition_tab)
    else:
        return pt_weight_formula(coefficients, m, partition_tab)
    
def pt_weighted_sum(coefficients, powers, shape, prec, base = None):
    """Calculate the weighted sum over weight*q^size"""
    if base is None:
        base = ZZ
    R = PowerSeriesRing(base, 'q')
    q = R.gen()
    return sum(
        partition_weight(coefficients, powers, part, side='pt') * q ** size
        for size in range(prec)
        for part in ReversePlanePartitions(shape, size)) + R([], prec)

def pt_weighted_sum_abc(powers, shape, prec):
    """Calculate the sum, but leave in symbols a, b, c"""
    base = PolynomialRing(ZZ, 'a,b,c')
    return pt_weighted_sum(base.gens(), powers, shape, prec, base)

def dt_weight_formula(coefficients, power, partition):
    return sum((a * i) ** power + (b * j) ** power + (c * k) ** power
               for i, j in partition.cells()
               for k in range(partition[i][j]))

def pt_weight_formula(coefficients, power, partition):
    return sum((a * i) ** power + (b * j) ** power + (c * (-1-k)) ** power
               for i, j in partition.cells()
               for k in range(partition[i][j]))

def dt_weighted_sum(coefficients, powers, shape, prec, base=None):
    if base is None:
        base = ZZ
    R = PowerSeriesRing(base, 'q')
    q = R.gen()
    return sum(
        negatives = [-value for value in coefficients]
        partition_weight(negatives, powers, part, side='dt') * q**size
        for size in range(prec)
        for part in SkewPlanePartitions(shape, size)) + R([], prec)

# Temp Stuff for Skew Plane Partitions


def SkewPlanePartitions(outer_shape, inner_shape, n):
    """ Generate a list of all skew plane partitions"""
    # The shape of the plane partitions to be generated
    plane_shape = SkewPartition(outer_shape, inner_shape)
    cells = plane_shape.cells()

def funcToPlanePartition(func):
    """use cumulative sum to convert functions from shape to non-negative integers
       into skew plane partitions.

       Args:
         func: A list of lists containing non-negative integers.
    """
    shape = [len(row) for row in func]
    current = [row[-1] for row in func]
