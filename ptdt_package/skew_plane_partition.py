# -*- mode: sage -*-
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux
from sage.combinat.partition import Partition, Partitions
from sage.sets.family import Family
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.rings.all import NN
from sage.categories.infinite_enumerated_sets import \
    InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.structure.parent import Parent

class SkewPlanePartition(SkewTableau):

    r"""
    A class to model Skew Plane Partitions

    INPUT:

    - ``t`` -- a SkewTableau

    OUTPUT:

    - A SkewPlanePartition object constructed from ``t``.

    A skew plane partition is a skew tableau whose entries are non-negative
    integers which are weakly decreasing in rows and columns.
    Has a bijection with skew plane partitions.


    EXAMPLES::

    sage: from ptdt_package import *
    sage: SPP = SkewPlanePartition([[None, None, 2],[1,1,1],[0]]); SPP
    [[None, None, 2], [1, 1, 1]]
    sage: SPP.shape()
    [3, 3] / [2]
    sage: SPP.pp()
    .  .  2
    1  1  1
    sage: SPP.partition_size()
    5

    sage: SkewPlanePartition([])
    []

    TESTS::

    """

    @staticmethod
    def __classcall_private__(cls, t):
        if isinstance(t, SkewPlanePartition):
            return t
        if t not in SkewTableaux():
            raise ValueError("%s is not a skew tableau" % t)
        tableau = SkewTableau(t)
        shape = tableau.inner_shape()
        SPP = SkewPlanePartitions_all(shape)
        return SPP.element_class(SPP, t)

    def __init__(self, parent, t):
        if t not in parent:
            raise ValueError("%s is not an element of %s" % (t, parent))
        t = [[x for x in row if x != 0] for row in t
             if any(x != 0 for x in row)]
        super(SkewPlanePartition, self).__init__(parent, t)

    def partition_size(self):
        return sum(x for r in self for x in r if x is not None)

class SkewPlanePartitions(SkewTableaux):

    r"""
    A factory class for the various classes of Skew Plane Partitions.

    INPUT:

    - ``shape`` -- the shape of the skew tableaux
    - ``size`` -- the size of the skew tableaux

    OUTPUT:

    - The appropriate class, after checking basic consistency tests.

    A skew plane partition is a skew tableau whose entries are non-negative
    integers which are weakly decreasing in rows and columns.

    EXAMPLES::

    sage: from ptdt_package import *
    sage: SPP = SkewPlanePartitions([2,1]);SPP.cardinality()
    +Infinity

    sage: SPP = SkewPlanePartitions([2,1], 1);SPP.cardinality()
    3

    sage: SPP = SkewPlanePartitions([3,2],1)
    sage: SPP.random_element()   #random
    [[None, None, None], [None, None], [1]]

    sage: SPP = SkewPlanePartitions([5,2],4);SPP
    Skew plane partitions with inner shape [5, 2] and size 4

    sage: SPP = SkewPlanePartitions([2,1],1);SPP.list()
    [[[None, None, 1], [None]],
     [[None, None], [None, 1]],
     [[None, None], [None], [1]]]


    sage: ([[None, None, 2], [None, None]]) in SkewPlanePartitions([2,2], 2)
    True
    sage: Tableau([[None, None, 2], [None, None]]) in SkewPlanePartitions([2,2], 2)
    True
    sage: ([[None, None, -2], [None, None]]) in SkewPlanePartitions([2,2], 2)
    False
    sage: ([[None, None, 2], [None, None]]) in SkewPlanePartitions([3,2], 2)
    False
    sage: ([1,1]) in SkewPlanePartitions([2],2)
    False

    sage: SPP = SkewPlanePartitions([5,1]);SPP.subset()
    Skew plane partitions with inner shape [5, 1]

    sage: SPP = SkewPlanePartitions([5,1]);SPP.subset(2)
    Skew plane partitions with inner shape [5, 1] and size 2

    TESTS::

    sage: SkewPlanePartition([[0,0,-1],[0,1]])
    Traceback (most recent call last):
    ...
    ValueError: [[0, 0, -1], [0, 1]] is not an element of Skew plane partitions with inner shape []

    sage: SkewPlanePartition([[None,1,1],[2,2]])
    Traceback (most recent call last):
    ...
    ValueError: [[None, 1, 1], [2, 2]] is not an element of Skew plane partitions with inner shape [1]

    """

    Element = SkewPlanePartition

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return SkewPlanePartitions_all(shape_part)
        elif size not in NN:
            raise ValueError("Size must be a non-negative integer")

        return SkewPlanePartitions_size(shape_part, size)

class SkewPlanePartitions_size(SkewPlanePartitions):
    def __init__(self, shape, size):
        self._shape = shape
        self._size = size
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _reverse_plane_partitions(self):
        from .reverse_plane_partition import ReversePlanePartitions
        extra = max(self._size, 1)
        width = (self._shape[0] if self._shape else 0) + extra
        shape = ([width] * extra +
                 [width - s for s in reversed(self._shape)])
        return ReversePlanePartitions(shape, self._size)

    def __iter__(self):
        return (self(rpp.to_SkewPlanePartition())
                for rpp in self._reverse_plane_partitions())

    def cardinality(self):
        return self._reverse_plane_partitions().cardinality()

    def random_element(self):
        rpp = self._reverse_plane_partitions().random_element()
        return self(rpp.to_SkewPlanePartition())

    def __contains__(self, x):
        return (x in SkewPlanePartitions_all(self._shape) and
                SkewPlanePartition(x).partition_size() == self._size)

    def _repr_(self):
        return ("Skew plane partitions with inner shape %s and size %d"
                % (self._shape, self._size))

class SkewPlanePartitions_all(SkewPlanePartitions,
                              DisjointUnionEnumeratedSets):
    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NN, partial(SkewPlanePartitions_size, shape))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True,
                                             keepkey=False,
                                             category=cat)

    def __contains__(self, x):
        if x not in SkewTableaux():
            return False
        if any(c not in NN and c is not None for r in x for c in r):
            return False
        tab = SkewTableau(x)
        lst = list(tab) + list(tab.conjugate())
        return (all(a is None or a >= b for r in lst
                    for a, b in zip(r, r[1:]))
                and tab.inner_shape() == self._shape)

    def subset(self, size=None):
        if size is None:
            return self
        return self._family[size]

    def _repr_(self):
        return "Skew plane partitions with inner shape %s" % self._shape
