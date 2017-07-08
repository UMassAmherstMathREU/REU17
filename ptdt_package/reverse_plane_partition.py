# -*- mode: sage -*-

from sage.combinat.tableau import Tableau, Tableaux
from sage.combinat.partition import Partition, Partitions
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.rings.all import NN
from sage.categories.infinite_enumerated_sets import \
    InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.structure.parent import Parent
from .hillman_grassl import Hillman_Grassl

class ReversePlanePartition(Tableau):
    @staticmethod
    def __classcall_private__(cls, t):
        if isinstance(t, ReversePlanePartition):
            return t
        if t not in Tableaux():
            raise ValueError("%s is not a tableau" % t)
        tableau = Tableau(t)
        shape = tableau.shape()
        if t not in ReversePlanePartitions_all(shape):
            raise ValueError(
                "%s is not a reverse plane partition" % t)

        RPP = ReversePlanePartitions(shape)
        return RPP.element_class(RPP, tableau)

    def partition_size(self):
        return sum(c for r in c for r in self)

    def to_HilmanGrasslTableau(self):
        from hillman_grassl_tableau import HillmanGrasslTableau
        return HillmanGrasslTableau(Hillman_Grassl(self))

    def to_SkewPlanePartition(self):
        from .skew_plane_partition import SkewPlanePartition
        if not self:
            return SkewPlanePartition([])
        width = len(self[0])
        return SkewPlanePartition([
            [None] * (width - len(row)) + list(reversed(row))
            for row in reversed(self)])

class ReversePlanePartitions(Tableaux):
    Element = ReversePlanePartition

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return ReversePlanePartitions_all(shape_part)
        elif size not in NN:
            raise ValueError("Size must be a non-negative integer")

        return ReversePlanePartitions_size(shape_part, size)

class ReversePlanePartitions_size(ReversePlanePartitions):
    def __init__(self, shape, size):
        self._size = size
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _hillman_grassl(self):
        from hillman_grassl_tableau import HillmanGrasslTableaux
        return HillmanGrasslTableaux(self._shape, self._size)

    def __iter__(self):
        return (self(hg.to_ReversePlanePartition())
                for hg in self._hillman_grassl())

    def cardinality(self):
        return self._hillman_grassl().cardinality()

    def random_element(self):
        hg = self._hillman_grassl().random_element()
        return self(hg.to_ReversePlanePartition())

    def __contains__(self, x):
        return (x in ReversePlanePartitions_all(self._shape) and
                sum(c for r in x for c in r) == self._size)

    def _repr_(self):
        return ("Reverse Plane Partitions of shape %s and size %d"
                % (self._shape, self._size))

class ReversePlanePartitions_all(ReversePlanePartitions,
                                 DisjointUnionEnumeratedSets):
    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NN, partial(ReversePlanePartitions, shape))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True,
                                             keepkey=False,
                                             category=cat)

    def __contains__(self, x):
        if x not in Tableaux():
            return False
        if any(c not in NN for r in x for c in r):
            return False
        tab = Tableau(x)
        if tab.shape() != self._shape:
            return False
        rows_cols = list(tab) + list(tab.conjugate())
        return all(a <= b for row in tab for a, b in zip(row, row[1:]))

    def subset(self, size=None):
        if size is None:
            return self
        return self._family[size]

    def _repr_(self):
        return "Reverse Plane Partitions of shape %s" % self._shape
