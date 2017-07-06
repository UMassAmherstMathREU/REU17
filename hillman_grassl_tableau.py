# -*- mode: sage -*-

from sage.combinat.tableau import Tableau, Tableaux
from sage.combinat.partition import Partition
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.rings.all import NN
from sage.categories.infinite_enumerated_sets import \
    InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading

class HillmanGrasslTableau(Tableau):
    @staticmethod
    def __classcall_private__(cls, t=None):
        if isinstance(t, HillmanGrasslTableau):
            return t
        if t not in Tableaux():
            raise ValueError("%s is not a tableau" % t)
        t = Tableau(t)
        if t not in HillmanGrasslTableaux_all(t.shape()):
            raise ValueError(
                "%s is not a Hillman-Grassl tableau" % t)
        HG = HillmanGrasslTableaux_all(t.shape())
        return HG.element_class(HG, t)

    def hg_size(self):
        shape = self.shape()
        return sum(t[i][j] * shape.hook_length(i, j)
                   for i, j in shape.cells())

    def to_ReversePlanePartition(self):
        from reverse_plane_partition import ReversePlanePartition
        from hillman_grassl import Inverse_HG
        return ReversePlanePartition(Inverse_HG(self))

class HillmanGrasslTableaux(Tableaux):
    Element = HillmanGrasslTableau

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return HillmanGrasslTableaux_all(shape_part)
        elif size not in NonNegativeIntegers():
            raise ValueError("Size must be a non-negative integer")

        return HillmanGrasslTableaux_size(shape_part, size)

    def from_integer_vector(self, vec):
        if vec not in IntegerVectors():
            raise ValueError("%s is not an integer vector" % vec)
        if self._shape.size() != len(vec):
            raise ValueError("%s is not length %d" %
                             (vec, self._shape.size()))
        tableau = []
        pos = 0
        for rowlen in self._shape:
            tableau += [vec[pos:(pos + rowlen)]]
            pos += rowlen

        if tableau not in self:
            raise ValueError("%s is not in %s" % (tableau, self))

        return self.element_class(self, tableau)

class HillmanGrasslTableaux_size(HillmanGrasslTableaux):
    def __init__(self, shape, size):
        self._size = size
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        return ("Hillman-Grassl tableaux of shape %s and size %d"
                % (self._shape, self._size))

    def __contains__(self, x):
        return (x in HillmanGrasslTableaux_all(self._shape) and
                HillmanGrasslTableau(x).hg_size() == self._size)

    def _weighted_integer_vectors(self):
        hooks = sum(self._shape.hook_lengths(), [])
        return WeightedIntegerVectors(self._size, hooks)

    def __iter__(self):
        return (self.from_integer_vector(vec)
                for vec in self._weighted_integer_vectors())

    def cardinality(self):
        return self._weighted_integer_vectors().cardinality()

    def random_element(self):
        vec = self._weighted_integer_vectors().random_element()
        return self.from_integer_vector(vec)
 
class HillmanGrasslTableaux_all(HillmanGrasslTableaux,
                                DisjointUnionEnumeratedSets):
    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NN, partial(HillmanGrasslTableaux, shape))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True,
                                             keepkey=False,
                                             category=cat)

    def __contains__(self, x):
        return (x in Tableaux() and
                all(c in NN for r in x for c in r) and
                Tableau(x).shape() == self._shape)

    def subset(self, size=None):
        if size is None:
            return self
        return self._family[size]

    def _repr_(self):
        return "Hillman-Grassl tableaux of shape %s" % self._shape
