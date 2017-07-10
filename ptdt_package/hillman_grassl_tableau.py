# -*- mode: sage -*-
from sage.combinat.tableau import Tableau, Tableaux
from sage.combinat.partition import Partition, Partitions
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.sets.family import Family
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.rings.all import NN
from sage.categories.infinite_enumerated_sets import \
    InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.structure.parent import Parent

class HillmanGrasslTableau(Tableau):
    r"""
    A class to model Hillman-Grassl Tableau

    INPUT:

    - ``t`` -- a Tableau

    OUTPUT:

    - A HillmanGrasslTableau object constructed from ``t``.

    A Hillman-Grassl tableau is a tableau whose entries are non-negative
    integers and whose size is a weighted sum over hook lengths.

    EXAMPLES::

        sage: from ptdt_package import *
        sage: hg = HillmanGrasslTableau([[1,3,2],[2,1]]); hg
        [[1, 3, 2], [2, 1]]
        sage: hg.shape()
        [3, 2]
        sage: hg.pp()
        1 3 2
        2 1

        sage: HillmanGrasslTableau([])
        []

        sage: HG = HillmanGrasslTableau([[1,2,0],[1,0,1],[1]])
        sage: HG.to_ReversePlanePartition()
        [[0, 1, 3], [2, 4, 4], [3]]
        sage: HG.hg_size()
        17

    TESTS::

        sage: HillmanGrasslTableau([[-1,1,2]])
        Traceback (most recent call last):
        ...
        ValueError: [[-1, 1, 2]] is not a Hillman-Grassl tableau
    """
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
        return sum(self[i][j] * shape.hook_length(i, j)
                   for i, j in shape.cells())

    def to_ReversePlanePartition(self):
        from reverse_plane_partition import ReversePlanePartition
        from hillman_grassl import Inverse_HG
        return ReversePlanePartition(Inverse_HG(self))

class HillmanGrasslTableaux(Tableaux):

    r"""
    A factory class for the various classes of Hillman-Grassl tableaux.

    INPUT:

    - ``shape`` -- the shape of the tableaux
    - ``size`` -- the size of the tableaux

    OUTPUT:

    - The appropriate class, after checking basic consistency tests.

    A Hillman-Grassl tableau is a tableau whose entries are non-negative
    integers and whose size is a weighted sum over hook lengths.

    EXAMPLES::

    sage: from ptdt_package import *
    sage: HG = HillmanGrasslTableaux([2,1]);HG.cardinality()
    +Infinity

    sage: HG = HillmanGrasslTableaux([4,2], 2);HG.cardinality()
    5

    sage: HG = HillmanGrasslTableaux([4,3],3)
    sage: HG.random_element()   #random
    [[0, 0, 0, 1], [1, 0, 1]]

    sage: HG = HillmanGrasslTableaux([5,2],4);HG
    Hillman-Grassl tableaux of shape [5, 2] and size 4

    sage: HG = HillmanGrasslTableaux([3,2],4);HG.list()
    [[[1, 0, 0], [0, 0]],
     [[0, 1, 0], [0, 1]],
     [[0, 1, 1], [0, 0]],
     [[0, 0, 0], [2, 0]],
     [[0, 0, 0], [1, 2]],
     [[0, 0, 1], [1, 1]],
     [[0, 0, 2], [1, 0]],
     [[0, 0, 0], [0, 4]],
     [[0, 0, 1], [0, 3]],
     [[0, 0, 2], [0, 2]],
     [[0, 0, 3], [0, 1]],
     [[0, 0, 4], [0, 0]]]

     sage: ([[0,0,1],[0,1]]) in HillmanGrasslTableaux([3,2], 2)
     True
     sage: Tableau([[0,0,1],[0,1]]) in HillmanGrasslTableaux([3,2], 2)
     True
     sage: ([[0,0,-1],[0,1]]) in HillmanGrasslTableaux([3,2], 2)
     False
     sage: ([[0,0,1],[0,1]]) in HillmanGrasslTableaux([4,2], 2)
     False

     sage: HG = HillmanGrasslTableaux([5,1]);HG.subset()
     Hillman-Grassl tableaux of shape [5, 1]

     sage: HG = HillmanGrasslTableaux([5,1]);HG.subset(2)
     Hillman-Grassl tableaux of shape [5, 1] and size 2

    TESTS::

    """

    Element = HillmanGrasslTableau

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return HillmanGrasslTableaux_all(shape_part)
        elif size not in NN:
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
        F = Family(NN, partial(HillmanGrasslTableaux_size, shape))
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
