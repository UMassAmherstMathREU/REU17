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
from hillman_grassl import Hillman_Grassl

class ReversePlanePartition(Tableau):
    r"""
    A class to model Reverse Plane Partitions

    INPUT:

    - ``t`` -- a Tableau

    OUTPUT:

    - A ReversePlanePartition object constructed from ``t``.

    description of ReversePlanePartition

    EXAMPLES::

        sage: RPP = ReversePlanePartition([[0,1,3],[2,4,4],[3]]); t
        [[0, 1, 3], [2, 4, 4], [3]]
        sage: RPP.shape()
        [3, 3, 2]
        sage: RPP.pp() # pretty print
        0 1 3
        2 4 4
        3
        sage: RPP.partition_size()
        17
        sage: RPP.to_HilmanGrasslTableau()
        [[1, 2, 0], [1, 0, 1], [1]]

        sage: ReversePlanePartition([]) # The empty tableau
        []

    TESTS::

    """
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

class ReversePlanePartitions(Tableaux):

    r"""
    A factory class for the various classes of Reverse Plane Partitions.

    INPUT:

    - ``shape`` -- the shape of the tableaux
    - ``size`` -- the size of the tableaux

    OUTPUT:

    - The appropriate class, after checking basic consistency tests.

    description of ReversePlanePartitions

    EXAMPLES::

    sage: RPP = ReversePlanePartitions([2,1]);RPP.cardinality()
    +Infinity
#Do these next
    sage: RPP = ReversePlanePartitions([4,2], 2);RPP.cardinality()
    5

    sage: RPP = ReversePlanePartitions([4,3],3)
    sage: RPP.random_element()   #random
    [[0, 0, 0, 1], [1, 0, 1]]

    sage: RPP = ReversePlanePartitions([5,2],4);RPP
    Reverse Plane Partitions of shape [5, 2] and size 4

    sage: RPP = ReversePlanePartitions([3,2],4);HG.list()
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

     sage: ([[0,0,1],[0,1]]) in ReversePlanePartitions([3,2], 2)
     True
     sage: Tableau([[0,0,1],[0,1]]) in ReversePlanePartitions([3,2], 2)
     True
     sage: ([[0,0,-1],[0,1]]) in ReversePlanePartitions([3,2], 2)
     False
     sage: ([[0,0,-1],[0,1]]) in ReversePlanePartitions([4,2], 2)
     False

     sage: RPP = ReversePlanePartitions([5,1]);RPP.subset()
     Reverse Plane Partitions of shape [5, 1]

     sage: RPP = ReversePlanePartitions([5,1]);RPP.subset(2)
     Reverse Plane Partitions of shape [5, 1] and size 2

    TESTS::

    """

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
        return (hg.to_ReversePlanePartition()
                for hg in self._hillman_grassl())

    def cardinality(self):
        return self._hillman_grassl().cardinality()

    def random_element(self):
        hg = self._hillman_grassl().random_element()
        rpp = hg.to_ReversePlanePartition()
        return self.element_class(self, rpp)

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
