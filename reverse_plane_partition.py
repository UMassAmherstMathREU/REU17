from sage.combinat.tableau import Tableau
from hillman_grassl.py import Hillman_Grassl

class ReversePlanePartition(Tableau):
    '''
    
    '''

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
                "%s is not a Reverse Plane Partition" % t)

        RPP = ReversePlanePartitions(shape)
        return RPP.element_class(RPP, tableau)

    def to_HilmanGrasslTableau(self):
        return Hillman_Grassl(self)

class ReversePlanePartitions(Tableaux):
    Element = ReversePlanePartition

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return ReversePlanePartitions_all(shape_part)
        elif size not in NonNegativeIntegers():
            raise ValueError("Size must be a non-negative integer")

        return (super(ReversePlanePartitions, cls)
                .__classcall__(cls, shape_part, size))

    def __init__(self, shape, size):
        self._size = size
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        return (ReversePlanePartition(vec=v, shape=self._shape)
                for v in self._weighted_integer_vectors())

    def cardinality(self):
        return self._weighted_integer_vectors().cardinality()

    def random_element(self):
        v = self._weighted_integer_vectors().random_element()
        return ReversePlanePartition(vec=v, shape=self._shape)

    def __contains__(self, x):
        return (x in ReversePlanePartitions_all(self._shape) and
                sum(c for r in x for c in r) == self._size)

    def _repr_(self):
        return ("Reverse Plane Partitions of shape %s and size %d"
                % (self._shape, self._size))

class ReversePlanePartitions_all(DisjointUnionEnumeratedSets):

    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NonNegativeIntegers(),
                   partial(ReversePlanePartitions, shape))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True,
                                             keepkey=False,
                                             category=cat)

    def __contains__(self, x):
        return (x in Tableaux() and
                all(c in NonNegativeIntegers()
                    for r in x for c in r) and
                Tableau(x).shape() == self._shape)

    def subset(self, size=None):
        if size is None:
            return self
        return self._family[size]

    def _repr_(self):
        return "Reverse Plane Partitions of shape %s" % self._shape
