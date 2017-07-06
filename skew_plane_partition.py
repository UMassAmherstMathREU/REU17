from sage.combinat.skew_tableau import SkewTableau
from hillman_grassl.py import Hillman_Grassl

class SkewPlanePartition(SkewTableau):

    @staticmethod
    def __classcall_private__(cls, t):
        if isinstance(t, SkewPlanePartition):
            return t
        if t not in SkewTableaux():
            raise ValueError("%s is not a skew tableau" % t)
        tableau = SkewTableau(t)
        shape = tableau.shape()
        if t not in SkewPlanePartitions_all(shape):
            raise ValueError(
                "%s is not a Skew Plane Partition" % t)

class SkewPlanePartitions(SkewTableaux):
    Element = ReversePlanePartition

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        pass

    def __init__(self, shape, size):
        pass

    def __iter__(self):
        pass

    def cardinality(self):
        pass

    def random_element(self):
        pass

    def __contains__(self, x):
        pass

    def _repr_(self):
        return ("Skew Plane Partitions of shape %s and size %d"
                % (self._shape, self._size))

class SkewPlanePartitions_all(DisjointUnionEnumeratedSets):

    def __init__(self, shape):
        pass

    def __contains__(self, x):
        pass

    def _repr_(self):
        return "Skew Plane Partitions of shape %s" % self._shape
