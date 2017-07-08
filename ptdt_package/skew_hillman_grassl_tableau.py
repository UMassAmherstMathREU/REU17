# -*- mode: sage -*-

# This module might not be necessary...
# Keep it around just in case

from sage.combinat.skew_tableau import SkewTableau, SkewTableaux

class SkewHillmanGrasslTableau(SkewTableau):
    @staticmethod
    def __classcall_private__(cls, t):
        if isinstance(t, SkewHillmanGrasslTableau):
            return t
        if t not in SkewTableaux():
            raise ValueError("%s is not a skew tableau" % t)
        t = SkewTableau(t)
        shape = t.inner_shape()
        SHG = SkewHillmanGrasslTableaux_all(shape)
        if t not in SHG:
            raise ValueError(
                "%s is not a skew Hillman-Grassl tableau" % t)
        return SHG.element_class(SHG, t)

    def __init__(self, parent, t):
        height = max(i + 1 for i, row in enumerate(t)
                     if any(x != 0 for x in row))
        width = max(i + 1 for row in t for i, x in enumerate(row)
                    if x != 0)
        t = [list(row[:width]) + [0] * max(width - len(row), 0)
             for row in t[:height]]
        super(SkewHillmanGrasslTableau, self).__init__(parent, t)

    def hg_size(self):
        shape = self.inner_shape()
        return sum(self[i][j] * outer_hook(shape, i, j)
                   for i, j in self.cells())

    def to_HillmanGrasslTableau(self, nrows=None, ncols=None):
        if nrows is None:
            nrows = len(self) + 1
        if ncols is None:
            ncols = (len(self[0]) + 1) if self else 1
        if nrows < len(self):
            raise ValueError("%s has more than %d rows" % (self, nrows))
        if self and ncols < len(self[0]):
            raise ValueError("%s has more than %d columns"
                             % (self, nrows))
        return HillmanGrasslTableau(
            [[0] * ncols for i in range(nrows - len(self))]
            + [[0] * (ncols - len(row))
               + [x for x in reversed(row) if x is not None]
               for row in reversed(self)])

def outer_hook(part, i, j):
    conj = part.conjugate()
    ii = conj[j] if j < len(conj) else 0
    jj = part[i] if i < len(part) else 0
    return i - ii + j - jj + 1

class SkewHillmanGrasslTableaux(SkewTableaux):
    Element = SkewHillmanGrasslTableau

    @staticmethod
    def __classcall_private__(cls, shape, size=None):
        if shape not in Partitions():
            raise ValueError("Shape must be a partition")
        shape_part = Partition(shape)

        if size is None:
            return SkewHillmanGrasslTableaux_all(shape_part)
        elif size not in NN:
            raise ValueError("Size must be a non-negative integer")

        return SkewHillmanGrasslTableaux_size(shape_part, size)

class SkewHillmanGrasslTableaux_size(SkewHillmanGrasslTableaux):
    def __init__(self, shape, size):
        self._shape = shape
        self._size = size
        Parent.__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        return ("Skew Hillman-Grassl tableaux of shape %s and size %s"
                % (self._shape, self._size))

    def __contains__(self, x):
        return (x in SkewHillmanGrasslTableaux_all(self._shape) and
                SkewHillmanGrasslTableau(x).hg_size() == self._size)

    def _hook_tableu(self):
        arr = ([[None] * s for s in self._shape]
               + [[] for i in range(self._size)])
        for i, row in enumerate(arr):
            for j in range(len(row), len(row) + self._size):
                h = outer_hook(self._shape, i, j)
                if h > self._size:
                    break
                row.append(h)
        return SkewTableau(arr)

    def _weighted_integer_vectors(self):
        hooks = self._hook_tableu()
        weights = [h for r in hooks for h in r if h]
        return WeightedIntegerVectors(self._size, weights)

    def _from_integer_vector(self, vec):
        row_lens = self._hook_tableu().shape().row_lengths()
        arr = []
        pos = 0
        for i, l in enumerate(row_lens):
            nones = self._shape[i] if i < len(self._shape) else 0
            arr += [[None] * nones + vec[pos:(pos+l)]]
            pos += l
        return self.element_class(self, arr)
        
    def __iter__(self):
        return (self._from_integer_vector(vec)
                for vec in self._weighted_integer_vectors())
        
    def cardinality(self):
        return self._weighted_integer_vectors().cardinality()

    def random_element(self):
        vec = self._weighted_integer_vectors().random_element()
        return self._from_integer_vector(vec)

class SkewHillmanGrasslTableaux_all(SkewHillmanGrasslTableaux,
                                    DisjointUnionEnumeratedSets):
    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NN, partial(SkewHillmanGrasslTableaux_size, shape))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True,
                                             keepkey=False,
                                             category=cat)

    def _repr_(self):
        return "Skew Hillman-Grassl tableaux of shape %s" % self._shape

    def __contains__(self, x):
        return (x in SkewTableaux() and
                all(c in NN or c is None for r in x for c in r) and
                SkewTableau(x).inner_shape() == self._shape)

    def subset(self, size=None):
        if size is None:
            return self
        return self._family[size]
