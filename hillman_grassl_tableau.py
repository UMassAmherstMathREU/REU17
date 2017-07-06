# -*- mode: sage -*-

class HillmanGrasslTableau(Tableau):
    @staticmethod
    def __classcall_private__(cls, t=None, vec=None, shape=None):
        if ((t is None) == (vec is None) or
            (vec is None) != (shape is None)):
            raise ValueError(
                "Incorrect syntax for HillmanGrasslTableau")
        elif t is not None:
            if isinstance(t, HillmanGrasslTableau):
                return t
            if t not in Tableaux():
                raise ValueError("%s is not a tableau" % t)
            tableau = Tableau(t)
            shape = tableau.shape()
            if t not in HillmanGrasslTableaux_all(shape):
                raise ValueError(
                    "%s is not a Hillman-Grassl tableau" % t)
        else:
            if vec not in IntegerVectors():
                raise ValueError("%s is not an integer vector" % vec)
            if shape not in Partitions(len(vec)):
                raise ValueError("%s is not a partition of %d" %
                                 (shape, len(vec)))
            tableau = []
            pos = 0
            for rowlen in shape:
                tableau += [vec[pos:(pos + rowlen)]]
                pos += rowlen
                
        size = _hg_size(tableau)
        HG = HillmanGrasslTableaux(shape, size)
        return HG.element_class(HG, tableau)

def _hg_size(t):
    shape = Partition([len(r) for r in t])
    return sum(t[i][j] * shape.hook_length(i, j)
               for i, j in shape.cells())
        
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

        return (super(HillmanGrasslTableaux, cls)
                .__classcall__(cls, shape_part, size))

    def __init__(self, shape, size):
        self._size = size
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        return ("Hillman-Grassl tableaux of shape %s and size %d"
                % (self._shape, self._size))

    def __contains__(self, x):
        return (x in HillmanGrasslTableaux_all(self._shape) and
                _hg_size(Tableau(x)) == self._size)

    def _weighted_integer_vectors(self):
        hooks = sum(self._shape.hook_lengths(), [])
        return WeightedIntegerVectors(self._size, hooks)
    
    def __iter__(self):
        return (HillmanGrasslTableau(vec=v, shape=self._shape)
                for v in self._weighted_integer_vectors())

    def cardinality(self):
        return self._weighted_integer_vectors().cardinality()

    def random_element(self):
        v = self._weighted_integer_vectors().random_element()
        return HillmanGrasslTableau(vec=v, shape=self._shape)

class HillmanGrasslTableaux_all(DisjointUnionEnumeratedSets):
    def __init__(self, shape):
        from functools import partial
        self._shape = shape
        F = Family(NonNegativeIntegers(),
                   partial(HillmanGrasslTableaux, shape))
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
        return "Hillman-Grassl tableaux of shape %s" % self._shape
