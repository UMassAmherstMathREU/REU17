from sage.combinat.tableau import Tableau
from hillman_grassl.py import Hillman_Grassl

class ReversePlanePartition(Tableau):

    @staticmethod
    def __classcall_private__():
        pass

    def __init__(self):
        pass

    def to_HilmanGrasslTableau(self):
        return Hillman_Grassl(self)

class ReversePlanePartitions(Tableaux):

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

    Element = ReversePlanePartition
