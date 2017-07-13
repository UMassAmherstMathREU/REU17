from ptdt_package import *
from ptdt_package.weights import partition_weight

R.<a,b> = QQ[]
c = - a - b

shape = Partition([3, 2, 1])
conj = shape.conjugate()

def T(i):
    return i * (i + 1) / 2

def hook_weight(i, j):
    return ((T(shape[i] - 1) - T(j) +conj[j] * j) * b +
            (T(conj[j] - 1) - T(i) + shape[i] * i) * a)

def hg_weight(hg):
    return sum(hg[i][j] * hook_weight(i, j)
               for i, j in hg.cells())

for hg in HillmanGrasslTableaux(shape, 10):
    rpp = hg.to_ReversePlanePartition()
    pw = partition_weight((a, b, c), 1, rpp, True)
    hw = hg_weight(hg)
    assert c.divides(pw - hw)
