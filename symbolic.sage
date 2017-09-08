''' Implementation of PT/DT invariants using symbolic rings

SymbolicRings can keep track of exponentials, so they capture
information about all insertions using finite information.
'''

R1.<b,c,z1,z2> = QQ[]
R2.<q> = SR[[]]
a = -b-c

def column_sum(value, i, j, z, invert):
    if value is None:
        # Infinite column of DT
        T = b * i + c * j
        return (1 - e^(b*z)) * (1 - e^(c*z)) * e^(T * z)
    elif not invert:
        # Finite part of DT
        T = lambda k: a * k + b * i + c * j
        P = prod(1 - e^(t * z) for t in [a,b,c])
        return P * sum(e^(T(k) * z) for k in range(value))
    else:
        # Infinite part of PT
        T = b * i + c * j - (value) * a
        return (1 - e^(b*z)) * (1 - e^(c*z)) * e^(T * z)

def chern_character(part, z, invert=False):
    return 1 - sum(column_sum(val, i, j, z, invert)
                   for i, row in enumerate(part)
                   for j, val in enumerate(row))

def chern_product(part, zs, invert=False):
    return expand(prod(chern_character(part, z, invert) for z in zs))

def all_pps(P, size):
    return tuple(P.graded_component(size))

def new_weighted_sum(P, zs, invert=False, prec=8):
    return sum(chern_product(pi, zs, invert) * q^size
               for size in range(prec)
               for pi in all_pps(P, size)) + O(q^prec)

def PT(shape, zs=(), prec=8):
    return new_weighted_sum(ReversePlanePartitions(shape),
                            zs, True, prec)

def DT(shape, zs=(), prec=8):
    return new_weighted_sum(SkewPlanePartitions(shape),
                            zs, False, prec)

def DTn(shape, zs=(), prec=8):
    """ Compute the normalized DT series, or DT' """
    return DT(shape, zs, prec) / DT([], (), prec)
