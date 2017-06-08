### Define SSYT to plane partition algorithm
def ssyt_to_pp(P, Q):
    newCols = [merge_cols(pcol, qcol) for pcol, qcol
               in zip(P.conjugate(), Q.conjugate())]
    pi = PlanePartition(newCols).transpose()
    conjRows = [Partition(r).conjugate() for r in pi]
    return PlanePartition(conjRows)

def merge_cols(pcol, qcol):
    if len(pcol) != len(qcol):
        raise ValueError("The two columns have different sizes")
    diagLen = len(pcol)
    pslant = [pcol[i] + i for i in range(diagLen)]
    qslant = [qcol[i] + i for i in range(diagLen)]
    qconj = Partition(qslant).conjugate()
    return [pslant[i] if i < diagLen else qconj[i]
            for i in range(len(qconj))]

### Define reverse algorithm
def recover_upper(part):
    rank = part.frobenius_rank()
    return [part[i] - i for i in range(rank)]

def recover_lower(part):
    return recover_upper(part.conjugate())

def pp_to_ssyt(part):
    conjRows = [Partition(r).conjugate() for r in part]
    pi = PlanePartition(conjRows)
    pcols = [recover_upper(Partition(c)) for c in pi.transpose()]
    qcols = [recover_lower(Partition(c)) for c in pi.transpose()]
    P = Tableau(pcols).conjugate()
    Q = Tableau(qcols).conjugate()
    return P, Q
