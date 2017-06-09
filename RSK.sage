from subprocess import call

### Implement row insertion
def row_insert_rownum(ssyt, k):
    """
    Insert an element into an SSYT using row insertion.

    The SSYT is not modified, instead, the resulting SSYT is
    returned.  The row number where the insertion happened is
    also returned, for convenience.
    
    Args:
        ssyt: The SemistandardTableau to insert into
        k: The element to insert

    Returns:
        A pair of the resulting SSYT, and the index of the row
        inserted (starting from 0).
    """
    row = 0
    while row < len(ssyt.shape()):
        # find the first col with entry >k
        col = 0
        while col < ssyt.shape()[row] and ssyt[row][col] <= k:
            col += 1

        # If we're at the end of the row, we're done
        # otherwise, bump the current entry and move to the next row
        if col == ssyt.shape()[row]:
            return ssyt.add_entry((row, col), k), row
        else:
            ssyt, k = ssyt.add_entry((row, col), k), ssyt[row][col]
            row += 1

    # We got to the last row, so just add a new row
    return ssyt.add_entry((row, 0), k), row

def row_insert(ssyt, k):
    """
    Insert an element into an SSYT using row insertion.

    The SSYT is not modified, instead, the resulting SSYT is
    returned.
    
    Args:
        ssyt: The SemistandardTableau to insert into
        k: The element to insert

    Returns:
        The resulting SSYT
    """
    return row_insert_rownum(ssyt, k)[0]

def ssyt_delete_entry(ssyt, row):
    """
    Remove the last entry in a given row from the SSYT.

    Args:
        ssyt: The SemistandardTableau to remove from
        row: The index of the row to remove from

    Returns:
        The resulting SSYT
    """
    l = [list(r) for r in ssyt]
    del l[row][-1]
    if not l[row]:
        del l[row]
    return SemistandardTableau(l)

def inverse_row_insert(ssyt, row):
    """
    Reverse the row insertion process.

    Reverse a row insertion, starting with the last entry on
    the given row.

    Args:
        ssyt: The SemistandardTableau to start with
        row: The index (starting from 0) of the row which was
            added to.

    Returns:
        A pair (newSSYT, k), with the resulting SSYT and the
        element which was bumped out.

    Raises:
        ValueError: if removing from row would leave an
            invalid tableau
    """
    shape = ssyt.shape()
    if row < len(ssyt) - 1 and shape[row] == shape[row + 1]:
        raise ValueError("removing from row would leave an "
                         "invalid tableau")
    # pull out last entry of row
    k = ssyt[row][-1]
    ssyt = ssyt_delete_entry(ssyt, row)
    # convert to mutable so we can mess with it
    ssytList = [list(r) for r in ssyt]
    # move up the table and reverse the bumping
    while row > 0:
        row -= 1
        col = max(i for i in range(ssyt.shape()[row])
                  if ssyt[row][i] < k)
        k, ssytList[row][col] = ssytList[row][col], k
    return SemistandardTableau(ssytList), k

### Check row insertion matches with builtin method
def test_row_insert():
    # edge case
    ssyt = SemistandardTableau([])
    assert row_insert(ssyt, 1) == SemistandardTableau([[1]])
    # random input
    for i in range(20):
        ssyt = SemistandardTableaux(10).random_element()
        for k in range(1, 11):
            newssyt, row = row_insert_rownum(ssyt, k)
            assert newssyt == ssyt.bump(k)
            assert ssyt, k == reverse_row_insert(ssyt, row)
            assert ssyt, k == ssyt.reverse_bump(row)
test_row_insert()

### Define RSK algorithm
def rsk(mat, reverse=False):
    """Apply the RSK algorithm to a matrix.

    Given a matrix of non-negative integers, produce two SSYTs of the
    same shape.

    Args:
        mat: A matrix of non-negative integers

    Returns:
        A pair (P, Q) of SSYTs of the same shape
    """
    # Initialize P and Q to be empty SSYTs
    P = SemistandardTableau([])
    Q = SemistandardTableau([])

    # The "two-line array" can be represented as a dict indexed by
    # pairs, where the (i, j) entries indicates the number of times
    # (i, j) appears.  Turning this into a list and sorting puts the
    # array in the order we want, where we first sort by i and then
    # sort by j
    omega = [((i + 1, j + 1), mat[i,j]) for i, j in mat.dict()]
    if reverse and omega:
        M = max(max(a, b) + 1 for (a, b), c in omega)
        omega = [((M - a, M - b), c) for (a, b), c in omega]
    omega.sort()
    for (i, j), a in omega:
        for _ in range(a):
            # Note that (i, j) start counting from 0 in python,
            # but most texts start counting from 1
            P, row = row_insert_rownum(P, j)
            col = Q.shape()[row] if row < len(Q.shape()) else 0
            Q = Q.add_entry((row, col), i)
    if reverse:
        P = Tableau([M - x for x in r] for r in P)
        Q = Tableau([M - x for x in r] for r in Q)
    return P, Q

def inverse_rsk(P, Q, nrow=None, ncol=None, reverse=False):
    """
    Apply the inverse RSK algorithm

    Args:
        P: The insertion tableau
        Q: The recording tableau
        nrow: The number of rows to make in the new matrix,
            or None to determine automatically
        ncol: The number of columns to make in the new matrix,
            or None to determine automatically
        reverse: Whether to use SSYT's or reverse SSYT's

    Returns:
        A matrix with non-negative integer entries 

    Raises:
        ValueError: If P and Q are not the same shape
    """
    if P.shape() != Q.shape():
        raise ValueError("P and Q are not the same shape")

    if reverse:
        M = max(x + 1 for x in P.entries() + Q.entries())
        P = SemistandardTableau([[M - x for x in r] for r in P])
        Q = SemistandardTableau([[M - x for x in r] for r in Q])

    omega = {}
    while P:
        # find the largest, rightmost entry
        i, _, r = max((Q[r][-1], len(Q[r]), r) for r in range(len(Q)))
        Q = ssyt_delete_entry(Q, r)
        P, j = inverse_row_insert(P, r)
        if (i, j) in omega:
            omega[i, j] += 1
        else:
            omega[i, j] = 1

    if reverse:
        omega = {(M - i - 1, M - j - 1): k for (i, j), k in omega.items()}
    else:
        omega = {(i - 1, j - 1): k for (i, j), k in omega.items()}

    if nrow is None:
        nrow = max(i + 1 for i, j in omega)
    if ncol is None:
        ncol = max(j + 1 for i, j in omega)

    return matrix(ZZ, nrow, ncol, omega)

### Test the RSK algorthim using the example in the text
def test_rsk():
    A = matrix(ZZ, 3, 3, [1,0,2,0,2,0,1,1,0])
    P = SemistandardTableau([[1,1,2,2],[2,3],[3]])
    Q = SemistandardTableau([[1,1,1,3],[2,2],[3]])
    assert rsk(A) == (P, Q)
    assert inverse_rsk(P, Q, 3, 3) == A
test_rsk()

### Define SSYT to plane partition algorithm
def ssyt_to_pp(P, Q):
    """
    Convert two reverse SSYTs to a plane partition

    Args:
        P: The insertion tableau
        Q: The recording tableau

    Returns:
        A plane partition resulting from merging P and Q
    """
    newCols = [merge_cols(pcol, qcol) for pcol, qcol
               in zip(P.conjugate(), Q.conjugate())]
    pi = PlanePartition(newCols).transpose()
    conjRows = [Partition(r).conjugate() for r in pi]
    return PlanePartition(conjRows)

def merge_cols(pcol, qcol):
    """
    Merge two columns by flipping and identifying diagonals.

    Args:
        pcol, qcol: Lists of integers representing partitions

    Returns:
        The result of the merge, as an integer list
    """
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

### Go the whole way from mat -> pp
def mat_to_pp(mat):
    P, Q = rsk(mat, reverse=True)
    return ssyt_to_pp(P, Q)

def pp_to_mat(pp, nrow=None, ncol=None):
    P, Q = pp_to_ssyt(pp)
    return inverse_rsk(P, Q, nrow, ncol, reverse=True)

### Write output into a table

def mat_to_table_row(mat):
    P, Q = rsk(mat, reverse=True)
    pp = ssyt_to_pp(P, Q)
    return mat, P, Q, pp

def pp_to_table_row(pp):
    P, Q = pp_to_ssyt(pp)
    mat = inverse_rsk(P, Q, reverse=True)
    return mat, P, Q, pp

class MatrixPartitionTable(object):
    @staticmethod
    def from_mats(mats):
        table = MatrixPartitionTable()
        table.entries = [mat_to_table_row(mat) for mat in mats]
        return table

    @staticmethod
    def from_parts(parts):
        table = MatrixPartitionTable()
        table.entries = [pp_to_table_row(pp) for pp in parts]
        return table

    def _latex_(self):
        header = """
        \\begin{tabular}{c | c | c | c | c}
        Matrix & $P$ & $Q$ & Partition (table) & Partition (diagram) \\\\
        \\hline
        """
        footer = "\\end{tabular}"
        rows = ["${}$ & {} & {} & {} & {}".format( latex(mat), latex(P),
                                                   latex(Q),
                                                   latex(pp.to_tableau()),
                                                   latex(pp))
                for mat, P, Q, pp in self.entries]
        body = "\\\\ \\hline".join(rows)
        return header + body + footer

    def __repr__(self):
        return "Matrix/Partition Table " + repr(self.entries)

    def __add__(self, other):
        table = MatrixPartitionTable()
        table.entries = self.entries + other.entries
        return table

    
mats = [matrix(ZZ, 2, 2, mat) for mat in [[1,0,0,0],
                                          [0,1,0,0],
                                          [0,0,1,0],
                                          [0,0,0,1],
                                          [1,0,0,1],
                                          [0,1,1,0],
                                          [1,1,1,1]]]

parts = [PlanePartition(part) for part in
         [[[2, 1], [1]],
          [[1, 1]],
          [[1]]]]

table = MatrixPartitionTable.from_mats(mats) + MatrixPartitionTable.from_parts(parts)
