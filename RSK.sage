from collections import Counter

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
         [[[1]],
          [[1, 1], [1]],
          [[2, 1], [1]],
          [[1, 1], [1, 1]],
          [[1, 1], [1], [1]],
          [[1, 1, 1], [1]]]]

#table = MatrixPartitionTable.from_mats(mats) + MatrixPartitionTable.from_parts(parts)
table = MatrixPartitionTable.from_parts(
    pp for pp in PlanePartitions((3,3,2))
    if pp.cells() and len(pp.cells()) <= 6)

### List all skew PP's of given shape and size
def is_skew_pp(pp):
    """ 
    Verify pp is a valid skew plane partition.
    
    Args:
        pp: A SkewTableau representing a plane partition.
    Returns:
        True if pp is weakly decreasing down rows and columns.
    """
    rows_and_cols = list(pp) + list(pp.conjugate())
    return all(a is None or a >= b
               for row in rows_and_cols
               for a, b in zip(row, row[1:]))

def skew_pps_with_shape(skew_part, size):
    """
    Find all skew pps with exactly the given shape and size.

    Args:
        skew_part: A SkewPartition representing the shape of the
            tableau.
        size: The sum of all the entries.
    Returns:
        A list of all possible skew plane partitions.
    """
    for comp in Compositions(size, length = skew_part.size()):
        i = 0
        rows = []
        for (outer, length) in zip(skew_part.outer(),
                                   skew_part.row_lengths()):
            rows.append([None] * (outer - length) + comp[i:i+length])
            i += length
        tab = SkewTableau(row for row in rows if row)
        if is_skew_pp(tab):
            yield tab

def all_pps(size, inner_shape = None, outer_shape = None):
    """
    Find all skew pps satisfying the given conditions.

    Args:
        size: The sum of all entries in the partition
        inner_shape: The shape at infinity.  If None, then
            there is no shape at infinity (and this gives a
            non-skew plane partition).
        outer_shape: A bound that the partition must fit in.
            It does not have to fit exactly.  If None, then
            there is no boundary.
    Returns:
        A list of all plane partitions satisfying the
        conditions.
    """
    if inner_shape is None:
        inner_shape = Partition([])
    return [pp
            for n in range(size + inner_shape.size() + 1)
            for shape in Partitions(n)
            if shape.contains(inner_shape)
            if outer_shape is None or outer_shape.contains(shape)
            for pp in skew_pps_with_shape(
                    SkewPartition([shape, inner_shape]), size)]

def pp_pairs(shape, size):
    """
    Find all plane partition pairs with the given shape and size.

    A plane partition pair is a pair of non-skew plane partitions
    where the first partition fits in shape, and the second has
    no restrictions, and the sum of the sizes is constant.

    Returns:
        A list of all pairs of plane partitions.
    """
    return [(pp1, pp2)
            for inner_size in range(size + 1)
            for pp1 in all_pps(inner_size, outer_shape = shape)
            for pp2 in all_pps(size - inner_size)]

def pp_both_lists(shape, size):
    return all_pps(size, inner_shape = shape), pp_pairs(shape, size)

def pp_both_lengths(shape, size):
    x, y = pp_both_lists(shape, size)
    return len(x), len(y)

def skew_pp_type(pp, max_entry):
    """
    Find the type of a given skew pp.

    The type is a tuple listing how many occurences of each number
    there are in pp.

    Args:
        pp: The skew plane partition to work with.
        max_entry: The length of the tuple to return.
    Returns:
        A tuple representing the type.
    """
    count = Counter(x for row in pp for x in row if x is not None)
    return tuple(count[i] for i in range(1, max_entry + 1))

def skew_pp_trace(pp):
    """
    Find the "trace" of a skew plane partition.

    This is a made up operation we thought would be presesrved.
    Turns out it's not.  It's essentially the sum of the two
    diagonals when the shape at infinity is a box.
    """
    startrow = [j for i, j in pp.cells() if i == 0]
    startcol = [i for i, j in pp.cells() if j == 0]
    trace = 0
    if startrow:
        s = min(startrow)
        trace += sum(pp[i][j] for i, j in pp.cells() if i + s == j)
    if startcol:
        s = min(startcol)
        trace += sum(pp[i][j] for i, j in pp.cells() if i == j + s)
    return trace

def examine_pp_problem(shape, size):
    """
    Find any useful information about this particular case.

    This method changes a lot based on what we think is
    "interesting".

    Args:
        shape: A partition describing the shape at infinity.
        size: The sum of all entries in the skew partition,
            and in the partition pair.
    Returns:
        Whether or not this case is "interseting" in some way
    """
    skews, pairs = pp_both_lists(shape, size)
    if len(skews) != len(pairs):
        raise AssertionError("Lengths do not match!")
    
    skew_size_count = Counter(pp.size() for pp in skews)
    inner_size_count = Counter(pp1.size() for pp1, pp2 in pairs)
    outer_size_count = Counter(pp2.size() for pp1, pp2 in pairs)
    pair_size_count = Counter(pp1.size() + pp2.size()
                              for pp1, pp2 in pairs)
    skew_type_count = Counter(skew_pp_type(pp, size) for pp in skews)
    inner_type_count = Counter(skew_pp_type(pp1, size)
                               for pp1, pp2 in pairs)
    outer_type_count = Counter(skew_pp_type(pp2, size)
                               for pp1, pp2 in pairs)
    combined_type_count = Counter(skew_pp_type(pp1, size)
                                  + skew_pp_type(pp2, size)
                                 for pp1, pp2 in pairs)
    summed_type_count = Counter(
        tuple(x + y for x, y in zip(skew_pp_type(pp1, size),
                                    skew_pp_type(pp2, size)))
        for pp1, pp2 in pairs)
    split_count = Counter(
        sum(x for row in pp1 for x in row)
        for pp1, pp2 in pairs)
    shape1_count = Counter(pp1.shape() for pp1, pp2 in pairs)
    shape2_count = Counter(pp2.shape() for pp1, pp2 in pairs)
    shape3_count = Counter((pp1.shape(), pp2.shape())
                           for pp1, pp2 in pairs)
    skew_shape_count = Counter(pp.shape() for pp in skews)
    skew_trace_count = Counter(skew_pp_trace(pp) for pp in skews)
    pair_trace_count = Counter(
        sum(pp1[i][j] for i, j in pp1.cells() if i == j) +
        sum(pp2[i][j] for i, j in pp2.cells() if i == j)
        for pp1, pp2 in pairs)
    
    print "----------------------------------"
    print "Skew sizes:", sorted(skew_size_count.values())
    print "Pair sizes:", sorted(pair_size_count.values())
    print "Skew types:", dict(skew_type_count)
    print "Summed types:", dict(summed_type_count)
    print "Splits counts:", sorted(split_count.values())
    print "Shape 1:", sorted(shape1_count.values())
    print "Shape 2:", sorted(shape2_count.values())
    print "Shape 3:", sorted(shape3_count.values())
    print "Skew Shape:", sorted(skew_shape_count.values())
    print "Skew Trace:", dict(skew_trace_count)
    print "Pair Trace:", dict(pair_trace_count)
    shapes_match = (sorted(shape3_count.values()) ==
                    sorted(skew_shape_count.values()))
    types_match = (skew_type_count == summed_type_count)
    # return true if this is somehow an excpetion
    #return (skew_size_count != pair_size_count or
    #        skew_type_count != summed_type_count)
    return skew_trace_count != pair_trace_count

"""
non_matches = [(i, j, k)
               for i in range(1,4)
               for j in range(1,4)
               for k in range(1,6)
               if examine_pp_problem(Partition([i] * j), k)]
"""

### Find all the cases we don't know how to handle
def unknown_examples(width, height, size):
    """
    Find all skew partitions and partition pairs which
    we don't know how to handle.

    Find all skew partitions which cannot be split up into
    a pair of partitions with one fitting inside a box in the
    most obvious way.  Also find all pairs of plane partitions
    in which the second partition exceeds both the width and
    height of the box.  These two lists should have the same
    length.

    Args:
        width, height: Dimensions of the shape at infinity.
        size: The sum of all entries in either case.
    Returns:
        A pair of lists, the first being the skew partitions
        and the second being pairs of partitions.
    """
    skew, pairs = pp_both_lists(Partition([width] * height), size)
    bad_pairs = [(p, q) for p, q in pairs
                 if (len(q.outer_shape()) > height
                     and q.outer_shape()[0] > width)]
    bad_skews = [s for s in skew
                 if any(i >= height and j >= width
                        for i, j in s.cells())
                 or (max(i + 1 for i, j in s.cells()) > height * 2 and
                     max(j + 1 for i, j in s.cells()) > width * 2)]
    assert len(bad_pairs) == len(bad_skews), \
        "There are {} leftover skews, but {} leftover pairs!".format(
            len(bad_skews), len(bad_pairs))
    return bad_skews, bad_pairs

def print_unknown_examples(width, height, size):
    """
    Print out the unknown examples in a nice way.

    Args:
        width, height: The shape at infinity.
        size: The sum of all entries.
    """
    skews, pairs = unknown_examples(width, height, size)
    print "-" * 20
    print "Bad Skew Partitions ({}):".format(len(skews))
    for s in skews:
        s.pp()
        print
    print "Bad Pairs of Parititons ({}):".format(len(pairs))
    for p, q in pairs:
        print "First:"
        if p:
            p.pp()
        else:
            print "  (null)"
        print "Second:"
        q.pp()
        print
