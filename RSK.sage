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
def rsk(mat):
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
    omega = sorted(list(mat.dict().items()))
    for (i, j), a in omega:
        for _ in range(a):
            # Note that (i, j) start counting from 0 in python,
            # but most texts start counting from 1
            P, row = row_insert_rownum(P, j + 1)
            col = Q.shape()[row] if row < len(Q.shape()) else 0
            Q = Q.add_entry((row, col), i + 1)
    return P, Q

def inverse_rsk(P, Q):
    """
    Apply the inverse RSK algorithm

    Args:
        P: The insertion tableau
        Q: The recording tableau

    Returns:
        An list of tuples representing a two-line array

    Raises:
        ValueError: If P and Q are not the same shape
    """
    if P.shape() != Q.shape():
        raise ValueError("P and Q are not the same shape")

    omega = []
    while P:
        # find the largest, rightmost entry
        i, _, r = max((Q[r][-1], len(Q[r]), r) for r in range(len(Q)))
        Q = ssyt_delete_entry(Q, r)
        P, j = inverse_row_insert(P, r)
        omega.append((i, j))
    # We built omega up backwards, flip it
    omega.reverse()
    return omega

### Test the RSK algorthim using the example in the text
def test_rsk():
    A = matrix(ZZ, 3, 3, [1,0,2,0,2,0,1,1,0])
    P = SemistandardTableau([[1,1,2,2],[2,3],[3]])
    Q = SemistandardTableau([[1,1,1,3],[2,2],[3]])
    omega = [(1,1),(1,3),(1,3),(2,2),(2,2),(3,1),(3,2)]
    assert rsk(A) == (P, Q)
    assert inverse_rsk(P, Q) == omega
test_rsk()
