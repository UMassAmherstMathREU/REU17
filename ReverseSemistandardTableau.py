class ReverseSemistandardTableau(Tableau):
    """
    A class to model a reversed semistandard tableau.

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    OUTPUT:

    - A ReverseSemistandardTableau object constructed from ``t``.

    A reverse semistandard tableau is a tableau whose entries are positive integers,
    which are decreasing increasing in rows and strictly decreasing down columns.

    EXAMPLES::
    TBD

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a SemistandardTableau from the appropriate
    :class:`Parent` object::

        sage: SST = SemistandardTableaux()
        sage: SST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::
        TBD
        sage: SemistandardTableau([[1,2,3],[1]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 3], [1]] is not a column strict tableau

        sage: SemistandardTableau([[1,2,1]])
        Traceback (most recent call last):
        ...
        ValueError: The rows of [[1, 2, 1]] are not weakly increasing

        sage: SemistandardTableau([[0,1]])
        Traceback (most recent call last):
        ...
        ValueError: entries must be positive integers
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::
            TBD
            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, SemistandardTableau):
            return t
        elif t in SemistandardTableaux():
            return SemistandardTableaux_all().element_class(SemistandardTableaux_all(), t)

        # t is not a reverse semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c,(int,Integer)) and c>0 for row in t for c in row):
            raise ValueError("entries must be positive integers"%t)

        if any(row[c]<row[c+1] for row in t for c in range(len(row)-1)):
            raise ValueError("The rows of %s are not weakly decreasing"%t)

        # If we're still here ``t`` cannot be column strict
        raise ValueError('%s is not a column strict tableau' % t)


    def __init__(self, parent, t):
        r"""
        Initialize a reverse semistandard tableau.

        TESTS::

            TBD
            sage: t = Tableaux()([[1,1],[2]])
            sage: s = SemistandardTableaux(3)([[1,1],[2]])
            sage: s==t
            True
            sage: s.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: r = SemistandardTableaux(3)(t); r.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: isinstance(r, Tableau)
            True
            sage: s2 = SemistandardTableaux(3)([(1,1),(2,)])
            sage: s2 == s
            True
            sage: s2.parent()
            Semistandard tableaux of size 3 and maximum entry 3
        """
        super(SemistandardTableau, self).__init__(parent, t)

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()

        for row in t:
            if any(c not in PI for c in row):
                raise ValueError("the entries of a semistandard tableau must be non-negative integers")
            if any(row[c] < row[c+1] for c in range(len(row)-1)):
                raise ValueError("the entries in each row of a semistandard tableau must be weakly decreasing")

        # and strictly increasing down columns
        if t:
            for row, next in zip(t, t[1:]):
                if not all(row[c] > next[c] for c in range(len(next))):
                    raise ValueError("the entries of each column of a semistandard tableau must be strictly decreasing")
