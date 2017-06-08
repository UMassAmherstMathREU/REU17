##########################
# Reverse Semi-standard tableaux #
##########################
class ReverseSemistandardTableaux(Tableaux):
    """
    A factory class for the various classes of reverse semistandard tableaux.

    INPUT:

    Keyword arguments:
    TBD
    - ``size`` -- The size of the tableaux
    - ``shape`` -- The shape of the tableaux
    - ``eval`` -- The weight (also called content or evaluation) of
      the tableaux
    - ``max_entry`` -- A maximum entry for the tableaux.  This can be a
      positive integer or infinity (``oo``). If ``size`` or ``shape`` are
      specified, ``max_entry`` defaults to be ``size`` or the size of
      ``shape``.

    Positional arguments:

    - The first argument is interpreted as either ``size`` or ``shape``
      according to whether it is an integer or a partition
    - The second keyword argument will always be interpreted as ``eval``

    OUTPUT:

    - The appropriate class, after checking basic consistency tests. (For
      example, specifying ``eval`` implies a value for `max_entry`).

    A reverse semistandard tableau is a tableau whose entries are positive integers,
    which are weakly decreasing in rows and strictly decreasing down columns.
    Note that Sage uses the English convention for partitions and tableaux;
    the longer rows are displayed on top.

    Classes of semistandard tableaux can be iterated over if and only if there
    is some restriction.

    EXAMPLES::
        TBD
        sage: SST = SemistandardTableaux([2,1]); SST
        Semistandard tableaux of shape [2, 1] and maximum entry 3
        sage: SST.list()
        [[[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]]]

        sage: SST = SemistandardTableaux(3); SST
        Semistandard tableaux of size 3 and maximum entry 3
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 1, 3]],
         [[1, 2, 2]],
         [[1, 2, 3]],
         [[1, 3, 3]],
         [[2, 2, 2]],
         [[2, 2, 3]],
         [[2, 3, 3]],
         [[3, 3, 3]],
         [[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

        sage: SST = SemistandardTableaux(3, max_entry=2); SST
        Semistandard tableaux of size 3 and maximum entry 2
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 2, 2]],
         [[2, 2, 2]],
         [[1, 1], [2]],
         [[1, 2], [2]]]

        sage: SST = SemistandardTableaux(3, max_entry=oo); SST
        Semistandard tableaux of size 3
        sage: SST[123]
        [[3, 4], [6]]

        sage: SemistandardTableaux(max_entry=2)[11]
        [[1, 1], [2]]

        sage: SemistandardTableaux()[0]
        []

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`SemistandardTableaux`
        for more information.

        TESTS::

            sage: SemistandardTableaux()
            Semistandard tableaux
            sage: SemistandardTableaux(3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(size=3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(0)
            Semistandard tableaux of size 0 and maximum entry 0
            sage: SemistandardTableaux([2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux([])
            Semistandard tableaux of shape [] and maximum entry 0
            sage: SemistandardTableaux(eval=[2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
            sage: SemistandardTableaux(3, [2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(3, shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(3, [2,1], shape=[2,1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, max_entry=4)
            Semistandard tableaux of size 3 and maximum entry 4
            sage: SemistandardTableaux(3, max_entry=oo)
            Semistandard tableaux of size 3
            sage: SemistandardTableaux([2, 1], max_entry=oo)
            Semistandard tableaux of shape [2, 1]
            sage: SemistandardTableaux([2, 1], [2, 1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: mu = Partition([2,1]); SemistandardTableaux(mu, mu)
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, [2, 1], max_entry=2)
            Semistandard tableaux of size 3 and weight [2, 1]

            sage: SemistandardTableaux(3, shape=[2])
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes

            sage: SemistandardTableaux(3, [2])
            Traceback (most recent call last):
            ...
            ValueError: size and eval are different sizes

            sage: SemistandardTableaux([2],[3])
            Traceback (most recent call last):
            ...
            ValueError: shape and eval are different sizes

            sage: SemistandardTableaux(2,[2], max_entry=4)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux(eval=[2], max_entry=oo)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape must be a (skew) partition
        """
        from sage.combinat.partition import Partition, _Partitions
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        mu = kwargs.get('eval', None)
        mu = kwargs.get("mu", mu)

        max_entry = kwargs.get('max_entry', None)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError( "size was specified more than once" )
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError( "the shape was specified more than once" )
                shape = args[0] # we check it's a partition later

        if len(args) == 2:
            # The second non-keyword argument is the weight
            if mu is not None:
                raise ValueError( "the weight was specified more than once" )
            else:
                mu = args[1]

        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError( "size must be an integer" )
            elif size < 0:
                raise ValueError( "size must be non-negative" )

        if shape is not None:
            from sage.combinat.skew_partition import SkewPartitions
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape in _Partitions:
                shape = Partition(shape)
            elif shape in SkewPartitions():
                from sage.combinat.skew_tableau import SemistandardSkewTableaux
                return SemistandardSkewTableaux(shape, mu)
            else:
                raise ValueError( "shape must be a (skew) partition" )

        if mu is not None:
            if (not mu in Compositions()) and\
                    (not mu in _Partitions):
                raise ValueError( "mu must be a composition" )
            mu = Composition(mu)

        is_inf = max_entry is PlusInfinity()

        if max_entry is not None:
            if not is_inf and not isinstance(max_entry, (int, Integer)):
                raise ValueError( "max_entry must be an integer or PlusInfinity" )
            elif max_entry <= 0:
                raise ValueError( "max_entry must be positive" )

        if (mu is not None) and (max_entry is not None):
            if max_entry != len(mu):
                raise ValueError( "the maximum entry must match the weight" )

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError( "size and shape are different sizes" )

        if (size is not None) and (mu is not None):
            if sum(mu) != size:
                # This could return an empty class instead of an error
                raise ValueError( "size and eval are different sizes" )

        # Dispatch appropriately
        if (shape is not None) and (mu is not None):
            if sum(shape) != sum(mu):
                # This could return an empty class instead of an error
                raise ValueError( "shape and eval are different sizes" )
            else:
                return SemistandardTableaux_shape_weight(shape, mu)

        if (shape is not None):
            if is_inf:
                return SemistandardTableaux_shape_inf(shape)
            return SemistandardTableaux_shape(shape, max_entry)

        if (mu is not None):
            return SemistandardTableaux_size_weight(sum(mu), mu)

        if (size is not None):
            if is_inf:
                return SemistandardTableaux_size_inf(size)
            return SemistandardTableaux_size(size, max_entry)

        return SemistandardTableaux_all(max_entry)

    Element = SemistandardTableau

    def __init__(self, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardTableaux()
            sage: TestSuite(S).run()
        """
        if 'max_entry' in kwds:
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        Tableaux.__init__(self, **kwds)

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__`` for enumerated sets
        does not allow slices so we override it.

        EXAMPLES::

            sage: StandardTableaux([4,3,3,2])[10:20]     # indirect doctest
            [[[1, 3, 9, 12], [2, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 3, 9, 12], [2, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 5, 8, 12], [2, 6, 10], [3, 7, 11], [4, 9]],
             [[1, 4, 8, 12], [2, 6, 10], [3, 7, 11], [5, 9]],
             [[1, 3, 8, 12], [2, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 2, 8, 12], [3, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 4, 8, 12], [2, 5, 10], [3, 7, 11], [6, 9]],
             [[1, 3, 8, 12], [2, 5, 10], [4, 7, 11], [6, 9]]]

            sage: SemistandardTableaux(size=2, max_entry=oo)[5]
            [[2, 3]]

            sage: SemistandardTableaux([2,1], max_entry=oo)[3]
            [[1, 2], [3]]

            sage: SemistandardTableaux(3, max_entry=2)[0:5]    # indirect doctest
            [[[1, 1, 1]],
            [[1, 1, 2]],
            [[1, 2, 2]],
            [[2, 2, 2]],
            [[1, 1], [2]]]

            sage: SemistandardTableaux([2,2], [2, 1, 1])[0]    # indirect doctest
            [[1, 1], [2, 3]]

            sage: SemistandardTableaux([1,1,1], max_entry=4)[0:4]
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: SemistandardTableaux(3, [2,1])[1]    # indirect doctest
            [[1, 1], [2]]

            sage: StandardTableaux(3)[:]  # indirect doctest
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

            sage: StandardTableaux([2,2])[1]   # indirect doctest
            [[1, 2], [3, 4]]

        TESTS::

            sage: SemistandardTableaux()[5]
            [[1], [2]]

            sage: SemistandardTableaux(max_entry=2)[5]
            [[2, 2]]

            sage: SemistandardTableaux()[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: SemistandardTableaux(size=2, max_entry=oo)[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None and not self.is_finite():
                raise ValueError( 'infinite set' )
        else:
            raise ValueError( 'r must be an integer or a slice' )
        count=0
        tabs=[]
        for t in self:
            if count==stop:
                break
            if count>=start:
                tabs.append(t)
            count+=1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count==stop or stop is None:
            return tabs
        raise IndexError('value out of range')

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`SemistandardTableau`.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardTableau([[1]]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True
            sage: [[1,3,2]] in T
            False

        Check that :trac:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.SemistandardTableaux()
            False
        """
        if isinstance(t, SemistandardTableau):
            return self.max_entry is None or \
                    len(t) == 0 or \
                    max(max(row) for row in t) <= self.max_entry
        elif not t:
            return True
        elif Tableaux.__contains__(self, t):
            for row in t:
                if not all(c > 0 for c in row):
                    return False
                if not all(row[i] >= row[i+1] for i in range(len(row)-1)):
                    return False
            for row, next in zip(t, t[1:]):
                if not all(row[c] > next[c] for c in range(len(next))):
                    return False
            return self.max_entry is None or max(max(row) for row in t) <= self.max_entry
        else:
            return False
