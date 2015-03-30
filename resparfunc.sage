"""
resparfunc

A restricted partition function counts the number of possible partitions
of a nonnegative integer, subject to various constraints.
This module is about restricted partition functions of a special kind,
also known as Sylvester's denumerant.
"""
#***********************************************************************
#       Copyright (C) 2015 Anne Huschitt <Anne.Huschitt@gmail.com>,
#                          Alfred Seiler
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

class RestrictedPartitionFunction(object):
    r"""
    A RestrictedPartitionFunction allows to get the number
    of partitions with parts in a given list of positive integers
    for arbitrary nonnegative integers very efficiently.

    Once constructed with ``RestrictedPartitionFunction(parts_in=A)``
    for a specific restriction list ``A``, a
    ``RestrictedPartitionFunction`` object can be used to calculate the
    number of partitions for an integer ``t`` using the object's method
    ``number_of_partitions(t)``.
    Note that this gives the same result as
    ``Partitions(t, parts_in=A).cardinality()``, but the once
    constructed object can be reused for different values of ``t`` and
    the values of ``t`` can be considerably larger.

    Constructing a ``RestrictedPartitionFunction`` object for a specific
    restriction list ``A`` can take a long time depending on the numbers
    in  ``A``. The needed time in seconds is printed when the
    calculation is finished, but this can be suppressed by setting the
    class variable ``RestrictedPartitionFunction.print_time=False``.
    The construction of the object is an iterative process with a number
    of levels equal to the size of the restriction list. To watch the
    progress of construction, one can set the class variable
    ``RestrictedPartitionFunction.print_info=True``. Before starting
    the iterative process, the numbers of the restriction list are
    first rearranged with the purpose to get a shorter construction
    time. For testing purposes, it is possible to suppress this
    rearrangement by setting the class variable
    ``RestrictedPartitionFunction.preserve_order=True``.

    The representation of a ``RestrictedPartitionFunction`` is based on
    the fact that the restricted partition function `d(t)` for
    `A=(a_{1},a_{2},...,a_{n})` can be written as
    `d(t)=\sum_{m=0}^{n-1}\sum_{i=1}^{u_{m}}f_{m,i}(t)t^{m}` with
    appropriate natural numbers `u_{m}` and periodical functions
    `f_{m,i}`, where `m=0,...,n-1` and `i=1,...,u_{m}`. The periods of
    the functions `f_{m,i}` are correlated with common factors of the
    numbers `a_{i}`.
    We represent a periodical function `f_{m,i}` by an object of type
    ``PeriodObject``  containing its period `P_{m,i}` and a list of
    function values `f_{m,i}(0),...,f_{m,i}(P_{m,i}-1)`.
    Remark that such a representation of a ``RestrictedPartitionFunction``
    is not unique.

    To explore the content of a ``RestrictedPartitionFunction`` object
    ``d``, one can use the object method ``get_content_str()``. With the
    argument ``output_format`` it is possible to select the extent and
    presentation of the returned content. As an example, with
    ``output_format=summary`` one gets the number of period objects
    and the respective periods for each power. The full content is
    returned if ``output_format=full`` is given. Be aware that the
    amount of information can be very large in this case. Another
    string representation of the full content which may be better
    suited for a subsequent evaluation can be obtained with
    ``output_format=list``.
    The full content in ``full`` or ``list`` format can also be written
    to a file with the object method ``to_file``. In contrast to the
    builtin ``save``, the content is saved in printable form, so it can
    be easily explored with a suitable text editor.
    With the class method ``from_file``, an object with the same
    content can be reconstructed from a file written with ``to_file``.

    If a ``RestrictedPartitionFunction`` object for a given restriction
    list has been already constructed, one can obtain a new
    ``RestrictedPartitionFunction`` for a restriction list with an
    additional positive integer ``a`` with the class method
    ``build_next()``.


    EXAMPLES:

    Construct some partition functions for specific restriction lists
    and use them to get the number of partitions for several numbers::

        sage: RestrictedPartitionFunction.print_time=False  # suppress timing for doctest
        sage: d=RestrictedPartitionFunction([1,2,3])
        sage: d.number_of_partitions(3)
        3
        sage: d.number_of_partitions(10000)
        8338334
        sage: d.number_of_partitions(10**20)
        833333333333333333383333333333333333334
        sage: d.number_of_partitions(10**10**6).ndigits()
        1999999
        sage: A=[5,10,10,2,8,20,15,2,9,9,7,4,12,13,19]
        sage: d=RestrictedPartitionFunction(parts_in=A)
        sage: d.number_of_partitions(9)
        8
        sage: d.number_of_partitions(10000)
        39327093201690271931461735037880

    Comparing ``RestrictedPartitionFunction().number_of_partitions()`` with
    ``Partitions().cardinality()``. The result must be the same::

        sage: RestrictedPartitionFunction.print_time=False  # suppress timing for doctest
        sage: A=[10,12,13]
        sage: d=RestrictedPartitionFunction(parts_in=A)
        sage: t=500
        sage: d.number_of_partitions(t) == Partitions(t, parts_in=A).cardinality()
        True
        sage: len([t for t in [100..200] if d.number_of_partitions(t) != Partitions(t, parts_in=A).cardinality()])
        0

    Showing full content of period objects contained in the representation
    of an object of type ``RestrictedPartitionFunction``::

        sage: A = [10, 5, 6]
        sage: d = RestrictedPartitionFunction(A)
        sage: print d.get_content_str(output_format='full')
        RestrictedPartitionFunction for [5, 6, 10], consisting of 3 lists with a total of 5 period objects:
        power=0, 2 period object(s):
        P=10, values=[2/3, 57/200, -7/75, -7/200, -28/75, 5/8, 17/75, 13/200, -22/75, -3/200]
        P=6, values=[1/3, -1/3, 0, 0, 1/3, 0]
        power=1, 2 period object(s):
        P=5, values=[29/300, 23/300, 17/300, 11/300, 1/60]
        P=2, values=[-1/75, -3/100]
        power=2, 1 period object(s):
        P=1, values=[1/600]

    Showing summary of period objects contained in the
    ``RestrictedPartitionFunction`` object constructed above::

        sage: print d.get_content_str(output_format='summary')
        RestrictedPartitionFunction for [5, 6, 10], consisting of 3 lists with a total of 5 period objects:
        power=0, 2 PO(s) with period: 10,6
        power=1, 2 PO(s) with period: 5,2
        power=2, 1 PO(s) with period: 1

    Showing summary of period objects contained in the representation of an
    object of type ``RestrictedPartitionFunction`` may be suited especially
    to get a quick overview if the full representation of the period
    objects is too large for a screen output::

        sage: A=[5, 10, 10, 2, 8, 20, 15, 2, 9, 9, 7, 4, 12, 13, 19]
        sage: d=RestrictedPartitionFunction(parts_in=A)
        sage: print d.get_content_str(output_format='summary')
        RestrictedPartitionFunction for [2, 5, 7, 9, 13, 19, 4, 2, 12, 8, 15, 10, 9, 20, 10], consisting of 15 lists with a total of 29 period objects:
        power=0, 8 PO(s) with period: 12,15,7,9,13,19,8,20
        power=1, 3 PO(s) with period: 10,4,9
        power=2, 3 PO(s) with period: 10,4,3
        power=3, 3 PO(s) with period: 4,5,3
        power=4, 2 PO(s) with period: 5,2
        power=5, 1 PO(s) with period: 2
        power=6, 1 PO(s) with period: 2
        power=7, 1 PO(s) with period: 2
        power=8, 1 PO(s) with period: 1
        power=9, 1 PO(s) with period: 1
        power=10, 1 PO(s) with period: 1
        power=11, 1 PO(s) with period: 1
        power=12, 1 PO(s) with period: 1
        power=13, 1 PO(s) with period: 1
        power=14, 1 PO(s) with period: 1

    Writing the full content in string form to a file where it can be
    viewed with a suitable text editor::

        sage: d.to_file(filename='fullout1.txt', output_format='full')

    Construct partition functions for some well-known restriction lists,
    showing the needed time in seconds to get an impression about the
    performance of the algorithm. Additionally, save the full content
    of each partition function to a file and show the file size in bytes
    to get an impression about the space needed for the representation
    of the restricted partition function (Note: the comments random and/or
    long time are for automatic doctesting as the times are not exactly
    reproducible and some of them exceed one second)::

        sage: RestrictedPartitionFunction.print_time=True # we want to see the timing
        sage: pathname = os.path.join(SAGE_TMP, 'fullout')
        sage: A=[1, 2, 3, 4, 5, 6]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random time
        time: 0.016
        sage: d.to_file(pathname, 'full')
        sage: print os.path.getsize(pathname)
        598
        sage: A=[12223, 12224, 36674, 61119, 85569]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random long time
        time: 14.232
        sage: d.to_file(pathname, 'full')       # depends on long time line
        sage: print os.path.getsize(pathname)   # depends on long time line
        7952568
        sage: A=[12137, 24269, 36405, 36407, 48545, 60683]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random long time
        time: 17.108
        sage: d.to_file(pathname, 'full')       # depends on long time line
        sage: print os.path.getsize(pathname)   # depends on long time line
        7357507
        sage: A=[20601,40429,40429,45415,53725,61919,64470,69340,78539,95043]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random long time
        time: 93.444
        sage: d.to_file(pathname, 'full')       # depends on long time line
        sage: print os.path.getsize(pathname)   # depends on long time line
        42745525
        sage: A=[5,10,10,2,8,20,15,2,9,9,7,4,12,13,19]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random time
        time: 0.144
        sage: d.to_file(pathname, 'full')
        sage: print os.path.getsize(pathname)
        6295
        sage: A = [i for i in range(1,61)]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random long time
        time: 16.564
        sage: d.to_file(pathname, 'full')       # depends on long time line
        sage: print os.path.getsize(pathname)   # depends on long time line
        379127
        sage: A = [i for i in range(1,101)]
        sage: d=RestrictedPartitionFunction(parts_in=A) # random long time
        time: 130.536
        sage: d.to_file(pathname, 'full')       # depends on long time line
        sage: print os.path.getsize(pathname)   # depends on long time line
        2146284

    """
    # some class variables that influence verbosity during initialization
    # of a RestrictedPartitionFunction object
    print_time = True
    print_info = False
    print_detail = False
    # class variable which allows to switch off initial sorting of
    # parts_in list
    preserve_order = False

    def __init__(self, parts_in=None, poly_coeffs=None, empty=False):
        """
        Initialize ``self``.

        INPUT:

        - ``parts_in`` -- a list of positive integers
        - ``poly_coeffs`` -- for internal use only
        - ``empty`` -- for internal use only

        See documentation of RestrictedPartitionFunction for details.
        """
        if poly_coeffs is None and not empty:
            cls = self.__class__
            if not isinstance(parts_in, (list, tuple)):
                raise ValueError('%s is not a list or tuple'%parts_in)
            for act_a in parts_in:
                if not act_a in ZZ or not act_a > 0:
                    raise ValueError('%s is not a positive integer'%act_a)
            act_d = cls._build(parts_in)
            self._poly_coeffs = act_d._poly_coeffs
            self._parts_in = act_d._parts_in
            self._level = act_d._level
        else:
            self._parts_in = parts_in
            self._level = len(parts_in)
            if empty:
                self._poly_coeffs = []
                for i in range(self._level):
                    self._poly_coeffs.append(PeriodObjectList())
            else:
                # poly_coeffs given
                self._poly_coeffs = poly_coeffs

    def __iter__(self):
        return iter(self._poly_coeffs)

    def __getitem__(self, pos):
        return self._poly_coeffs[pos]

    def number_of_partitions(self, t):
        r"""
        Get the number of partitions of a nonnegative integer ``t`` with
        parts in a specific restriction list as contained in the current
        object of type ``RestrictedPartitionFunction``.

        Note that this gives the same result as
        ``Partitions(t, parts_in=A).cardinality()``, but the once
        constructed object of type ``RestrictedPartitionFunction`` can
        be used to call ``number_of_partitions`` consecutively for
        different values of ``t`` and the values of ``t`` can be
        considerably larger.

        INPUT:

        - ``t`` -- a nonnegative integer

        OUTPUT:

        - the number of partitions of ``t`` with parts in the restriction
          list as contained in the current object

        EXAMPLES:

        Construct a partition function for a specific restriction list
        and use it to get the number of partitions for several numbers::

            sage: RestrictedPartitionFunction.print_time=False  # suppress timing for doctest
            sage: A=[5,10,10,2,8,20,15,2,9,9,7,4,12,13,19]
            sage: d=RestrictedPartitionFunction(parts_in=A)
            sage: d.number_of_partitions(9)
            8
            sage: d.number_of_partitions(100000)
            3591162466582613140951155095382336216938321372
            sage: d.number_of_partitions(10**6)
            355852439929171688810360576684162614436054571631359229484001
            sage: d.number_of_partitions(10**7)
            35552753497121464994061233188681496764739760087879328128319410451219615958

        For increasing powers of 10**10, get the number of partitions for the
        same same restriction list and show the number of digits of the
        resulting numbers::

            sage: print [d.number_of_partitions(10**10**i).ndigits() for i in range(1,8)] # long time (1 min)
            [116, 1376, 13976, 139976, 1399976, 13999976, 139999976]

        """
        if not t in ZZ or not t > 0:
            raise ValueError('%s is not a nonnnegative integer'%t)
        retval = self._number_of_partitions(t)
        return ZZ(retval)

    def _number_of_partitions(self, t):
        """Internally used by number_of_partitons"""
        sval = 0
        for power in range(self._level-1, -1, -1):
            for pf in self._poly_coeffs[power]:
                sval += pf[t % pf.period]
            if power == 0: break
            sval *= t
        return sval

    def __repr__(self):
        number_pos = sum(len(pobl._list) for pobl in self._poly_coeffs)
        return (
            "RestrictedPartitionFunction for {0!r}, consisting of "
            "{1!r} lists with a total of {2!r} period objects".
            format(self._parts_in, len(self._poly_coeffs), number_pos)
        )

    def get_content_str(self, output_format='summary'):
        r"""
        Return a string representation of the content of an object of
        type ``RestrictedPartitionFunction`` where the amount and
        representation of the information is selected by the argument
        ``output_format``.

        INPUT:

        - ``output_format`` - (default: ``'summary'``)

          - ``'summary'`` -- show number of period objects and their period
          - ``'full'`` -- show full content, including values of period objects
          - ``'list'`` -- show full content in list form

        OUTPUT:

        - a string with the representation of the content of the object in
          the requested ``output_format``.

        EXAMPLES:

        Let us have a look at the various representations for a very small
        ``RestrictedParttionFunction`` object::

            sage: RestrictedPartitionFunction.print_time=False  # suppress timing for doctest
            sage: d=RestrictedPartitionFunction(parts_in=[1,2,3])

        The usual string representation does not reveal the content::

            sage: print d
            RestrictedPartitionFunction for [1, 2, 3], consisting of 3 lists with a total of 4 period objects

        With ``get_content_str(output_format='summary')`` or, equivalently,
        ``get_content_str()`` one gets the number of period objects and
        their periods for each power::

            sage: print d.get_content_str()
            RestrictedPartitionFunction for [1, 2, 3], consisting of 3 lists with a total of 4 period objects:
            power=0, 2 PO(s) with period: 2,3
            power=1, 1 PO(s) with period: 1
            power=2, 1 PO(s) with period: 1

        With ``get_content_str(output_format='full')``, the values stored in
        the period objects are shown additionally::

            sage: print d.get_content_str('full')
            RestrictedPartitionFunction for [1, 2, 3], consisting of 3 lists with a total of 4 period objects:
            power=0, 2 period object(s):
            P=2, values=[1/4, 0]
            P=3, values=[3/4, 5/12, 5/12]
            power=1, 1 period object(s):
            P=1, values=[1/2]
            power=2, 1 period object(s):
            P=1, values=[1/12]

        With ``output_format='list'`` the full information is formatted in
        list form, which may be better suited for a subsequent
        evaluation::

            sage: print d.get_content_str('list')
            [[1, 2, 3], [[[1/4, 0], [3/4, 5/12, 5/12]], [[1/2]], [[1/12]]]]

        We now create a ``RestrictedPartitionFunction`` for a larger
        restriction list and show the summary of its content::

            sage: A=[5, 10, 10, 2, 8, 20, 15, 2, 9, 9, 7, 4, 12, 13, 19]
            sage: d=RestrictedPartitionFunction(parts_in=A)
            sage: print d.get_content_str('summary')
            RestrictedPartitionFunction for [2, 5, 7, 9, 13, 19, 4, 2, 12, 8, 15, 10, 9, 20, 10], consisting of 15 lists with a total of 29 period objects:
            power=0, 8 PO(s) with period: 12,15,7,9,13,19,8,20
            power=1, 3 PO(s) with period: 10,4,9
            power=2, 3 PO(s) with period: 10,4,3
            power=3, 3 PO(s) with period: 4,5,3
            power=4, 2 PO(s) with period: 5,2
            power=5, 1 PO(s) with period: 2
            power=6, 1 PO(s) with period: 2
            power=7, 1 PO(s) with period: 2
            power=8, 1 PO(s) with period: 1
            power=9, 1 PO(s) with period: 1
            power=10, 1 PO(s) with period: 1
            power=11, 1 PO(s) with period: 1
            power=12, 1 PO(s) with period: 1
            power=13, 1 PO(s) with period: 1
            power=14, 1 PO(s) with period: 1

        Then we use the string created with ``output_format='list'`` to
        analyze the content of the partition function.
        First, we want to see the for example the second period object
        from power 4. The next example extracts the greatest value from
        the period objects of each power and shows them::

            sage: d_list_str = d.get_content_str(output_format='list')
            sage: d_list = sage_eval(d_list_str)
            sage: d_list[1][4][1]
            [145092185570084783/624488216002560000000,
             138207696216866033/624488216002560000000]
            sage: [(i, max([v for p in pl for v in p])) for i, pl in enumerate(d_list[1])]
            [(0, 3076481423136499854771839/506341290417035673600000),
             (1, 683728757821142046661/3091216669212672000000),
             (2, 33978463148176350767/1545608334606336000000),
             (3, 1591376310823161673/429335648501760000000),
             (4, 145092185570084783/624488216002560000000),
             (5, 10889579327/1069915392000000),
             (6, 4339575160381/14636442562560000000),
             (7, 761837855243/122946117525504000000),
             (8, 32978781481/351274621501440000000),
             (9, 1673677/1621267483852800000),
             (10, 11454959/1405098486005760000000),
             (11, 278429/6182433338425344000000),
             (12, 15313/92736500076380160000000),
             (13, 29/80371633399529472000000),
             (14, 1/2813007168983531520000000)]

        """
        if output_format == 'full':
            return (
                "{0}:\n".format(self) +
                "\n".join("power={0}, {1}".format(i, str(pobl))
                     for i, pobl in enumerate(self._poly_coeffs)
                )
            )
        elif output_format == 'list':
            return (
                "[{0}, {1}]".format(
                    "[{0}]".format(", ".join(str(a) for a in self._parts_in)),
                    "[{0}]".format(", ".
                        join("[{0}]".format(", ".
                            join("[{0}]".format(", ".
                                join(str(v) for v in po._vals))
                            for po in pobl._list))
                        for pobl in self._poly_coeffs))
                )
            )
        elif output_format == 'summary':
            return (
                "{0}:\n".format(self) +
                "\n".join("power={0}, {1}".format(i, pobl.get_summary())
                     for i, pobl in enumerate(self._poly_coeffs)
                )
            )
        else:
            raise ValueError(
                 "Invalid output_format specified: {0}. Possible values are: "
                 "{1}, {2}, {3}".format(output_format, 'summary', 'full', 'list')
            )

    @classmethod
    def build_next(cls, a, dl):
        r"""
        Build a new object of type ``RestrictedPartitionFunction`` from an
        existing one, with an additional positive integer ``a`` in the
        restriction list.

        Using ``build_next``, one can take advantage of an already
        constructed ``RestrictedPartitionFunction`` object that took a
        long time to construct if additions to the restriction list are
        needed later. However, remark that if the complete restriction list
        is known from the beginning, a considerably shorter overall
        construction time can be often achieved as compared to successive
        calls to ``build_next``, because the numbers from ``parts_in``
        are first put into an order that is presumably better suited
        for the calculations.

        Also remark that the content of a ``RestrictedPartitionFunction``
        object constructed by successive calls to ``build_next``
        will not necessarily be the same as when the list is given at the
        beginning, but it will always lead to the same results when used
        to get the ``number_of_partitions`` for the same number ``t``.

        INPUT:

        - ``a`` -- a positive integer
        - ``dl`` -- an object of type ``RestrictedPartitionFunction``

        OUTPUT:

        - an object of type ``RestrictedPartitionFunction``

        EXAMPLES:

        Construct partition function ``d1`` for a restriction list. Starting
        from ``d1``, construct ``d2`` with an additional number in the
        restriction list. Show summary of both restricted partition functions
        and use them to get number of partitions::

            sage: RestrictedPartitionFunction.print_time=False  # suppress timing for doctest
            sage: d1=RestrictedPartitionFunction(parts_in=[5, 6, 10])
            sage: d2=RestrictedPartitionFunction.build_next(a=2, dl=d1)
            sage: print d1.get_content_str('summary')
            RestrictedPartitionFunction for [5, 6, 10], consisting of 3 lists with a total of 5 period objects:
            power=0, 2 PO(s) with period: 10,6
            power=1, 2 PO(s) with period: 5,2
            power=2, 1 PO(s) with period: 1
            sage: print d2.get_content_str('summary')
            RestrictedPartitionFunction for [5, 6, 10, 2], consisting of 4 lists with a total of 6 period objects:
            power=0, 2 PO(s) with period: 10,6
            power=1, 2 PO(s) with period: 2,5
            power=2, 1 PO(s) with period: 2
            power=3, 1 PO(s) with period: 1
            sage: d1.number_of_partitions(11)
            1
            sage: d2.number_of_partitions(11)
            2

        Now construct ``d3`` for the same restriction list as the resulting
        list for ``d2``.
        Verify that ``d2`` and ``d3`` give the same
        ``number_of_partitions`` for several input numbers.
        Note that the content of ``d2`` and ``d3`` is not the same::

            sage: d3=RestrictedPartitionFunction(parts_in=[5, 6, 10, 2])
            sage: len([t for t in [100..200] if d2.number_of_partitions(t) != d3.number_of_partitions(t)])
            0
            sage: print d3.get_content_str('summary')
            RestrictedPartitionFunction for [2, 5, 6, 10], consisting of 4 lists with a total of 6 period objects:
            power=0, 2 PO(s) with period: 6,10
            power=1, 2 PO(s) with period: 2,5
            power=2, 1 PO(s) with period: 2
            power=3, 1 PO(s) with period: 1
            sage: print d2[1]
            2 period object(s):
            P=2, values=[71/600, 9/400]
            P=5, values=[29/600, 17/600, 29/600, 1/120, 1/120]
            sage: print d3[1]
            2 period object(s):
            P=2, values=[19/150, 37/1200]
            P=5, values=[1/25, 1/50, 1/25, 0, 0]

        """
        if not a in ZZ or not a > 0:
            raise ValueError('%s is not a positive integer'%a)
        if not isinstance(dl, cls):
            raise ValueError('Parameter dl is not of type %s'%cls)
        return cls._build_next(a, dl)

    @classmethod
    def _build_next(cls, a, dl):
        r"""
        Build a new object of type ``RestrictedPartitionFunction`` from an
        existing one. For interal use.
        """
        l = dl._level
        if cls.print_info:
            print "Starting calculation for level {0}, a = {1}".format(l+1, a)
        dlplus1 = cls(parts_in=dl._parts_in + [a], empty=True)
        for m in range(l):
            for fmi in dl[m]:
                P = fmi.period
                if m == 0:
                    Fmi_mplus1, Fmis = _find_Fs_m0(fmi, m, a, cls.print_detail)
                elif (m*m*P + 2*a) * 3 > m*P*P:
                    Fmi_mplus1, Fmis = _find_Fs_std(fmi, m, a, cls.print_detail)
                else:
                    Fmi_mplus1, Fmis = _find_Fs_qp(fmi, m, a, cls.print_detail)
                dlplus1[m+1].merge_or_append(Fmi_mplus1)
                for s in range(m+1):
                    dlplus1[s].merge_or_append(Fmis[s])
        if cls.print_detail:
            print "Level finalization",
        t1 = cputime()
        correction_po = PeriodObject(
            [dl._number_of_partitions(i) - dlplus1._number_of_partitions(i)
            for i in range(a)]
        )
        if cls.print_detail:
            print " time: " + str(cputime(t1))
        dlplus1[0].merge_or_append(correction_po)
        return dlplus1

    @classmethod
    def _build(cls, parts_in):
        r"""
        Build an object of type ``RestrictedPartitionFunction`` for a given
        restriction list. For interal use.

        INPUT:

        - `parts_in`` -- a list of positive integers

        Output:

        - an object of type ``RestrictedPartitionFunction``

        """
        if cls.print_info:
            print ("Building restricted partition function for the "
                   "{0} parts: {1}").format(len(parts_in), parts_in)
        t = cputime()
        if not cls.preserve_order:
            if cls.print_info:
                print "Sorting Input list"
            parts_in = build_wrk(parts_in)
            if cls.print_info:
                print "Sorted list: " + str(parts_in)
        a = parts_in[0]
        act_d = cls(
            parts_in=[a],
            poly_coeffs=[
                PeriodObjectList(list_=[PeriodObject(values=(1,), period=a)]            )
            ]
        )
        for act_a in parts_in[1:]:
            act_d = cls._build_next(act_a, act_d)
        if cls.print_time:
            print "time: " + str(cputime(t))
        return act_d

    def to_file(self, filename, output_format='full'):
        r"""
        Write the string representation of the full content of an
        object of type ``RestrictedPartitionFunction`` to a file.

        INPUT:

        - ``filename`` -- name of file
        - ``output_format`` -- (default: ``'full'``)

          - ``'full'`` -- write full content, including values of period objects
          - ``'list'`` -- write full content in list form

        OUTPUT:

        A string representation of the full content of the object is
        written to the specified file in the requested ``output_format``

        EXAMPLES:

        Create a restricted partition function object and save its full
        content to a file::

            sage: d=RestrictedPartitionFunction(parts_in=[6,10,31])
            sage: d
            RestrictedPartitionFunction for [6, 31, 10], consisting of 3 lists with a total of 5 period objects
            sage: pathname = os.path.join(SAGE_TMP, 'outfile1')
            sage: d.to_file(pathname)

        """
        if output_format not in ('full', 'list'):
            raise ValueError(
                 "Invalid output_format specified: {0}. Possible values are: "
                 "{1}, {2}".format(output_format, 'full', 'list')
            )
        fh = open(filename, 'w')
        fh.write(self.get_content_str(output_format=output_format))
        fh.close()


    @classmethod
    def from_file(cls, filename):
        r"""
        Create a new object of type ``RestrictedPartitionFunction``
        from a file that has been created before by means of
        the object method ``to_file``.

        INPUT:

        - ``filename`` -- name of file

        OUTPUT:

        - an object of type ``RestrictedPartitionFunction``

        EXAMPLES:

        Suppose the example from the description of ``to_file`` has been
        run in a previous session. We simulate this by assigning another
        value to ``d``. Then we reconstruct ``d`` from the file:: 

            sage: d=1
            sage: d
            1
            sage: pathname = os.path.join(SAGE_TMP, 'outfile1')
            sage: d=RestrictedPartitionFunction.from_file(pathname)
            sage: d
            RestrictedPartitionFunction for [6, 31, 10], consisting of 3 lists with a total of 5 period objects

        """
        parts_in = []
        poly_coeffs = []
        fh = open(filename, 'r')
        line1 = fh.readline()
        if line1.startswith('[['):
            # list form
            content = sage_eval(line1)
            parts_in = content[0]
            poly_coeffs = [PeriodObjectList(list_=[PeriodObject(po)
                for po in pobl]) for pobl in content[1]]

        elif line1.startswith('RestrictedPartitionFunction'):
            parts_in = [ZZ(vs) for vs in
                ((line1.partition('[')[2]).partition(']')[0]).split(', ')]
            pobs = []
            actp = 0
            for line in fh:
                if line.startswith('power'):
                    if actp > 0:
                        poly_coeffs.append(
                            PeriodObjectList(list_=pobs[:]))
                    actp += 1
                    pobs = []
                elif line.startswith('P='):
                    postr = (line.partition('[')[2]).partition(']')[0]
                    vals = [QQ(vs) for vs in postr.split(', ')]
                    pobs.append(PeriodObject(vals))
            poly_coeffs.append(PeriodObjectList(list_=pobs[:]))
        fh.close()
        return cls(parts_in=parts_in, poly_coeffs=poly_coeffs)


def _find_Fs_std(f, m, a, print_detail):
    P = f.period
    if print_detail:
        print "algo=std (m={0}, P={1})".format(m, P),
    t1 = cputime()
    gcdaP = gcd(a, P)
    lcmaP = a / gcdaP * P
    steps = P / gcdaP
    Fma_mplus1 = PeriodObject(
        values=[sum([f[t+j*gcdaP] for j in range(steps)]) / (m+1) / lcmaP
                   for t in range(gcdaP)]
    )
    gks = find_gks_0_to_m(P, m, a)
    Fmas = PeriodObjectS(m+1, P)
    for s in range(m+1):
        for t in range(P):
            x = 0
            for z in range(t, t+P, gcdaP):
                if f[z % P] != 0:
                    x += gks[s][(z-t)%P] * f[z%P]
            Fmas[s][t] = x
    if print_detail:
        print " time: " + str(cputime(t1))
    return Fma_mplus1, Fmas

def _find_Fs_qp(f, m, a, print_detail):
    P = f.period
    if print_detail:
        print "algo=qp  (m={0}, P={1})".format(m, P),
    t1 = cputime()
    gcdaP = gcd(a, P)
    lcmaP = a / gcdaP * P
    steps = P / gcdaP
    Fma_mplus1 = PeriodObject(period=gcdaP)
    for t in range(gcdaP):
        val = 0
        for j in range(steps):
            val += f[t+j*gcdaP]
        Fma_mplus1[t] = (val / (m+1)) / lcmaP
    Fvals = [0 for t in range((m+1)*P + a)]
    for l in range(a):
        val = 0
        for pos in range(l, (m+1)*P + a, a):
            val += f[pos % P] * pos**m
            Fvals[pos] = val
    FDvals = [Fvals[t] - Fma_mplus1[t % gcdaP] * t**(m+1)
                 for t in range((m+1)*P + a)
             ]
    FDDeltaVals = [FDvals[t+a] - FDvals[t]
                      for t in range((m+1)*P)
                  ]
    Hcl = find_qp_coeffs(m, P, FDDeltaVals)
    Fmas = find_hks(P, m, a, Hcl)
    if print_detail:
        print " time: " + str(cputime(t1))
    return Fma_mplus1, Fmas


def _find_Fs_m0(f, m, a, print_detail):
    P = f.period
    if print_detail:
        print "algo=m0  (m={0}, P={1})".format(m, P),
    t1 = cputime()
    gcdaP = gcd(a, P)
    steps = P/gcdaP
    mfa = [sum([f[(t+j*a)%P] for j in range(steps)]) / P * gcdaP
              for t in range(gcdaP)
          ]
    g0 = PeriodObjectS(1, P)
    if a % P != 0:
        for l in range(gcdaP):
            val = 0
            pos = l
            for step in range(steps):
                val += f[pos] - mfa[l]
                g0[0][pos] = val
                pos = (pos + a) % P
    g1 = PeriodObject([mfa[t] / a for t in range(gcdaP)])
    if print_detail:
        print " time: " + str(cputime(t1))
    return g1, g0


@cached_function
def _cached_binomial(m, n):
    return binomial(m, n)

@cached_function
def _cached_bernoulli(n):
    return bernoulli(n)

def find_gks_0_to_m(P, m, a):
    gks = PeriodObjectS(m+1, P)
    gcdaP = gcd(a, P)
    lcmaP = a / gcdaP * P
    steps = P / gcdaP
    lcmaP_pot = [1/lcmaP, 1]
    tmp = 1
    for i in range(m):
        tmp *= lcmaP
        lcmaP_pot.append(tmp)
    mp1 = m + 1
    mp1_is_even = mp1 % 2 == 0
    for r in range(steps):
        alfa = - r * a
        z = (lcmaP + alfa) % P
        mp1_pk_is_even = mp1_is_even
        for k in range(mp1):
            val = 0
            act_alfapot = 1
            is_even = mp1_pk_is_even
            for u in range(mp1-k+1):
                s = k + u
                if s != 0:
                    tmpval = (_cached_binomial(s, k) *
                        act_alfapot / (mp1) *
                        _cached_binomial(mp1, mp1-s) *
                        _cached_bernoulli(mp1-s) *
                        lcmaP_pot[(m-s) + 1])
                    if is_even:
                        val += tmpval
                    else:
                        val -= tmpval
                is_even = not is_even
                if u == mp1-k:
                    break
                act_alfapot *= alfa
            mp1_pk_is_even = not mp1_pk_is_even
            gks[k][z] = val
    return gks

def find_qp_coeffs(m, p, vals):
    if len(vals) != (m+1)*p:
        print ("Error, incorrect number of vals: " + str(len(vals)) +
               ", should be " + str((m+1)*p))
        return None
    coeff_list = []
    for t in range(p):
        coeff = [0 for r in range(m+1)]
        for r in range(m+1):
            coeff_prod = 1
            tmp_coeff = [1] + [0 for s in range(1, m+1)]
            for s in range(m+1):
                if s != r:
                    for k in range(m, 0, -1):
                        tmp_coeff[k] = tmp_coeff[k] * (-t-s*p) + tmp_coeff[k-1]
                    tmp_coeff[0] = tmp_coeff[0] * (-t-s*p)
                    coeff_prod *= (r-s)*p
            for j in range(m+1):
                coeff[j] += tmp_coeff[j] / coeff_prod * vals[t + r*p]
        coeff_list.append(coeff)
    return coeff_list

def find_hks(P, m, a, Hcl):
    gcdaP = gcd(a, P)
    lcmaP = a / gcdaP * P
    h = PeriodObjectS(m+1, P)
    TPA = [[0 for i in range(m+1)] for j in range(m+1)]
    for i in range(m+1):
        for k in range(i+1):
            TPA[k][i] = _cached_binomial(i, k) * (-a)**(i-k)
    TPD = [[0 for i in range(m+1)] for j in range(m+1)]
    for i in range(m+1):
        for k in range(i+1):
            TPD[k][i] = _cached_binomial(i, k) * (-lcmaP)**(i-k)

    steps = P / gcdaP
    for t in range(gcdaP):
        x = [0 for i in range(m+1)]
        pos = t
        for step in range(steps):
            y = [x[i] + Hcl[pos][i] for i in range(m+1)]
            x = [0 for k in range(m+1)]
            for k in range(m+1):
                for i in range(m+1):
                    x[k] += y[i] * TPA[k][i]
            pos = (pos + a) % P
        h[m][t] = -1/TPD[m-1][m] * x[m-1]
        for k in range(m-1, 0, -1):
            val = x[k-1]
            for i in range(k+1, m+1):
                val += TPD[k-1][i] * h[i][t]
            h[k][t] = -1/TPD[k-1][k] * val
        h[0][t] = 0
        pos = t
        while (pos + a) % P != t:
            y = [h[i][pos] + Hcl[pos][i] for i in range(m+1)]
            pos = (pos + a) % P
            for k in range(m+1):
                for i in range(m+1):
                    h[k][pos] += y[i] * TPA[k][i]
    return h

def build_wrk(input_list):
    l1 = [[gcd(a, b), i, j] for i, a in enumerate(input_list)
            for j, b in enumerate(input_list)
            if j > i
    ]
    gcds = sorted(list(set([el[0] for el in l1 if el[0] > 1])), reverse=True)
    last_part_is = []
    for gcd_val in gcds:
        lsg = [elem for elem in l1 if elem[0] == gcd_val]
        lsg_is = list(set([e[1] for e in lsg]) | set([e[2] for e in lsg]))
        found_fpi = False
        for i in [e[1] for e in sorted([(input_list[i], i) for i in lsg_is])]:
            if not found_fpi:
                if not i in last_part_is:
                    found_fpi = True
            else:
                if not i in last_part_is:
                    last_part_is.append(i)
    result = []
    for i in range(len(input_list)):
        if not i in last_part_is:
            result.append(input_list[i])
    result.sort()
    for i in reversed(last_part_is):
        result.append(input_list[i])
    return result

class PeriodObject(object):
    """Holds a period object."""
    def __init__(self, values=None, period=None):
        if period is None:
            # no period - get it from length of values
            if isinstance(values, list):
                self._vals = values
            else:
                self._vals = list(values)
            self.period = len(self._vals)
        else:
            self.period = period
            if values is None:
                # no values - make them all zeroes
                self._vals = [0 for _ in range(period)]
            elif isinstance(values, list) and len(values) == period:
                self._vals = values
            else:
                # take values and fill with zeroes or truncate to match period
                vallist = list(values)
                numvals = len(vallist)
                if numvals <= period:
                    vallist.extend([0 for t in range(period-numvals)])
                    self._vals = vallist
                else:
                    self._vals = vallist[:period]
    def __str__(self):
        return "P={0.period}, values={0._vals}".format(self)
    def __repr__(self):
        return ("PeriodObject(period={0.period!r},"
                "values={0._vals!r})".format(self))
    def __iter__(self):
        return iter(self._vals)
    def __getitem__(self, pos):
        return self._vals[pos]
    def __setitem__(self, pos, val):
        self._vals[pos] = val

class PeriodObjectS(object):
    """Holds a list of period objects with same period."""
    def __init__(self, s, P):
        self._pfs = [PeriodObject(period=P) for k in range(s)]
    def __iter__(self):
        return iter(self._pfs)
    def __getitem__(self, pos):
        return self._pfs[pos]

class PeriodObjectList(object):
    """Holds a list of period objects for a power."""
    def __init__(self, list_=None):
        if list_ is None:
            self._list = []
        else:
            self._list = list_
    def __str__(self):
        return (
            "{0} period object(s):\n".format(len(self._list)) +
            "\n".join(str(po) for po in self._list)
        )
    def get_summary(self):
        return (
            "{0} PO(s) with period: ".format(len(self._list)) +
            ",".join(str(po.period) for po in self._list)
        )
    def __repr__(self):
        return ("PeriodObjectList(list_={0._list!r})").format(self)
    def __iter__(self):
        return iter(self._list)
    def __getitem__(self, pos):
        return self._list[pos]
    def append(self, new_po):
        self._list.append(new_po)
    def merge_or_append(self, new_po):
        found_contained = False
        contained_ixs = []
        new_p = new_po.period
        for i, ex_po in enumerate(self._list):
            ex_p = ex_po.period
            if ex_p == new_p:
                for j in range(ex_p):
                    ex_po[j] += new_po[j]
                return
            elif ex_p % new_p == 0:
                for j in range(ex_p):
                    ex_po[j] += new_po[j % new_p]
                return
            elif new_p % ex_p == 0:
                found_contained = True
                contained_ixs.append(i)
        if found_contained:
            if len(contained_ixs) == 1:
                ex_po = self._list[contained_ixs[0]]
                ex_p = ex_po.period
                for j in range(new_p):
                    new_po[j] += ex_po[j % ex_p]
                self._list[contained_ixs[0]] = new_po
                return
            new_list = []
            for i, ex_po in enumerate(self._list):
                if i in contained_ixs:
                    ex_p = ex_po.period
                    for j in range(new_p):
                        new_po[j] += ex_po[j % ex_p]
                else:
                    new_list.append(ex_po)
            new_list.append(new_po)
            self._list = new_list
            return
        self._list.append(new_po)

