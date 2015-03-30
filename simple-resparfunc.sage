"""
simple-resparfunc

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

    This class implements the simple version of the algorithm and
    is not optimized. Its only purpose is to illustrate the idea of 
    the algorithm. An optimized version is also available.

    Once constructed with ``RestrictedPartitionFunction(parts_in=A)``
    for a specific restriction list ``A``, a
    ``RestrictedPartitionFunction`` object can be used to calculate the
    number of partitions for an integer ``t`` using the object's method
    ``.number_of_partitions(t)``.
    Note that this gives the same result as
    ``Partitions(t, parts_in=A).cardinality()``, but the once
    constructed object can be reused for different values of ``t`` and
    the values of ``t`` can be considerably larger.

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

    If a ``RestrictedPartitionFunction`` object for a given restriction
    list has been already constructed, one can obtain a new
    ``RestrictedPartitionFunction`` for a restriction list with an
    additional positive integer ``a`` with the class method
    ``build_next()``.


    EXAMPLES:

    Construct some partition functions for specific restriction lists
    and use them to get the number of partitions for several numbers::

        sage: d=RestrictedPartitionFunction([1,2,3])
        sage: d.number_of_partitions(3)
        3
        sage: d.number_of_partitions(10000)
        8338334
        sage: d.number_of_partitions(10**20)
        833333333333333333383333333333333333334
        sage: d.number_of_partitions(10**10**6).ndigits()
        1999999

    Comparing ``RestrictedPartitionFunction().number_of_partitions()`` with
    ``Partitions().cardinality()``. The result must be the same::

        sage: A=[10,12,13]
        sage: d=RestrictedPartitionFunction(parts_in=A)
        sage: t=500
        sage: d.number_of_partitions(t) == Partitions(t, parts_in=A).cardinality()
        True
        sage: len([t for t in [100..200] if d.number_of_partitions(t) != Partitions(t, parts_in=A).cardinality()])
        0

    Showing the content of a ``RestrictedPartitionFunction``::

        sage: d = RestrictedPartitionFunction(parts_in=[10, 5, 6])
        sage: print d.get_content_str(output_format='full')
        RestrictedPartitionFunction for [10, 5, 6], consisting of 3 lists with a total of 5 period objects:
        power=0, 2 period object(s):
        P=10, values=[0, -17/50, -19/25, -33/50, -26/25, 0, -11/25, -14/25, -24/25, -16/25]
        P=6, values=[1, 7/24, 2/3, 5/8, 1, 5/8]
        power=1, 2 period object(s):
        P=2, values=[1/30, 1/60]
        P=5, values=[1/20, 3/100, 1/100, -1/100, -3/100]
        power=2, 1 period object(s):
        P=1, values=[1/600]

    Showing summary of period objects contained in the representation of an
    object of type ``RestrictedPartitionFunction``::

        sage: A = [i for i in range(1,12)]
        sage: d=RestrictedPartitionFunction(parts_in=A)
        sage: print d.get_content_str(output_format='summary')
        RestrictedPartitionFunction for [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], consisting of 11 lists with a total of 19 period objects:
        power=0, 6 PO(s) with period: 8,11,10,9,7,6
        power=1, 3 PO(s) with period: 4,5,3
        power=2, 2 PO(s) with period: 2,3
        power=3, 1 PO(s) with period: 2
        power=4, 1 PO(s) with period: 2
        power=5, 1 PO(s) with period: 1
        power=6, 1 PO(s) with period: 1
        power=7, 1 PO(s) with period: 1
        power=8, 1 PO(s) with period: 1
        power=9, 1 PO(s) with period: 1
        power=10, 1 PO(s) with period: 1

   """

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
                # poly_coeffs is not None:
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

        INPUT:

        - ``t`` -- a nonnegative integer

        OUTPUT:

        - the number of partitions of ``t`` with parts in the restriction
          list as contained in the current object
          
        """
        if not t in ZZ or not t > 0:
            raise ValueError('%s is not a nonnnegative integer'%t)
        sval = 0
        for power in range(self._level-1, -1, -1):
            for pf in self._poly_coeffs[power]:
                sval += pf[t % pf.period]
            if power == 0: break
            sval *= t
        return ZZ(sval)

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
        needed later.
        
        INPUT:

        - ``a`` -- a positive integer
        - ``dl`` -- an object of type ``RestrictedPartitionFunction``

        OUTPUT:

        - an object of type ``RestrictedPartitionFunction``

        """
        dlplus1 = cls(parts_in=dl._parts_in+[a], empty=True)
        for k in [1..dl._level]:
            for m in [k-1..dl._level-1]:
                for fmi in dl[m]:
                    for z in range(fmi.period):
                        gk = calc_g(k, fmi.period, z, m, a)
                        for t in range(gk.period):
                            gk[t] *= fmi[z]
                        dlplus1[k].merge_or_append(gk)
        for m in [0..dl._level-1]:
            for fmi in dl[m]:
                for z in range(fmi.period):
                    g0 = calc_g(0, fmi.period, z, m, a)
                    for t in range(g0.period):
                        g0[t] *= fmi[z]
                    dlplus1[0].merge_or_append(g0)
                    gd0 = calc_gd0(fmi.period, z, m, a)
                    for t in range(gd0.period):
                        gd0[t] *= fmi[z]
                    dlplus1[0].merge_or_append(gd0)
        return dlplus1

    @classmethod
    def _build(cls, parts_in):
        r"""
        Build an object of type ``RestrictedPartitionFunction`` for a given
        restriction list. For interal use.

        INPUT:

        - `parts_in`` -- a list of positive integers

        Output:

        An object of type ``RestrictedPartitionFunction``.
        """
        a = parts_in[0]
        act_d = cls([a], poly_coeffs=[PeriodObjectList(
                list_=[PeriodObject(values=(1,), period=a)])])
        for act_a in parts_in[1:]:
            act_d = cls.build_next(act_a, act_d)
        return act_d

def cms(m, s, a, P):
    if s == 0: return 0
    return (-1)**(m+1-s) / (m+1) * binomial(m+1, m+1-s) \
           * bernoulli(m+1-s) * lcm(a, P)**(m-s)

def find_alfa(t, P, z, a):
    if (t - z) % gcd(a, P) != 0: return 0
    else:
        result = t
        while (result % P) != z: result -= a
        result -= t
        return result

def find_beta(t, P, z, a):
    if (t - z) % gcd(a, P) != 0: return 0
    else:
        result = t % a
        while (result % P) != z: result += a
        return result

def calc_g(k, P, z, m, a):
    if k == m + 1:
        gk = PeriodObject(period=gcd(a, P))
        val = 1 / (m+1) / lcm(a, P)
        for t in range(gcd(a, P)):
            if (t - z) % gcd(a, P) == 0: gk[t] = val
        return gk
    gk = PeriodObject(period=P)
    for t in range(P):
        if (t - z) % gcd(a, P) == 0:
            alfa = find_alfa(t, P, z, a)
            gk[t] = sum(binomial(k+u, k) * cms(m, k+u, a, P) *
                        alfa**u for u in [0..m+1-k])
    return gk

def calc_gd0(P, z, m, a):
    gd0 = PeriodObject(period=a)
    for t in range(a):
        if (t - z) % gcd(a, P) == 0:
            beta = find_beta(t, P, z, a)
            gd0[t] = -sum(cms(m, s, a, P) * beta**s
                          for s in [1..m+1]) + beta**m
    return gd0


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
        return ("P={0.period}, values={0._vals}".format(self))
    def __repr__(self):
        return ("PeriodObject(period={0.period!r},"
                "values={0._vals!r})".format(self))
    def __iter__(self):
        return iter(self._vals)
    def __getitem__(self, pos):
        return self._vals[pos]
    def __setitem__(self, pos, val):
        self._vals[pos] = val


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

