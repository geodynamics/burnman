# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import operator
import bisect
import pkgutil
from collections import Counter, OrderedDict
import numpy as np
from .math import round_to_n, linear_interpol


class OrderedCounter(Counter, OrderedDict):
    """
    Counter that remembers the order elements are first encountered
    """

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)


def copy_documentation(copy_from):
    """
    Decorator @copy_documentation(another_function) will copy the documentation found in a different
    function (for example from a base class). The docstring applied to some function a() will be ::

        (copied from BaseClass.some_function):
        <documentation from BaseClass.some_function>
        <optionally the documentation found in a()>

    """

    def mydecorator(func):
        def wrapper(*args):
            return func(*args)

        old = ""
        if func.__doc__:
            old = "\n" + func.__doc__

        copied_from = ""
        if hasattr(copy_from, "__name__"):
            copied_from = "(copied from " + copy_from.__name__ + "):\n"
        wrapper.__doc__ = copied_from + copy_from.__doc__ + old
        wrapper.__name__ = func.__name__
        return wrapper

    return mydecorator


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def flatten(arr):
    return (
        flatten(arr[0]) + (flatten(arr[1:]) if len(arr) > 1 else [])
        if type(arr) is list or type(arr) is np.ndarray
        else [arr]
    )


def pretty_print_values(popt, pcov, params, extra_decimal_places=0):
    """
    Takes a numpy array of parameters, the corresponding covariance matrix
    and a set of parameter names and prints the parameters and
    principal 1-s.d.uncertainties (np.sqrt(pcov[i][i]))
    in a nice text based format.
    """
    for i, p in enumerate(params):
        p_rnd = round_to_n(popt[i], np.sqrt(pcov[i][i]), 1 + extra_decimal_places)
        c_rnd = round_to_n(
            np.sqrt(pcov[i][i]), np.sqrt(pcov[i][i]), 1 + extra_decimal_places
        )

        if p_rnd != 0.0:
            p_expnt = np.floor(np.log10(np.abs(p_rnd)))
        else:
            p_expnt = 0.0

        scale = np.power(10.0, p_expnt)
        nd = p_expnt - np.floor(np.log10(np.abs(c_rnd))) + extra_decimal_places
        print(
            "{0:s}: ({1:{4}{5}f} +/- {2:{4}{5}f}) x {3:.0e}".format(
                p, p_rnd / scale, c_rnd / scale, scale, 0, (nd) / 10.0
            )
        )


def pretty_print_table(table, use_tabs=False):
    """
    Takes a 2d table and prints it in a nice text based format. If
    use_tabs=True then only \t is used as a separator. This is useful for
    importing the data into other apps (Excel, ...). The default is to pad
    the columns with spaces to make them look neat. The first column is
    left aligned, while the remainder is right aligned.
    """
    if use_tabs:
        for r in table:
            print("\t".join(r).replace("_", "\\_"))
        return

    def col_width(table, colidx):
        return max([len(str(row[colidx])) for row in table])

    # create a format string with the first column left aligned, the others right
    # example:   {:<27}{:>11}{:>6}{:>8}
    frmt = "".join(
        [
            ("{:<" if i == 0 else "{:>") + str(1 + col_width(table, i)) + "}"
            for i in range(len(table[0]))
        ]
    )
    for r in table:
        print(frmt.format(*r))


def sort_table(table, col=0):
    """
    Sort the table according to the column number
    """
    return sorted(table, key=operator.itemgetter(col))


def read_table(filename):
    datastream = pkgutil.get_data("burnman", "data/" + filename)
    datalines = [
        line.strip() for line in datastream.decode("ascii").split("\n") if line.strip()
    ]
    table = []

    for line in datalines:
        if line[0] != "#":
            numbers = np.fromstring(line, sep=" ")
            table.append(numbers)
    return np.array(table)


def array_from_file(filename):
    """
    Generic function to read a file containing floats and commented lines
    into a 2D numpy array.

    Commented lines are prefixed by the characters # or %.
    """
    f = open(filename, "r")
    data = []
    datastream = f.read()
    f.close()
    datalines = [
        line.strip().split() for line in datastream.split("\n") if line.strip()
    ]
    for line in datalines:
        if line[0] != "#" and line[0] != "%":
            data.append(map(float, line))

    data = np.array(zip(*data))
    return data


def cut_table(table, min_value, max_value):
    tablen = []
    for i in range(min_value, max_value, 1):
        tablen.append(table[i, :])
    return tablen


def lookup_and_interpolate(table_x, table_y, x_value):
    idx = bisect.bisect_left(table_x, x_value) - 1
    if idx < 0:
        return table_y[0]
    elif idx < len(table_x) - 1:
        return linear_interpol(
            x_value, table_x[idx], table_x[idx + 1], table_y[idx], table_y[idx + 1]
        )
    else:
        return table_y[idx]


def attribute_function(m, attributes, powers=[]):
    """
    Function which returns a function which can be used to
    evaluate material properties at a point. This function
    allows the user to define the property returned
    as a string. The function can itself be passed to another
    function
    (such as nonlinear_fitting.confidence_prediction_bands()).

    Properties can either be simple attributes (e.g. K_T) or
    a product of attributes, each raised to some power.

    :param m: The material instance evaluated by the output function.
    :type m: :class:`burnman.Material`

    :param attributes: The list of material attributes / properties to
        be evaluated in the product.
    :type attributes: list of str

    :param powers: The powers to which each attribute should be raised
        during evaluation.
    :type powers: list of floats

    :returns: Function which returns the value of product(a_i**p_i)
        as a function of condition (x = [P, T, V]).
    :rtype: function
    """
    if type(attributes) is str:
        attributes = [attributes]
    if len(powers) == 0:
        powers = [1.0 for a in attributes]

    def f(x):
        P, T, V = x
        m.set_state(P, T)
        value = 1.0
        for a, p in zip(*[attributes, powers]):
            value *= np.power(getattr(m, a), p)
        return value

    return f
