# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import subprocess
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


def pretty_string_values(
    popt, pcov, extra_decimal_places=0, combine_value_and_sigma=False
):
    """
    Takes a numpy array of parameters, the corresponding covariance matrix
    and a set of parameter names and returns the scaled variables and
    principal 1-s.d.uncertainties (np.sqrt(pcov[i][i])) and scaling factor
    as three separate lists of strings.

    :param popt: Parameter values
    :type popt: numpy array
    :param pcov: Variance-covariance matrix
    :type pcov: 2D numpy array
    :param extra_decimal_places: extra precision for values, defaults to 0
    :type extra_decimal_places: int, optional
    :param combine_value_and_sigma: give values in value(sigma) format, defaults to False
    :type combine_value_and_sigma: bool, optional
    :return: values, uncertainties and the scaling factors
    :rtype: tuple of 3 lists
    """
    pval = []
    psig = []
    pscale = []
    for i, p in enumerate(popt):
        p_rnd = round_to_n(p, np.sqrt(pcov[i][i]), 1 + extra_decimal_places)
        c_rnd = round_to_n(
            np.sqrt(pcov[i][i]), np.sqrt(pcov[i][i]), 1 + extra_decimal_places
        )

        if p_rnd != 0.0:
            p_expnt = np.floor(np.log10(np.abs(p_rnd)))
        else:
            p_expnt = 0.0

        scale = np.power(10.0, p_expnt)
        nd = p_expnt - np.floor(np.log10(np.abs(c_rnd))) + extra_decimal_places
        pval.append(f"{p_rnd / scale:0{nd / 10.0}f}")
        if combine_value_and_sigma:
            pval[-1] = f"{pval[-1]}({int(c_rnd / scale * np.power(10, nd))})"
        psig.append(f"{c_rnd / scale:0{nd / 10.0}f}")
        pscale.append(f"{scale:.0e}")
    return (pval, psig, pscale)


def pretty_print_values(popt, pcov, params, extra_decimal_places=0):
    """
    Takes a numpy array of parameters, the corresponding covariance matrix
    and a set of parameter names and prints the scaled variables and
    principal 1-s.d.uncertainties (np.sqrt(pcov[i][i])) and scaling factor
    in an easy to read format.

    :param popt: Parameter values
    :type popt: numpy array
    :param pcov: Variance-covariance matrix
    :type pcov: 2D numpy array
    :param extra_decimal_places: extra precision for values, defaults to 0
    :type extra_decimal_places: int, optional
    """
    pval, psig, pscale = pretty_string_values(popt, pcov, extra_decimal_places)
    for i, p in enumerate(params):
        print(f"{p:s}: ({pval[i]} +/- {psig[i]}) x {pscale[i]}")


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


def attribute_function(m, attributes, powers=[], using_volume=False):
    """
    Returns a callable function that evaluates a derived material property
    defined as a product of material attributes (optionally raised to given
    powers) at specified thermodynamic conditions.

    The returned function is designed to be compatible with
    :func:`burnman.optimize.nonlinear_fitting.confidence_prediction_bands`,
    which expects functions of the form ``f(x)`` where `x` contains
    composition and thermodynamic state variables.

    The expected order of arguments in `x` depends on the value of
    `using_volume`:
    - If ``using_volume=False``, the function expects ``x = [X1, X2, ..., P, T, _]`` and the state of the material is set using `set_state`
    - If ``using_volume=True``, the function expects ``x = [X1, X2, ..., V, T, _]`` and the state of the material is set using `set_state_with_volume`

    The function sets the state of the material using the appropriate method,
    evaluates the specified attributes, raises them to the corresponding
    powers, and returns the product.

    :param m: The material instance whose attributes are evaluated.
    :type m: :class:`burnman.Material`

    :param attributes: List of material attributes to include in the product
        expression (e.g., ['K_T', 'rho']).
    :type attributes: list of str

    :param powers: List of exponents to raise each corresponding attribute to.
        If empty, defaults to 1.0 for each attribute.
    :type powers: list of float

    :param using_volume: If True, the input function defines the
        thermodynamic state using volume and temperature (V, T).
        If False, pressure and temperature (P, T) are used.
        This argument also determines the expected order of the
        final three elements in the input vector `x`.
    :type using_volume: bool

    :returns: A function ``f(x)`` that evaluates the product of the
        specified attributes, each raised to the corresponding power,
        at the composition and thermodynamic conditions specified
        by `x`. The format of `x` is determined by `using_volume`,
        and is compatible with the `x_array` argument of
        :func:`burnman.optimize.nonlinear_fitting.confidence_prediction_bands`.
    :rtype: function
    """
    if type(attributes) is str:
        attributes = [attributes]
    if len(powers) == 0:
        powers = [1.0 for a in attributes]

    def f(x):
        if len(x) > 3:
            m.set_composition(x[:-3])
        if using_volume:
            m.set_state_with_volume(x[-3], x[-2])
        else:
            m.set_state(x[-3], x[-2])
        value = 1.0
        for a, p in zip(*[attributes, powers]):
            value *= np.power(getattr(m, a), p)
        return value

    return f


def extract_lines_between_markers(file, start_string, end_string, inclusive=False):
    """
    Extract lines from a file between two marker strings.

    :param file: Path to the input text file.
    :type file: str
    :param start_string: The marker string indicating where to start collecting lines.
    :type start_string: str
    :param end_string: The marker string indicating where to stop collecting lines.
    :type end_string: str
    :param inclusive: Whether to include the lines containing start_string and end_string, defaults to False.
    :type inclusive: bool, optional
    :return: A list of lines between the two markers.
    :rtype: list[str]
    """
    lines = []
    with open(file, encoding="latin-1") as f:
        inside = False
        for line in f:
            if start_string in line and not inside:
                inside = True
                if inclusive:
                    lines.append(line.strip())
                continue  # Don't reprocess the start line if not inclusive

            if inside:
                if end_string in line:
                    if inclusive:
                        lines.append(line.strip())
                    break
                lines.append(line.strip())

    return lines


def run_cli_program_with_input(program_path, stdin, verbose=True):
    """
    Run a command-line program with provided input and capture its output.

    :param program_path: Path to the CLI executable or a list of command-line arguments.
    :type program_path: str | list[str]
    :param stdin: The input string to pass to the program via standard input.
    :type stdin: str
    :param verbose: If True, prints the program's output to stdout, defaults to True.
    :type verbose: bool, optional
    :return: The standard output produced by the program.
    :rtype: str
    """
    process = subprocess.Popen(
        program_path,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdoutput, stderr = process.communicate(stdin)

    if verbose:
        print(stdoutput)
        if stderr:
            print("Error output:", stderr)

    return stdoutput
