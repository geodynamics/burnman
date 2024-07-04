# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
from sympy import Matrix, Rational
import scipy.integrate as integrate
from collections import Counter
import itertools

from .reductions import row_reduce


def round_to_n(x, xerr, n):
    return round(x, -int(np.floor(np.log10(np.abs(xerr)))) + (n - 1))


def unit_normalize(a, order=2, axis=-1):
    """
    Calculates the L2 normalized array of numpy array a
    of a given order and along a given axis.
    """
    l2 = np.atleast_1d(np.apply_along_axis(np.linalg.norm, axis, a, order))

    l2[l2 == 0] = 1
    return a / np.expand_dims(l2, axis)[0][0]


def float_eq(a, b):
    """
    Test if two floats are almost equal to each other
    """
    return abs(a - b) < 1e-10 * max(1e-5, abs(a), abs(b))


def linear_interpol(x, x1, x2, y1, y2):
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    assert x1 <= x
    assert x2 >= x
    assert x1 <= x2

    alpha = (x - x1) / (x2 - x1)
    return (1.0 - alpha) * y1 + alpha * y2


def bracket(fn, x0, dx, args=(), ratio=1.618, maxiter=100):
    """
    Given a function and a starting guess, find two
    inputs for the function that bracket a root.

    :param fn: The function to bracket.
    :type fn: function
    :param x0: The starting guess.
    :type x0: float
    :param dx: Small step for starting the search.
    :type dx: float
    :param args: Additional arguments to give to fn.
    :type args: tuple
    :param ratio: The step size increases by this ratio
        every step in the search. Defaults to
        the golden ratio.
    :type ratio: float
    :param maxiter: The maximum number of steps before giving up.
    :type maxiter: int

    :returns: xa, xb, fa, fb. xa and xb are the inputs which
        bracket a root of fn. fa and fb are the values of the
        function at those points.
        If the bracket function takes more than maxiter steps,
        it raises a ValueError.
    :rtype: tuple of floats

    """
    niter = 0
    dx = np.abs(dx)
    assert ratio > 1.0

    # Get the starting positions
    f0 = fn(x0, *args)
    x_left = x0 - dx
    x_right = x0 + dx
    f_left = fn(x_left, *args)
    f_right = fn(x_right, *args)

    # Overshot zero, try making dx smaller
    if (f0 - f_left) * (f_right - f0) < 0.0:
        while (
            (f0 - f_left) * (f_right - f0) < 0.0
            and dx > np.finfo("float").eps
            and niter < maxiter
        ):
            dx /= ratio
            x_left = x0 - dx
            x_right = x0 + dx
            f_left = fn(x_left, *args)
            f_right = fn(x_right, *args)
            niter += 1
        if (
            niter == maxiter
        ):  # Couldn't find something with same slope in both directions
            raise ValueError("Cannot find zero.")

    niter = 0
    slope = f_right - f0
    if slope > 0.0 and f0 > 0.0:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left
    elif slope > 0.0 and f0 < 0.0:  # Walk right
        x1 = x_right
        f1 = f_right
    elif slope < 0.0 and f0 > 0:  # Walk right
        x1 = x_right
        f1 = f_right
    else:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left

    # Do the walking
    while f0 * f1 > 0.0 and niter < maxiter:
        dx *= ratio
        xnew = x1 + dx
        fnew = fn(xnew, *args)
        x0 = x1
        f0 = f1
        x1 = xnew
        f1 = fnew
        niter += 1

    if f0 * f1 > 0.0:
        raise ValueError("Cannot find zero.")
    else:
        return x0, x1, f0, f1


def _pad_ndarray_inverse_mirror(array, padding):
    """
    Pads an ndarray according to an inverse mirror
    scheme. For example, for a 1D array
    [2, 4, 6, 7, 8] padded by 3 cells, we have:

     padding  |  original array |  padding

    -3 -2  0  |  2  4  6  7  8  |  9 10 12

    :param array: The array to be padded
    :type array: numpy.ndarray

    :param padding: The number of elements with which to pad the
        array in each dimension.
    :type padding: tuple

    :returns: The padded array.
    :type: numpy.ndarray
    """
    padded_shape = [n + 2 * padding[i] for i, n in enumerate(array.shape)]
    padded_array = np.zeros(padded_shape)

    slices = tuple(
        [slice(padding[i], padding[i] + l) for i, l in enumerate(array.shape)]
    )
    padded_array[slices] = array

    padded_array_indices = list(
        itertools.product(
            *[range(n + 2 * padding[i]) for i, n in enumerate(array.shape)]
        )
    )
    inserted_indices = list(
        itertools.product(
            *[range(padding[i], padding[i] + l) for i, l in enumerate(array.shape)]
        )
    )
    padded_array_indices.extend(inserted_indices)

    counter = Counter(padded_array_indices)
    keys = list(counter.keys())
    padded_indices = [keys[i] for i, value in enumerate(counter.values()) if value == 1]
    edge_indices = tuple(
        [
            tuple(
                [
                    np.min(
                        [
                            np.max([axis_idx, padding[dimension]]),
                            padded_array.shape[dimension] - padding[dimension] - 1,
                        ]
                    )
                    for dimension, axis_idx in enumerate(idx)
                ]
            )
            for idx in padded_indices
        ]
    )
    mirror_indices = tuple(
        [
            tuple(
                [
                    2 * edge_indices[i][j] - padded_indices[i][j]
                    for j in range(len(array.shape))
                ]
            )
            for i in range(len(padded_indices))
        ]
    )

    for i, idx in enumerate(padded_indices):
        padded_array[idx] = (
            2.0 * padded_array[edge_indices[i]] - padded_array[mirror_indices[i]]
        )

    return padded_array


def smooth_array(
    array, grid_spacing, gaussian_rms_widths, truncate=4.0, mode="inverse_mirror"
):
    """
    Creates a smoothed array by convolving it with a gaussian filter.
    Grid resolutions and gaussian RMS widths are required for each of
    the axes of the numpy array. The smoothing is truncated at a
    user-defined number of standard deviations. The edges of the array
    can be padded in a number of different ways given by the
    'mode' parameter.

    :param array: The array to smooth.
    :type array: numpy.ndarray
    :param grid_spacing: The spacing of points along each axis.
    :type grid_spacing: numpy.array of floats
    :param gaussian_rms_widths: The Gaussian RMS widths/standard deviations for the
        Gaussian convolution.
    :type gaussian_rms_widths: numpy.array of floats
    :param truncate: The number of standard deviations at which to truncate
        the smoothing.
    :type truncate: float
    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap', 'inverse_mirror'}
        The mode parameter determines how the array borders are handled
        either by scipy.ndimage.filters.gaussian_filter.
        Default is 'inverse_mirror', which uses
        :func:`burnman.tools.math._pad_ndarray_inverse_mirror`.
    :type mode: str

    :returns: The smoothed array
    :rtype: numpy.ndarray
    """

    # gaussian_filter works with standard deviations normalised to
    # the grid spacing.
    sigma = tuple(np.array(gaussian_rms_widths) / np.array(grid_spacing))

    if mode == "inverse_mirror":
        padding = tuple([int(np.ceil(truncate * s)) for s in sigma])
        padded_array = _pad_ndarray_inverse_mirror(array, padding)
        smoothed_padded_array = gaussian_filter(padded_array, sigma=sigma)
        slices = tuple(
            [slice(padding[i], padding[i] + l) for i, l in enumerate(array.shape)]
        )
        smoothed_array = smoothed_padded_array[slices]
    else:
        smoothed_array = gaussian_filter(array, sigma=sigma, mode=mode)

    return smoothed_array


def interp_smoothed_array_and_derivatives(
    array,
    x_values,
    y_values,
    x_stdev=0.0,
    y_stdev=0.0,
    truncate=4.0,
    mode="inverse_mirror",
    indexing="xy",
):
    """
    Creates a smoothed array on a regular 2D grid. Smoothing
    is achieved using :func:`burnman.tools.math.smooth_array`.
    Outputs scipy.interpolate.RegularGridInterpolator() interpolators
    which can be used to query the array, or its derivatives in the
    x- and y- directions.

    :param array: The array to smooth. Each element array[i][j]
        corresponds to the position x_values[i], y_values[j]
    :type array: numpy.array (2D)
    :param x_values: The gridded x values over which to create the smoothed grid
    :type x_values: numpy.array (1D)
    :param y_values: The gridded y_values over which to create the smoothed grid
    :type y_values: numpy.array (1D)
    :param x_stdev: The standard deviation for the Gaussian filter along the x axis
    :type x_stdev: float
    :param y_stdev: The standard deviation for the Gaussian filter along the y axis
    :type y_stdev: float
    :param truncate: The number of standard deviations at which to truncate
        the smoothing.
    :type truncate: float
    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap', 'inverse_mirror'}
        The mode parameter determines how the array borders are handled
        either by scipy.ndimage.filters.gaussian_filter.
        Default is 'inverse_mirror', which uses
        :func:`burnman.tools.math._pad_ndarray_inverse_mirror`.
    :type mode: str
    :param indexing: {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') indexing of output.
        See numpy.meshgrid for more details.
    :type indexing: str

    :returns: Three RegularGridInterpolator functors
        interpolation functions for the smoothed property and
        the first derivatives with respect to x and y.
    :rtype: tuple
    """

    dx = x_values[1] - x_values[0]
    dy = y_values[1] - y_values[0]

    if indexing == "xy":
        smoothed_array = smooth_array(
            array=array,
            grid_spacing=np.array([dy, dx]),
            gaussian_rms_widths=np.array([y_stdev, x_stdev]),
            truncate=truncate,
            mode=mode,
        ).T

    elif indexing == "ij":
        smoothed_array = smooth_array(
            array=array,
            grid_spacing=np.array([dx, dy]),
            gaussian_rms_widths=np.array([x_stdev, y_stdev]),
            truncate=truncate,
            mode=mode,
        )

    else:
        raise Exception("Indexing scheme not recognised. Should be ij or xy.")

    dSAdxdx, dSAdydy = np.gradient(smoothed_array)

    interps = (
        RegularGridInterpolator((x_values, y_values), smoothed_array, method="linear"),
        RegularGridInterpolator((x_values, y_values), dSAdxdx / dx, method="linear"),
        RegularGridInterpolator((x_values, y_values), dSAdydy / dy, method="linear"),
    )

    return interps


def compare_l2(depth, calc, obs):
    """
    Computes the L2 norm for N profiles at a time (assumed to be linear between points).

    :type depths: array of float
    :param depths: depths. :math:`[m]`
    :type calc: list of arrays of float
    :param calc: N arrays calculated values, e.g. [mat_vs,mat_vphi]
    :type obs: list of arrays of float
    :param obs: N arrays of values (observed or calculated) to compare to , e.g. [seis_vs, seis_vphi]

    :returns: array of L2 norms of length N
    :rtype: array of floats
    """
    err = []
    for i in range(len(calc)):
        err.append(l2(depth, calc[i], obs[i]))

    return err


def compare_chifactor(calc, obs):
    """
    Computes the chi factor for N profiles at a time. Assumes a 1% a priori uncertainty on the seismic model.


    :type calc: list of arrays of float
    :param calc: N arrays calculated values, e.g. [mat_vs,mat_vphi]
    :type obs: list of arrays of float
    :param obs: N arrays of values (observed or calculated) to compare to , e.g. [seis_vs, seis_vphi]

    :returns: error array of length N
    :rtype: array of floats
    """
    err = []
    for i in range(len(calc)):
        err.append(chi_factor(calc[i], obs[i]))

    return err


def l2(x, funca, funcb):
    """
    Computes the L2 norm for one profile(assumed to be linear between points).

    :type x: array of float
    :param x: depths :math:`[m]`.
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: L2 norm
    :rtype: array of floats
    """
    diff = np.array(funca - funcb)
    diff = diff * diff

    return integrate.trapezoid(diff, x)


def nrmse(x, funca, funcb):
    """
    Normalized root mean square error for one profile
    :type x: array of float
    :param x: depths in m.
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: RMS error
    :rtype: array of floats
    """
    diff = np.array(funca - funcb)
    diff = diff * diff
    rmse = np.sqrt(np.sum(diff) / x)
    nrmse = rmse / (np.max(funca) - np.min(funca))

    return nrmse


def chi_factor(calc, obs):
    """
    :math:`\\chi` factor for one profile assuming 1% uncertainty on the reference model (obs)
    :type calc: list of arrays of float
    :param calc: array calculated values
    :type obs: list of arrays of float
    :param obs: array of reference values to compare to

    :returns: :math:`\\chi` factor
    :rtype: array of floats
    """

    err = np.empty_like(calc)
    for i in range(len(calc)):
        err[i] = pow((calc[i] - obs[i]) / (0.01 * np.mean(obs)), 2.0)

    err_tot = np.sum(err) / len(err)

    return err_tot


def independent_row_indices(array):
    """
    Returns the indices corresponding to an independent set of rows
    for a given array. The independent rows are determined from the pivots
    used during row reduction/Gaussian elimination.

    :param array: The input array.
    :type array: 2D numpy.array of floats

    :returns: The indices corresponding to a set of independent rows
        of the input array.
    :rtype: 1D numpy array of integers
    """
    m = Matrix(
        array.shape[0],
        array.shape[1],
        lambda i, j: Rational(array[i, j]).limit_denominator(1000),
    )
    _, pivots, swaps = row_reduce(
        m,
        iszerofunc=lambda x: x.is_zero,
        simpfunc=lambda x: Rational(x).limit_denominator(1000),
    )
    indices = np.array(range(len(array)))
    for swap in np.array(swaps):
        indices[swap] = indices[swap[::-1]]
    return sorted(indices[: len(pivots)])


def array_to_rational_matrix(array):
    """
    Converts a numpy array into a sympy matrix
    filled with rationals
    """
    return Matrix([[Rational(v).limit_denominator(1000) for v in row] for row in array])


def complete_basis(basis, flip_columns=True):
    """
    Creates a full basis by filling remaining rows with
    selected rows of the identity matrix that have indices
    not in the column pivot list of the basis RREF

    :param basis: A partial basis
    :type basis: numpy array
    :param flip_columns: whether to calculate the RREF
        from a flipped basis. This can be useful, for example,
        when dependent columns are known to be further to the
        right in the basis. defaults to True
    :type flip_columns: bool, optional
    :return: basis
    :rtype: numpy array
    """

    n, m = basis.shape
    if n < m:
        if flip_columns:
            pivots = list(m - 1 - np.array(Matrix(np.flip(basis, axis=1)).rref()[1]))
        else:
            pivots = list(Matrix(basis)).rref()[1]

        return np.concatenate(
            (basis, np.identity(m)[[i for i in range(m) if i not in pivots], :]), axis=0
        )
    else:
        return basis


def generate_complete_basis(incomplete_basis, array):
    """
    Given a 2D array with independent rows and a second 2D array that spans a
    larger space, creates a complete basis for the combined array using all
    the rows of the first array, followed by any required rows of the
    second array. So, for example, if the first array is:
    [[1, 0, 0], [1, 1, 0]] and the second array is:
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]], the complete basis will be:
    [[1, 0, 0], [1, 1, 0], [0, 0, 1]].

    :param incomplete_basis: An array containing the basis to be completed.
    :type incomplete_basis: 2D numpy.array

    :param array: An array spanning the full space for which a basis is required.
    :type array: 2D numpy.array

    :returns: An array containing the basis vectors spanning both of the
        input arrays.
    :rtype: 2D numpy array
    """

    incomplete_rank = array_to_rational_matrix(incomplete_basis).rank()
    if incomplete_rank < len(incomplete_basis):
        raise Exception(
            "The incomplete basis is rank-deficient. " "Remove one or more endmembers."
        )

    a = np.concatenate((incomplete_basis, array))
    complete_basis = np.array(a[independent_row_indices(a)], dtype=float)

    # Store the rank of the matrix for later comparison
    len_basis = array_to_rational_matrix(complete_basis).rank()

    # This next step ensures that all of the original
    # rows are contained in the new basis in their original order
    c = np.linalg.lstsq(
        np.array(complete_basis).astype(float).T,
        np.array(incomplete_basis).astype(float).T,
        rcond=None,
    )[0].T

    for row in np.eye(len(c[0])):
        old_rank = array_to_rational_matrix(c).rank()
        c2 = np.concatenate((c, [row]))
        new_rank = array_to_rational_matrix(c2).rank()

        if new_rank > old_rank:
            c = c2

    complete_basis = c.dot(complete_basis)

    # Check that the matrix rank has not changed
    if len_basis != array_to_rational_matrix(complete_basis).rank():
        raise Exception(
            "Basis length changed during conversion. " "Report this bug to developers."
        )

    return complete_basis.round(decimals=12) + 0.0


def is_positive_definite(matrix):
    """
    Checks if a matrix is positive definite

    :param matrix: Input matrix
    :type matrix: 2D numpy array

    :return: Whether or not the matrix is positive definite
    :rtype: bool
    """
    return np.all(np.linalg.eigvals(matrix) > 0)
