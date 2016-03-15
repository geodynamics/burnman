# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import operator
import bisect
import os
import pkgutil
import numpy as np
from scipy.optimize import fsolve, curve_fit
from . import constants


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
            print("\t".join(r).replace("_", "\_"))
        return

    def col_width(table, colidx):
        return max([len(str(row[colidx])) for row in table])

    # create a format string with the first column left aligned, the others right
    # example:   {:<27}{:>11}{:>6}{:>8}
    frmt = "".join(
        [('{:<' if i == 0 else '{:>') + str(1 + col_width(table, i)) + '}' for i in range(len(table[0]))])
    for r in table:
        print(frmt.format(*r))


def sort_table(table, col=0):
    """
    Sort the table according to the column number
    """
    return sorted(table, key=operator.itemgetter(col))


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
    assert(x1 <= x)
    assert(x2 >= x)
    assert(x1 <= x2)

    alpha = (x - x1) / (x2 - x1)
    return (1. - alpha) * y1 + alpha * y2


def read_table(filename):
    datastream = pkgutil.get_data('burnman', 'data/' + filename)
    datalines = [line.strip()
                 for line in datastream.decode('ascii').split('\n') if line.strip()]
    table = []

    for line in datalines:
        if (line[0] != '#'):
            numbers = np.fromstring(line, sep=' ')
            table.append(numbers)
    return np.array(table)


def array_from_file(filename):
    """
    Generic function to read a file containing floats and commented lines
    into a 2D numpy array.

    Commented lines are prefixed by the characters # or %.
    """
    f = open(filename, 'r')
    data = []
    datastream = f.read()
    f.close()
    datalines = [line.strip().split()
                 for line in datastream.split('\n') if line.strip()]
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
    if (idx < 0):
        return table_y[0]
    elif (idx < len(table_x) - 1):
        return linear_interpol(x_value, table_x[idx], table_x[idx + 1],
                               table_y[idx], table_y[idx + 1])
    else:
        return table_y[idx]


def molar_volume_from_unit_cell_volume(unit_cell_v, z):
    """
    Converts a unit cell volume from Angstroms^3 per unitcell,
    to m^3/mol.

    Parameters
    ----------
    unit_cell_v : float
        Unit cell volumes [A^3/unit cell]

    z : float
        Number of formula units per unit cell


    Returns
    -------
    V : float
        Volume [m^3/mol]
    """
    V = unit_cell_v * constants.Avogadro / 1.e30 / z
    return V


def fit_PVT_data(mineral, fit_params, PT, V, V_sigma=None):
    """
    Given a mineral of any type, a list of fit parameters
    and a set of PTV points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    Parameters
    ----------
    mineral : mineral
        Mineral for which the parameters should be optimized

    fit_params : list of strings
        List of dictionary keys contained in mineral.params
        corresponding to the variables to be optimized
        during fitting. Initial guesses are taken from the existing
        values for the parameters

    PT : list of two lists of floats (or numpy arrays)
        The two lists contain a set of pressures [Pa]
        and temperatures [K] of the volume data points

    V : list (or numpy array) of floats
        Volumes [m^3/mol] of the mineral at the P,T condition
        given by parameter PT

    V_sigma : list (or numpy array) of floats
        Optional uncertainties on the volumes [m^3/mol]
        of the mineral at the P,T condition
        given by parameter PT

    Returns
    -------
    popt : numpy array of floats
        A list of optimized parameters

    pcov : 2D numpy array of floats
        The covariance matrix of the optimized parameters
    """
    def fit_data(PT, *params):
        for i, param in enumerate(fit_params):
            mineral.params[param] = params[i]
        volumes = []

        for P, T in zip(*PT):
            mineral.set_state(P, T)
            volumes.append(mineral.V)
        return volumes

    guesses = [mineral.params[param] for param in fit_params]
    popt, pcov = curve_fit(fit_data, PT, V, guesses, V_sigma)

    return popt, pcov


def equilibrium_pressure(minerals, stoichiometry, temperature, pressure_initial_guess=1.e5):
    """
    Given a list of minerals, their reaction stoichiometries
    and a temperature of interest, compute the
    equilibrium pressure of the reaction.

    Parameters
    ----------
    minerals : list of minerals
        List of minerals involved in the reaction.

    stoichiometry : list of floats
        Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol]

    temperature : float
        Temperature of interest [K]

    pressure_initial_guess : optional float
        Initial pressure guess [Pa]

    Returns
    -------
    pressure : float
        The equilibrium pressure of the reaction [Pa]
    """
    def eqm(P, T):
        gibbs = 0.
        for i, mineral in enumerate(minerals):
            mineral.set_state(P[0], T)
            gibbs = gibbs + mineral.gibbs * stoichiometry[i]
        return gibbs

    pressure = fsolve(eqm, [pressure_initial_guess], args=(temperature))[0]

    return pressure


def equilibrium_temperature(minerals, stoichiometry, pressure, temperature_initial_guess=1000.):
    """
    Given a list of minerals, their reaction stoichiometries
    and a pressure of interest, compute the
    equilibrium temperature of the reaction.

    Parameters
    ----------
    minerals : list of minerals
        List of minerals involved in the reaction.

    stoichiometry : list of floats
        Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol]

    pressure : float
        Pressure of interest [Pa]

    temperature_initial_guess : optional float
        Initial temperature guess [K]

    Returns
    -------
    temperature : float
        The equilibrium temperature of the reaction [K]
    """
    def eqm(T, P):
        gibbs = 0.
        for i, mineral in enumerate(minerals):
            mineral.set_state(P, T[0])
            gibbs = gibbs + mineral.gibbs * stoichiometry[i]
        return gibbs

    temperature = fsolve(eqm, [temperature_initial_guess], args=(pressure))[0]

    return temperature


def invariant_point(minerals_r1, stoichiometry_r1,
                    minerals_r2, stoichiometry_r2,
                    pressure_temperature_initial_guess=[1.e9, 1000.]):
    """
    Given a list of minerals, their reaction stoichiometries
    and a pressure of interest, compute the
    equilibrium temperature of the reaction.

    Parameters
    ----------
    minerals : list of minerals
        List of minerals involved in the reaction.

    stoichiometry : list of floats
        Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol]

    pressure : float
        Pressure of interest [Pa]

    temperature_initial_guess : optional float
        Initial temperature guess [K]

    Returns
    -------
    temperature : float
        The equilibrium temperature of the reaction [K]
    """
    def eqm(PT):
        P, T = PT
        gibbs_r1 = 0.
        for i, mineral in enumerate(minerals_r1):
            mineral.set_state(P, T)
            gibbs_r1 = gibbs_r1 + mineral.gibbs * stoichiometry_r1[i]
        gibbs_r2 = 0.
        for i, mineral in enumerate(minerals_r2):
            mineral.set_state(P, T)
            gibbs_r2 = gibbs_r2 + mineral.gibbs * stoichiometry_r2[i]
        return [gibbs_r1, gibbs_r2]

    pressure, temperature = fsolve(eqm, pressure_temperature_initial_guess)
    return pressure, temperature


def hugoniot(mineral, P_ref, T_ref, pressures, reference_mineral=None):
    """
    Calculates the temperatures (and volumes) along a Hugoniot
    as a function of pressure according to the Hugoniot equation
    U2-U1 = 0.5*(p2 - p1)(V1 - V2) where U and V are the
    internal energies and volumes (mass or molar) and U = F + TS


    Parameters
    ----------
    mineral : mineral
        Mineral for which the Hugoniot is to be calculated.

    P_ref : float
        Reference pressure [Pa]

    T_ref : float
        Reference temperature [K]

    pressures : numpy array of floats
        Set of pressures [Pa] for which the Hugoniot temperature
        and volume should be calculated

    reference_mineral : mineral
        Mineral which is stable at the reference conditions
        Provides an alternative U_0 and V_0 when the reference
        mineral transforms to the mineral of interest at some
        (unspecified) pressure.

    Returns
    -------
    temperatures : numpy array of floats
        The Hugoniot temperatures at pressure

    volumes : numpy array of floats
        The Hugoniot volumes at pressure
    """

    def Ediff(T, mineral, P, P_ref, U_ref, V_ref):
        mineral.set_state(P, T[0])
        U = mineral.helmholtz + T[0] * mineral.S
        V = mineral.V

        return (U - U_ref) - 0.5 * (P - P_ref) * (V_ref - V)

    if reference_mineral is None:
        reference_mineral = mineral

    reference_mineral.set_state(P_ref, T_ref)
    U_ref = reference_mineral.helmholtz + T_ref * reference_mineral.S
    V_ref = reference_mineral.V

    temperatures = np.empty_like(pressures)
    volumes = np.empty_like(pressures)

    for i, P in enumerate(pressures):
        temperatures[i] = fsolve(
            Ediff, [T_ref], args=(mineral, P, P_ref, U_ref, V_ref))[0]
        volumes[i] = mineral.V

    return temperatures, volumes


def convert_fractions(composite, phase_fractions, input_type, output_type):
    """
    Takes a composite with a set of user defined molar, volume
    or mass fractions (which do not have to be the fractions
    currently associated with the composite) and
    converts the fractions to molar, mass or volume.

    Conversions to and from mass require a molar mass to be
    defined for all phases. Conversions to and from volume
    require set_state to have been called for the composite.

    Parameters
    ----------
    composite : Composite
        Composite for which fractions are to be defined.

    phase_fractions : list of floats
        List of input phase fractions (of type input_type)

    input_type : string
        Input fraction type: 'molar', 'mass' or 'volume'

    output_type : string
        Output fraction type: 'molar', 'mass' or 'volume'

    Returns
    -------
    output_fractions : list of floats
        List of output phase fractions (of type output_type)
    """
    if input_type == 'volume' or output_type == 'volume':
        if composite.temperature == None:
            raise Exception(
                composite.to_string() + ".set_state(P, T) has not been called, so volume fractions are currently undefined. Exiting.")

    if input_type == 'molar':
        molar_fractions = phase_fractions
    if input_type == 'volume':
        total_moles = sum(
            volume_fraction / phase.molar_volume for volume_fraction,
            phase in zip(phase_fractions, composite.phases))
        molar_fractions = [volume_fraction / (phase.molar_volume * total_moles)
                           for volume_fraction, phase in zip(phase_fractions, composite.phases)]
    if input_type == 'mass':
        total_moles = sum(mass_fraction / phase.molar_mass for mass_fraction,
                          phase in zip(phase_fractions, composite.phases))
        molar_fractions = [mass_fraction / (phase.molar_mass * total_moles)
                           for mass_fraction, phase in zip(phase_fractions, composite.phases)]

    if output_type == 'volume':
        total_volume = sum(
            molar_fraction * phase.molar_volume for molar_fraction,
            phase in zip(molar_fractions, composite.phases))
        output_fractions = [molar_fraction * phase.molar_volume /
                            total_volume for molar_fraction, phase in zip(molar_fractions, composite.phases)]
    elif output_type == 'mass':
        total_mass = sum(molar_fraction * phase.molar_mass for molar_fraction,
                         phase in zip(molar_fractions, composite.phases))
        output_fractions = [molar_fraction * phase.molar_mass /
                            total_mass for molar_fraction, phase in zip(molar_fractions, composite.phases)]
    elif output_type == 'molar':
        output_fractions = molar_fractions

    return output_fractions


def bracket(fn, x0, dx, args=(), ratio=1.618, maxiter=100):
    """
    Given a function and a starting guess, find two
    inputs for the function that bracket a root.

    Parameters
    ----------
    fn : function
        The function to bracket
    x0 : float
        The starting guess
    dx : float
        Small step for starting the search
    args : parameter list
        Additional arguments to give to fn
    ratio :
        The step size increases by this ratio
        every step in the search. Defaults to
        the golden ratio.
    maxiter : int
        The maximum number of steps before giving up.

    Returns
    -------
    xa, xb, fa, fb: floats
        xa and xb are the inputs which bracket a root of fn.
        fa and fb are the values of the function at those points.
        If the bracket function takes more than maxiter steps,
        it raises a ValueError.
    """
    niter = 0
    dx = np.abs(dx)
    assert(ratio > 1.0)

    # Get the starting positions
    f0 = fn(x0, *args)
    x_left = x0 - dx
    x_right = x0 + dx
    f_left = fn(x_left, *args)
    f_right = fn(x_right, *args)

    # Overshot zero, try making dx smaller
    if (f0 - f_left) * (f_right - f0) < 0.:
        while (f0 - f_left) * (f_right - f0) < 0. and dx > np.finfo('float').eps and niter < maxiter:
            dx /= ratio
            x_left = x0 - dx
            x_right = x0 + dx
            f_left = fn(x_left, *args)
            f_right = fn(x_right, *args)
            niter += 1
        if niter == maxiter:  # Couldn't find something with same slope in both directions
            raise ValueError('Cannot find zero.')

    niter = 0
    slope = f_right - f0
    if slope > 0. and f0 > 0.:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left
    elif slope > 0. and f0 < 0.:  # Walk right
        x1 = x_right
        f1 = f_right
    elif slope < 0. and f0 > 0:  # Walk right
        x1 = x_right
        f1 = f_right
    else:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left

    # Do the walking
    while f0 * f1 > 0. and niter < maxiter:
        dx *= ratio
        xnew = x1 + dx
        fnew = fn(xnew, *args)
        x0 = x1
        f0 = f1
        x1 = xnew
        f1 = fnew
        niter += 1

    if f0 * f1 > 0.:
        raise ValueError('Cannot find zero.')
    else:
        return x0, x1, f0, f1
