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
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp2d
from collections import Counter
import itertools

from . import constants
from . import nonlinear_fitting
import itertools

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

def flatten(l): return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list or type(l) is np.ndarray else [l]

def round_to_n(x, xerr, n):
    return round(x, -int(np.floor(np.log10(np.abs(xerr)))) + (n - 1))

def pretty_print_values(popt, pcov, params):
    """
    Takes a numpy array of parameters, the corresponding covariance matrix 
    and a set of parameter names and prints the parameters and 
    principal 1-s.d.uncertainties (np.sqrt(pcov[i][i])) 
    in a nice text based format. 
    """
    for i, p in enumerate(params):
        p_rnd = round_to_n(popt[i], np.sqrt(pcov[i][i]), 1)
        c_rnd = round_to_n(np.sqrt(pcov[i][i]), np.sqrt(pcov[i][i]), 1)

        if p_rnd != 0.:
            p_expnt = np.floor(np.log10(np.abs(p_rnd)))
        else:
            p_expnt = 0.
            
        scale = np.power(10., p_expnt)
        nd = p_expnt - np.floor(np.log10(np.abs(c_rnd)))
        print ('{0:s}: ({1:{4}{5}f} +/- {2:{4}{5}f}) x {3:.0e}'.format(p, p_rnd/scale, c_rnd/scale, scale, 0, (nd)/10.))

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
        
def pretty_plot():
    """
    Makes pretty plots. Overwrites the matplotlib default settings to allow
    for better fonts. Slows down plotting
    """
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\\usepackage{relsize}'
    plt.rc('font', family='sanserif')

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


def fit_PTp_data(mineral, fit_params, p_flags, PTp, PTp_covariances=[], mle_tolerances=[], max_lm_iterations=50, verbose=True):
    """
    Given a mineral of any type, a list of fit parameters
    and a set of P-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    Parameters
    ----------
    mineral : mineral
        Mineral for which the parameters should be optimized

    observed_property : string
        Attribute name for the property to be fit (e.g. 'V')

    fit_params : list of strings
        List of dictionary keys contained in mineral.params
        corresponding to the variables to be optimized
        during fitting. Initial guesses are taken from the existing
        values for the parameters

    PTp : numpy array of observed P-T-property values

    PTp_covariances : numpy array of P-T-property covariances (optional)
        If not given, all covariance matrices are chosen 
        such that C00 = 1, otherwise Cij = 0. 
        In other words, all data points have equal weight, 
        with all error in the pressure.

    Returns
    -------
    model : instance of fitted model
        Fitting-related attributes are as follows:
            n_dof : integer
                Degrees of freedom of the system
            data_mle : 2D numpy array
                Maximum likelihood estimates of the observed data points 
                on the best-fit curve
            jacobian : 2D numpy array
                d(weighted_residuals)/d(parameter)
            weighted_residuals : numpy array
                Weighted residuals
            weights : numpy array
                1/(data variances normal to the best fit curve)
            WSS : float
                Weighted sum of squares residuals
            popt : numpy array
                Optimized parameters
            pcov : 2D numpy array
                Covariance matrix of optimized parameters
            noise_variance : float
                Estimate of the variance of the data normal to the curve
    """

    class Model(object):
        def __init__(self, mineral, data, data_covariance, flags, fit_params, guessed_params, delta_params, mle_tolerances):
            self.m = mineral
            self.data = data
            self.data_covariance = data_covariance
            self.flags = flags
            self.fit_params = fit_params
            self.set_params(guessed_params)
            self.delta_params = delta_params
            self.mle_tolerances = mle_tolerances
            
        def set_params(self, param_values):
            i=0
            for param in self.fit_params:
                if isinstance(self.m.params[param], float):
                    self.m.params[param] = param_values[i]
                    i += 1
                else:
                    for j in range(len(self.m.params[param])):
                        self.m.params[param][j] = param_values[i]
                        i += 1

        def get_params(self):
            params = []
            for i, param in enumerate(self.fit_params):
                params.append(self.m.params[param])
            return np.array(flatten([mineral.params[prm] for prm in fit_params]))

        def function(self, x, p_flag):
            P, T, p = x
            self.m.set_state(P, T)
            return np.array([P, T, getattr(self.m, p_flag)])

        def normal(self, x, p_flag):
            P, T, p = x
            
            if p_flag == 'V':
                self.m.set_state(P, T)
                n = np.array([-1./self.m.K_T, self.m.alpha, -1./self.m.V])
            else:
                dP = 1.e5
                dT = 1.
                
                dPdp = (2.*dP)/(self.function([P+dP, T, 0.], p_flag)[2] - self.function([P-dP, T, 0.], p_flag)[2])
                dpdT = (self.function([P, T+dT, 0.], p_flag)[2] - self.function([P, T-dT, 0.], p_flag)[2])/(2.*dT)
                dPdT = -dPdp*dpdT
                n = np.array([-1., dPdT, dPdp])
                
            return n/np.linalg.norm(n)

        
    # If only one property flag is given, assume it applies to all data
    if type(p_flags) is str:
        p_flags = np.array([p_flags] * len(PTp[:,0]))

    # Apply mle tolerances if they dont exist
    if mle_tolerances == []:
        mineral.set_state(1.e5, 300.)
        mle_tolerance_factor = 1.e-5
        mle_tolerances = np.empty(len(p_flags))
        for i, p_flag in enumerate(p_flags):
            if p_flag in ['gibbs', 'enthalpy', 'H', 'helmholtz']:
                mle_tolerances[i] = 1. # 1 J
            else:
                mle_tolerances[i] = mle_tolerance_factor*getattr(mineral, p_flag)
        
    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and p)
    if PTp_covariances == []:
        PTp_covariances = np.zeros((len(PTp[:,0]), len(PTp[0]), len(PTp[0])))
        for i in range(len(PTp_covariances)):
            PTp_covariances[i][0][0] = 1.
    
    guessed_params = np.array(flatten([mineral.params[prm] for prm in fit_params]))
    model = Model(mineral = mineral,
                  data = PTp,
                  data_covariance = PTp_covariances,
                  flags = p_flags,
                  fit_params = fit_params,
                  guessed_params = guessed_params,
                  delta_params = guessed_params*1.e-5,
                  mle_tolerances = mle_tolerances)

    nonlinear_fitting.nonlinear_least_squares_fit(model, max_lm_iterations = max_lm_iterations, verbose=verbose)

    return model


def fit_PTV_data(mineral, fit_params, PTV, PTV_covariances=[], max_lm_iterations=50, verbose=True):
    """
    A simple alias for the fit_PTp_data for when all the data is volume data
    """
        
    return fit_PTp_data(mineral=mineral, p_flags='V',
                        PTp=PTV, PTp_covariances=PTV_covariances, 
                        fit_params=fit_params, max_lm_iterations=max_lm_iterations, verbose=verbose)


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

def check_eos_consistency(m, P=1.e9, T=300., tol=0.01, verbose=False):
    """
    Compute numerical derivatives of the gibbs free energy of a mineral
    under given conditions, and check these values against those provided 
    analytically by the equation of state

    Parameters
    ----------
    m : mineral
        The mineral for which the equation of state 
        is to be checked for consistency
    P : float
        The pressure at which to check consistency
    T : float
        The temperature at which to check consistency
    tol : float
        The fractional tolerance for each of the checks
    verbose : boolean
        Decide whether to print information about each
        check

    Returns
    -------
    consistency: boolean
        If all checks pass, returns True

    """
    dT = 1.
    dP = 1.
    
    m.set_state(P, T)
    G0 = m.gibbs
    S0 = m.S
    V0 = m.V

    expr = ['G = F + PV', 'G = H - TS', 'G = E - TS + PV']
    eq = [[m.gibbs, (m.helmholtz + P*m.V)],
          [m.gibbs, (m.H - T*m.S)],
          [m.gibbs, (m.internal_energy - T*m.S + P*m.V)]]
    
    m.set_state(P, T + dT)
    G1 = m.gibbs
    S1 = m.S
    V1 = m.V
    
    
    m.set_state(P + dP, T)
    G2 = m.gibbs
    V2 = m.V
    
    # T derivatives
    m.set_state(P, T + 0.5*dT)
    expr.extend(['S = dG/dT', 'alpha = 1/V dV/dT', 'C_p = T dS/dT'])
    eq.extend([[m.S, -(G1 - G0)/dT],
               [m.alpha, (V1 - V0)/dT/m.V],
               [m.heat_capacity_p, (T + 0.5*dT)*(S1 - S0)/dT]])
    
    # P derivatives
    m.set_state(P + 0.5*dP, T)
    expr.extend(['V = dG/dP', 'K_T = -V dP/dV'])
    eq.extend([[m.V, (G2 - G0)/dP],
               [m.K_T, -0.5*(V2 + V0)*dP/(V2 - V0)]])

    consistencies = [np.abs(e[0]/e[1] - 1) < tol for e in eq]
    consistency = np.all(consistencies)
    
    if verbose == True:
        print('Checking EoS consistency for {0:s}'.format(m.to_string()))
        print('Expressions within tolerance of {0:2f}'.format(tol))
        for i, c in enumerate(consistencies):
            print('{0:10s} : {1:5s}'.format(expr[i], str(c)))
        if consistency == True:
            print('All EoS consistency constraints satisfied for {0:s}'.format(m.to_string()))
        else:
            print('Not satisfied all EoS consistency constraints for {0:s}'.format(m.to_string()))
            
    return consistency

def _pad_ndarray_inverse_mirror(array, padding):
    """
    Pads an ndarray according to an inverse mirror
    scheme. For example, for a 1D array 
    [2, 4, 6, 7, 8] padded by 3 cells, we have:

     padding  |  original array |  padding

    -3 -2  0  |  2  4  6  7  8  |  9 10 12

    Parameters
    ----------
    array : numpy ndarray
        The array to be padded
    padding : tuple
        The number of elements with which to pad the
        array in each dimension.

    Returns
    -------
    padded_array: numpy ndarray
        The padded array

    """
    padded_shape = [n + 2*padding[i] for i, n in enumerate(array.shape)]
    padded_array = np.zeros(padded_shape)

    slices = tuple([ slice(padding[i], padding[i] + l) for i, l in enumerate(array.shape)])
    padded_array[slices] = array

    padded_array_indices = list(itertools.product(*[range(n + 2*padding[i]) for i, n in enumerate(array.shape)]))
    inserted_indices = list(itertools.product(*[range(padding[i], padding[i] + l) for i, l in enumerate(array.shape)]))
    padded_array_indices.extend(inserted_indices)

    counter = Counter(padded_array_indices)
    keys = list(counter.keys())
    padded_indices = [keys[i] for i, value in enumerate(counter.values()) if value == 1]
    edge_indices = tuple([tuple([np.min([np.max([axis_idx, padding[dimension]]), padded_array.shape[dimension] - padding[dimension] - 1])
                                 for dimension, axis_idx in enumerate(idx)]) for idx in padded_indices])
    mirror_indices = tuple([tuple([2*edge_indices[i][j] - padded_indices[i][j] for j in range(len(array.shape))]) for i in range(len(padded_indices))])

    for i, idx in enumerate(padded_indices):
        padded_array[idx] = 2.*padded_array[edge_indices[i]] - padded_array[mirror_indices[i]]

    return padded_array


def smooth_array(array, grid_spacing,
                 gaussian_rms_widths, truncate=4.0,
                 mode='inverse_mirror'):
    """
    Creates a smoothed array by convolving it with a gaussian filter. 
    Grid resolutions and gaussian RMS widths are required for each of
    the axes of the numpy array. The smoothing is truncated at a 
    user-defined number of standard deviations. The edges of the array
    can be padded in a number of different ways given by the
    'mode' parameter.

    Parameters
    ----------
    array : numpy ndarray
        The array to smooth
    grid_spacing : numpy array of floats
        The spacing of points along each axis
    gaussian_rms_widths : numpy array of floats
        The Gaussian RMS widths/standard deviations for the 
        Gaussian convolution.
    truncate : float (default=4.) 
        The number of standard deviations at which to truncate 
        the smoothing.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap', 'inverse_mirror'}
        The mode parameter determines how the array borders are handled
        either by scipy.ndimage.filters.gaussian_filter.
        Default is 'inverse_mirror', which uses
        burnman.tools._pad_ndarray_inverse_mirror().

    Returns
    -------
    smoothed_array: numpy ndarray
       The smoothed array

    """

    # gaussian_filter works with standard deviations normalised to
    # the grid spacing.
    sigma = tuple(np.array(gaussian_rms_widths)/np.array(grid_spacing))
    
    if mode == 'inverse_mirror':
        padding = tuple([int(np.ceil(truncate*s)) for s in sigma])
        padded_array = _pad_ndarray_inverse_mirror(array, padding)
        smoothed_padded_array = gaussian_filter(padded_array,
                                                sigma=sigma)
        slices = tuple([ slice(padding[i], padding[i] + l) for i, l in enumerate(array.shape)])
        smoothed_array = smoothed_padded_array[slices]
    else:
        smoothed_array = gaussian_filter(array, sigma=sigma, mode=mode)
        
    return smoothed_array


def interp_smoothed_array_and_derivatives(array,
                                          x_values, y_values,
                                          x_stdev=0., y_stdev=0.,
                                          truncate=4.,
                                          mode='inverse_mirror',
                                          indexing='xy'):
    """
    Creates a smoothed array on a regular 2D grid. Smoothing 
    is achieved using burnman.tools.smooth_array(). 
    Outputs scipy.interpolate.interp2d() interpolators 
    which can be used to query the array, or its derivatives in the 
    x- and y- directions.

    Parameters
    ----------
    array : 2D numpy array
        The array to smooth. Each element array[i][j]
        corresponds to the position x_values[i], y_values[j]
    x_values : 1D numpy array
        The gridded x values over which to create the smoothed grid
    y_values : 1D numpy array
        The gridded y_values over which to create the smoothed grid
    x_stdev : float
        The standard deviation for the Gaussian filter along the x axis
    y_stdev : float
        The standard deviation for the Gaussian filter along the x axis
    truncate : float (optional) 
        The number of standard deviations at which to truncate 
        the smoothing (default = 4.).
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap', 'inverse_mirror'}
        The mode parameter determines how the array borders are handled
        either by scipy.ndimage.filters.gaussian_filter.
        Default is 'inverse_mirror', which uses
        burnman.tools._pad_ndarray_inverse_mirror().
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') indexing of output.
        See numpy.meshgrid for more details.

    Returns
    -------
    interps: tuple of three interp2d functors
        interpolation functions for the smoothed property and 
        the first derivatives with respect to x and y.

    """

    
    dx = x_values[1] - x_values[0]
    dy = y_values[1] - y_values[0]

    if indexing == 'xy':
        smoothed_array = smooth_array(array = array,
                                      grid_spacing = np.array([dy, dx]),
                                      gaussian_rms_widths = np.array([y_stdev, x_stdev]),
                                      truncate=truncate,
                                      mode=mode)

    elif indexing == 'ij':
        smoothed_array = smooth_array(array = array,
                                      grid_spacing = np.array([dx, dy]),
                                      gaussian_rms_widths = np.array([x_stdev, y_stdev]),
                                      truncate=truncate,
                                      mode=mode).T

    else:
        raise Exception('Indexing scheme not recognised. Should be ij or xy.')
    
    dSAdydy, dSAdxdx = np.gradient(smoothed_array)

    interps = (interp2d(x_values, y_values, smoothed_array, kind='linear'),
               interp2d(x_values, y_values, dSAdxdx/dx, kind='linear'),
               interp2d(x_values, y_values, dSAdydy/dy, kind='linear'))

    return interps


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

    Parameters
    ----------
    m : Material
        The material instance evaluated by the output function.
    attributes : list of strings
        The list of material attributes / properties to
        be evaluated in the product
    powers : list of floats
        The powers to which each attribute should be raised 
        during evaluation
    Returns
    -------
    f : function(x)
        Function which returns the value of product(a_i**p_i)
        as a function of condition (x = [P, T, V])
    """
    if type(attributes) is str:
        attributes = [attributes]
    if powers == []:
        powers = [1. for a in attributes]
    def f(x):
        P, T, V = x
        m.set_state(P, T)
        value = 1.
        for a, p in zip(*[attributes, powers]):
            value *= np.power(getattr(m, a), p)
        return value
    return f
