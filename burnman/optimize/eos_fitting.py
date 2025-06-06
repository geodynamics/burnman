# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np

from . import nonlinear_fitting
from ..utils.misc import flatten
from ..utils.math import unit_normalize
from .nonlinear_fitting import NonLinearModel, nonlinear_least_squares_fit


def default_mle_tolerances(material, flags):
    """
    Computes default tolerances for maximum likelihood estimation (MLE).

    The tolerances are computed per property flag, using either:
    - A fixed absolute tolerance (e.g., 1.0 J) for energy-related properties
    - A fixed absolute tolerance (e.g., 100 Pa) for pressure
    - A relative tolerance (1e-5 x standard state value) for all other properties

    All computed tolerances must be > 0.

    :param material: The material instance to evaluate properties on.
    :type material: :class:`burnman.Material`

    :param flags: List of property names (strings) for which MLE tolerances
        should be generated.
    :type flags: list of str

    :returns: A NumPy array of tolerances, one per property flag.
    :rtype: np.ndarray

    :raises Exception: If any computed tolerance is not strictly positive.
    """
    material.set_state(1.0e5, 300.0)
    mle_tolerance_factor = 1.0e-5
    mle_tolerances = np.empty(len(flags))

    energy_flags = {
        "H",
        "energy",
        "molar_internal_energy",
        "helmholtz",
        "molar_helmholtz",
        "gibbs",
        "molar_gibbs",
        "enthalpy",
        "molar_enthalpy",
    }
    pressure_flags = {"P", "pressure"}

    for i, flag in enumerate(flags):
        if flag in energy_flags:
            mle_tolerances[i] = 1.0  # Absolute tolerance in joules
        elif flag in pressure_flags:
            mle_tolerances[i] = 100.0  # Absolute tolerance in pascals
        else:
            value = getattr(material, flag)
            mle_tolerances[i] = mle_tolerance_factor * value

    if not np.all(mle_tolerances > 0.0):
        raise Exception(
            f"All MLE tolerances should be > 0. They are currently: {mle_tolerances}"
        )

    return mle_tolerances


class MineralFit(NonLinearModel):
    """
    Class for fitting mineral parameters to experimental data.
    Instances of this class are passed to
    :func:`burnman.nonlinear_least_squares_fit`.

    For attributes added to this model when fitting is done,
    please see the documentation for that function.
    """

    def __init__(
        self,
        mineral,
        data,
        data_covariances,
        flags,
        fit_params,
        mle_tolerances,
        delta_params=None,
        bounds=None,
    ):
        self.m = mineral
        self.data = data
        self.data_covariances = data_covariances
        self.flags = flags
        self.fit_params = fit_params
        self.mle_tolerances = mle_tolerances
        if delta_params is None:
            self.delta_params = self.get_params() * 1.0e-5 + 1.0e-10
        else:
            self.delta_params = delta_params
        self.bounds = bounds

    def set_params(self, param_values):
        i = 0

        if self.bounds is not None:
            param_values = np.clip(param_values, self.bounds[:, 0], self.bounds[:, 1])

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
        return np.array(flatten([self.m.params[prm] for prm in self.fit_params]))

    def function(self, x, flag):
        P, T, p = x
        self.m.set_state(P, T)
        return np.array([P, T, getattr(self.m, flag)])

    def normal(self, x, flag):
        P, T, p = x

        if flag == "V":
            self.m.set_state(P, T)
            dPdp = -self.m.isothermal_bulk_modulus_reuss / self.m.V
            dpdT = self.m.alpha * self.m.V
        elif flag == "H":
            self.m.set_state(P, T)
            dPdp = 1.0 / ((1.0 - T * self.m.alpha) * self.m.V)
            dpdT = self.m.molar_heat_capacity_p
        elif flag == "S":
            self.m.set_state(P, T)
            dPdp = -1.0 / (self.m.alpha * self.m.V)
            dpdT = self.m.molar_heat_capacity_p / T
        elif flag == "gibbs":
            self.m.set_state(P, T)
            dPdp = 1.0 / self.m.V
            dpdT = -self.m.S
        else:
            dP = 1.0e5
            dT = 1.0
            dPdp = (2.0 * dP) / (
                self.function([P + dP, T, 0.0], flag)[2]
                - self.function([P - dP, T, 0.0], flag)[2]
            )
            dpdT = (
                self.function([P, T + dT, 0.0], flag)[2]
                - self.function([P, T - dT, 0.0], flag)[2]
            ) / (2.0 * dT)
        dPdT = -dPdp * dpdT
        n = np.array([-1.0, dPdT, dPdp])
        return unit_normalize(n)


def fit_PTp_data(
    mineral,
    fit_params,
    flags,
    data,
    data_covariances=[],
    mle_tolerances=[],
    param_tolerance=1.0e-5,
    delta_params=None,
    bounds=None,
    max_lm_iterations=50,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=True,
):
    """
    Given a mineral of any type, a list of fit parameters
    and a set of P-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    :param mineral: Mineral for which the parameters should be optimized.
    :type mineral: :class:`burnman.Mineral`

    :param fit_params: List of dictionary keys contained in mineral.params
        corresponding to the variables to be optimized
        during fitting. Initial guesses are taken from the existing
        values for the parameters
    :type fit_params: list of str

    :param flags: Attribute names for the property to be fit for the whole
        dataset or each datum individually (e.g. 'V')
    :type flags: string or list of strings

    :param data: Observed X-P-T-property values
    :type data: 2D numpy.array

    :param data_covariances: X-P-T-property covariances (optional)
        If not given, all covariance matrices are chosen
        such that all data points have equal weight,
        with all error in the pressure.
    :type data_covariances: 3D numpy.array

    :param mle_tolerances: Tolerances for termination of the
        maximum likelihood iterations (optional).
    :type mle_tolerances: numpy.array

    :param param_tolerance: Fractional tolerance for termination
        of the nonlinear optimization (optional).
    :type param_tolerance: float

    :param delta_params: Initial values for the change in parameters (optional).
    :type delta_params: numpy.array

    :param bounds: Minimum and maximum bounds for the parameters (optional).
        The shape must be (n_parameters, 2).
    :type bounds: 2D numpy.array

    :param max_lm_iterations: Maximum number of Levenberg-Marquardt iterations.
    :type max_lm_iterations: int

    :param verbose: Whether to print detailed information about the
        optimization to screen.
    :type verbose: bool

    :returns: Model with optimized parameters.
    :rtype: :class:`burnman.optimize.eos_fitting.MineralFit`
    """

    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:, 0]))

    if len(flags) != len(data):
        raise Exception(
            f"The number of flags (n = {len(flags)}) must be equal "
            f"to the number of data (n = {len(data)})."
        )

    # Apply mle tolerances if they dont exist
    if len(mle_tolerances) == 0:
        mle_tolerances = default_mle_tolerances(mineral, flags)

    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and p)
    covariances_defined = True
    if len(data_covariances) == 0:
        covariances_defined = False
        data_covariances = np.zeros((len(data[:, 0]), len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][0][0] = 1.0

    model = MineralFit(
        mineral=mineral,
        data=data,
        data_covariances=data_covariances,
        flags=flags,
        fit_params=fit_params,
        delta_params=delta_params,
        mle_tolerances=mle_tolerances,
        bounds=bounds,
    )

    nonlinear_least_squares_fit(
        model,
        max_lm_iterations=max_lm_iterations,
        param_tolerance=param_tolerance,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=verbose,
    )

    if verbose is True and covariances_defined is True:
        confidence_interval = 0.9
        d = nonlinear_fitting.extreme_values(
            model.weighted_residuals, confidence_interval
        )
        confidence_bound, indices, probabilities = d
        if indices != []:
            print(
                "The function nonlinear_fitting.extreme_values"
                "(model.weighted_residuals, confidence_interval) "
                f"has determined that there are {len(indices):d} data points"
                " which have residuals which are not expected at the "
                f"{confidence_interval*100.:.1f}% confidence level "
                f"(> {confidence_bound:.1f} s.d. away from the model fit).\n"
                "Their indices and the probabilities of finding "
                "such extreme values are:"
            )
            for i, idx in enumerate(indices):
                print(
                    f"[{idx:d}]: {probabilities[i]:.4f} "
                    f"({np.abs(model.weighted_residuals[idx]):.1f} s.d. "
                    "from the model)"
                )
            print(
                "You might consider removing them from your fit, "
                "or increasing the uncertainties in their "
                "measured values.\n"
            )

    return model


def fit_PTV_data(
    mineral,
    fit_params,
    data,
    data_covariances=[],
    delta_params=None,
    bounds=None,
    param_tolerance=1.0e-5,
    max_lm_iterations=50,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=True,
):
    """
    A simple alias for the fit_PTp_data for when all the data is volume data
    """

    return fit_PTp_data(
        mineral=mineral,
        flags="V",
        data=data,
        data_covariances=data_covariances,
        fit_params=fit_params,
        param_tolerance=param_tolerance,
        delta_params=delta_params,
        bounds=bounds,
        max_lm_iterations=max_lm_iterations,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=verbose,
    )


class MineralFitV(NonLinearModel):
    """
    Class for fitting mineral parameters to experimental data
    using volume as an independent variable.
    Instances of this class are passed to
    :func:`burnman.nonlinear_least_squares_fit`.

    For attributes added to this model when fitting is done,
    please see the documentation for that function.
    """

    def __init__(
        self,
        mineral,
        data,
        data_covariances,
        flags,
        fit_params,
        mle_tolerances,
        delta_params=None,
        bounds=None,
    ):
        self.m = mineral
        self.data = data
        self.data_covariances = data_covariances
        self.flags = flags
        self.fit_params = fit_params
        self.mle_tolerances = mle_tolerances
        if delta_params is None:
            self.delta_params = self.get_params() * 1.0e-5 + 1.0e-10
        else:
            self.delta_params = delta_params
        self.bounds = bounds

    def set_params(self, param_values):
        i = 0

        if self.bounds is not None:
            param_values = np.clip(param_values, self.bounds[:, 0], self.bounds[:, 1])

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
        return np.array(flatten([self.m.params[prm] for prm in self.fit_params]))

    def function(self, x, flag):
        V, T, p = x
        self.m.set_state_with_volume(V, T)
        return np.array([V, T, getattr(self.m, flag)])

    def normal(self, x, flag):
        V, T, p = x

        # See Stacey and Hodgkinson (2019), Table A2.
        if flag == "P" or flag == "pressure":
            self.m.set_state_with_volume(V, T)
            dVdp = -self.m.V / self.m.isothermal_bulk_modulus_reuss  # dVdP|T
            dpdT = self.m.alpha * self.m.isothermal_bulk_modulus_reuss  # dPdT|V
        elif flag == "S" or flag == "molar_entropy":
            self.m.set_state_with_volume(V, T)
            dVdp = 1.0 / (self.m.alpha * self.m.isothermal_bulk_modulus_reuss)  # dVdS|T
            dpdT = self.m.molar_heat_capacity_v / T  # dSdT|V
        elif flag == "helmholtz" or flag == "molar_helmholtz":
            self.m.set_state_with_volume(V, T)
            dVdp = -1.0 / self.m.pressure  # dVdF|T
            dpdT = -self.m.S  # dFdT|V
        else:
            dV = V * 1.0e-5
            dT = 1.0
            dVdp = (2.0 * dV) / (
                self.function([V + dV, T, 0.0], flag)[2]
                - self.function([V - dV, T, 0.0], flag)[2]
            )
            dpdT = (
                self.function([V, T + dT, 0.0], flag)[2]
                - self.function([V, T - dT, 0.0], flag)[2]
            ) / (2.0 * dT)
        dVdT = -dVdp * dpdT
        n = np.array([-1.0, dVdT, dVdp])
        return unit_normalize(n)


def fit_VTp_data(
    mineral,
    fit_params,
    flags,
    data,
    data_covariances=[],
    mle_tolerances=[],
    param_tolerance=1.0e-5,
    delta_params=None,
    bounds=None,
    max_lm_iterations=50,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=True,
):
    """
    Given a mineral of any type, a list of fit parameters
    and a set of V-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    :param mineral: Mineral for which the parameters should be optimized.
    :type mineral: :class:`burnman.Mineral`

    :param fit_params: List of dictionary keys contained in mineral.params
        corresponding to the variables to be optimized
        during fitting. Initial guesses are taken from the existing
        values for the parameters
    :type fit_params: list of str

    :param flags: Attribute names for the property to be fit for the whole
        dataset or each datum individually (e.g. 'P')
    :type flags: string or list of strings

    :param data: Observed V-T-property values
    :type data: 2D numpy.array

    :param data_covariances: V-T-property covariances (optional)
        If not given, all covariance matrices are chosen
        such that all data points have equal weight,
        with all error in the properties.
    :type data_covariances: 3D numpy.array

    :param mle_tolerances: Tolerances for termination of the
        maximum likelihood iterations (optional).
    :type mle_tolerances: numpy.array

    :param param_tolerance: Fractional tolerance for termination
        of the nonlinear optimization (optional).
    :type param_tolerance: float

    :param delta_params: Initial values for the change in parameters (optional).
    :type delta_params: numpy.array

    :param bounds: Minimum and maximum bounds for the parameters (optional).
        The shape must be (n_parameters, 2).
    :type bounds: 2D numpy.array

    :param max_lm_iterations: Maximum number of Levenberg-Marquardt iterations.
    :type max_lm_iterations: int

    :param verbose: Whether to print detailed information about the
        optimization to screen.
    :type verbose: bool

    :returns: Model with optimized parameters.
    :rtype: :class:`burnman.optimize.eos_fitting.MineralFitV`
    """

    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:, 0]))

    if len(flags) != len(data):
        raise Exception(
            f"The number of flags (n = {len(flags)}) must be equal "
            f"to the number of data (n = {len(data)})."
        )

    # Apply mle tolerances if they dont exist
    if len(mle_tolerances) == 0:
        mle_tolerances = default_mle_tolerances(mineral, flags)

    # If covariance matrix is not given, apply unit weighting to all properties
    # (with zero errors on V and T)
    covariances_defined = True
    if len(data_covariances) == 0:
        covariances_defined = False
        data_covariances = np.zeros((len(data[:, 0]), len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][2][2] = 1.0

    model = MineralFitV(
        mineral=mineral,
        data=data,
        data_covariances=data_covariances,
        flags=flags,
        fit_params=fit_params,
        delta_params=delta_params,
        mle_tolerances=mle_tolerances,
        bounds=bounds,
    )

    nonlinear_least_squares_fit(
        model,
        max_lm_iterations=max_lm_iterations,
        param_tolerance=param_tolerance,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=verbose,
    )

    if verbose is True and covariances_defined is True:
        confidence_interval = 0.9
        d = nonlinear_fitting.extreme_values(
            model.weighted_residuals, confidence_interval
        )
        confidence_bound, indices, probabilities = d
        if indices != []:
            print(
                "The function nonlinear_fitting.extreme_values"
                "(model.weighted_residuals, confidence_interval) "
                f"has determined that there are {len(indices):d} data points"
                " which have residuals which are not expected at the "
                f"{confidence_interval*100.:.1f}% confidence level "
                f"(> {confidence_bound:.1f} s.d. away from the model fit).\n"
                "Their indices and the probabilities of finding "
                "such extreme values are:"
            )
            for i, idx in enumerate(indices):
                print(
                    f"[{idx:d}]: {probabilities[i]:.4f} "
                    f"({np.abs(model.weighted_residuals[idx]):.1f} s.d. "
                    "from the model)"
                )
            print(
                "You might consider removing them from your fit, "
                "or increasing the uncertainties in their "
                "measured values.\n"
            )

    return model


def fit_VTP_data(
    mineral,
    fit_params,
    data,
    data_covariances=[],
    delta_params=None,
    bounds=None,
    param_tolerance=1.0e-5,
    max_lm_iterations=50,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=True,
):
    """
    A simple alias for the fit_VTp_data for when all the data is pressure data
    """

    return fit_VTp_data(
        mineral=mineral,
        flags="P",
        data=data,
        data_covariances=data_covariances,
        fit_params=fit_params,
        param_tolerance=param_tolerance,
        delta_params=delta_params,
        bounds=bounds,
        max_lm_iterations=max_lm_iterations,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=verbose,
    )


class SolutionFit(NonLinearModel):
    """
    Class for fitting mineral parameters to experimental data.
    Instances of this class are passed to
    :func:`burnman.nonlinear_least_squares_fit`.

    For attributes added to this model when fitting is done,
    please see the documentation for that function.
    """

    def __init__(
        self,
        solution,
        data,
        data_covariances,
        flags,
        fit_params,
        mle_tolerances,
        delta_params=None,
        bounds=None,
    ):
        self.m = solution
        self.data = data
        self.data_covariances = data_covariances
        self.flags = flags
        self.fit_params = fit_params
        self.fit_params_strings = []
        for p in fit_params:
            if isinstance(p, list):
                csv_list_mbrs = ",".join([str(i) for i in p[1:]])
                self.fit_params_strings.append(f"{p[0]} ({csv_list_mbrs})")
            else:
                self.fit_params_strings.append(p)

        self.mle_tolerances = mle_tolerances
        if delta_params is None:
            self.delta_params = self.get_params() * 1.0e-5 + 1.0e-10
        else:
            self.delta_params = delta_params
        self.bounds = bounds

    def set_params(self, param_values):
        # fit_params is a list of lists
        # if the list has length 2, the first item should be an integer
        # indicating the endmember number in the solution
        # if the list has length 3, the first two items should be endmember
        # numbers, and the third should be the interaction parameter type
        # (E, S or V).
        i = 0

        if self.bounds is not None:
            param_values = np.clip(param_values, self.bounds[:, 0], self.bounds[:, 1])

        for param in self.fit_params:
            value = param_values[i]
            if len(param) == 2:
                key, imbr = param
                if isinstance(self.m.endmembers[imbr][0].params[key], float):
                    self.m.endmembers[imbr][0].params[key] = value
                    i += 1
                else:
                    n_values = len(self.m.endmembers[imbr][0].params[key])
                    for j in range(n_values):
                        self.m.endmembers[imbr][0].params[key][j] = value
                        i += 1
            elif len(param) == 3:
                key, imbr, jmbr = param
                ai = self.m.solution_model.alphas[imbr]
                aj = self.m.solution_model.alphas[jmbr]
                if key == "E":
                    self.m.solution_model.We[imbr, jmbr] = 2.0 * value / (ai * aj)
                if key == "S":
                    self.m.solution_model.Ws[imbr, jmbr] = 2.0 * value / (ai * aj)
                if key == "V":
                    self.m.solution_model.Wv[imbr, jmbr] = 2.0 * value / (ai * aj)

                i += 1
            else:
                raise Exception("param length must be two or three")

    def get_params(self):
        params = []
        for param in self.fit_params:
            if len(param) == 2:
                key, imbr = param
                value = self.m.endmembers[imbr][0].params[key]
                if isinstance(value, float):
                    params.append(value)
                else:
                    params.extend(list(value))

            elif len(param) == 3:
                key, imbr, jmbr = param
                ai = self.m.solution_model.alphas[imbr]
                aj = self.m.solution_model.alphas[jmbr]
                if key == "E":
                    params.append(
                        self.m.solution_model.We[imbr, jmbr] * (ai * aj) / 2.0
                    )
                if key == "S":
                    params.append(
                        self.m.solution_model.Ws[imbr, jmbr] * (ai * aj) / 2.0
                    )
                if key == "V":
                    params.append(
                        self.m.solution_model.Wv[imbr, jmbr] * (ai * aj) / 2.0
                    )
            else:
                raise Exception("param length must be two or three")
        return np.array(params)

    def function(self, x, flag):
        self.m.set_composition(x[: self.m.n_endmembers])
        P, T, p = x[self.m.n_endmembers :]
        self.m.set_state(P, T)

        f = np.copy(x)
        f[-1] = getattr(self.m, flag)
        return f

    def normal(self, x, flag):
        self.m.set_composition(x[: self.m.n_endmembers])
        P, T, p = x[self.m.n_endmembers :]

        if flag == "V":
            self.m.set_state(P, T)
            dPdp = -self.m.isothermal_bulk_modulus_reuss / self.m.V
            dpdT = self.m.alpha * self.m.V
        elif flag == "H":
            self.m.set_state(P, T)
            dPdp = 1.0 / ((1.0 - T * self.m.alpha) * self.m.V)
            dpdT = self.m.molar_heat_capacity_p
        elif flag == "S":
            self.m.set_state(P, T)
            dPdp = -1.0 / (self.m.alpha * self.m.V)
            dpdT = self.m.molar_heat_capacity_p / T
        elif flag == "gibbs":
            self.m.set_state(P, T)
            dPdp = 1.0 / self.m.V
            dpdT = -self.m.S
        else:
            dP = 1.0e5
            dT = 1.0
            xP0 = np.copy(x)
            xP1 = np.copy(x)
            xT0 = np.copy(x)
            xT1 = np.copy(x)
            xP0[-3] = xP1[-3] - dP
            xP1[-3] = xP1[-3] + dP
            xT0[-2] = xP1[-2] - dT
            xT1[-2] = xP1[-2] + dT

            dPdp = (2.0 * dP) / (
                self.function(xP1, flag)[2] - self.function(xP0, flag)[2]
            )
            dpdT = (self.function(xT1, flag)[2] - self.function(xT0, flag)[2]) / (
                2.0 * dT
            )
        dPdT = -dPdp * dpdT
        n = np.zeros(len(x))
        n[-3:] = np.array([-1.0, dPdT, dPdp])
        return unit_normalize(n)


def fit_XPTp_data(
    solution,
    fit_params,
    flags,
    data,
    data_covariances=[],
    mle_tolerances=[],
    param_tolerance=1.0e-5,
    delta_params=None,
    bounds=None,
    max_lm_iterations=50,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=True,
):
    """
    Given a symmetric solution, a list of fit parameters
    and a set of P-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    :param solution: Solution for which the parameters should be optimized.
    :type solution: :class:`burnman.Solution`

    :param fit_params: Variables to be optimized
        during fitting. Each list is either of length two or three.
        The first item of length-2 lists should be a
        dictionary key contained in one of the endmember
        mineral.params, and the second item should be the index of
        the endmember in the solution (indexing starts from 0).
        The first item of length-3 lists should be one of 'E', 'S' or
        'V' (the excess energies, entropies or volumes in each binary).
        The second two items should be the indices of the pair of
        endmembers bounding the binary, in ascending order
        (indexing starts from 0). Initial guesses are taken from the existing
        values for the parameters.
    :type fit_params: list of lists

    :param flags: Attribute names for the property to be fit for the whole
        dataset or each datum individually (e.g. 'V')
    :type flags: string or list of strings

    :param data: Observed X-P-T-property values
    :type data: 2D numpy.array

    :param data_covariances: X-P-T-property covariances (optional).
        If not given, all covariance matrices are chosen
        such that all data points have equal weight,
        with all error in the pressure.
    :type data_covariances: 3D numpy.array

    :param mle_tolerances: Tolerances for termination of the
        maximum likelihood iterations (optional).
    :type mle_tolerances: numpy.array

    :param param_tolerance: Fractional tolerance for termination
        of the nonlinear optimization (optional).
    :type param_tolerance: float

    :param delta_params: Initial values for the change in parameters (optional).
    :type delta_params: numpy.array

    :param bounds: Minimum and maximum bounds for the parameters (optional).
        The shape must be (n_parameters, 2).
    :type bounds: 2D numpy.array

    :param max_lm_iterations: Maximum number of Levenberg-Marquardt iterations.
    :type max_lm_iterations: int

    :param verbose: Whether to print detailed information about the
        optimization to screen.
    :type verbose: bool

    :returns: Model with optimized parameters.
    :rtype: :class:`burnman.optimize.eos_fitting.SolutionFit`
    """
    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:, 0]))

    if len(flags) != len(data):
        raise Exception(
            f"The number of flags (n = {len(flags)}) must be equal "
            f"to the number of data (n = {len(data)})."
        )

    # Apply mle tolerances if they dont exist
    if len(mle_tolerances) == 0:
        mle_tolerances = default_mle_tolerances(solution, flags)

    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and property)
    covariances_defined = True
    if len(data_covariances) == 0:
        covariances_defined = False
        nX = solution.n_endmembers
        data_covariances = np.zeros((len(data[:, 0]), len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][nX][nX] = 1.0

    model = SolutionFit(
        solution=solution,
        data=data,
        data_covariances=data_covariances,
        flags=flags,
        fit_params=fit_params,
        mle_tolerances=mle_tolerances,
        delta_params=delta_params,
        bounds=bounds,
    )

    nonlinear_least_squares_fit(
        model,
        max_lm_iterations=max_lm_iterations,
        param_tolerance=param_tolerance,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=verbose,
    )

    if verbose is True and covariances_defined is True:
        confidence_interval = 0.9
        v_extreme = nonlinear_fitting.extreme_values(
            model.weighted_residuals, confidence_interval
        )
        confidence_bound, indices, probabilities = v_extreme
        if indices != []:
            print(
                "The function nonlinear_fitting.extreme_values"
                "(model.weighted_residuals, confidence_interval) "
                f"has determined that there are {len(indices):d} "
                "data points which have residuals which are not "
                f"expected at the {confidence_interval*100.:.1f}% "
                "confidence level "
                f"(> {confidence_bound:.1f} s.d. away from the model fit).\n"
                "Their indices and the probabilities of "
                "finding such extreme values are:"
            )
            for i, idx in enumerate(indices):
                print(
                    f"[{idx:d}]: {probabilities[i]:.4f} "
                    f"({np.abs(model.weighted_residuals[idx]):.1f} s.d. "
                    "from the model)"
                )
            print(
                "You might consider removing them from your fit, "
                "or increasing the uncertainties in "
                "their measured values.\n"
            )

    return model
