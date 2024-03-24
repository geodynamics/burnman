# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
from sympy import Matrix, nsimplify
import scipy.optimize as opt

from burnman.classes.solutionmodel import IdealSolution

from ..constants import gas_constant
from .material import material_property, cached_property
from .mineral import Mineral
from .elasticsolutionmodel import ElasticMechanicalSolution
from .averaging_schemes import reuss_average_function

from ..utils.math import bracket
from ..utils.reductions import independent_row_indices
from ..utils.chemistry import sum_formulae, sort_element_list_to_IUPAC_order


class ElasticSolution(Mineral):
    """
    This is the base class for all Elastic solutions.
    Site occupancies, endmember activities and the constant
    and volume and temperature dependencies of the excess
    properties can be queried after using set_composition().
    States of the solution can only be queried after setting
    the pressure, temperature and composition using set_state()
    and set_composition.

    This class is available as :class:`burnman.ElasticSolution`.
    It uses an instance of :class:`burnman.ElasticSolutionModel` to
    calculate interaction terms between endmembers.

    All the solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T
    and V derivatives in J/K/mol and Pa/mol.

    The parameters are relevant to all Elastic solution models. Please
    see the documentation for individual models for details about
    other parameters.

    :param name: Name of the solution.
    :type name: string
    :param solution_model: The ElasticSolutionModel object defining the
        properties of the solution.
    :type solution_model: :class:`burnman.ElasticSolutionModel`
    :param molar_fractions: The molar fractions of each endmember in the solution.
        Can be reset using the set_composition() method.
    :type molar_fractions: numpy.array
    """

    def __init__(self, name=None, solution_model=None, molar_fractions=None):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.
        """
        Mineral.__init__(self)

        # Solution needs a method attribute to call Mineral.set_state().
        # Note that set_method() below will not change self.method
        self.method = "ElasticSolutionMethod"

        if name is not None:
            self.name = name
        if solution_model is not None:
            self.solution_model = solution_model

        if isinstance(solution_model, ElasticMechanicalSolution):
            self.solution_type = "mechanical"
        else:
            self.solution_type = "chemical"

        # Starting guess and delta for pressure iteration
        self.min_V0 = min(
            [mbr[0].params["V_0"] for mbr in self.solution_model.endmembers]
        )
        self.dV = 0.01 * self.min_V0

        # Equation of state
        for i in range(self.n_endmembers):
            self.solution_model.endmembers[i][0].set_method(
                self.solution_model.endmembers[i][0].params["equation_of_state"]
            )

        # Molar fractions
        if molar_fractions is not None:
            self.set_composition(molar_fractions)

    @cached_property
    def endmembers(self):
        return self.solution_model.endmembers

    def set_composition(self, molar_fractions):
        """
        Set the composition for this solution.
        Resets cached properties.

        :param molar_fractions: Molar abundance for each endmember,
            needs to sum to one.
        :type molar_fractions: list of float
        """
        assert len(self.solution_model.endmembers) == len(molar_fractions)

        if self.solution_type != "mechanical":
            assert sum(molar_fractions) > 0.9999
            assert sum(molar_fractions) < 1.0001

        self.reset()
        self.molar_fractions = np.array(molar_fractions)

        if self.temperature is not None:
            _ = self.molar_volume

    def set_method(self, method):
        for i in range(self.n_endmembers):
            self.solution_model.endmembers[i][0].set_method(method)
        # note: do not set self.method here!
        self.reset()

    def set_state(self, pressure, temperature):
        Mineral.set_state(self, pressure, temperature)

        try:
            _ = self.molar_volume
        except AttributeError:
            pass

    @material_property
    def formula(self):
        """
        Returns molar chemical formula of the solution.
        """
        return sum_formulae(self.endmember_formulae, self.molar_fractions)

    @material_property
    def activities(self):
        """
        Returns a list of endmember activities [unitless].
        """
        volumes = [
            self.solution_model.endmembers[i][0].method.volume(
                self.pressure,
                self.temperature,
                self.solution_model.endmembers[i][0].params,
            )
            for i in range(self.n_endmembers)
        ]

        gibbs_pure = [
            self.solution_model.endmembers[i][0].method.gibbs_free_energy(
                self.pressure,
                self.temperature,
                volumes[i],
                self.solution_model.endmembers[i][0].params,
            )
            for i in range(self.n_endmembers)
        ]

        acts = np.exp(
            (self.partial_gibbs - np.array(gibbs_pure))
            / (gas_constant * self.temperature)
        )
        return acts

    @material_property
    def activity_coefficients(self):
        """
        Returns a list of endmember activity coefficients
        (gamma = activity / ideal activity) [unitless].
        """
        return np.exp(
            np.log(self.activities)
            - IdealSolution._log_ideal_activities(
                self.solution_model, self.molar_fractions
            )
        )

    @material_property
    def molar_internal_energy(self):
        """
        Returns molar internal energy of the mineral [J/mol].
        Aliased with self.energy
        """
        return self.molar_helmholtz + self.temperature * self.molar_entropy

    @material_property
    def _excess_partial_helmholtz(self):
        """
        Returns excess partial molar helmholtz energy
        at constant volume [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_helmholtz_energies(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _excess_partial_pressures(self):
        """
        Returns excess partial pressures at constant volume [Pa].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_pressures(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _excess_partial_entropies(self):
        """
        Returns excess partial entropies at constant volume [J/K].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_entropies(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _partial_helmholtz(self):
        """
        Returns endmember partial molar Helmholtz energy at constant volume [J/mol].
        Property specific to solutions.
        """
        return (
            np.array(
                [
                    self.solution_model.endmembers[i][0].helmholtz
                    for i in range(self.n_endmembers)
                ]
            )
            + self._excess_partial_helmholtz
        )

    @material_property
    def _partial_pressures(self):
        """
        Returns endmember partial pressures at constant volume [Pa].
        Property specific to solutions.
        """
        return (
            np.array(
                [
                    self.solution_model.endmembers[i][0].pressure
                    for i in range(self.n_endmembers)
                ]
            )
            + self._excess_partial_pressures
        )

    @material_property
    def _partial_entropies(self):
        """
        Returns endmember partial entropies at constant volume [J/K].
        Property specific to solutions.
        """
        return (
            np.array(
                [
                    self.solution_model.endmembers[i][0].molar_entropy
                    for i in range(self.n_endmembers)
                ]
            )
            + self._excess_partial_entropies
        )

    @material_property
    def partial_gibbs(self):
        """
        Returns endmember partial molar Gibbs energy
        at constant pressure [J/mol].
        Property specific to solutions.
        """
        return self._partial_helmholtz + self.pressure * self.molar_volume

    @material_property
    def _dPdX(self):
        """
        Returns the change in pressure with amount of each endmember
        at constant volume.
        """
        sumX = np.sum(self.molar_fractions)
        sumXP = np.einsum("i,i->", self.molar_fractions, self._partial_pressures)
        return (self._partial_pressures * sumX - sumXP) / (sumX * sumX)

    @material_property
    def _dVdX(self):
        """
        Returns the change in pressure with amount of each endmember
        at constant pressure.
        """
        return self.molar_volume / self.isothermal_bulk_modulus_reuss * self._dPdX

    @material_property
    def _dSdX_mod(self):
        """
        Returns the additional change in entropy with
        amount of each endmember due to converting from constant volume
        to constant pressure
        """
        return self.alpha * self.molar_volume * self._dPdX

    @material_property
    def partial_volumes(self):
        """
        Returns endmember partial molar volumes [m^3/mol].
        Property specific to solutions.
        """
        A = np.eye(self.n_endmembers) - self.molar_fractions
        Vs = self.molar_volume + np.einsum("ij, j->i", A, self._dVdX)
        return Vs

    @material_property
    def partial_entropies(self):
        """
        Returns endmember partial molar entropies [J/K/mol].
        Property specific to solutions.
        """
        A = np.eye(self.n_endmembers) - self.molar_fractions
        Ss = self._partial_entropies + np.einsum("ij, j->i", A, self._dSdX_mod)
        return Ss

    @material_property
    def _excess_helmholtz(self):
        """
        Returns molar excess Helmholtz energy at constant volume [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_helmholtz_energy(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _helmholtz_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the Helmholtz energy at constant volume [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.helmholtz_hessian(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _entropy_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the entropy at constant volume [J/K].
        Property specific to solutions.
        """
        return self.solution_model.entropy_hessian(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def _pressure_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the pressure at constant volume [Pa].
        Property specific to solutions.
        """
        return self.solution_model.pressure_hessian(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def gibbs_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the Gibbs energy at constant pressure [J/mol].
        Property specific to solutions.
        """
        raise NotImplementedError

    @material_property
    def molar_helmholtz(self):
        """
        Returns molar Helmholtz energy of the solution [J/mol].
        Aliased with self.helmholtz.
        """
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].molar_helmholtz
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self._excess_helmholtz
        )

    @material_property
    def molar_gibbs(self):
        """
        Returns molar Gibbs free energy of the solution [J/mol].
        Aliased with self.gibbs.
        """
        return self.molar_helmholtz + self.pressure * self.molar_volume

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the solution [kg/mol].
        """
        return sum(
            [
                self.solution_model.endmembers[i][0].molar_mass
                * self.molar_fractions[i]
                for i in range(self.n_endmembers)
            ]
        )

    @material_property
    def excess_pressure(self):
        """
        Returns excess pressure of the solution [Pa].
        Specific property for solutions.
        """
        return self.solution_model.excess_pressure(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the solution [m^3/mol].
        Aliased with self.V.
        """

        def _delta_pressure(volume):
            self._ptmp = [
                self.solution_model.endmembers[i][0].method.pressure(
                    self.temperature,
                    volume,
                    self.solution_model.endmembers[i][0].params,
                )
                for i in range(self.n_endmembers)
            ]

            pressure_try = sum(
                [
                    self._ptmp[i] * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            ) + self.solution_model.excess_pressure(
                volume, self.temperature, self.molar_fractions
            )

            return pressure_try - self.pressure

        def _K_T(volume):
            # Note, this only works when the excess pressure is not a function
            # of V or T.

            return sum(
                [
                    self.solution_model.endmembers[i][
                        0
                    ].method.isothermal_bulk_modulus_reuss(
                        0.0,
                        self.temperature,
                        volume,
                        self.solution_model.endmembers[i][0].params,
                    )
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )

        try:
            # The first attempt to find a bracket for
            # root finding uses V_0 as a starting point
            sol = bracket(_delta_pressure, self.min_V0, self.dV)
        except Exception:
            # At high temperature, the naive bracketing above may
            # try a volume guess that exceeds the point at which the
            # bulk modulus goes negative at that temperature.
            # In this case, we try a more nuanced approach by
            # first finding the volume at which the bulk modulus goes
            # negative, and then either (a) raising an exception if the
            # desired pressure is less than the pressure at that volume,
            # or (b) using that pressure to create a better bracket for
            # brentq.

            sol_K_T = bracket(_K_T, self.min_V0, self.dV)
            V_crit = opt.brentq(_K_T, sol_K_T[0], sol_K_T[1])
            P_min = self.pressure + _delta_pressure(V_crit)
            if P_min > self.pressure:
                raise Exception(
                    "The desired pressure is not achievable "
                    "at this temperature. The minimum pressure "
                    f"achievable is {P_min:.2e} Pa."
                )
            else:
                try:
                    sol = bracket(_delta_pressure, V_crit - self.dV, self.dV)
                except Exception:
                    raise Exception(
                        "Cannot find a volume, perhaps you are "
                        "outside of the range of validity for "
                        "the equation of state?"
                    )

        V = opt.brentq(_delta_pressure, sol[0], sol[1])

        _delta_pressure(V)
        for i in range(self.n_endmembers):
            self.solution_model.endmembers[i][0].set_state(
                self._ptmp[i], self.temperature
            )

        return V

    @material_property
    def density(self):
        """
        Returns density of the solution [kg/m^3].
        Aliased with self.rho.
        """
        return self.molar_mass / self.molar_volume

    @material_property
    def excess_entropy(self):
        """
        Returns excess molar entropy [J/K/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_entropy(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropy of the solution [J/K/mol].
        Aliased with self.S.
        """
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].S * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_entropy
        )

    @material_property
    def excess_enthalpy(self):
        """
        Returns excess molar enthalpy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_enthalpy(
            self.molar_volume, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the solution [J/mol].
        Aliased with self.H.
        """
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].H * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_enthalpy
        )

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        Returns isothermal bulk modulus of the solution [Pa].
        Aliased with self.K_T.
        """
        return sum(
            [
                self.solution_model.endmembers[i][0].isothermal_bulk_modulus_reuss
                * self.molar_fractions[i]
                for i in range(self.n_endmembers)
            ]
        )

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Returns adiabatic bulk modulus of the solution [Pa].
        Aliased with self.K_S.
        """
        if self.temperature < 1e-10:
            return self.isothermal_bulk_modulus_reuss
        else:
            return (
                self.isothermal_bulk_modulus_reuss
                * self.molar_heat_capacity_p
                / self.molar_heat_capacity_v
            )

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns isothermal compressibility of the solution.
        (or inverse isothermal bulk modulus) [1/Pa].
        Aliased with self.K_T.
        """
        return 1.0 / self.isothermal_bulk_modulus_reuss

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        Returns adiabatic compressibility of the solution.
        (or inverse adiabatic bulk modulus) [1/Pa].
        Aliased with self.K_S.
        """
        return 1.0 / self.isentropic_bulk_modulus_reuss

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the solution [Pa].
        Aliased with self.G.
        """
        G_list = np.fromiter(
            (e[0].G for e in self.solution_model.endmembers),
            dtype=float,
            count=self.n_endmembers,
        )
        return reuss_average_function(self.molar_fractions, G_list)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the solution [m/s].
        Aliased with self.v_p.
        """
        return np.sqrt(
            (self.isentropic_bulk_modulus_reuss + 4.0 / 3.0 * self.shear_modulus)
            / self.density
        )

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the solution [m/s].
        Aliased with self.v_phi.
        """
        return np.sqrt(self.isentropic_bulk_modulus_reuss / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the solution [m/s].
        Aliased with self.v_s.
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the solution [unitless].
        Aliased with self.gr.
        """
        if self.temperature < 1e-10:
            return float("nan")
        else:
            return (
                self.thermal_expansivity
                * self.isothermal_bulk_modulus_reuss
                * self.molar_volume
                / self.molar_heat_capacity_v
            )

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha)
        of the solution [1/K].
        Aliased with self.alpha.
        """
        alphaKT = sum(
            [
                self.solution_model.endmembers[i][0].isothermal_bulk_modulus_reuss
                * self.solution_model.endmembers[i][0].alpha
                * self.molar_fractions[i]
                for i in range(self.n_endmembers)
            ]
        )
        return alphaKT / self.isothermal_bulk_modulus_reuss

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the
        solution [J/K/mol].
        Aliased with self.C_v.
        """
        return sum(
            [
                self.solution_model.endmembers[i][0].molar_heat_capacity_v
                * self.molar_fractions[i]
                for i in range(self.n_endmembers)
            ]
        )

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure
        of the solution [J/K/mol].
        Aliased with self.C_p.
        """
        return (
            self.molar_heat_capacity_v
            + self.molar_volume
            * self.temperature
            * self.thermal_expansivity
            * self.thermal_expansivity
            * self.isothermal_bulk_modulus_reuss
        )

    @cached_property
    def stoichiometric_matrix(self):
        """
        A sympy Matrix where each element M[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """

        def f(i, j):
            e = self.elements[j]
            if e in self.endmember_formulae[i]:
                return nsimplify(self.endmember_formulae[i][e])
            else:
                return 0

        return Matrix(len(self.endmember_formulae), len(self.elements), f)

    @cached_property
    def stoichiometric_array(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """
        return np.array(self.stoichiometric_matrix)

    @cached_property
    def reaction_basis(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of moles of endmember[j] involved in reaction[i].
        """
        reaction_basis = np.array(
            [v[:] for v in self.stoichiometric_matrix.T.nullspace()]
        )

        if len(reaction_basis) == 0:
            reaction_basis = np.empty((0, len(self.endmember_names)))

        return reaction_basis

    @cached_property
    def n_reactions(self):
        """
        The number of reactions in reaction_basis.
        """
        return len(self.reaction_basis[:, 0])

    @cached_property
    def independent_element_indices(self):
        """
        A list of an independent set of element indices. If the amounts of
        these elements are known (element_amounts),
        the amounts of the other elements can be inferred by
        -compositional_null_basis[independent_element_indices].dot(element_amounts).
        """
        return sorted(independent_row_indices(self.stoichiometric_matrix.T))

    @cached_property
    def dependent_element_indices(self):
        """
        The element indices not included in the independent list.
        """
        return [
            i
            for i in range(len(self.elements))
            if i not in self.independent_element_indices
        ]

    @cached_property
    def compositional_null_basis(self):
        """
        An array N such that N.b = 0 for all bulk compositions that can
        be produced with a linear sum of the endmembers in the solution.
        """
        null_basis = np.array([v[:] for v in self.stoichiometric_matrix.nullspace()])

        M = null_basis[:, self.dependent_element_indices]
        assert (M.shape[0] == M.shape[1]) and (M == np.eye(M.shape[0])).all()

        return null_basis

    @cached_property
    def endmember_formulae(self):
        """
        A list of formulae for all the endmember in the solution.
        """
        return [mbr[0].params["formula"] for mbr in self.solution_model.endmembers]

    @cached_property
    def endmember_names(self):
        """
        A list of names for all the endmember in the solution.
        """
        return [mbr[0].name for mbr in self.solution_model.endmembers]

    @cached_property
    def n_endmembers(self):
        """
        The number of endmembers in the solution.
        """
        return len(self.solution_model.endmembers)

    @cached_property
    def elements(self):
        """
        A list of the elements which could be contained in the solution,
        returned in the IUPAC element order.
        """
        keys = []
        for f in self.endmember_formulae:
            keys.extend(f.keys())

        return sort_element_list_to_IUPAC_order(set(keys))


ElasticSolidSolution = ElasticSolution
