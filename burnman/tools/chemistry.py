# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

# This module provides higher level chemistry-related functions.

from __future__ import absolute_import
import numpy as np
from scipy.optimize import fsolve
import itertools
from sympy import Rational

from .. import constants

# Import common lower level functions for backwards compatibility
from ..classes.polytope import MaterialPolytope
from ..utils.chemistry import dictionarize_formula, formula_mass
from ..utils.chemistry import formula_to_string, site_occupancies_to_strings
from ..utils.chemistry import compositional_array
from ..utils.chemistry import reaction_matrix_as_strings


def fugacity(standard_material, assemblage):
    """
    Calculates the fugacity of a standard material in another assemblage.

    .. note:: set_method and set_state should already have been
        used on both assemblages.

    :param standard_material: Standard material for which to calculate the fugacity.
        The material must have a formula as a dictionary parameter.
    :type assemblage: :class:`burnman.Material`

    :param assemblage: Assemblage for which to calculate the fugacity.
    :type assemblage: :class:`burnman.Composite`

    :returns: Value of the fugacity of the component with respect to
        the standard material.
    :rtype: float
    """
    component_formula = standard_material.params["formula"]
    chemical_potential = assemblage.chemical_potential([component_formula])[0]

    fugacity = np.exp(
        (chemical_potential - standard_material.gibbs)
        / (constants.gas_constant * assemblage.temperature)
    )
    return fugacity


def relative_fugacity(component_formula, assemblage, reference_assemblage):
    """
    Calculates the fugacity of a chemical component in one assemblage
    relative to another one.

    .. note:: set_method and set_state should already have been
        used on both assemblages.

    :param component_formula: Chemical formula for which to compute the
        relative fugacity.
    :type component_formula: dictionary

    :param assemblage: Assemblage for which to calculate the fugacity.
    :type assemblage: :class:`burnman.Composite`

    :param reference_assemblage: Reference assemblage against which to
        measure the fugacity.
    :type reference_assemblage: :class:`burnman.Composite`

    :returns: Value of the fugacity of the component in the assemblage
        with respect to the reference_assemblage.
    :rtype: float
    """
    chemical_potential = assemblage.chemical_potential([component_formula])[0]
    reference_chemical_potential = reference_assemblage.chemical_potential(
        [component_formula]
    )[0]

    relative_fugacity = np.exp(
        (chemical_potential - reference_chemical_potential)
        / (constants.gas_constant * assemblage.temperature)
    )
    return relative_fugacity


def equilibrium_pressure(
    minerals, stoichiometry, temperature, pressure_initial_guess=1.0e5
):
    """
    Given a list of minerals, their reaction stoichiometries
    and a temperature of interest, compute the
    equilibrium pressure of the reaction.

    :param minerals: List of minerals involved in the reaction.
    :type minerals: list of :class:`burnman.Mineral`

    :param stoichiometry: Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol].
    :type stoichiometry: list of floats

    :param temperature: Temperature of interest [K].
    :type temperature: float

    :param pressure_initial_guess: Initial pressure guess [Pa].
    :type pressure_initial_guess: float

    :returns: The equilibrium pressure of the reaction [Pa].
    :rtype: float
    """

    def eqm(P, T):
        gibbs = 0.0
        for i, mineral in enumerate(minerals):
            mineral.set_state(P[0], T)
            gibbs = gibbs + mineral.gibbs * stoichiometry[i]
        return gibbs

    pressure = fsolve(eqm, [pressure_initial_guess], args=(temperature))[0]

    return pressure


def equilibrium_temperature(
    minerals, stoichiometry, pressure, temperature_initial_guess=1000.0
):
    """
    Given a list of minerals, their reaction stoichiometries
    and a pressure of interest, compute the
    equilibrium temperature of the reaction.

    :param minerals: List of minerals involved in the reaction.
    :type minerals: list of :class:`burnman.Mineral`

    :param stoichiometry: Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol].
    :type stoichiometry: list of floats

    :param pressure: Pressure of interest [Pa].
    :type pressure: float

    :param temperature_initial_guess: Initial temperature guess [K].
    :type temperature_initial_guess: float

    :returns: The equilibrium temperature of the reaction [K].
    :rtype: float
    """

    def eqm(T, P):
        gibbs = 0.0
        for i, mineral in enumerate(minerals):
            mineral.set_state(P, T[0])
            gibbs = gibbs + mineral.gibbs * stoichiometry[i]
        return gibbs

    temperature = fsolve(eqm, [temperature_initial_guess], args=(pressure))[0]

    return temperature


def invariant_point(
    minerals_r1,
    stoichiometry_r1,
    minerals_r2,
    stoichiometry_r2,
    pressure_temperature_initial_guess=[1.0e9, 1000.0],
):
    """
    Given a list of minerals, their reaction stoichiometries
    and a pressure of interest, compute the
    equilibrium temperature of the reaction.

    :param minerals: List of minerals involved in the reaction.
    :type minerals: list of :class:`burnman.Mineral`

    :param stoichiometry: Reaction stoichiometry for the minerals provided.
        Reactants and products should have the opposite signs [mol].
    :type stoichiometry: list of floats

    :param pressure: Pressure of interest [Pa].
    :type pressure: float

    :param temperature_initial_guess: Initial temperature guess [K].
    :type temperature_initial_guess: float

    :returns: The equilibrium temperature of the reaction [K].
    :rtype: float
    """

    def eqm(PT):
        P, T = PT
        gibbs_r1 = 0.0
        for i, mineral in enumerate(minerals_r1):
            mineral.set_state(P, T)
            gibbs_r1 = gibbs_r1 + mineral.gibbs * stoichiometry_r1[i]
        gibbs_r2 = 0.0
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

    :param mineral:  Mineral for which the Hugoniot is to be calculated.
    :type mineral: :class:`burnman.Mineral`

    :param P_ref: Reference pressure [Pa]
    :type P_ref: float

    :param T_ref:  Reference temperature [K]
    :type T_ref: float

    :param pressures: Set of pressures [Pa] for which the Hugoniot temperature
        and volume should be calculated.
    :type pressures: numpy.array of floats

    :param reference_mineral: Mineral which is stable at the reference conditions
        Provides an alternative U_0 and V_0 when the reference
        mineral transforms to the mineral of interest at some
        (unspecified) pressure.
    :type reference_mineral: :class:`burnman.Mineral`

    :returns: The Hugoniot temperatures and volumes at the given pressures.
    :rtype: tuple of numpy.arrays
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
            Ediff, [T_ref], args=(mineral, P, P_ref, U_ref, V_ref)
        )[0]
        volumes[i] = mineral.V

    return temperatures, volumes


def reactions_from_stoichiometric_matrix(stoichiometric_matrix):
    """
    Returns a list of all the balanced reactions between compounds
    of fixed chemical composition. Includes both the forward and
    reverse reactions
    (so there will always be an even number of reactions).

    :param stoichiometric_matrix: An array of the stoichiometric
        (molar) amounts of component j in compound i.
    :type stoichiometric_matrix: 2D numpy array

    :returns: An array of the stoichiometric (molar) amounts of
        compound j in reaction i.
    :rtype: 2D numpy array
    """
    n_components = stoichiometric_matrix.shape[1]

    equalities = np.concatenate(
        ([np.zeros(n_components)], np.array(stoichiometric_matrix).astype(float))
    ).T

    polys = [
        MaterialPolytope(equalities, np.diag(v).astype(float))
        for v in itertools.product(*[[-1, 1]] * len(equalities[0]))
    ]
    reactions = []
    for p in polys:
        v = np.array(
            [
                [Rational(value).limit_denominator(1000000) for value in v]
                for v in p.raw_vertices
            ],
            dtype=float,
        )

        if v is not []:
            reactions.extend(v)

    reactions = np.unique(np.array(reactions), axis=0)
    reactions = np.array(
        [[Rational(value).limit_denominator(1000000) for value in v] for v in reactions]
    )

    assert np.max(reactions[:-1, 0]) == 0
    assert np.max(reactions[-1, 1:]) == 0
    reactions = reactions[:-1, 1:]
    return reactions


def reactions_from_formulae(formulae, compound_names, return_strings=True):
    """
    Returns a list of all the balanced reactions between compounds
    of fixed chemical composition. Includes both the forward and
    reverse reactions
    (so there will always be an even number of reactions).

    :param formulae: List of the chemical formulae, either as strings or
        as a list of dictionaries of elements.
    :type formulae: list of dictionaries or list of strings

    :param compound_names: List of the compound names in the formula list.
    :type compound_names: list of strings

    :param return_strings: Whether to return the reactions as strings or array.
    :type return_strings: bool

    :returns: Either a 2D array of the stoichiometric (molar) amounts of
        compound j in reaction i, or a list of strings.
        The parameter compound_names is only used if strings
        are requested.
    :rtype: 2D numpy array or list of strings
    """
    if isinstance(formulae[0], str):
        dict_formulae = [dictionarize_formula(f) for f in formulae]
    else:
        dict_formulae = formulae
    stoichiometric_matrix, elements = compositional_array(dict_formulae)
    R = reactions_from_stoichiometric_matrix(stoichiometric_matrix)
    if return_strings:
        return reaction_matrix_as_strings(R, compound_names)
    else:
        return R
