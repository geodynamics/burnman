# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from sympy import Matrix, nsimplify
import warnings

from .material import Material, material_property, cached_property
from .mineral import Mineral
from .solution import Solution
from . import averaging_schemes

from ..utils.reductions import independent_row_indices
from ..utils.chemistry import sum_formulae, sort_element_list_to_IUPAC_order
from ..utils.chemistry import reaction_matrix_as_strings


def check_pairs(phases, fractions):
    if len(fractions) < 1:
        raise Exception("ERROR: we need at least one phase")

    if len(phases) != len(fractions):
        raise Exception("ERROR: different array lengths for phases and fractions")

    total = sum(fractions)
    if abs(total - 1.0) > 1e-10:
        raise Exception("ERROR: list of molar fractions does not add up to one")
    for p in phases:
        if not isinstance(p, Mineral):
            raise Exception(
                "ERROR: object of type " "%s" " is not of type Mineral" % (type(p))
            )


# static composite of minerals/composites
class Composite(Material):
    """
    Base class for a composite material.
    The static phases can be minerals or materials,
    meaning composite can be nested arbitrarily.

    The fractions of the phases can be input
    as either 'molar' or 'mass' during instantiation,
    and modified (or initialised) after this point by
    using set_fractions.

    This class is available as ``burnman.Composite``.
    """

    def __init__(
        self, phases, fractions=None, fraction_type="molar", name="Unnamed composite"
    ):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        :param phases: List of phases.
        :type phases: list of :class:`burnman.Material`
        :param fractions: molar or mass fraction for each phase.
        :type fractions: list of floats
        :param fraction_type: 'molar' or 'mass' (optional, 'molar' as standard)
            specify whether molar or mass fractions are specified.
        :type fraction_type: str
        """

        Material.__init__(self)

        assert len(phases) > 0
        self.phases = phases

        if fractions is not None:
            self.set_fractions(fractions, fraction_type)
        else:
            self.molar_fractions = None

        self.set_averaging_scheme("VoigtReussHill")
        self.name = name
        self.equilibrium_tolerance = 1.0e-3  # J/reaction
        self.print_precision = 4  # number of significant figures used by self.__str__

    def __str__(self):
        string = "Composite: {0}".format(self.name)
        try:
            string += "\n  P, T: {0:.{sf}g} Pa, {1:.{sf}g} K".format(
                self.pressure, self.temperature, sf=self.print_precision
            )
        except:
            pass
        string += "\nPhase and endmember fractions:"
        for phase, fraction in zip(*self.unroll()):
            string += "\n  {0}: {1:0.{sf}f}".format(
                phase.name, fraction, sf=self.print_precision
            )
            if isinstance(phase, Solution):
                for i in range(phase.n_endmembers):
                    string += "\n    {0}: {1:0.{sf}f}".format(
                        phase.endmember_names[i],
                        phase.molar_fractions[i],
                        sf=self.print_precision,
                    )
        return string

    def set_fractions(self, fractions, fraction_type="molar"):
        """
        Change the fractions of the phases of this Composite.
        Resets cached properties

        :param fractions: list or numpy array of floats
            molar or mass fraction for each phase.
        :param fraction_type: 'molar' or 'mass'
            specify whether molar or mass fractions are specified.
        """
        assert len(self.phases) == len(fractions)

        if isinstance(fractions, list):
            fractions = np.array(fractions)

        try:
            total = sum(fractions)
        except TypeError:
            raise Exception(
                "Since v0.8, burnman.Composite takes an array of Materials, "
                "then an array of fractions"
            )

        assert np.all(fractions >= -1e-12)

        self.reset()

        if abs(total - 1.0) > 1e-12:
            warnings.warn(
                "Warning: list of fractions does not add "
                f"up to one but {total:g}. Normalizing."
            )
            fractions /= total

        if fraction_type == "molar":
            molar_fractions = fractions
        elif fraction_type == "mass":
            molar_fractions = self._mass_to_molar_fractions(self.phases, fractions)
        else:
            raise Exception(
                "Fraction type not recognised. " "Please use 'molar' or mass"
            )

        # Set minimum value of a molar fraction at 0.0 (rather than -1.e-12)
        self.molar_fractions = molar_fractions.clip(0.0)

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for phase in self.phases:
            phase.set_method(method)
        # Clear the cache on resetting method
        self.reset()

    def set_averaging_scheme(self, averaging_scheme):
        """
        Set the averaging scheme for the moduli in the composite.
        Default is set to VoigtReussHill, when Composite is initialized.
        """

        if type(averaging_scheme) == str:
            self.averaging_scheme = getattr(averaging_schemes, averaging_scheme)()
        else:
            self.averaging_scheme = averaging_scheme
        # Clear the cache on resetting averaging scheme
        self.reset()

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        Material.set_state(self, pressure, temperature)
        for phase in self.phases:
            phase.set_state(pressure, temperature)

    def debug_print(self, indent=""):
        print("{0}Composite: {1}".format(indent, self.name))
        indent += "  "
        if self.molar_fractions is None:
            for i, phase in enumerate(self.phases):
                phase.debug_print(indent + "  ")
        else:
            for i, phase in enumerate(self.phases):
                print("%s%g of" % (indent, self.molar_fractions[i]))
                phase.debug_print(indent + "  ")

    def unroll(self):
        if self.molar_fractions is None:
            raise Exception("Unroll only works if the composite has defined fractions.")
        phases = []
        fractions = []
        for i, phase in enumerate(self.phases):
            p_mineral, p_fraction = phase.unroll()
            check_pairs(p_mineral, p_fraction)
            fractions.extend([f * self.molar_fractions[i] for f in p_fraction])
            phases.extend(p_mineral)
        return phases, fractions

    def to_string(self):
        """
        return the name of the composite
        """
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    @material_property
    def formula(self):
        """
        Returns molar chemical formula of the composite
        """
        return sum_formulae([ph.formula for ph in self.phases], self.molar_fractions)

    @material_property
    def molar_internal_energy(self):
        """
        Returns molar internal energy of the mineral [J/mol]
        Aliased with self.energy
        """
        U = sum(
            phase.molar_internal_energy * molar_fraction
            for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
        )
        return U

    @material_property
    def molar_gibbs(self):
        """
        Returns molar Gibbs free energy of the composite [J/mol]
        Aliased with self.gibbs
        """
        G = sum(
            phase.molar_gibbs * molar_fraction
            for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
        )
        return G

    @material_property
    def molar_helmholtz(self):
        """
        Returns molar Helmholtz free energy of the mineral [J/mol]
        Aliased with self.helmholtz
        """
        F = sum(
            phase.molar_helmholtz * molar_fraction
            for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
        )
        return F

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the composite [m^3/mol]
        Aliased with self.V
        """
        volumes = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        return np.sum(volumes)

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the composite [kg/mol]
        """
        return sum(
            [
                phase.molar_mass * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )

    @material_property
    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        Aliased with self.rho
        """
        densities = np.array([phase.density for phase in self.phases])
        volumes = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        return self.averaging_scheme.average_density(volumes, densities)

    @material_property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.S
        """
        S = sum(
            phase.molar_entropy * molar_fraction
            for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
        )
        return S

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.H
        """
        H = sum(
            phase.molar_enthalpy * molar_fraction
            for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
        )
        return H

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        Returns isothermal bulk modulus of the composite [Pa]
        Aliased with self.K_T
        """
        V_frac = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        K_ph = np.array([phase.isothermal_bulk_modulus_reuss for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        """
        V_frac = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        K_ph = np.array([phase.isentropic_bulk_modulus_reuss for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns isothermal compressibility of the composite
        (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_T
        """
        return 1.0 / self.isothermal_bulk_modulus_reuss

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        Returns isothermal compressibility of the composite
        (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_S
        """
        return 1.0 / self.isentropic_bulk_modulus_reuss

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        """
        V_frac = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        K_ph = np.array([phase.isentropic_bulk_modulus_reuss for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_shear_moduli(V_frac, K_ph, G_ph)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the composite [m/s]
        Aliased with self.v_p
        """
        return np.sqrt(
            (self.isentropic_bulk_modulus_reuss + 4.0 / 3.0 * self.shear_modulus)
            / self.density
        )

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the composite [m/s]
        Aliased with self.v_phi
        """
        return np.sqrt(self.isentropic_bulk_modulus_reuss / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the composite [m/s]
        Aliased with self.v_s
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the composite [unitless]
        Aliased with self.gr
        """
        return (
            self.thermal_expansivity
            * self.isothermal_bulk_modulus_reuss
            * self.molar_volume
            / self.molar_heat_capacity_v
        )

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the composite [1/K]
        Aliased with self.alpha
        """
        volumes = np.array(
            [
                phase.molar_volume * molar_fraction
                for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)
            ]
        )
        alphas = np.array([phase.thermal_expansivity for phase in self.phases])
        return self.averaging_scheme.average_thermal_expansivity(volumes, alphas)

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar_heat capacity at constant volume of the composite [J/K/mol]
        Aliased with self.C_v
        """
        c_v = np.array([phase.molar_heat_capacity_v for phase in self.phases])
        return self.averaging_scheme.average_heat_capacity_v(self.molar_fractions, c_v)

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure of the composite [J/K/mol]
        Aliased with self.C_p
        """
        c_p = np.array([phase.molar_heat_capacity_p for phase in self.phases])
        return self.averaging_scheme.average_heat_capacity_p(self.molar_fractions, c_p)

    @material_property
    def endmember_partial_gibbs(self):
        """
        Returns the partial Gibbs energies for all
        the endmember minerals in the Composite
        """
        partial_gibbs = np.empty(self.n_endmembers)
        j = 0
        for i, n_endmembers in enumerate(self.endmembers_per_phase):
            if n_endmembers == 1:
                partial_gibbs[j] = self.phases[i].gibbs
            else:
                partial_gibbs[j : j + n_endmembers] = self.phases[i].partial_gibbs
            j += n_endmembers
        return partial_gibbs

    @material_property
    def reaction_affinities(self):
        """
        Returns the affinities corresponding to each reaction in reaction_basis
        """
        return self.reaction_basis.dot(self.endmember_partial_gibbs)

    @material_property
    def equilibrated(self):
        """
        Returns True if the reaction affinities are all zero
        within a given tolerance given by self.equilibrium_tolerance.
        """
        return np.all(np.abs(self.reaction_affinities) < self.equilibrium_tolerance)

    def set_components(self, components):
        """
        Sets the components and components_array attributes of the
        composite material. The components attribute is a list of dictionaries
        containing the chemical formulae of the components.
        The components_array attribute is a 2D numpy array describing the
        linear transformation between endmember amounts and component amounts.

        The components do not need to be linearly independent, not do they need
        to form a complete basis set for the composite.
        However, it must be possible to obtain the composition of each
        component from a linear sum of the endmember compositions of
        the composite. For example, if the composite was composed of
        MgSiO3 and Mg2SiO4, SiO2 would be a valid component, but Si would not.
        The method raises an exception if any of the chemical potentials are
        not defined by the assemblage.

        :param components: List of formulae of the components.
        :type components: list of dictionaries
        """
        # Convert components into array form
        b = np.array(
            [
                [component[el] if el in component else 0.0 for component in components]
                for el in self.elements
            ]
        )

        # Solve to find a set of endmember proportions that
        # satisfy each of the component formulae
        p = np.linalg.lstsq(self.stoichiometric_array.T, b, rcond=None)

        res = np.abs((self.stoichiometric_array.T.dot(p[0]) - b).T)
        res = np.sum(res, axis=1)
        # Check that all components can be described by linear sums of
        # the endmembers
        if not np.all(res < 1.0e-12):
            bad_indices = np.argwhere(res > 1.0e-12)

            raise Exception(
                f"Components {bad_indices} not defined by " "prescribed assemblage"
            )

        self.components = components
        self.component_array = p[0]

    def chemical_potential(self, components=None):
        """
        Returns the chemical potentials of the currently defined components
        in the composite. Raises an exception if
        the assemblage is not equilibrated.

        :param components: List of formulae of the desired components.
            If not specified, the method uses the components specified
            by a previous call to set_components.
        :type components: list of dictionaries

        :returns: The chemical potentials of the desired components in the
            equilibrium composite.
        :rtype: numpy.array of floats
        """
        if not self.equilibrated:
            raise Exception(
                "This composite is not equilibrated, so "
                "it cannot have a defined chemical potential."
            )

        if components is not None:
            self.set_components(components)

        # Return the chemical potential of each component
        return np.dot(self.component_array.T, self.endmember_partial_gibbs)

    def _mass_to_molar_fractions(self, phases, mass_fractions):
        """
        Converts a set of mass fractions for phases into a set of molar fractions.

        :param phases: The list of phases for which fractions should be converted.
        :type phases: list of :class:`burnman.Material`

        :param mass_fractions: An array of mass fractions of the input phases.
        :type mass_fractions: numpy.array of floats

        :returns: An array of molar fractions corresponding to the
            input molar fractions.
        :rtype: numpy.array of floats
        """
        molar_masses = np.array([phase.molar_mass for phase in phases])
        moles = mass_fractions / molar_masses
        return moles / sum(moles)

    @cached_property
    def stoichiometric_matrix(self):
        """
        An sympy Matrix where each element M[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """

        def f(i, j):
            e = self.elements[j]
            if e in self.endmember_formulae[i]:
                return nsimplify(self.endmember_formulae[i][e])
            else:
                return 0

        return Matrix(self.n_endmembers, self.n_elements, f)

    @cached_property
    def stoichiometric_array(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """
        return np.array(self.stoichiometric_matrix).astype(float)

    @cached_property
    def reaction_basis(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of moles of endmember[j] involved in reaction[i].
        """
        reaction_basis = np.array(
            [v[:] for v in self.stoichiometric_matrix.T.nullspace()], dtype=float
        )

        if len(reaction_basis) == 0:
            reaction_basis = np.empty((0, self.n_endmembers))

        return reaction_basis

    @cached_property
    def reaction_basis_as_strings(self):
        """
        Returns a list of string representations of all the reactions in
        reaction_basis.
        """
        return reaction_matrix_as_strings(self.reaction_basis, self.endmember_names)

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
        the amounts of the other elements can be
        inferred by
        -compositional_null_basis[independent_element_indices].dot(element_amounts)
        """
        return sorted(independent_row_indices(self.stoichiometric_matrix.T))

    @cached_property
    def dependent_element_indices(self):
        """
        The element indices not included in the independent list.
        """
        return [
            i
            for i in range(self.n_elements)
            if i not in self.independent_element_indices
        ]

    @cached_property
    def reduced_stoichiometric_array(self):
        """
        The stoichiometric array including only the independent elements
        """
        return self.stoichiometric_array[:, self.independent_element_indices]

    @cached_property
    def compositional_null_basis(self):
        """
        An array N such that N.b = 0 for all bulk compositions that can
        be produced with a linear sum of the endmembers in the composite.
        """
        null_basis = np.array(
            [v[:] for v in self.stoichiometric_matrix.nullspace()], dtype=float
        )
        if null_basis.shape[0] != 0:
            M = null_basis[:, self.dependent_element_indices]
            assert (M.shape[0] == M.shape[1]) and (M == np.eye(M.shape[0])).all()

        return null_basis

    @cached_property
    def endmember_formulae(self):
        """
        A list of the formulae in the composite.
        """
        self._set_endmember_properties()
        return self.__dict__["endmember_formulae"]

    @cached_property
    def endmember_names(self):
        """
        A list of the endmember names contained in the composite.
        Mineral names are returned as given in Mineral.name.
        Solution endmember names are given in the format
        `Mineral.name in Solution.name`.
        """
        self._set_endmember_properties()
        return self.__dict__["endmember_names"]

    @cached_property
    def endmembers_per_phase(self):
        """
        A list of integers corresponding to the number of endmembers
        stored within each phase.
        """
        self._set_endmember_properties()
        return self.__dict__["endmembers_per_phase"]

    @cached_property
    def elements(self):
        """
        A list of the elements which could be contained in the composite,
        returned in the IUPAC element order.
        """
        self._set_endmember_properties()
        return self.__dict__["elements"]

    @cached_property
    def n_endmembers(self):
        """
        Returns the number of endmembers in the composite.
        """
        return len(self.endmember_names)

    @cached_property
    def n_elements(self):
        """
        Returns the total number of distinct elements
        which might be in the composite.
        """
        return len(self.elements)

    def _set_endmember_properties(self):
        """
        Sets endmember_formulae, endmember_names, endmembers_per_phase and
        elements as properties of the Composite. This helper function
        is used to set all properties at the same time while still allowing
        the properties to be stored and documented as individual
        cached_properties.
        """
        endmember_formulae = []
        endmember_names = []
        endmembers_per_phase = []
        for ph_idx, ph in enumerate(self.phases):
            if isinstance(ph, Solution):
                endmember_formulae.extend(ph.endmember_formulae)
                endmember_names.extend(
                    [name + " in " + ph.name for name in ph.endmember_names]
                )
                endmembers_per_phase.append(ph.n_endmembers)

            elif isinstance(ph, Mineral):
                endmember_formulae.append(ph.formula)
                endmember_names.append(ph.name)
                endmembers_per_phase.append(1)

            else:
                raise Exception(
                    "Unsupported Material type, can only read"
                    "burnman.Mineral or burnman.Solution"
                )

        # Populate the stoichiometric matrix
        keys = []
        for f in endmember_formulae:
            keys.extend(f.keys())

        # Save to dict so that we only need to do this once
        self.__dict__["endmember_formulae"] = endmember_formulae
        self.__dict__["endmember_names"] = endmember_names
        self.__dict__["endmembers_per_phase"] = endmembers_per_phase
        self.__dict__["elements"] = sort_element_list_to_IUPAC_order(set(keys))
