# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


class EquationOfState(object):
    """
    This class defines the interface for an equation of state
    that a mineral uses to determine its properties at a
    given :math:`P, T`.  In order define a new equation of state, you
    should define these functions.

    All functions should accept and return values in SI units.

    In general these functions are functions of pressure,
    temperature, and volume, as well as a "params" object,
    which is a Python dictionary that stores the material
    parameters of the mineral, such as reference volume,
    Debye temperature, reference moduli, etc.

    The functions for volume and density are just functions
    of temperature, pressure, and "params"; after all, it
    does not make sense for them to be functions of volume or density.
    """

    def volume(self, pressure, temperature, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Molar volume of the mineral :math:`[m^3]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def pressure(self, temperature, volume, params):
        """
        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Pressure of the mineral, including cold and thermal parts
            :math:`[m^3]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def density(self, volume, params):
        """
        Calculate the density of the mineral :math:`[kg/m^3]`.
        The params object must include a "molar_mass" field.

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Density of the mineral. :math:`[kg/m^3]`
        :rtype: float
        """
        return params["molar_mass"] / volume

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Grueneisen parameter of the mineral. :math:`[unitless]`
        :rtype: float
        """
        raise NotImplementedError("")

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Isothermal bulk modulus of the mineral. :math:`[Pa]`
        :rtype: float
        """
        raise NotImplementedError("")

    def isentropic_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Adiabatic bulk modulus of the mineral. :math:`[Pa]`
        :rtype: float
        """
        raise NotImplementedError("")

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Shear modulus of the mineral. :math:`[Pa]`
        :rtype: float
        """
        raise NotImplementedError("")

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Heat capacity at constant volume of the mineral. :math:`[J/K/mol]`
        :rtype: float
        """
        raise NotImplementedError("")

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Heat capacity at constant pressure of the mineral. :math:`[J/K/mol]`
        :rtype: float
        """
        raise NotImplementedError("")

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Thermal expansivity of the mineral. :math:`[1/K]`
        :rtype: float
        """
        raise NotImplementedError("")

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Gibbs energy of the mineral :math:`[J/mol]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Helmholtz energy of the mineral :math:`[J/mol]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def entropy(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Entropy of the mineral :math:`[J/K/mol]`.
        :rtype: float
        """

        raise NotImplementedError("")

    def enthalpy(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Enthalpy of the mineral :math:`[J/mol]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        :param pressure: Pressure at which to evaluate the equation of state
            :math:`[Pa]`.
        :type pressure: float

        :param temperature: Temperature at which to evaluate the equation of state
            :math:`[K]`.
        :type temperature: float

        :param volume: Molar volume of the mineral. For consistency this should be
            calculated using :func:`volume` :math:`[m^3]`.
        :type volume: float

        :param params: Dictionary containing material parameters required by
            the equation of state.
        :type params: dict

        :returns: Internal energy of the mineral :math:`[J/mol]`.
        :rtype: float
        """
        raise NotImplementedError("")

    def validate_parameters(self, params):
        """
        The params object is just a dictionary associating mineral physics parameters
        for the equation of state.
        Different equation of states can have different parameters,
        and the parameters may have ranges of validity.  The intent of this function is
        twofold. First, it can check for the existence of the parameters that the
        equation of state needs, and second, it can check whether the parameters
        have reasonable values.  Unreasonable values will frequently be due
        to unit issues (e.g., supplying bulk moduli in GPa instead of Pa).
        In the base class this function does nothing, and an equation of state
        is not required to implement it.  This function will
        not return anything, though it may raise warnings or errors.

        :param params: Dictionary containing material parameters required by the
            equation of state.
        :type params: dict
        """
        pass
