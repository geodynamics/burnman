# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
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
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        volume : float
            Molar volume of the mineral. :math:`[m^3]`
        """
        raise NotImplementedError("")

    def pressure(self, temperature, volume, params):
        """
        Parameters
        ----------
        volume : float
            Molar volume at which to evaluate the equation of state. [m^3]
        temperature : float
            Temperature at which to evaluate the equation of state. [K]
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        pressure : float
            Pressure of the mineral, including cold and thermal parts. [m^3]
        """
        raise NotImplementedError("")

    def density(self, volume, params):
        """
        Calculate the density of the mineral :math:`[kg/m^3]`.
        The params object must include a "molar_mass" field.

        Parameters
        ----------
        volume : float
        Molar volume of the mineral.  For consistency this should be calculated
        using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        density : float
            Density of the mineral. :math:`[kg/m^3]`
        """
        return params["molar_mass"] / volume

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        gamma : float
            Grueneisen parameter of the mineral. :math:`[unitless]`
        """
        raise NotImplementedError("")

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        K_T : float
            Isothermal bulk modulus of the mineral. :math:`[Pa]`
        """
        raise NotImplementedError("")

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        K_S : float
            Adiabatic bulk modulus of the mineral. :math:`[Pa]`
        """
        raise NotImplementedError("")

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        G : float
            Shear modulus of the mineral. :math:`[Pa]`
        """
        raise NotImplementedError("")

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        C_V : float
            Heat capacity at constant volume of the mineral. :math:`[J/K/mol]`
        """
        raise NotImplementedError("")

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        C_P : float
            Heat capacity at constant pressure of the mineral. :math:`[J/K/mol]`
        """
        raise NotImplementedError("")

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. :math:`[Pa]`
        temperature : float
            Temperature at which to evaluate the equation of state. :math:`[K]`
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. :math:`[m^3]`
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        alpha : float
            Thermal expansivity of the mineral. :math:`[1/K]`
        """
        raise NotImplementedError("")

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. [Pa]
        temperature : float
            Temperature at which to evaluate the equation of state. [K]
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. [m^3]
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        G : float
            Gibbs free energy of the mineral
        """
        raise NotImplementedError("")

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        temperature : float
            Temperature at which to evaluate the equation of state. [K]
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. [m^3]
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        F : float
            Helmholtz free energy of the mineral
        """
        raise NotImplementedError("")

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """

        raise NotImplementedError("")

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. [Pa]
        temperature : float
            Temperature at which to evaluate the equation of state. [K]
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        H : float
            Enthalpy of the mineral
        """
        raise NotImplementedError("")

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the equation of state. [Pa]
        temperature : float
            Temperature at which to evaluate the equation of state. [K]
        volume : float
            Molar volume of the mineral.  For consistency this should be calculated
            using :func:`volume`. [m^3]
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        U : float
            Internal energy of the mineral
        """
        raise NotImplementedError("")

    def validate_parameters(self, params):
        """
        The params object is just a dictionary associating mineral physics parameters
        for the equation of state.  Different equation of states can have different parameters,
        and the parameters may have ranges of validity.  The intent of this function is
        twofold. First, it can check for the existence of the parameters that the
        equation of state needs, and second, it can check whether the parameters have reasonable
        values.  Unreasonable values will frequently be due to unit issues (e.g., supplying
        bulk moduli in GPa instead of Pa). In the base class this function does nothing,
        and an equation of state is not required to implement it.  This function will
        not return anything, though it may raise warnings or errors.

        Parameters
        ----------
        params : dictionary
            Dictionary containing material parameters required by the equation of state.
        """
        pass
