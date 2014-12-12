
# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

class EquationOfState(object):
    """
    This class defines the interface for an equation of state
    that a mineral uses to determine its properties at a
    given P,T.  In order define a new equation of state, you
    should define these functions.

    All functions should accept and return values in SI units.

    In general these functions are functions of pressure,
    temperature, and volume, as well as a "params" object,
    which is a Python dictionary that stores the material 
    parameters of the mineral, such as reference volume,
    Debye temperature, reference moduli, etc.

    The functions for volume and density are just functions 
    of temperature, pressure, and "params"; after all, it 
    does not make sense for them to be functions of volume/density.
    """

    def volume(self, pressure, temperature, params):
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
        volume : float
            Molar volume of the mineral. [m^3]
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

    def density(self, pressure, temperature, params):
        """
        Calculate the density of the mineral.  
        The params object must include a "molar_mass" field.

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
        density : float
            Density of the mineral. [kg/m^3]
        """
        return params["molar_mass"] / self.volume(pressure, temperature, params)

    def grueneisen_parameter(self, pressure, temperature, volume, params):
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
        gamma : float
            Grueneisen parameter of the mineral. [unitless]
        """
        raise NotImplementedError("")

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
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
        K_T : float
            Isothermal bulk modulus of the mineral. [Pa]
        """
        raise NotImplementedError("")

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
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
        K_S : float
            Adiabatic bulk modulus of the mineral. [Pa]
        """
        raise NotImplementedError("")

    def shear_modulus(self, pressure, temperature, volume, params):
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
            Shear modulus of the mineral. [Pa]
        """
        raise NotImplementedError("")

    def heat_capacity_v(self, pressure, temperature, volume, params):
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
        C_V : float
            Heat capacity at constant volume of the mineral. [J/K/mol]
        """
        raise NotImplementedError("")

    def heat_capacity_p(self, pressure, temperature, volume, params):
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
        C_P : float
            Heat capacity at constant pressure of the mineral. [J/K/mol]
        """
        raise NotImplementedError("")

    def thermal_expansivity(self, pressure, temperature, volume, params):
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
        alpha : float
            Thermal expansivity of the mineral. [1/K]
        """
        raise NotImplementedError("")

    def reference_temperature( self, params ):
        """
        Parameters
        ----------
        params : dictionary
            Dictionary containing material parameters required by the equation of state.

        Returns
        -------
        T_0 : float
            If params contains a "T_0" entry, return that reference temperature, otherwise
            return 300.0 Kelvin
        """
        if 'T_0' in params:
            return params['T_0']
        else:
            return 300.0


    def gibbs_free_energy( self, pressure, temperature, params ):
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

    def helmholtz_free_energy( self, temperature, volume, params ):
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

    def entropy( self, pressure, temperature, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """

        raise NotImplementedError("")

    def enthalpy( self, pressure, temperature, params ):
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

    def internal_energy( self, pressure, temperature, volume, params ):
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
