# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

class EquationOfState:
    """
    This class defines the interface for an equation of state
    that a mineral uses to determine its properties at a 
    given P,T.  In order define a new equation of state, you
    should define these functions.

    All functions should accept and return values in SI units.
   
    In general these functions are functions of pressure,
    temperature, and volume, as well as a "params" object,
    which stores the material parameters of the stuff,
    such as reference volume, Debye temperature, etc. 
    The exceptions are volume and density, which are
    just assumed to be functions of pressure and temperature.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume at the pressure and temperature [m^3]
        """
        raise NotImplementedError("")

    def density(self, pressure, temperature, params):
        """
        Returns density at the pressure and temperature [kg/m^3]
        """
        return params["molar_mass"] / self.volume(pressure, temperature, params)

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter at the pressure, temperature, and volume
        """
        raise NotImplementedError("")

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus at the pressure, temperature, and volume [Pa]
        """
        raise NotImplementedError("")

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus at the pressure, temperature, and volume [Pa]
        """
        raise NotImplementedError("")

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus at the pressure, temperature, and volume [Pa]
        """
        raise NotImplementedError("")

    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature, and volume [J/K/mol]
        """
        raise NotImplementedError("")

    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        raise NotImplementedError("")

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        """
        raise NotImplementedError("")

