# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

class equation_of_state:
    """
    This class defines the interface for an equation of state
    that a mineral uses to determine its properties at a 
    given P,T.  In order define a new equation of state, you
    should define these functions"

    All functions should accept and return values in SI units.
   
    In general these functions are functions of pressure,
    temperature, and volume, as well as a "params" object,
    which stores the material parameters of the stuff,
    such as reference volume, debye temperature, etc...
    The exceptions are volume and density, which are
    just assumed to be functions of pressure and temperature
    """

    def volume(self, pressure, temperature, params):
        raise NotImplementedError("")

    def density(self, pressure, temperature, params):
        raise NotImplementedError("")

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        raise NotImplementedError("")

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        raise NotImplementedError("")

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        raise NotImplementedError("")

    def shear_modulus(self, pressure, temperature, volume, params):
        raise NotImplementedError("")

    def heat_capacity_v(self, pressure, temperature, volume, params):
        """ heat capacity at constant volume """
        raise NotImplementedError("")

    def heat_capacity_p(self, pressure, temperature, volume, params):
        """ heat capacity at constant pressure """
        raise NotImplementedError("")

    def thermal_expansivity(self, pressure, temperature, volume, params):
        raise NotImplementedError("")

