# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
This module provides several helper minerals/materials.

"""
from __future__ import absolute_import
from __future__ import print_function

from .material import Material
from .composite import Composite
from .material import Material, material_property


class HelperRockSwitcher(Material):
    """
    A Helper that represents a Material that switches between different rocks based on a user specified
    select_rock() function based on current temperature and pressure. This class can be used in several
    ways:
    1. By creating an instance and setting select_rock to a lambda that returns a rock
    2. By deriving from this class and implementing select_rock.
    """
    def __init__(self):
        self.current_rock = None
        Material.__init__(self)

    def select_rock(self):
        raise NotImplementedError("Need to implement select_rock() in derived class!")

    def set_method(self, method):
        raise NotImplementedError("Need to implement select_rock() in derived class!")

    def debug_print(self, indent=""):
        print("%sHelperRockSwitcher" % (indent))

    def set_state(self, pressure, temperature):
        Material.set_state(self, pressure, temperature)

        self.current_rock = self.select_rock()
        self.current_rock.set_state(pressure, temperature)

    def unroll(self):
        return self.current_rock.unroll()

    @material_property
    def internal_energy(self):
        return self.current_rock.internal_energy

    @material_property
    def molar_gibbs(self):
        return self.current_rock.molar_gibbs

    @material_property
    def molar_helmholtz(self):
        return self.current_rock.molar_helmholtz

    @material_property
    def molar_mass(self):
        return self.current_rock.molar_mass

    @material_property
    def molar_volume(self):
        return self.current_rock.molar_volume

    @material_property
    def density(self):
        return self.current_rock.density

    @material_property
    def molar_entropy(self):
        return self.current_rock.molar_entropy

    @material_property
    def molar_enthalpy(self):
        return self.current_rock.molar_enthalpy

    @material_property
    def isothermal_bulk_modulus(self):
        return self.current_rock.isothermal_bulk_modulus

    @material_property
    def adiabatic_bulk_modulus(self):
        return self.current_rock.adiabatic_bulk_modulus

    @material_property
    def isothermal_compressibility(self):
        return self.current_rock.isothermal_compressibility

    @material_property
    def adiabatic_compressibility(self):
        return self.current_rock.adiabatic_compressibility

    @material_property
    def shear_modulus(self):
        return self.current_rock.shear_modulus

    @material_property
    def p_wave_velocity(self):
        return self.current_rock.p_wave_velocity

    @material_property
    def bulk_sound_velocity(self):
        return self.current_rock.bulk_sound_velocity

    @material_property
    def shear_wave_velocity(self):
        return self.current_rock.shear_wave_velocity

    @material_property
    def grueneisen_parameter(self):
        return self.current_rock.grueneisen_parameter

    @material_property
    def thermal_expansivity(self):
        return self.current_rock.thermal_expansivity

    @material_property
    def heat_capacity_v(self):
        return self.current_rock.heat_capacity_v

    @material_property
    def heat_capacity_p(self):
        return self.current_rock.heat_capacity_p


class HelperLowHighPressureRockTransition(HelperRockSwitcher):
    """
    A Helper that represents a Material that switches between two given rocks based on a given transition pressure.
    """
    def __init__(self, transition_pressure, low_pressure_rock, high_pressure_rock):
        self.transition_pressure = transition_pressure
        self.rocks = [low_pressure_rock, high_pressure_rock]
        HelperRockSwitcher.__init__(self)
        self._name = "HelperLowHighPressureRockTransition("+str(self.transition_pressure)+ " GPa, " + self.rocks[0].name +  ", " + self.rocks[1].name + ")"

    def select_rock(self):
        if self._pressure < self.transition_pressure:
            return self.rocks[0]
        else:
            return self.rocks[1]

    def set_method(self, method):
        for r in self.rocks:
            r.set_method(method)

    def debug_print(self, indent=""):
        print("%sHelperLowHighPressureRockTransition (%f GPa):" % (indent, self.transition_pressure))
        indent += "  "
        for r in self.rocks:
            r.debug_print(indent)



class HelperSpinTransition(Composite):

    """
    Helper class that makes a mineral that switches between two materials
    (for low and high spin) based on some transition pressure [Pa]
    """

    def __init__(self, transition_pressure, ls_mat, hs_mat):
        """
        Takes a transition pressure, and two minerals.  Use the
        thermoelastic parameters for ls_mat below the transition
        pressure, and the thermoelastic parameters for hs_mat
        above the transition pressure
        """
        Material.__init__(self)
        self.transition_pressure = transition_pressure
        self.ls_mat = ls_mat
        self.hs_mat = hs_mat
        Composite.__init__(self, [ls_mat, hs_mat])

    def debug_print(self, indent=""):
        print("%sHelperSpinTransition:" % indent)
        self.ls_mat.debug_print(indent + "  ")
        self.hs_mat.debug_print(indent + "  ")

    def set_state(self, pressure, temperature):
        if (pressure >= self.transition_pressure):
            Composite.set_fractions(self, [1.0, 0.0])
        else:
            Composite.set_fractions(self, [0.0, 1.0])

        Composite.set_state(self, pressure, temperature)
