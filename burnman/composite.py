# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import warnings

from .material import Material, material_property
from .mineral import Mineral
from . import averaging_schemes
from . import chemicalpotentials


def check_pairs(phases, fractions):
        if len(fractions) < 1:
            raise Exception('ERROR: we need at least one phase')

        if len(phases) != len(fractions):
            raise Exception('ERROR: different array lengths for phases and fractions')

        total = sum(fractions)
        if abs(total-1.0)>1e-10:
            raise Exception('ERROR: list of molar fractions does not add up to one')
        for p in phases:
            if not isinstance(p, Mineral):
                raise Exception('ERROR: object of type ''%s'' is not of type Mineral' % (type(p)))


# static composite of minerals/composites
class Composite(Material):
    """
    Base class for a static composite material with fixed molar fractions. The
    elements can be minerals or materials, meaning composite can be nested
    arbitrarily.

    This class is available as ``burnman.Composite``.
    """
    def __init__(self, phases, fractions=None, fraction_type='molar'):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        Parameters
        ----------
        phases: list of :class:`burnman.Material`
            list of phases.
        fractions: list of floats
            molar fraction for each phase.
        """

        
        assert(len(phases)>0)
        self.phases = phases
        
        if fractions is not None:
            self.set_fractions(fractions, fraction_type)
        else:
            self.molar_fractions=None

        self.set_averaging_scheme('VoigtReussHill')

        Material.__init__(self)


    def set_fractions(self, fractions, fraction_type='molar'):
        assert(len(self.phases)==len(fractions))

        try:
            total = sum(fractions)
        except TypeError:
            raise Exception("Since v0.8, burnman.Composite takes an array of Materials, then an array of fractions")

        for f in fractions:
            assert (f >= -1e-12)

        if abs(total - 1.0) > 1e-12:
            warnings.warn("Warning: list of fractions does not add up to one but %g. Normalizing." % total)
            corrected_fractions = [fr / total for fr in fractions]
            fractions = corrected_fractions


        if fraction_type == 'molar':
            molar_fractions = fractions
        elif fraction_type == 'mass':
            molar_fractions = self._mass_to_molar_fractions(self.phases, fractions)
        else:
            raise Exception("Fraction type not recognised. Please use 'molar' or mass")
           
        # Set minimum value of a molar fraction at 0.0 (rather than -1.e-12)
        self.molar_fractions = [max(0.0, fraction) for fraction in molar_fractions]  


    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for phase in self.phases:
            
            phase.set_method(method)

    def set_averaging_scheme(self,averaging_scheme):
        """
        Set the averaging scheme for the moduli in the composite.
        Default is set to VoigtReussHill, when Composite is initialized.
        """
        
        if type(averaging_scheme)==str:
            self.averaging_scheme = getattr(averaging_schemes,averaging_scheme)()
        else:
            self.averaging_scheme = averaging_scheme


    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        self.reset()
        self._pressure = pressure
        self._temperature = temperature
        for phase in self.phases:
                phase.set_state(pressure, temperature)

            
    def debug_print(self, indent=""):
        print("%sComposite:" % indent)
        indent += "  "
        try: # if fractions have been defined
            for i, phase in enumerate(self.phases):
                print("%s%g of" % (indent, self.molar_fractions[i]))
                phase.debug_print(indent + "  ")
        except:
            for i, phase in enumerate(self.phases):
                phase.debug_print(indent + "  ")


    def unroll(self):
        phases = []
        fractions = []
        try:
            for i, phase in enumerate(self.phases):
                p_mineral, p_fraction = phase.unroll()
                check_pairs(p_mineral, p_fraction)
                fractions.extend([f*self.molar_fractions[i] for f in p_fraction])
                phases.extend(p_mineral)
        except:
            raise Exception("Unroll only works if the composite has defined fractions.")
        return (phases, fractions)

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    @material_property
    def internal_energy(self):
        """
        Returns internal energy of the mineral [J]
        Aliased with self.energy
        """
        U = sum(phase.internal_energy*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return U
    
    @material_property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the composite [J]
        Aliased with self.gibbs
        """
        G = sum(phase.molar_gibbs*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return G

    @material_property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        Aliased with self.helmholtz
        """
        F = sum(phase.molar_helmholtz*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return F

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the composite [m^3/mol]
        Aliased with self.V
        """
        volumes = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        return  np.sum(volumes)
    
    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the composite [kg/mol]
        """
        return sum([ phase.molar_mass*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])


    @material_property
    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        Aliased with self.rho
        """
        densities = np.array([phase.density for phase in self.phases])
        volumes = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        return  self.averaging_scheme.average_density(volumes,densities)

            
    @material_property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.S
        """
        S = sum(phase.molar_entropy*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return S

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.H
        """
        H = sum(phase.molar_enthalpy*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return H


    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the composite [Pa]
        Aliased with self.K_T
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.isothermal_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
   
        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)


    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
        
        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)


    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the composite (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_T
        """
        return 1./self.isothermal_bulk_modulus()


    @material_property
    def adiabatic_compressibility(self):
        """
        Returns isothermal compressibility of the composite (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_S
        """
        return 1./self.adiabatic_bulk_modulus()


    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
        
        return self.averaging_scheme.average_shear_moduli(V_frac, K_ph, G_ph)
                                                        


    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the composite [m/s]
        Aliased with self.v_p
        """
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
                            self.shear_modulus) / self.density)


    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the composite [m/s]
        Aliased with self.v_phi
        """
        return  np.sqrt(self.adiabatic_bulk_modulus / self.density)
  

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
        return self.thermal_expansivity * self.isothermal_bulk_modulus / (self.density * self.heat_capacity_v)



    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the composite [1/K]
        Aliased with self.alpha
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        alphas = np.array([phase.thermal_expansivity for phase in self.phases])
        return self.averaging_scheme.average_thermal_expansivity(volumes,alphas)


    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the composite [J/K/mol]
        Aliased with self.C_v
        """
        c_v = np.array([phase.heat_capacity_v for phase in self.phases])
        return self.averaging_scheme.averaging_heat_capacity_v(self.molar_fractions,c_v)


    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the composite [J/K/mol]
        Aliased with self.C_p
        """
        c_p = np.array([phase.heat_capacity_v for phase in self.phases])
        return self.averaging_scheme.averaging_heat_capacity_p(self.molar_fractions,c_p)


    def _mass_to_molar_fractions(self, minerals, mass_fractions):
        total_moles=0.
        try:
            for i, mineral in enumerate(minerals):
                total_moles+=mass_fractions[i]/mineral.molar_mass
            molar_fractions=[]
            for i, mineral in enumerate(minerals):
                molar_fractions.append(mass_fractions[i]/(mineral.molar_mass*total_moles))
        except AttributeError:
            raise Exception("Mass fractions cannot be set before composition has been set for all the phases in the composite.")
    
        return molar_fractions

