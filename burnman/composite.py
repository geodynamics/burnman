# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import warnings

from .material import Material
from .mineral import Mineral


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
        if fractions is not None:
            self.set_fractions(fractions, fraction_type)
        else:
            self.molar_fractions=None

        self.set_averaging_scheme('VoigtReussHill')
        self.set_state(1.e5,300.)

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
        elif fraction_type == 'volume':
            molar_fractions = self._volume_to_molar_fractions(self.phases, fractions)
        else:
            raise Exception("Fraction type not recognised. Please use 'molar', 'mass' or 'volume'")

           
        # Set minimum value of a molar fraction at 0.0 (rather than -1.e-12)
        self.molar_fractions = [max(0.0, fraction) for fraction in molar_fractions]  

    def set_composition(self, phase_compositions):
        assert(len(self.phases)==len(phase_compositions))
        for i, composition in enumerate(phase_compositions):
            if composition != []:
                self.phases[i].set_composition(composition)
            

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
        self.pressure = pressure
        self.temperature = temperature
        for phase in self.phases:
                phase.set_state(pressure, temperature)
        #self.rho = self.Vp = self.Vs = self.Vphi = self.K = self.G = None #Do we need this?
            
    def composition(self):
        self.phase_compositions=[self.phases[i].composition() for i in range(len(self.phases))]
        bulk_composition=dict()
        for i, composition in enumerate(self.phase_compositions):
            for element in composition:
                if element not in bulk_composition:
                    bulk_composition[element] = self.molar_fractions[i]*composition[element]
                else:
                    bulk_composition[element] += self.molar_fractions[i]*composition[element]
        return bulk_composition
            
>>>>>>> Added examples and core functionality. Elastic/thermodynamic functionality not yet tested.
    def debug_print(self, indent=""):
        print("%sComposite:" % indent)
        indent += "  "
        try: # if fractions have been defined
            for i, phase in enumerate(self.phases):
                print "%s%g of" % (indent, self.molar_fractions[i])
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

    def evaluate_phase_properties(self,vars_list):
        values = np.empty((len(vars_list),len(self.phases)))
        for v,var in enumerate(vars_list):
            for p,phase in enumerate(self.phases):
                values[v,p] = getattr(phase.var)()



    def calculate_elastic_properties(self, averaging_scheme=averaging_schemes.VoigtReussHill()):
        """
        Calculates the seismic velocities of the Assemblage, using
        an averaging scheme for the velocities of individual phases

        :type averaging_scheme: :class:`burnman.averaging_schemes.averaging_scheme`
        :param averaging_scheme: Averaging scheme to use.

        :returns: :math:`\\rho` :math:`[kg/m^3]` , :math:`V_p, V_s,` and :math:`V_{\phi}` :math:`[m/s]`, bulk modulus :math:`K` :math:`[Pa]`,shear modulus :math:`G` :math:`[Pa]`
        :rtype: lists of floats

        """
        moduli = [self._calculate_moduli()]
        moduli = burnman.average_moduli(moduli, averaging_scheme)[0]
        self.Vp, self.Vs, self.Vphi = burnman.compute_velocity(moduli)
        self.rho = moduli.rho
        self.K = moduli.K
        self.G = moduli.G

    def chemical_potentials(self, component_formulae):
        component_formulae_dict=[chemicalpotentials.dictionarize_formula(f) for f in component_formulae]
        return chemicalpotentials.chemical_potentials(self.phases, component_formulae_dict)

    def fugacity(self, standard_material):
        return chemicalpotentials.fugacity(standard_material, self.phases)

    @property
    def molar_mass(self):
        """
        Returns molar mass of the composite [kg/mol]
        """
        self._molar_mass = sum([ phase.molar_mass*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        return self._molar_mass

    @property
    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        densities = np.array([phase.density for phase in self.phases])
        volumes = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        self._density = self.averaging_scheme.average_density(volumes,densities)
        return self._density

    @property
    def molar_volume(self):
        """
        Returns molar volume of the composite [m^3/mol]
        """
        volumes = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        self._molar_volume = np.sum(volumes)
        return self._molar_volume

    @property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the composite [unitless]
        """
        self._grueneisen_paramteer = self.thermal_expansivity * self.isothermal_bulk_modulus / (self.density * self.heat_capacity_v)
        return self._grueneisen_parameter

    @property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the composite [Pa]
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.isothermal_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
   
        self._isothermal_bulk_modulus = self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)
        return self._isothermal_bulk_modulus

    @property
    def compressibility(self):
        """
        Returns compressibility of the composite (or inverse isothermal bulk modulus) [1/Pa]
        """
        self._compressibility = 1./self.isothermal_bulk_modulus()
        return self._compressibility

    @property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
        
        self._adiabatic_bulk_modulus = self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)
        
        return self._adiabatic_bulk_modulus

    @property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array([phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])
        
        self._shear_modulus = self.averaging_scheme.average_shear_moduli(V_frac, K_ph, G_ph)
                                                        
        return self._shear_modulus

    @property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the composite [1/K]
        """
        V_frac = np.array([phase.molar_volume*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        alphas = np.array([phase.thermal_expansivity for phase in self.phases])
        
        self._thermal_expansivity = self.averaging_scheme.average_thermal_expansivity(volumes,alphas)
        return self._thermal_expansivity

    @property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the composite [J/K/mol]
        """
        c_v = np.array([phase.heat_capacity_v for phase in self.phases])
        self._heat_capacity_v = self.averaging_scheme.averaging_heat_capacity_v(self.molar_fractions,c_v)
        return self._heat_capacity_v

    @property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the composite [J/K/mol]
        """
        c_p = np.array([phase.heat_capacity_v for phase in self.phases])
        self._heat_capacity_p = self.averaging_scheme.averaging_heat_capacity_p(self.molar_fractions,c_p)
        return self._heat_capacity_p

    @property
    def v_s(self):
        """
        Returns shear wave speed of the composite [m/s]
        """
        self._v_s = np.sqrt(self.shear_modulus / self.density)
        return self._v_s

    @property
    def v_p(self):
        """
        Returns P wave speed of the composite [m/s]
        """
        self._v_p =  np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
                            self.shear_modulus) / self.density)
        return self._v_p

    @property
    def v_phi(self):
        """
        Returns bulk sound speed of the composite [m/s]
        """
        self._v_phi =  np.sqrt(self.adiabatic_bulk_modulus / self.density)
        return self._v_phi

    @property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the composite [J]
        """
        raise NotImplementedError("need to implement molar_gibbs() for Composite()!")
        return None

    @property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        """
        raise NotImplementedError("need to implement molar_helmholtz() for Composite()!")
        return None

    @property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        raise NotImplementedError("need to implement molar_enthalpy() for Composite!")
        return None

    @property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        raise NotImplementedError("need to implement molar_entropy() in derived class!")
        return None




    def _calculate_moduli(self):
        """
        Calculate the elastic moduli and densities of the individual phases in the assemblage.

        :returns:
        answer -- an array of (n_evaluation_points by n_phases) of
        elastic_properties(), so the result is of the form
        answer[pressure_idx][phase_idx].V
        :rtype: list of list of :class:`burnman.elastic_properties`
        """
        elastic_properties_of_phases = []
        (minerals, molar_fractions) = self.unroll()
        volume_fractions = self._molar_to_volume_fractions(minerals, molar_fractions)

        for (mineral, volume_fraction) in zip(minerals, volume_fractions):

            e = burnman.ElasticProperties()
            e.V = volume_fraction * mineral.molar_volume
            e.K = mineral.adiabatic_bulk_modulus()
            e.G = mineral.shear_modulus()
            e.rho = mineral.molar_mass / mineral.molar_volume
            e.fraction = volume_fraction
            elastic_properties_of_phases.append(e)

        return elastic_properties_of_phases

    def _molar_to_volume_fractions(self, minerals, molar_fractions):
        total_volume=0.
        try:
            for i, mineral in enumerate(minerals):
                total_volume+=molar_fractions[i]*mineral.V
            volume_fractions=[]
            for i, mineral in enumerate(minerals):
                volume_fractions.append(molar_fractions[i]*mineral.V/total_volume)
        except AttributeError:
            raise Exception("Volume fractions cannot be set before set_state.")
    
        return volume_fractions

    def _volume_to_molar_fractions(self, minerals, volume_fractions):
        total_moles=0.
        for i, mineral in enumerate(minerals):
            total_moles+=volume_fractions[i]/mineral.V
        molar_fractions=[]
        for i, mineral in enumerate(minerals):
            molar_fractions.append(volume_fractions[i]/(mineral.V*total_moles))
            
        return molar_fractions

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


