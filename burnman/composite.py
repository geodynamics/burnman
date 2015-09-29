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
            assert(len(phases)==len(fractions))

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
                
        else:
                self.molar_fractions=None



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

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for phase in self.phases:
            phase.set_method(method)

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

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        self.pressure = pressure
        self.temperature = temperature
        for phase in self.phases:
            phase.set_state(pressure, temperature)

    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        try:
            densities = np.array([phase.density() for phase in self.phases])
            volumes = np.array([phase.molar_volume()*molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        except:
            raise Exception("Density cannot be computed without phase fractions being defined.")
        return np.sum(densities*volumes)/np.sum(volumes)

    def calculate_moduli(self):
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
            e.V = volume_fraction * mineral.molar_volume()
            e.K = mineral.adiabatic_bulk_modulus()
            e.G = mineral.shear_modulus()
            e.rho = mineral.molar_mass() / mineral.molar_volume()
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
                total_moles+=mass_fractions[i]/mineral.molar_mass()
            molar_fractions=[]
            for i, mineral in enumerate(minerals):
                molar_fractions.append(mass_fractions[i]/(mineral.molar_mass()*total_moles))
        except AttributeError:
            raise Exception("Mass fractions cannot be set before composition has been set for all the phases in the composite.")
    
        return molar_fractions

    def elastic_properties(self, averaging_scheme=averaging_schemes.VoigtReussHill()):
        """
        Calculates the seismic velocities of the Assemblage, using
        an averaging scheme for the velocities of individual phases

        :type averaging_scheme: :class:`burnman.averaging_schemes.averaging_scheme`
        :param averaging_scheme: Averaging scheme to use.

        :returns: :math:`\\rho` :math:`[kg/m^3]` , :math:`V_p, V_s,` and :math:`V_{\phi}` :math:`[m/s]`, bulk modulus :math:`K` :math:`[Pa]`,shear modulus :math:`G` :math:`[Pa]`
        :rtype: lists of floats

        """
        moduli = [self.calculate_moduli()]
        moduli = burnman.average_moduli(moduli, averaging_scheme)[0]
        self.Vp, self.Vs, self.Vphi = burnman.compute_velocity(moduli)
        self.rho = moduli.rho
        self.K = moduli.K
        self.G = moduli.G

        return self.rho, self.Vp, self.Vs, self.Vphi, self.K, self.G

    def chemical_potentials(self, component_formulae):
        component_formulae_dict=[burnman.chemicalpotentials.dictionarize_formula(f) for f in component_formulae]
        return burnman.chemicalpotentials.chemical_potentials(self.phases, component_formulae_dict)

    def fugacity(self, standard_material):
        return burnman.chemicalpotentials.fugacity(standard_material, self.phases)

    def molar_mass(self):
        """
        Returns molar mass of the composite [kg/mol]
        """
        molar_mass = sum([ self.phases[i][0].molar_mass()*self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return molar_mass
