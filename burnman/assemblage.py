# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings

import burnman
from burnman import Material
from burnman import Mineral
import burnman.averaging_schemes

def check_pairs(fractions, minerals):
        if len(fractions) < 1:
            raise Exception('ERROR: we need at least one mineral')

        if len(fractions) != len(minerals):
            raise Exception('ERROR: different array lengths')

        total = sum(fractions)
        if abs(total-1.0)>1e-10:
            raise Exception('ERROR: list of molar fractions does not add up to one')
        for p in minerals:
            if not isinstance(p, Mineral):
                raise Exception('ERROR: object of type ''%s'' is not of type material' % (type(p)))

class Assemblage(Material):
    """
    Base class for a composite material. The elements can be 
    minerals or materials, meaning composite can be nested
    arbitrarily.

    Assemblage differs from Composite in that material fractions 
    and compositions are defined with calls to set_phase_fractions,
    set_composition and set_state.

    This class is available as ``burnman.Assemblage``.
    """
    def __init__(self, phases=None):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        Parameters
        ----------
        phases: list of :class:`burnman.Material`
            list of phases.
        """
        if len(phases)==0:
            warnings.warn("Assemblage must be populated with phases before you can do anything with it!", stacklevel=2)
        self.phases=phases


    def debug_print(self, indent=""):
        print "%sComposite:" % indent
        indent += "  "
        for (fraction, phase) in self.children:
            print "%s%g of" % (indent, fraction)
            phase.debug_print(indent + "  ")    

    def set_phase_fractions(self, fractions, fraction_type='molar'):
        assert(len(self.phases)==len(fractions))

        for f in fractions:
            # we would like to check for >= 0, but this creates nasty behavior due to
            # floating point rules: 1.0-0.8-0.1-0.1 is not 0.0 but -1e-14.
            assert (f >= -1e-12)

        if fraction_type == 'molar':
            self.molar_fractions=fractions
            total = sum(fractions)
            if abs(total - 1.0) > 1e-12:
                warnings.warn("List of molar fractions sums to %g. Normalizing." % total, stacklevel=2)
                self.molar_fractions = [fr / total for fr in fractions]
            try:
                del self.mass_fractions
            except AttributeError:
                pass
            try:
                del self.volume_fractions
            except AttributeError:
                pass
        elif fraction_type == 'mass': # only if set_composition has been called
            try:
                self.mass_fractions=fractions
                self.molar_fractions = self.mass_to_molar_fractions(fractions, self.phases)
                try:
                    del self.volume_fractions
                except AttributeError:
                    pass
            except AttributeError:
                raise Exception("Mass fractions cannot be set before set_composition for all the phases in the assemblage.")
        elif fraction_type == 'volume': # only if set_composition and set_state have been called
            try:
                self.volume_fractions=fractions
                self.molar_fractions = self.volume_to_molar_fractions(fractions, self.phases)
                try:
                    del self.mass_fractions
                except AttributeError:
                    pass
            except AttributeError:
                raise Exception("Volume fractions cannot be set before set_state.") 
        else:
            raise Exception("Fraction type not recognised. Please use 'molar', 'mass' or 'volume'")

        # Make sure all fractions >= 0.0
        self.molar_fractions = [max(0.0, fr) for fr in self.molar_fractions]

        try:
            self.composition=self.assemblage_composition() # if set_composition has been called
        except AttributeError:
            pass

    def set_composition(self, phase_compositions):
        assert(len(self.phases)==len(phase_compositions))
        for i, composition in enumerate(phase_compositions):
            if composition != []:
                self.phases[i].set_composition(composition)

	if hasattr(self, 'molar_fractions'):
            self.composition=self.assemblage_composition()

    def assemblage_composition(self):
        self.phase_compositions=[self.phases[i].composition for i in range(len(self.phases))]
        self.composition=dict()
        for i, composition in enumerate(self.phase_compositions):
            for element in composition:
                if element not in self.composition:
                    self.composition[element] = self.molar_fractions[i]*composition[element]
                else:
                    self.composition[element] += self.molar_fractions[i]*composition[element]
        return self.composition

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for phase in self.phases:
            phase.set_method(method)

    def unroll(self):
        fractions = []
        minerals = []

        for idx, phase in enumerate(self.phases):
            fraction=self.molar_fractions[idx]
            p_fr,p_min = phase.unroll()
            check_pairs(p_fr,p_min)
            fractions.extend([i*fraction for i in p_fr])
            minerals.extend(p_min)
        return (fractions, minerals)

    def to_string(self):
        """
        Return the name of the assemblage
        """
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature, set_type="all", averaging_scheme=burnman.averaging_schemes.VoigtReussHill()):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        self.pressure = pressure
        self.temperature = temperature
        for phase in self.phases:
            phase.set_state(pressure, temperature)

        if set_type=="none":
            pass
        elif set_type=="gibbs_only":
            self.assemblage_gibbs()
        elif set_type=="elastic":
            self.assemblage_elastic_properties(averaging_scheme)
        elif set_type=="thermodynamic":
            self.assemblage_thermodynamic_properties()
        elif set_type=="all":
            self.assemblage_thermodynamic_properties()
            self.assemblage_elastic_properties(averaging_scheme)
	else:
            raise Exception("Type for set_state not recognised.")
    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        densities = np.array([ph.density() for (_,ph) in self.children])
        volumes = np.array([ph.molar_volume()*fraction for (fraction, ph) in self.children])
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
        (molar_fractions,minerals) = self.unroll()
        volume_fractions = self.molar_to_volume_fractions(molar_fractions, minerals)

        for (volume_fraction,mineral) in zip(volume_fractions,minerals):

            e = burnman.ElasticProperties()
            e.V = volume_fraction * mineral.molar_volume()
            e.K = mineral.adiabatic_bulk_modulus()
            e.G = mineral.shear_modulus()
            e.rho = mineral.molar_mass() / mineral.molar_volume()
            e.fraction = volume_fraction
            elastic_properties_of_phases.append(e)

        return elastic_properties_of_phases

    def molar_to_volume_fractions(self, molar_fractions, minerals):
        total_volume=0.
        for i, mineral in enumerate(minerals):
            total_volume+=molar_fractions[i]*mineral.V
        volume_fractions=[]
        for i, mineral in enumerate(minerals):
            volume_fractions.append(molar_fractions[i]*mineral.V/total_volume)
        return volume_fractions

    def volume_to_molar_fractions(self, volume_fractions, minerals):
        total_moles=0.
        for i, mineral in enumerate(minerals):
            total_moles+=volume_fractions[i]/mineral.V
        molar_fractions=[]
        for i, mineral in enumerate(minerals):
            molar_fractions.append(volume_fractions[i]/(mineral.V*total_moles))
        return molar_fractions

    def mass_to_molar_fractions(self, mass_fractions, minerals):
        total_moles=0.
        for i, mineral in enumerate(minerals):
            total_moles+=mass_fractions[i]/mineral.molar_mass()
        molar_fractions=[]
        for i, mineral in enumerate(minerals):
            molar_fractions.append(mass_fractions[i]/(mineral.molar_mass()*total_moles))
        return molar_fractions

    def assemblage_thermodynamic_properties(self):
        """
        Calculates the thermodynamic properties of the assemblage
        :returns: :math:`\mathcal{G}` :math:`[J/mol]`, :math:`\mathcal{H}` :math:`[J/mol]` , :math:`\mathcal{S}` :math:`[J/K/mol]`, :math:`V` :math:`[m^3/mol]`, :math:`C_p` :math:`[J/K/mol]`, :math:`C_v` :math:`[J/K/mol]`, :math:`alpha` :math:`[1/K]`, :math:`\gamma` :math:`[]`,
        :rtype: lists of floats
        """
        thermodynamic_properties_of_phases = []
        (fractions,minerals) = self.unroll()
        n_minerals=len(minerals)
        self.gibbs = sum([ minerals[i].gibbs * fractions[i] for i in range(n_minerals) ])
        self.H = sum([ minerals[i].H * fractions[i] for i in range(n_minerals) ])
        self.S = sum([ minerals[i].S * fractions[i] for i in range(n_minerals) ])
        self.V = sum([ minerals[i].V * fractions[i] for i in range(n_minerals) ])
        self.C_p = 0. #sum([ minerals[i].C_p * fractions[i] for i in range(n_minerals) ])
        self.C_v = 0. #sum([ minerals[i].C_v * fractions[i] for i in range(n_minerals) ])
        self.alpha = 0. #sum([ minerals[i].alpha * fractions[i] for i in range(n_minerals) ])
        self.gr = 0. #sum([ minerals[i].gr * fractions[i] for i in range(n_minerals) ])
        return self.gibbs, self.H, self.S, self.V, self.C_p, self.C_v, self.alpha, self.gr

    def assemblage_gibbs(self):
        """
        Calculates the molar gibbs free energy of the assemblage
        :returns: :math:`\mathcal{G}` :math:`[J/mol]`
        """
        thermodynamic_properties_of_phases = []
        (fractions,minerals) = self.unroll()
        n_minerals=len(minerals)
        self.gibbs = sum([ minerals[i].gibbs * fractions[i] for i in range(n_minerals) ])
        return self.gibbs


    def assemblage_elastic_properties(self, averaging_scheme=burnman.averaging_schemes.VoigtReussHill()):
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
        Returns molar mass of the assemblage [kg/mol]
        """
        molar_mass = sum([ self.phases[i][0].molar_mass()*self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return molar_mass
