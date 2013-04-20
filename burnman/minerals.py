# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from composition import *
import numpy as np
import voigt_reuss_hill as vrh
import slb_finitestrain as slb
import mie_grueneisen_debye as mgd
import slb_thirdorder as slb3
import birch_murnaghan as bm
import warnings


class material:
    """
    This is the base class for all minerals. States of the mineral
    can only be queried after setting the pressure and temperature
    using set_state(). The method for computing properties of
    the material is set using set_method(), which should be done
    once after creating the material.
    
    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the material parameters are expected to be in plain SI units.  This
    means that the elastic moduli should be in Pascals and NOT Gigapascals,
    and the Debye temperature should be in K not C.  Additionally, the
    reference volume should be in m^3/(mol molecule) and not in unit cell
    volume and 'n' should be the number of atoms per molecule.  Frequently in
    the literature the reference volume is given in Angstrom^3 per unit cell.
    To convert this to m^3/(mol molecule) you should multiply by 10^(-30) *
    N_a / Z, where N_a is Avogadro's number and Z is the number of atoms per
    unit cell.  You can look up Z in many places, including www.mindat.org
    """

    def __init__(self):
        self.params = {    'name':'generic',
            'equation_of_state': 'slb3', #Equation of state used to fit the parameters
            'ref_V': 0., #Molar volume [m^3/(mole molecules)] at room pressure/temperature
            'ref_K': 0., #Reference bulk modulus [Pa] at room pressure/temperature
            'K_prime': 0., #pressure derivative of bulk modulus
            'ref_mu': 0., #reference shear modulus at room pressure/temperature
            'mu_prime': 0., #pressure derivative of shear modulus
            'molar_mass': 0., #molar mass in units of [kg/mol]
            'n': 0., #number of atoms per molecule
            'ref_Debye': 0., #Debye temperature for material. See Stixrude & Lithgow-Bertelloni, 2005 for values 
            'ref_grueneisen': 0., #Gruneisen parameter for material. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'q0': 0., #q value used in caluclations. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'eta_0s': 0.0} #eta value used in calculations. See Stixrude & Lithgow-Bertelloni, 2005 for values
#        self.pressure = 0.0
#        self.temperature = 300
        self.method = []

    def set_method(self, method):
        """ use "slb" or "mgd" or slb3 """
        if (method == "slb"):
            self.method = slb;
        elif (method == "mgd"):
            self.method = mgd;
        elif (method == "slb3"):
            self.method = slb3;
        elif (method == "bm"):
            self.method = bm;
        else:
            raise("unsupported material method " + method)

    def to_string(self):
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        """ Update the material to the given pressure [Pa] and temperature.
        
        This updates the other properties of this class (v_s, v_p, ...).
        """
        #in an effort to avoid additional work, don't do all the calculations if nothing has changed
        try:
            if self.pressure == pressure and self.temperature == temperature and self.old_params == self.params:
                return
        except AttributeError:
            pass  #do nothing

        self.pressure = pressure
        self.temperature = temperature
        self.old_params = self.params
                
        if (self.method == bm):
            self.V = self.method.volume(self.pressure, self.params)
            self.K_T = self.method.bulk_modulus(self.V, self.params)
            self.K_S = self.K_T  # since you are ignoring temperature in this method
            if (self.params.has_key('ref_mu') and self.params.has_key('mu_prime')):
                self.mu = self.method.shear_modulus(self.V, self.params)
            else:
                self.mu = -1
                warnings.warn(('Warning: mu and or mu_prime are undefined for ' + self.to_string()))
        else:
            self.V = self.method.volume(self.pressure, self.temperature, self.params)
            self.K_T = self.method.bulk_modulus(self.temperature, self.V, self.params)
            self.K_S = self.method.bulk_modulus_adiabatic(self.temperature, self.V, self.params)
            self.C_v = self.method.heat_capacity_v(self.temperature,self.V, self.params)
            self.C_p = self.method.heat_capacity_p(self.temperature,self.V, self.params)
            self.alpha = self.method.thermal_expansivity(self.temperature, self.V, self.params)
            self.gr = self.method.grueneisen_parameter(self.params['ref_V']/self.V, self.params)
            if (self.params.has_key('ref_mu') and self.params.has_key('mu_prime')):
                self.mu = self.method.shear_modulus(self.temperature, self.V, self.params)
            else:    
                self.mu = -1
                warnings.warn(('Warning: mu and or mu_prime are undefined for ' + self.to_string()))

    def molar_mass(self):
        return self.params['molar_mass']

    def density(self):
        return  self.params['molar_mass'] / self.V

    # gives volume in m^3/mol
    def molar_volume(self):
        return self.V
    def grueneisen_parameter(self):
        return self.gr
    def bulk_modulus(self):
        return self.K_T
    def adiabatic_bulk_modulus(self):
        return self.K_S
    def thermal_expansivity(self):
        return self.alpha
    def heat_capacity_v(self):
        return self.C_v
    def heat_capacity_p(self):
        return self.C_p
    def shear_modulus(self):
        return self.mu
    def v_s(self):
        return np.sqrt(self.shear_modulus() / \
            self.density()) / 1000.
    def v_p(self):
        return np.sqrt((self.bulk_modulus() + 4. / 3. * \
            self.shear_modulus()) / self.density()) / 1000.

# combines two or more materials given a fixed molar_fraction
class helper_volumetric_mixing(material):
    # base_materials: list of materials
    # molar_fraction: list of molar ratios (sum up to 1)
    def __init__(self, base_materials, molar_fraction):
        self.base_materials = base_materials
        self.molar_fraction = molar_fraction
        assert(len(base_materials) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.999)
        assert(sum(molar_fraction) < 1.001)

    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.method = self.method
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_materials))


        self.params = {}

        self.params['n'] = sum([self.base_materials[i].params['n'] for i in itrange ]) / len(self.base_materials)
        phase_volume = [self.base_materials[i].molar_volume() * self.molar_fraction[i] for i in itrange ]

        # some properties need weighted averaging
        for prop in ['ref_V', 'molar_mass', 'ref_Debye', 'ref_grueneisen', 'q0', 'eta_0s']:
            self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])


        # some need VRH averaging
        for prop in ['ref_K', 'K_prime', 'ref_mu', 'mu_prime']:

            X = [ mat.params[prop] for mat in self.base_materials ]
            self.params[prop] = vrh.vrh_average(phase_volume, X)

        material.set_state(self, pressure, temperature)

class helper_spin_transition(material):
    """ helper class that switches between two materials (for low and high spin) 
    based on pressure
    """
    
    def __init__(self, lowspeed_pressure, ls_mat, hs_mat):
        """ if pressure>=lowspeed_pressure [Pa], use ls_mat, else hs_mat """
        self.lowspeed_pressure = lowspeed_pressure
        self.ls_mat = ls_mat
        self.hs_mat = hs_mat
                
    def set_state(self, pressure, temperature):
        if (pressure >= self.lowspeed_pressure):
            mat = self.ls_mat
        else:
            mat = self.hs_mat
            
        mat.method = self.method
        mat.set_state(pressure, temperature)
        self.params = mat.params
        material.set_state(self, pressure, temperature)                



################ Built-in Minerals
class test_mineral (material):
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'ref_V': 6.844e-6,
            'ref_K': 135.19e9,
            'K_prime': 6.04,
            'ref_mu': 175.0e9,
            'mu_prime': 1.7,
            'molar_mass': .055845,
            'n': 1,
            'ref_Debye': 998.85,
            'ref_grueneisen': 1.368,
            'q0': 0.917}

class stishovite (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'ref_V': 14.02e-6,
            'ref_K': 314.0e9,
            'K_prime': 4.4,
            'ref_mu': 220.0e9,
            'mu_prime': 1.6,
            'molar_mass': .0601,
            'n': 3,
            'ref_Debye': 1044.,
            'ref_grueneisen': 134,
            'q0': 2.4,
            'eta_0s': 5.0 }

class periclase (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.24e-6,
            'ref_K': 161.0e9,
            'K_prime': 3.9,
            'ref_mu': 130.0e9,
            'mu_prime': 2.2,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 773.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 2.3 }

class wuestite (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 12.06e-6,
            'ref_K': 152.0e9,
            'K_prime': 4.9,
            'ref_mu': 47.0e9,
            'mu_prime': 0.7,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 455.,
            'ref_grueneisen': 1.28,
            'q0': 1.5,
            'eta_0s': 0.8 }



class ferropericlase(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [periclase(), wuestite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_volumetric_mixing.__init__(self, base_materials, molar_fraction)



class mg_fe_perovskite(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_volumetric_mixing.__init__(self, base_materials, molar_fraction)


class mg_perovskite(material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein  
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 24.45e-6,
            'ref_K': 251.0e9,   
            'K_prime': 4.1,     
            'ref_mu': 175.0e9,  
            'mu_prime': 1.7,  
            'molar_mass': .1020,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.6 }

class fe_perovskite(material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 25.48e-6,
            'ref_K': 281.0e9, 
            'K_prime': 4.1,  
            'ref_mu': 138.0e9,
            'mu_prime': 1.7,   
            'molar_mass': .1319, 
            'n': 5,
            'ref_Debye': 841.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.1 }

class Matas_mg_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 24.43e-6,
            'ref_K': 250.0e9,   
            'K_prime': 4.0,     
            'ref_mu': 175.0e9,  
            'mu_prime': 1.8,    
            'molar_mass': .1020,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4} 

class Matas_fe_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 25.34e-6,
            'ref_K': 250.0e9, 
            'K_prime': 4.0,  
            'ref_mu': 135.0e9, 
            'mu_prime': 1.3,  
            'molar_mass': .1319, 
            'n': 5,
            'ref_Debye': 841.,
            'ref_grueneisen': 1.48,
            'q0': 1.4} 

class Matas_periclase (material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 11.25e-6,
            'ref_K': 160.1e9,
            'K_prime': 3.83,
            'ref_mu': 130.0e9,
            'mu_prime': 2.2,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 673.,
            'ref_grueneisen': 1.41,
            'q0': 1.3 }

class Matas_wuestite (material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 12.26e-6,
            'ref_K': 160.1e9,
            'K_prime': 3.83,
            'ref_mu': 46.0e9,
            'mu_prime':  0.6,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 673.,
            'ref_grueneisen': 1.41,
            'q0': 1.3 }


class Speziale_fe_periclase(helper_spin_transition):  
    def __init__(self):
        helper_spin_transition.__init__(self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
        self.cite = 'Speziale et al. 2007'

class Speziale_fe_periclase_HS(material):
    """
    Speziale et al. 2007, Mg#=83
    """ 
    def __init__(self):
            self.params = {
                        'equation_of_state': 'mgd3',
                        'ref_V': 22.9e-6,
                        'ref_K': 157.5e9,
                        'K_prime': 3.92,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2 }

class Speziale_fe_periclase_LS(material): 
    """
    Speziale et al. 2007, Mg#=83
    """
    def __init__(self):
        self.params = {
                        'equation_of_state': 'mgd3',
                        'ref_V': 21.49e-6,
                        'ref_K': 186.0e9,
                        'K_prime': 4.6,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587.,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2  }

    
    
        
class Speziale_fe_periclase(helper_spin_transition): 
    def __init__(self):        
        helper_spin_transition.__init__(self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
        self.cite = 'Speziale et al. 2007'

        


class Murakami_mg_perovskite(material):  
    """
    Murakami et al. (2012) supplementary table 5 and references therein, ref_V from Stixrude & Lithgow-Bertolloni 2005
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 24.45e-6,  # stixrud & L-B 2005
            'ref_K': 251.9e9,
            'K_prime': 4.01,
            'ref_mu': 164.7e9,
            'mu_prime': 1.58,
            'molar_mass': .102165,
            'n': 5,
            'ref_Debye': 1054.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 } 

class Murakami_fe_perovskite(material): 
    """
    Murakami et al. (2012), personal communication, Mg#=94, Al=4%
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 24.607e-6,
            'ref_K': 251.9e9, 
            'K_prime': 4.01,
            'ref_mu': 164.7e9,
            'mu_prime': 1.58,
            'molar_mass': .102165,
            'n': 5,
            'ref_Debye': 1054.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 1.48 } 
            
class Murakami_mg_periclase(material):
    """
    Murakami et al. (2012) supplementary table 5 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 24.607e-6,
            'ref_K': 251.9e9,
            'K_prime': 4.01,
            'ref_mu': 164.7e9,
            'mu_prime': 1.58,
            'molar_mass': .102165,
            'n': 5,
            'ref_Debye': 1054.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 }             

class Murakami_fe_periclase(helper_spin_transition):
    def __init__(self):
        helper_spin_transition.__init__(self, 63.0e9, Murakami_fe_periclase_LS(), Murakami_fe_periclase_HS())

class Murakami_fe_periclase_HS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 11.412e-6,
            'ref_K': 159.1e9,
            'K_prime': 4.11,
            'ref_mu': 105.43e9,
            'mu_prime': 1.773,
            'molar_mass': .0494,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54 }

class Murakami_fe_periclase_LS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 11.171e-6,
            'ref_K': 170.0e9,
            'K_prime': 4,
            'ref_mu': 116.34e9,
            'mu_prime': 1.668,
            'molar_mass': .0494,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54}




                
class helper_fe_dependent(material):
    def __init__(self, iron_number_with_pt, idx):
        self.iron_number_with_pt = iron_number_with_pt
        self.which_index = idx  # take input 0 or 1 from iron_number_with_pt()

    def create_inner_material(self, iron_number):
        return []  # needs to be overwritten in class deriving from this one

    def set_state(self, pressure, temperature):
        self.pressure = pressure
        self.temperature = temperature
        self.base_material = self.create_inner_material(self.iron_number())
        self.base_material.method = self.method
        self.base_material.set_state(pressure, temperature)
        self.params = self.base_material.params
        material.set_state(self, pressure, temperature)

    def iron_number(self):
        return self.iron_number_with_pt(self.pressure, self.temperature)[self.which_index]
    def molar_mass(self):
        return self.base_material.molar_mass()
    def density(self):
        return self.base_material.density()
    def molar_volume(self):
        return self.base_material.molar_volume()
    def bulk_modulus(self):
        return self.base_material.bulk_modulus()
    def v_s(self):
        return self.base_material.v_s()
    def v_p(self):
        return self.base_material.v_p()
    def geotherm(self):
        return self.base_material.v_s()


class mg_fe_perovskite_pt_dependent(helper_fe_dependent):
    def __init__(self, iron_number_with_pt, idx):
        helper_fe_dependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return mg_fe_perovskite(iron_number)

class ferropericlase_pt_dependent(helper_fe_dependent):
    def __init__(self, iron_number_with_pt, idx):
        helper_fe_dependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return ferropericlase(iron_number)
