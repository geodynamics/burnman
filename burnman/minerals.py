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
    """
    def __init__(self):
        self.params = {    'name':'generic',
            'ref_V': 0.,
            'ref_K': 0.,
            'K_prime': 0.,
            'ref_mu': 0.,
            'mu_prime': 0.,
            'molar_mass': 0.,
            'n': 0.,
            'ref_Debye': 0.,
            'ref_grueneisen': 0.,
            'q0': 0.}
        self.pressure = 0.0
        self.temperature = 300
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
        """ Update the material to the given pressure [GPa] and temperature.
        
        This updates the other properties of this class (v_s, v_p, ...).
        """
        
        self.pressure = pressure
        self.temperature = temperature
                
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
        return np.sqrt(self.shear_modulus() * 1.e9 / \
            self.density()) / 1000.
    def v_p(self):
        return np.sqrt((self.bulk_modulus() * 1.e9 + 4. / 3. * \
            self.shear_modulus()) / self.density()) / 1000.


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
            'ref_V': 6.844e-6,
            'ref_K': 135.19,
            'K_prime': 6.04,
            'ref_mu': 175.,
            'mu_prime': 1.7,
            'molar_mass': .055845,
            'n': 1,
            'ref_Debye': 998.85,
            'ref_grueneisen': 1.368,
            'q0': 0.917}

class stishovite (material):
    def __init__(self):
        self.params = {
            'ref_V': 14.02e-6,
            'ref_K': 314.,
            'K_prime': 4.4,
            'ref_mu': 220.,
            'mu_prime': 1.6,
            'molar_mass': .0601,
            'n': 3,
            'ref_Debye': 1044.,
            'ref_grueneisen': 1.6,
            'q0': 2.4,
            'eta_0s': 3.0 }

class periclase (material):
    """ 
    MgO
    """
    def __init__(self):
        self.params = {
            'ref_V': 11.24e-6,
            'ref_K': 161.,
            'K_prime': 3.9,
            'ref_mu': 130.,
            'mu_prime': 2.2,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 773.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 3.0 }
class wustite (material):
    """
    FeO
    """
    def __init__(self):
        self.params = {
            'ref_V': 12.06e-6,
            'ref_K': 152.,
            'K_prime': 4.9,
            'ref_mu': 47.,
            'mu_prime': 0.7,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 455.,
            'ref_grueneisen': 1.28,
            'q0': 1.5,
            'eta_0s': 3.0 }

# combines two or more materials given a fixed molar_abundance
# based on their volume at a certain T,p.
class helper_volumetric_mixing(material):
    # base_materials: list of materials
    # molar_abundances: list of molar ratios (sum up to 1)
    def __init__(self, base_materials, molar_abundances):
        self.base_materials = base_materials
        self.molar_abundances = molar_abundances
        assert(len(base_materials) == len(molar_abundances))
        assert(sum(molar_abundances) > 0.999)
        assert(sum(molar_abundances) < 1.001)

    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.method = self.method
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_materials))


        self.params = {}
        
        self.params['n'] = sum([self.base_materials[i].params['n'] for i in itrange ]) / len(self.base_materials)            
        phase_volume = [self.base_materials[i].molar_volume() * self.molar_abundances[i] for i in itrange ]

        # some properties need weighted averaging
        for prop in ['ref_V', 'molar_mass', 'ref_Debye', 'ref_grueneisen', 'q0', 'eta_0s']:
            self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_abundances[i] for i in itrange ])


        # some need VRH averaging
        for prop in ['ref_K', 'K_prime', 'ref_mu', 'mu_prime']:
            
            X = [ mat.params[prop] for mat in self.base_materials ]
            self.params[prop] = vrh.vhr_average(phase_volume, X)

        material.set_state(self, pressure, temperature)
        


class ferropericlase(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [periclase(), wustite()]
        molar_abundances = (1. - fe_num, fe_num)
        helper_volumetric_mixing.__init__(self, base_materials, molar_abundances)


class mg_fe_perovskite(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_abundances = (1. - fe_num, fe_num)
        helper_volumetric_mixing.__init__(self, base_materials, molar_abundances)


class ferropericlase_old_and_probably_wrong(material):
    def __init__(self, fe_num):
        self.mg = 1. - fe_num
        self.fe = fe_num
        self.pe = periclase()
        self.wu = wustite()
        self.params = {
            'ref_V': self.pe.params['ref_V'] * self.mg + self.wu.params['ref_V'] * self.fe,
            'ref_K': self.pe.params['ref_K'] * self.mg + self.wu.params['ref_K'] * self.fe,
            'K_prime': self.pe.params['K_prime'] * self.mg + self.wu.params['K_prime'] * self.fe,
            'ref_mu': self.pe.params['ref_mu'] * self.mg + self.wu.params['ref_mu'] * self.fe,
            'mu_prime': self.pe.params['mu_prime'] * self.mg + self.wu.params['mu_prime'] * self.fe,
            'molar_mass': self.pe.params['molar_mass'] * self.mg + self.wu.params['molar_mass'] * self.fe,
            'n': 5,
            'ref_Debye': self.pe.params['ref_Debye'] * self.mg + self.wu.params['ref_Debye'] * self.fe,
            'ref_grueneisen': self.pe.params['ref_grueneisen'] * self.mg + self.wu.params['ref_grueneisen'] * self.fe,
            'q0': self.pe.params['q0'] * self.mg + self.wu.params['q0'] * self.fe ,
            'eta_0s': self.pe.params['eta_0s'] * self.mg + self.wu.params['eta_0s'] * self.fe }



class mg_fe_perovskite_old_and_probably_wrong(material):
    def __init__(self, fe_num):
        self.mg = 1.0 - fe_num
        self.fe = fe_num
        self.mg_pv = mg_perovskite()
        self.fe_pv = fe_perovskite()
        self.params = {
            'ref_V': self.mg_pv.params['ref_V'] * self.mg + self.fe_pv.params['ref_V'] * self.fe,
            'ref_K': self.mg_pv.params['ref_K'] * self.mg + self.fe_pv.params['ref_K'] * self.fe,
            'K_prime': self.mg_pv.params['K_prime'] * self.mg + self.fe_pv.params['K_prime'] * self.fe,
            'ref_mu': self.mg_pv.params['ref_mu'] * self.mg + self.fe_pv.params['ref_mu'] * self.fe,
            'mu_prime': self.mg_pv.params['mu_prime'] * self.mg + self.fe_pv.params['mu_prime'] * self.fe,
            'molar_mass': self.mg_pv.params['molar_mass'] * self.mg + self.fe_pv.params['molar_mass'] * self.fe,
            'n': 5,
            'ref_Debye': self.mg_pv.params['ref_Debye'] * self.mg + self.fe_pv.params['ref_Debye'] * self.fe,
            'ref_grueneisen': self.mg_pv.params['ref_grueneisen'] * self.mg + self.fe_pv.params['ref_grueneisen'] * self.fe,
            'q0': self.mg_pv.params['q0'] * self.mg + self.fe_pv.params['q0'] * self.fe ,
            'eta_0s': self.mg_pv.params['eta_0s'] * self.mg + self.fe_pv.params['eta_0s'] * self.fe}

class mg_perovskite(material):
    def __init__(self):
        self.params = {
            'ref_V': 24.45e-6,
            'ref_K': 251.,
            'K_prime': 4.1,
            # 'ref_mu': 166.,
            # 'mu_prime': 1.57,
                        'ref_mu': 175.,  # from S & L.-B. 2005
                        'mu_prime': 1.8,
            'molar_mass': .1020,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 }

class fe_perovskite(material):
    def __init__(self):
        self.params = {
            'ref_V': 25.48e-6,
            'ref_K': 281.,
            'K_prime': 4.1,
            'ref_mu': 161.,
            'mu_prime': 1.57,
            'molar_mass': .1319,
            'n': 5,
            'ref_Debye': 1021.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 }
        
class Speziale_fe_periclase_LS(material):  # Speciale et al. 2007
    def __init__(self):
        self.params = {
                        'P_LS': 60,  # in GPa
                        'ref_V': 42.99e-5,
                        'ref_K': 186.,
                        'K_prime': 4.6,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587.,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2,
                        'eta_0s': 2.4 }  # not from paper

    
    
class Speziale_fe_periclase_HS(material):  # Speciale et al. 2007
    def __init__(self):    
            self.params = {
                        'ref_V': 45.83e-5,
                        'ref_K': 157.5,
                        'K_prime': 3.92,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2,
                        'eta_0s': 2.4 }  # not from paper
    
        
class Speziale_fe_periclase(helper_spin_transition):  # Speciale et al. 2007
    def __init__(self):        
        helper_spin_transition.__init__(self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
        self.cite = 'Speziale et al. 2007'

        
class Catalli_fe_perovskite_LS(material):  # Catalli et al 2009
    def __init__(self):        
        self.params = {
                        'ref_V': 165.78e-6,
                        'ref_K': 304.,
                        'K_prime': 4.5,
                        'ref_mu': 161.,
                        'mu_prime': 1.57,
                        'molar_mass': .1319,
                        'n': 5,
                        'ref_Debye': 1021.,
                        'ref_grueneisen': 1.48,
                        'q0': 1.4,
                        'eta_0s': 2.4 }

class Catalli_fe_perovskite_HS(material):  # Catalli et al 2009
    def __init__(self):    
        self.params = {
                        'ref_V': 165.78e-6,
                        'ref_K': 237.,
                        'K_prime': 4.5,
                        'ref_mu': 161.,
                        'mu_prime': 1.57,
                        'molar_mass': .1319,
                        'n': 5,
                        'ref_Debye': 1021.,
                        'ref_grueneisen': 1.48,
                        'q0': 1.4,
                        'eta_0s': 2.4 }    
        
class Catalli_fe_perovskite(helper_spin_transition):  # Catalli et al 2009
    def __init__(self):        
        helper_spin_transition.__init__(self, 63.0e9, Catalli_fe_perovskite_LS(), Catalli_fe_perovskite_HS())

        
class Murakami_fe_perovskite(material):  # From Murakami's emails, see Cayman for details, represents 4 wt% Al X_mg = .94
    def __init__(self):
        self.params = {
            'ref_V': 24.607e-6,
            'ref_K': 251.9,
            'K_prime': 4.01,
            'ref_mu': 164.7,
            # 'ref_mu': 157.39,  #refitted to second order  
                    'mu_prime': 1.58,
            # 'mu_prime': 2.08, #refitted to second order
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
    def __init__(self):
        self.params = {
            'ref_V': 11.412e-6,
            'ref_K': 159.,
            'K_prime': 4.11,
            'ref_mu': 105.43,
            'mu_prime': 1.773,
            'molar_mass': .0494,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 3.0 }

class Murakami_fe_periclase_LS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    def __init__(self):
        self.params = {
            'ref_V': 11.412e-6,  # 11.171e-6,modified by Sanne
            'ref_K': 159.,  # 170.,
            'K_prime': 4.11,  # 4,
            'ref_mu': 116.34,
            'mu_prime': 1.668,
            'molar_mass': .0494,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 3.0}



                
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
