# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import voigt_reuss_hill as vrh
import mie_grueneisen_debye as mgd
import birch_murnaghan as bm
import warnings
import slb
import equation_of_state as eos


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
        self.method = None

    def set_method(self, method):
        """ use "slb" or "mgd" or slb3 """
        if( isinstance(method, basestring)):
            if (method == "slb2"):
                self.method = slb.slb2()
            elif (method == "mgd2"):
                self.method = mgd.mgd2()
            elif (method == "mgd3"):
                self.method = mgd.mgd3()
            elif (method == "slb3"):
                self.method = slb.slb3()
            elif (method == "bm2"):
                self.method = bm.bm2()
            elif (method == "bm3"):
                self.method = bm.bm3()
            else:
                raise Exception("unsupported material method " + method)
        elif ( isinstance(method, eos.equation_of_state) ):
            self.method = method
        else:
            raise Exception("unsupported material method " + method.__class__.__name__ )

    def to_string(self):
        return "'" + self.__class__.__module__ + "." + self.__class__.__name__ + "'"

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
        
        self.V = self.method.volume(self.pressure, self.temperature, self.params)
        self.gr = self.method.grueneisen_parameter(self.pressure, self.temperature, self.V, self.params)
        self.K_T = self.method.isothermal_bulk_modulus(self.pressure, self.temperature, self.V, self.params)
        self.K_S = self.method.adiabatic_bulk_modulus(self.pressure, self.temperature, self.V, self.params)
        self.C_v = self.method.heat_capacity_v(self.pressure, self.temperature, self.V, self.params)
        self.C_p = self.method.heat_capacity_p(self.pressure, self.temperature, self.V, self.params)
        self.alpha = self.method.thermal_expansivity(self.pressure, self.temperature, self.V, self.params)
        
        if (self.params.has_key('ref_mu') and self.params.has_key('mu_prime')):
            self.mu = self.method.shear_modulus(self.pressure, self.temperature, self.V, self.params)
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
            self.density()) 
    def v_p(self):
        return np.sqrt((self.bulk_modulus() + 4. / 3. * \
            self.shear_modulus()) / self.density()) 

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


