# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
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
            'ref_G': 0., #reference shear modulus at room pressure/temperature
            'G_prime': 0., #pressure derivative of shear modulus
            'molar_mass': 0., #molar mass in units of [kg/mol]
            'n': 0., #number of atoms per molecule
            'ref_Debye': 0., #Debye temperature for material. See Stixrude & Lithgow-Bertelloni, 2005 for values 
            'ref_grueneisen': 0., #Gruneisen parameter for material. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'q0': 0., #q value used in caluclations. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'eta_0s': 0.0} #eta value used in calculations. See Stixrude & Lithgow-Bertelloni, 2005 for values
        self.method = None

    def set_method(self, method):
        """
        Set the equation of state to be used for this mineral.
        Takes a string corresponding to any of the predefined
        equations of state:  'bm2', 'bm3', 'mgd2', 'mgd3', 'slb2',
        or 'slb3'.  Alternatively, you can pass a user defined
        class which derives from the equation_of_state base class.
        """
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
        elif ( issubclass(method, eos.equation_of_state) ):
            self.method = method()
        else:
            raise Exception("unsupported material method " + method.__class__.__name__ )

    def to_string(self):
        """
        Returns the name of the mineral class
        """
        return "'" + self.__class__.__module__.replace(".minlib_",".") + "." + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        
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
        
        if (self.params.has_key('ref_G') and self.params.has_key('G_prime')):
            self.G = self.method.shear_modulus(self.pressure, self.temperature, self.V, self.params)
        else:    
            self.G = float('nan') #nan if there is no G, this should propagate through calculations to the end
            warnings.warn(('Warning: G and or G_prime are undefined for ' + self.to_string()))

    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        return self.params['molar_mass']

    def density(self):
        """
        Returns density of the mineral [kg/m^3]
        """
        return  self.params['molar_mass'] / self.V

    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        """
        return self.V
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
        """
        return self.gr
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        """
        return self.K_T
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        """
        return self.K_S
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        """
        return self.G
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the mineral [1/K]
        """
        return self.alpha
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        """
        return self.C_v
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        """
        return self.C_p
    def v_s(self):
        """
        Returns shear wave speed of the mineral [m/s]
        """
        return np.sqrt(self.shear_modulus() / \
            self.density()) 
    def v_p(self):
        """
        Returns P wave speed of the mineral [m/s]
        """
        return np.sqrt((self.adiabatic_bulk_modulus() + 4. / 3. * \
            self.shear_modulus()) / self.density()) 
    def v_phi(self):
        """
        Returns bulk sound speed of the mineral [m/s]
        """
        return np.sqrt(self.adiabatic_bulk_modulus() / self.density())

# combines two or more materials given a fixed molar_fraction
class helper_solid_solution(material):
    """
    Class for coming up with a new mineral based based on a solid
    solution between two or more end member minerals.  It is not
    completely clear how to do this, or how valid this approximation
    is, but here we just do a weighted arithmetic average of the 
    thermoelastic properties of the end members according to their molar fractions
    """
    def __init__(self, base_materials, molar_fraction):
        """
        Takes a list of end member minerals, and a matching list of
        molar fractions of those minerals for mixing them.  Simply
        comes up with a new mineral by doing a weighted arithmetic
        average of the end member minerals
        """
        self.base_materials = base_materials
        self.molar_fraction = molar_fraction
        assert(len(base_materials) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        #does not make sense to do a solid solution with different number of 
        #atoms per formula unit, at least not simply...
        for m in base_materials:
            if(base_materials[0].params.has_key('n')):
                assert(m.params['n'] == base_materials[0].params['n'])

    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.method = self.method
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_materials))

        self.params = {}

        # some do arithmetic averaging of the end members
        for prop in self.base_materials[0].params:
           try:
               self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
               #if there is a type error, it is probably a string.  Just go with the value of the first base_material.
               self.params[prop] = self.base_materials[0].params[prop]
        material.set_state(self, pressure, temperature)
        print self.params

class helper_spin_transition(material):
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
        self.transition_pressure = transition_pressure
        self.ls_mat = ls_mat
        self.hs_mat = hs_mat
                
    def set_state(self, pressure, temperature):
        if (pressure >= self.transition_pressure):
            mat = self.ls_mat
        else:
            mat = self.hs_mat
            
        mat.method = self.method
        mat.set_state(pressure, temperature)
        self.params = mat.params
        material.set_state(self, pressure, temperature)                


class helper_uncertainty(material):
   """
   incorporates uncertainties in mineral parameters
   """

   def __init__(self,perturbations):
       self.params['ref_K']=self.params['ref_K']+self.params['err_ref_K']*(1.0+perturbations[0])
       self.params['K_prime']=self.params['K_prime']+self.params['err_K_prime']*(1.0+perturbations[1])
       self.params['ref_G']=self.params['ref_G']+self.params['err_ref_G']*(1.0+perturbations[2])
       self.params['G_prime']=self.params['G_prime']+self.params['err_G_prime']*(1.0+perturbations[3])
       self.params['ref_Debye']=self.params['ref_Debye']+self.params['err_ref_Debye']*(1.0+perturbations[4])
       self.params['ref_grueneisen']=self.params['ref_grueneisen']+self.params['err_ref_grueneisen']*(1.0+perturbations[5])
       self.params['q0']=self.params['q0']+self.params['err_q0']*(1.0+perturbations[6])
       self.params['eta_0s']=self.params['eta_0s']+self.params['err_eta_0s']*(1.0+perturbations[7])

   def set_state(self,pressure, temperature):
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


