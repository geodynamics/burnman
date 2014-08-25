# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import warnings

import numpy as np

kd = lambda x,y : 1 if x==y else 0

class SolidSolution(Mineral):
    """
    This is the base class for all solid solutions. 
    States of the solid solution can only be queried after setting 
    the pressure and temperature using set_state(). 
    The method for computing properties of
    the solution is set using set_method(), which should be done
    once after creating the material.

    This class is available as ``burnman.SolidSolution``.

    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T 
    and P derivatives in J/K/mol and m^3/(mol molecule).
    """

    # init sets up matrices to speed up calculations when P, T, X is defined.
    def __init__(self, base_material, interaction_parameter):
        self.base_material = base_material
        self.molar_fraction = molar_fraction
        assert(len(base_material) == len(molar_fraction))
        assert(len(base_material) == len(van_laar_parameter))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        # make sure the number of atoms in each endmember is the same
        # (when we start adding minerals with vacancies
        # this will need to be changed).
        for m in base_material:
            if(base_material[0].params.has_key('n')):
                assert(m.params['n'] == base_material[0].params['n'])

    # Set the state of all the endmembers
    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.method = self.method
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_material))

        self.params = {}

        # NEEDS CHANGING FROM HERE!!!
        # Description of RTlngamma terms given by the asymmetric formalism
        # is included in Holland and Powell (2003)

        # the volume is calculated in the same way as Gex (p.493),
        # replacing Wij with Wvij
        Vex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*self.interaction_parameter[i][2] / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        V= sum([ self.base_materials[i].molar_volume() * self.molar_fraction[i] for i in itrange ]) + Vex

        # the volume is calculated in the same way as Gex (p.493),
        # replacing Wij with Wvij
        # NEED TO CORRECT INDICES FOR INT PARAM
        # W[i][j] == W[(itrange-1)*i+(j-1)]
        Vex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*self.interaction_parameter[(itrange-1)*i+(j-1)][2] / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        V= sum([ self.base_materials[i].molar_volume() * self.molar_fraction[i] for i in itrange ]) + Vex


        # Same for the Gibbs free energy
        # NEED TO CORRECT INDICES FOR INT PARAM
        Gex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*(pressure*self.interaction_parameter[(itrange-1)*i+(j-1)][0] + temperature*self.interaction_parameter[(itrange-1)*i+(j-1)][1] + pressure*self.interaction_parameter[(itrange-1)*i+(j-1)][2]) / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        G= sum([ self.base_materials[i].molar_gibbs() * self.molar_fraction[i] for i in itrange ]) + Gex

        

        

        for prop in self.base_materials[0].params:
           try:
               self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
               #if there is a type error, it is probably a string.  Just go with the value of the first base_material.
               self.params[prop] = self.base_materials[0].params[prop]
        Mineral.set_state(self, pressure, temperature)

    # BELOW IS FROM THE MINERAL CLASS.
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
                self.method = slb.SLB2()
            elif (method == "mgd2"):
                self.method = mgd.MGD2()
            elif (method == "mgd3"):
                self.method = mgd.MGD3()
            elif (method == "slb3"):
                self.method = slb.SLB3()
            elif (method == "bm2"):
                self.method = bm.BM2()
            elif (method == "bm3"):
                self.method = bm.BM3()
            else:
                raise Exception("unsupported material method " + method)
        elif ( issubclass(method, eos.EquationOfState) ):
            self.method = method()
        else:
            raise Exception("unsupported material method " + method.__class__.__name__ )

    def to_string(self):
        """
        Returns the name of the mineral class
        """
        return "'" + self.__class__.__module__.replace(".minlib_",".") + "." + self.__class__.__name__ + "'"

    def unroll(self):
        return ([1.0],[self])

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

        if (self.params.has_key('G_0') and self.params.has_key('Gprime_0')):
            self.G = self.method.shear_modulus(self.pressure, self.temperature, self.V, self.params)
        else:
            self.G = float('nan') #nan if there is no G, this should propagate through calculations to the end
            warnings.warn(('Warning: G_0 and or Gprime_0 are undefined for ' + self.to_string()))

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
