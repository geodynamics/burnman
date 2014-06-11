# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import warnings

import numpy as np

from burnman.material import Material
import burnman.equation_of_state as eos
import burnman.birch_murnaghan as bm
import burnman.slb as slb
import burnman.mie_grueneisen_debye as mgd

class Mineral(Material):
    """
    This is the base class for all minerals. States of the mineral
    can only be queried after setting the pressure and temperature
    using set_state(). The method for computing properties of
    the material is set using set_method(), which should be done
    once after creating the material.

    This class is available as ``burnman.Mineral``.

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
            'V_0': 0., #Molar volume [m^3/(mole molecules)] at room pressure/temperature
            'K_0': 0., #Reference bulk modulus [Pa] at room pressure/temperature
            'Kprime_0': 0., #pressure derivative of bulk modulus
            'G_0': 0., #reference shear modulus at room pressure/temperature
            'Gprime_0': 0., #pressure derivative of shear modulus
            'molar_mass': 0., #molar mass in units of [kg/mol]
            'n': 0., #number of atoms per molecule
            'Debye_0': 0., #Debye temperature for material. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'grueneisen_0': 0., #Gruneisen parameter for material. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'q_0': 0., #q value used in caluclations. See Stixrude & Lithgow-Bertelloni, 2005 for values
            'eta_s_0': 0.0} #eta value used in calculations. See Stixrude & Lithgow-Bertelloni, 2005 for values
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
