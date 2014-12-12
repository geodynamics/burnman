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
import inspect
import burnman.modified_tait as mt
import burnman.cork as cork

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
    To convert this to m^3/(mol of molecule) you should multiply by 10^(-30) *
    N_a / Z, where N_a is Avogadro's number and Z is the number of formula units per
    unit cell. You can look up Z in many places, including www.mindat.org
    """

    def __init__(self):
        if 'params' not in self.__dict__:
            self.params={}
        self.method = None
        if 'equation_of_state' in self.params:
            self.set_method(self.params['equation_of_state'])

    def set_method(self, method):
        """
        Set the equation of state to be used for this mineral.
        Takes a string corresponding to any of the predefined
        equations of state:  'bm2', 'bm3', 'mgd2', 'mgd3', 'slb2',
        'slb3' or 'mtait'.  Alternatively, you can pass a user defined
        class which derives from the equation_of_state base class.
        """

        if method is None:
            self.method = None
            return

        def to_method_obj(method):

            if isinstance(method, basestring):
                if (method == "slb2"):
                    return slb.SLB2()
                elif (method == "mgd2"):
                    return mgd.MGD2()
                elif (method == "mgd3"):
                    return mgd.MGD3()
                elif (method == "slb3"):
                    return slb.SLB3()
                elif (method == "bm2"):
                    return bm.BM2()
                elif (method == "bm3"):
                    return bm.BM3()
                elif (method == "mtait"):
                    return mt.MT()
                elif (method == "cork"):
                    return cork.CORK()
                else:
                    raise Exception("unsupported material method " + method)
            elif isinstance(method, eos.EquationOfState):
                return method
            elif inspect.isclass(method) and issubclass(method, eos.EquationOfState):
                return method()
            else:
                raise Exception("unsupported material method " + method.__class__.__name__ )

        new_method = to_method_obj(method)
        if self.method is not None and 'equation_of_state' in self.params:
            self.method = to_method_obj(self.params['equation_of_state'])

        if type(new_method).__name__ == 'instance':
            raise Exception("Please derive your method from object (see python old style classes)")

        if self.method is not None and type(new_method) is not type(self.method):
            warnings.warn('Overriding database equation of state. From ' + self.method.__class__.__name__ + ' to ' + new_method.__class__.__name__, stacklevel=2)

        self.method = new_method

    def to_string(self):
        """
        Returns the name of the mineral class
        """
        return "'" + self.__class__.__module__.replace(".minlib_",".") + "." + self.__class__.__name__ + "'"

    def debug_print(self, indent=""):
        print "%s%s" % (indent, self.to_string())

    def unroll(self):
        return ([1.0],[self])

    def eos_pressure(self, temperature, volume):
        return self.method.pressure(temperature, volume, self.params)

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

        if self.method is None:
            raise AttributeError, "no method set for mineral, or equation_of_state given in mineral.params"

        self.V = self.method.volume(self.pressure, self.temperature, self.params)
        self.gr = self.method.grueneisen_parameter(self.pressure, self.temperature, self.V, self.params)
        self.K_T = self.method.isothermal_bulk_modulus(self.pressure, self.temperature, self.V, self.params)
        self.K_S = self.method.adiabatic_bulk_modulus(self.pressure, self.temperature, self.V, self.params)
        self.C_v = self.method.heat_capacity_v(self.pressure, self.temperature, self.V, self.params)
        self.C_p = self.method.heat_capacity_p(self.pressure, self.temperature, self.V, self.params)
        self.alpha = self.method.thermal_expansivity(self.pressure, self.temperature, self.V, self.params)

        # Attempt to calculate the gibbs free energy and helmholtz free energy, but don't complain if the
        # equation of state does not calculate it, or if the mineral params do not have the requisite entries.
        try:
            self.gibbs = self.method.gibbs_free_energy(self.pressure, self.temperature, self.params)
        except (KeyError, NotImplementedError):
            self.gibbs = float('nan')
        try:
            self.helmholtz = self.method.helmholtz_free_energy(self.temperature, self.V, self.params)
        except (KeyError, NotImplementedError):
            self.helmholtz = float('nan')
        try:
            self.S = self.method.entropy(self.pressure, self.temperature, self.params)
        except (KeyError, NotImplementedError):
            self.S = float('nan')
        try:
            self.H = self.method.enthalpy(self.pressure, self.temperature, self.params)
        except (KeyError, NotImplementedError):
            self.H = float('nan')

        if (self.params.has_key('G_0') and self.params.has_key('Gprime_0')):
            self.G = self.method.shear_modulus(self.pressure, self.temperature, self.V, self.params)
        else:
            self.G = float('nan') #nan if there is no G, this should propagate through calculations to the end
            if self.params['equation_of_state'] != 'mtait':
                warnings.warn(('Warning: G_0 and or Gprime_0 are undefined for ' + self.to_string()))

    # The following gibbs function avoids having to calculate a bunch of unnecessary parameters over P-T space. This will be useful for gibbs minimisation.
    def calcgibbs(self, pressure, temperature):
        return self.method.gibbs_free_energy(pressure, temperature, self.params)

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
    def compressibility(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        """
        return 1./self.K_T
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

    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        """
        return self.gibbs

    def molar_helmholtz(self):
        """
        Returns Gibbs free energy of the mineral [J]
        """
        return self.helmholtz
