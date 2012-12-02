
import numpy as np
import matplotlib.pyplot as plt

#eos imports
import birch_murnaghan as bm
import mie_grueneisen_debye as mgd
import slb_finitestrain as slb

class material:
	def __init__(self):
		self.params = {	'name':'generic',
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
		""" use "slb" or "mgd" """
		if (method=="slb"):
			self.method = slb;
		elif (method=="mgd"):
			self.method = mgd;
		else:
			raise("unsupported material method " + method)

        def set_state(self, pressure, temperature):
                self.pressure = pressure
                self.temperature = temperature

                self.V = self.method.volume(self.pressure, self.temperature, self.params)
		self.K_T = self.method.bulk_modulus(self.pressure, self.temperature, self.V, self.params)
		self.K_S = self.method.bulk_modulus_adiabatic(self.pressure, self.temperature, self.V, self.params)
		self.mu = self.method.shear_modulus(self.pressure, self.temperature, self.V, self.params)

	def to_string(self):
		return "'" + self.__class__.__name__ + "'"
	def molar_mass(self):
		return self.params['molar_mass']
	def density(self):
		return  self.params['molar_mass']/self.V
        def molar_volume(self):
                return self.V
	def bulk_modulus(self):
		return self.K_T
        def adiabatic_bulk_modulus(self):
                return self.K_S
	def shear_modulus(self):
		return self.mu
	def v_s(self):
		return np.sqrt(self.shear_modulus()*1.e9 / \
			self.density())/1000.
	def v_p(self):
		return np.sqrt((self.bulk_modulus() *1.e9 +4./3. * \
			self.shear_modulus())/self.density())/1000.
