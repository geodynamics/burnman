#standard numpy, scipy imports
import numpy as np
import matplotlib.pyplot as plt

#eos imports
import birch_murnaghan as bm
import mie_grueneisen_debye as mgd
import slb_finitestrain as slb
import voigt_reuss_hill as vrh
import userinput

import sys
import imp
userinput=imp.load_source("userinput",sys.argv[1])
method = vars()[userinput.method]
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
	def molar_mass(self,pressure,temperature):
		return self.params['molar_mass']
	def density(self, pressure, temperature):
		V = bm.bm_volume(pressure, self.params)
		return  self.params['molar_mass']/V
	def molar_volume(self, pressure, temperature):
		V = method.volume(pressure, temperature, self.params)
		return V
	def bulk_modulus(self, pressure, temperature):
		V = bm.bm_volume(pressure,  self.params)
		K_T = method.bulk_modulus(pressure, temperature, V, self.params)
		return K_T
        def adiabatic_bulk_modulus(self, pressure, temperature):
                V = bm.bm_volume(pressure, self.params)
                K_S = method.bulk_modulus_adiabatic(pressure, temperature, V, self.params)
                return K_S
	def shear_modulus(self, pressure, temperature):
		V = bm.bm_volume(pressure, self.params)
		mu = method.shear_modulus(pressure, temperature, V, self.params)
		return mu
	def v_s(self, pressure, temperature):
		return np.sqrt(self.shear_modulus(pressure, temperature)*1.e9 / \
			self.density(pressure,temperature))/1000.
	def v_p(self, pressure, temperature):
		return np.sqrt((self.bulk_modulus(pressure,temperature) *1.e9 +4./3. * \
			self.shear_modulus(pressure,temperature))/self.density(pressure,temperature))/1000.

# is this used? Sanne		
	def geotherm(self, pressure):
		return bm.geotherm_brown_shankland(pressure, self.params) 



