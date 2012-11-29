import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt


# compute the bulk modulus as per the third order
# birch-murnaghan equation of state.  Returns bulk
# modulus in the same units as the reference bulk
# modulus.  Pressure must be in Pa
def bm_bulk_modulus(pressure, ref_rho, ref_K, K_prime):
	rho = bm_density(pressure, ref_rho, ref_K, K_prime)
	x = rho/ref_rho 
	xi = -3./4. * (K_prime - 4.)
	return ref_K/2. * ( (7.*pow(x, 7./3.) - 5.*pow(x, 5./3.))*(1. - xi*pow(x,2./3.)) + 2.*xi*(pow(x,3.) - pow(x,7./3.)))
  
# equation for the third order birch-murnaghan
# equation of  state, returns pressure in the same
# units that are supplied for the reference bulk
# modulus (ref_K)
def birch_murnaghan(rho, ref_rho, ref_K, K_prime):
	x = rho/ref_rho
	return 3.*ref_K/2. * (pow(x, 7./3.) - pow(x, 5./3.)) \
	* (1 + .75*(K_prime - 4)*(pow(x, 2./3.) - 1))

# get the birch-murnaghan density at a reference 
# temperature for a given pressure.  Give pressure
# in Pa, ref_rho in g/cc and ref_K in GPa
def bm_density(pressure, ref_rho, ref_K, K_prime):
	return opt.brentq(lambda x: birch_murnaghan(x, ref_rho, ref_K, K_prime)*1e9-pressure, 0, 12)

def bm_shear_modulus(pressure, ref_rho, ref_K, K_prime):
	x = rho/ref_rho
	rho = bm_density(pressure, ref_rho, ref_K, K_prime)
		
	

def geotherm(pressure, ref_rho, ref_K, K_prime, grueneisen):
        minP = 0
        T0 = 1700
        dP = 1.e9
      
	# integrate the adiabatic gradient equation
	lnT = integrate.quad( lambda x: (grueneisen/(1.e9*bm_bulk_modulus(x, ref_rho, ref_K, K_prime)))/(bm_density(x, ref_rho, ref_K, K_prime)), minP, pressure)
	T = T0*np.exp(lnT[0])
	print lnT[0]
	return T

def plot_geotherm():
        ref_rho = 3.37
        ref_K = 130.3
        K_prime = 3.9 
        pressures = np.arange(1.0e9, 135.0e9, 1.0e9)
        temperatures = np.arange(1.0e9, 135.0e9, 1.0e9)
        for i in range(len(pressures)):
		temperatures[i] = geotherm(pressures[i], ref_rho, ref_K, K_prime, 3.5)
        
	plt.plot(pressures, temperatures)
        plt.show()



if __name__ == "__main__":
	plot_geotherm()

