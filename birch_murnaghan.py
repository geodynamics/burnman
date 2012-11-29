import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt
import mie_grueneisen_debye as mgd
# compute the bulk modulus as per the third order
# birch-murnaghan equation of tate.  Returns bulk
# modulus in the same units as the reference bulk
# modulus.  Pressure must be in Pa
def bm_bulk_modulus(pressure, params):
	V = bm_volume(pressure, params)
	x = params['ref_V']/V
	
	xi = -3./4. * (params['K_prime'] - 4.)
##  derived by Ian
	test_K=params['ref_K']/2. * ( (7.*pow(x, 7./3.) - 5.*pow(x, 5./3.))*(1. - xi*pow(x,2./3.)) + 2.*xi*(pow(x,3.) - pow(x,7./3.)))
## Matas et al. equation A2
	A= -2.*(params['K_prime']-4.)+8./3.
	B=-4.*(params['K_prime']-4.)+8./3.
	C= B-A
	top=5./3.*A*pow(x,5./3.)-7./3.*B*pow(x,7./3.)+3.*C*pow(x,3.)
	bottom=A*pow(x,5./3.)-B*pow(x,7./3.)+C*pow(x,3.)
	test_K2=pressure*top/bottom/1e9
##
 	return test_K2
 
# equation for the third order birch-murnaghan
# equation of  state, returns pressure in the same
# units that are supplied for the reference bulk
# modulus (params['ref_K'])
def birch_murnaghan(x, params):
	return 3.*params['ref_K']/2. * (pow(x, 7./3.) - pow(x, 5./3.)) \
	* (1 - .75*(4-params['K_prime'] )*(pow(x, 2./3.) - 1))

# get the birch-murnaghan density at a reference 
# temperature for a given pressure.  Give pressure in
# Pa, returns density in kg/m^3
def bm_density(pressure, params):
	ratio = opt.brentq(lambda x: birch_murnaghan(x, params)*1e9-pressure, 0, 12)
        return ratio*params

# get the birch-murnaghan density at a reference 
# temperature for a given pressure (Pa).  Returns
# molar volume in m^3
def bm_volume(pressure, params):
	ratio = opt.brentq(lambda x: birch_murnaghan(x, params)*1e9-pressure, 0.1, 10.0)
        return 1./ratio * params['ref_V']

# get the birch murnaghan shear modulus at a reference temperature,
# for a given pressure.  Returns shear modulus in GPa (the same units
# as in params['ref_mu']
def bm_shear_modulus(pressure, params):
	V = bm_volume(pressure, params)
	x = params['ref_V']/V
	G=params['ref_mu'] * pow(x,5./3.)*(1.-0.5*(pow(x,2./3.)-1.)*(5.-3.*params['mu_prime']*params['ref_K']/params['ref_mu']))
	f=-1*1./2.*(1.-pow(x,2./3.))
	# exact same as stixrude up to first order, no thermal correction
    	G_slb=pow(1.+2.*f,5./2.)*(params['ref_mu']+(3.*(params['ref_K']*params['mu_prime']) - 5.*params['ref_mu'])*f)#+ \
         #((6.*params['ref_K']*params['mu_prime']) - 24.*params['ref_K'] - 14.*params['ref_mu']+ 9./2.*params['ref_K']*params['K_prime'])*pow(f,2.)) 
	return params['ref_mu'] * pow(x,5./3.)*(1.-0.5*(pow(x,2./3.)-1.)*(5.-3.*params['mu_prime']*params['ref_K']/params['ref_mu']))


#if __name__ == "__main__":
# not working 
#        ref_rho = 3.37
#        params = {'ref_K':130.3}
#        params['K_prime'] = 3.9 
#	params['ref_mu'] =  150.0
#	params['mu_prime'] = 2.0
#        pressures = np.arange(1.0e9, 135.0e9, 1.0e9)
#        temperatures = np.arange(1.0e9, 135.0e9, 1.0e9)
#        for i in range(len(pressures)):
#		temperatures[i] = bm_shear_modulus(pressures[i], params)
#        
#	plt.plot(pressures, temperatures)
#        plt.show()
