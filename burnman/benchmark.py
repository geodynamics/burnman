import numpy as np
import birch_murnaghan as bm
import mie_grueneisen_debye as mgd
import slb_finitestrain as slb
import matplotlib.pyplot as plt
from minerals import *


#Attempt to recreate Stixrude and Lithgow-Bertelloni (2005) Figure 1
def check_birch_murnaghan():
  test_mineral = material()
  test_mineral.params ={'name':'test',
			'ref_V': 6.844e-6,
			'ref_K': 259,
			'K_prime': 4.0,
			'ref_mu': 175.,
			'mu_prime': 1.7,
			'molar_mass': .0,
			'n': 0.,
			'ref_Debye': 0.,
			'ref_grueneisen': 0.,
			'q0': 0.}
 
  pressure = np.linspace(0., 140.e9, 100)
  volume = np.empty_like(pressure)
  bulk_modulus = np.empty_like(pressure)
  shear_modulus = np.empty_like(pressure)

  for i in range(len(pressure)):
    volume[i] = bm.bm_volume(pressure[i], test_mineral.params)
    bulk_modulus[i] = bm.bm_bulk_modulus(pressure[i], test_mineral.params)
    shear_modulus[i] = bm.bm_shear_modulus(pressure[i], test_mineral.params)

  plt.plot(pressure, bulk_modulus, pressure, shear_modulus)
  print bulk_modulus
  print shear_modulus
  plt.show()

def check_mgd_shim_duffy_kenichi():

  #Create gold material from Table 1
  gold = material()
  gold.params = {       'name': 'gold',
			'ref_V': 40.86e-6,
			'ref_K': 167,
			'K_prime': 5.0,
			'ref_mu': 0.,
			'mu_prime': 0.,
			'molar_mass': .196966,
			'n': 4.0,
			'ref_Debye': 170.,
			'ref_grueneisen': 2.97,
			'q0': 1.0}
  gold.method = mgd

  #Total pressures, pulled from Table 2  
  pressures = [np.array([0., 3.55, 7.55, 12.06, 17.16, 22.91, 29.42, 36.77, 45.11, 54.56, 65.29, 77.50, 91.42, 107.32, 125.51, 146.38, 170.38, 198.07])]
  pressures.append( np.array([4.99,8.53,12.53,17.04,22.13,27.88,34.38,41.73,50.06,59.50,70.22,82.43,96.33,112.22,130.40,151.25,175.24,202.90]))
  pressures.append( np.array([12.14,15.69,19.68,24.19,29.28,35.03,41.53,48.88,57.20,66.64,77.37,89.57,103.47,119.35,137.53,158.38,182.36,210.02]))
  pressures.append( np.array([19.30,22.84,26.84,31.35,36.44,42.19,48.68,56.03,64.35,73.80,84.52,96.72,110.62,126.50,144.68,165.53,189.51,217.17]))

  dv = np.empty_like(pressures)
  ref_dv = np.linspace(0.0, 0.34, len(pressures[0]))
  T= np.array([300., 1000.,2000.,3000.])
  for t in range(len(pressures)):
    for i in range(len(pressures[t])):
      gold.set_state(pressures[t][i]*1e9, T[t])
      V = gold.molar_volume()
      dv[t][i] = 1.0 - V/gold.params['ref_V']
    plt.plot(pressures[t],(ref_dv-dv[t]))
    print ref_dv-dv[t]
  plt.show()


#check_mgd_shim_duffy()
check_mgd_shim_duffy_kenichi()
#check_birch_murnaghan()

