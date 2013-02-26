# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import birch_murnaghan as bm
import mie_grueneisen_debye as mgd
import slb_finitestrain as slb
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from minerals import *


#Attempt to recreate Stixrude and Lithgow-Bertelloni (2005) Figure 1
def check_birch_murnaghan():
    plt.close()

    #make a test mineral
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
    test_mineral.method = slb
 
    pressure = np.linspace(0., 140.e9, 100)
    volume = np.empty_like(pressure)
    bulk_modulus = np.empty_like(pressure)
    shear_modulus = np.empty_like(pressure)

    #calculate its static properties
    for i in range(len(pressure)):
        volume[i] = bm.volume(pressure[i], test_mineral.params)
        bulk_modulus[i] = bm.bulk_modulus(volume[i], test_mineral.params)
        shear_modulus[i] = bm.shear_modulus(volume[i], test_mineral.params)

    #compare with figure 1
    plt.plot(pressure/1.e9, bulk_modulus, pressure/1.e9, shear_modulus)
    fig1 = mpimg.imread('../data/slb_fig1.png')
    plt.imshow(fig1, extent=[0,140,0,800], aspect='auto')
    plt.plot(pressure/1.e9, bulk_modulus, 'g+', pressure/1.e9, shear_modulus, 'g+')
    plt.ylim(0,800)
    plt.xlim(0,140)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Modulus (GPa)")
    plt.title("Comparing with Figure 1 of Stixrude and Lithgow-Bertelloni (2005)")

    plt.show()

def check_mgd_shim_duffy_kenichi():
    plt.close()
    #Create gold material from Table 1
    gold = material()
    gold.params = {'name': 'gold',
                   'ref_V': 40.86e-6,
                   'ref_K': 167,
                   'K_prime': 5.0,
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

    volumes = np.empty_like(pressures)
    ref_dv = np.linspace(0.0, 0.34, len(pressures[0]))
    ref_volumes = (1-ref_dv)*gold.params['ref_V']
    T= np.array([300., 1000.,2000.,3000.])
    for t in range(len(pressures)):
        for i in range(len(pressures[t])):
            gold.set_state(pressures[t][i]*1e9, T[t])
            volumes[t][i] = gold.molar_volume()
        plt.plot(pressures[t],(volumes[t]-ref_volumes)/ref_volumes*100)
    plt.ylim(-1,1)
    plt.ylabel("% difference in volume")
    plt.xlabel("Pressure (GPa)")
    plt.title("Comparing with Shim, Duffy, and Kenichi (2002)")
    plt.show()


def check_mgd_fei_mao_shu_hu():
    mgfeo = material() 
    mgfeo.params = {       'name': 'MgFeO',
                    'ref_V': 11.657e-6,
                    'ref_K': 157.0,
                    'K_prime': 4.0,
                    'molar_mass': .196966,
                    'n': 2.0,
                    'ref_Debye': 500.,
                    'ref_grueneisen': 1.50,
                    'q0': 1.1}
    mgfeo.method = mgd

    #pulled from table 1
    temperatures = np.array([300,300,483,483,483,590,593,593,593,700,600,500,650,600,600,650,700,737,727,673,600,543,565,585,600,628,654,745,768,747,726,700,676])
    volumes = np.array([77.418,72.327,74.427,73.655,72.595,74.1,73.834,73.101,70.845,73.024,72.630,68.644,72.969,72.324,71.857,72.128,73.283,73.337,72.963,71.969,69.894,67.430,67.607,67.737,68.204,68.518,68.955,70.777,72.921,72.476,72.152,71.858,71.473])
    #change from cubic angstroms per unit cell to cubic meters per mol of molecules.
    volumes = volumes/1.e30*6.022141e23/4.0
    ref_pressures = np.array([0.0, 12.23,7.77,9.69,12.54,9.21,9.90,11.83,18.35,12.68,13.15,25.16,12.53,14.01,15.34,14.86,11.99,12.08,13.03,15.46,21.44,29.98,29.41,29.05,27.36,26.38,24.97,19.49,13.39,14.48,15.27,15.95,16.94])
    ref_pressures = ref_pressures*1.e9
    pressures = np.empty_like(volumes)
  
    for i in range(len(temperatures)):
        pressures[i] = mgd.pressure(volumes[i], temperatures[i], mgfeo.params)
 
    plt.scatter(temperatures, (pressures-ref_pressures)/ref_pressures*100)
    plt.ylim(-1,1)
    plt.title("Comparing with Fei, Mao, Shu, and Hu (1991)")
    plt.xlabel("Temperature (K)")
    plt.ylabel("% difference in total pressure")
    plt.show()

if __name__ == "__main__":
    check_mgd_shim_duffy_kenichi()
    check_birch_murnaghan()
    check_mgd_fei_mao_shu_hu()
