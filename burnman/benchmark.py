# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 
sys.path.insert(1,os.path.abspath('.')) 
import burnman

import burnman.birch_murnaghan as bm
import burnman.mie_grueneisen_debye as mgd
import burnman.slb as slb
import matplotlib.image as mpimg



#Attempt to recreate Stixrude and Lithgow-Bertelloni (2005) Figure 1
def check_birch_murnaghan():
    plt.close()

    #make a test mineral
    test_mineral = burnman.material()
    test_mineral.params ={'name':'test',
                          'ref_V': 6.844e-6,
                          'ref_K': 259.0e9,
                          'K_prime': 4.0,
                          'ref_mu': 175.0e9,
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

    #calculate its static properties
    for i in range(len(pressure)):
        volume[i] = bm.volume(pressure[i], test_mineral.params)
        bulk_modulus[i] = bm.bulk_modulus(volume[i], test_mineral.params)
        shear_modulus[i] = bm.shear_modulus_third_order(volume[i], test_mineral.params) #third order is used for the plot we are comparing against

    #compare with figure 1
    plt.plot(pressure/1.e9, bulk_modulus/1.e9, pressure/1.e9, shear_modulus/1.e9)
    fig1 = mpimg.imread('input_figures/slb_fig1.png')
    plt.imshow(fig1, extent=[0,140,0,800], aspect='auto')
    plt.plot(pressure/1.e9, bulk_modulus/1.e9, 'g+', pressure/1.e9, shear_modulus/1.e9, 'g+')
    plt.ylim(0,800)
    plt.xlim(0,140)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Modulus (GPa)")
    plt.title("Comparing with Figure 1 of Stixrude and Lithgow-Bertelloni (2005)")

    plt.show()

def check_mgd_shim_duffy_kenichi():
    plt.close()
    #Create gold material from Table 1
    gold = burnman.material()
    gold.params = {'name': 'gold',
                   'ref_V': 10.22e-6,
                   'ref_K': 167.0e9,
                   'K_prime': 5.0,
                   'molar_mass': .196966,
                   'n': 1.0,
                   'ref_Debye': 170.,
                   'ref_grueneisen': 2.97, #this does better with gr = 2.93.  Why?
                   'q0': 1.0}
    gold.method = mgd

    #Total pressures, pulled from Table 2  
    ref_pressures = [np.array([0., 3.55, 7.55, 12.06, 17.16, 22.91, 29.42, 36.77, 45.11, 54.56, 65.29, 77.50, 91.42, 107.32, 125.51, 146.38, 170.38, 198.07])]
    ref_pressures.append( np.array([4.99,8.53,12.53,17.04,22.13,27.88,34.38,41.73,50.06,59.50,70.22,82.43,96.33,112.22,130.40,151.25,175.24,202.90]))
    ref_pressures.append( np.array([12.14,15.69,19.68,24.19,29.28,35.03,41.53,48.88,57.20,66.64,77.37,89.57,103.47,119.35,137.53,158.38,182.36,210.02]))
    ref_pressures.append( np.array([19.30,22.84,26.84,31.35,36.44,42.19,48.68,56.03,64.35,73.80,84.52,96.72,110.62,126.50,144.68,165.53,189.51,217.17]))
 
    eos = mgd.mgd3()
 
    pressures = np.empty_like(ref_pressures)
    ref_dv = np.linspace(0.0, 0.34, len(pressures[0]))
    ref_volumes = (1-ref_dv)*gold.params['ref_V']
    T= np.array([300., 1000.,2000.,3000.])
    for t in range(len(pressures)):
        for i in range(len(pressures[t])):
            pressures[t][i] = eos.pressure(T[t],ref_volumes[i], gold.params)
        plt.plot(ref_dv, (pressures[t]/1.e9-ref_pressures[t]))
    plt.ylim(-1,1)
    plt.ylabel("Difference in pressure (GPa)")
    plt.xlabel("1-dV/V")
    plt.title("Comparing with Shim, Duffy, and Kenichi (2002)")
    plt.show()


def check_mgd_fei_mao_shu_hu():
    mgfeo = burnman.material() 
    mgfeo.params = {       'name': 'MgFeO',
                    'ref_V': 11.657e-6,
                    'ref_K': 157.0e9,
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
    ref_pressures = ref_pressures
    pressures = np.empty_like(volumes)

    eos = mgd.mgd3()
  
    for i in range(len(temperatures)):
        pressures[i] = eos.pressure(temperatures[i],volumes[i], mgfeo.params)
 
    plt.scatter(temperatures, (pressures/1.e9-ref_pressures))
    plt.ylim(-1,1)
    plt.title("Comparing with Fei, Mao, Shu, and Hu (1991)")
    plt.xlabel("Temperature (K) at various volumes")
    plt.ylabel("Difference in total pressure (GPa)")
    plt.show()


def check_slb_fig3():
    perovskite= burnman.material() 
    perovskite.params = {       'name': 'perovksite',
                    'ref_V': burnman.tools.molar_volume_from_unit_cell_volume(168.27, 4.),
                    'ref_grueneisen': 1.63,
                    'q0': 1.7}

    volume = np.linspace(0.6, 1.0, 100)
    grueneisen_slb = np.empty_like(volume)
    grueneisen_mgd = np.empty_like(volume)
    q_slb = np.empty_like(volume)
    q_mgd = np.empty_like(volume)

    slb_eos = slb.slb2()
    mgd_eos = mgd.mgd2()
    

    #calculate its thermal properties
    for i in range(len(volume)):
        #call with dummy pressure and temperatures, they do not change it
        grueneisen_slb[i] = slb_eos.grueneisen_parameter(0., 0., volume[i]*perovskite.params['ref_V'], perovskite.params)
        grueneisen_mgd[i] = mgd_eos.grueneisen_parameter(0., 0., volume[i]*perovskite.params['ref_V'], perovskite.params)
        q_slb[i] = slb_eos.volume_dependent_q(1./volume[i], perovskite.params)
        q_mgd[i] = perovskite.params['q0']

    #compare with figure 7
    fig1 = mpimg.imread('input_figures/slb_fig3.png')
    plt.imshow(fig1, extent=[0.6, 1.0,0.35,2.0], aspect='auto')
    plt.plot(volume, grueneisen_slb, 'g+', volume, grueneisen_mgd, 'b+')
    plt.plot(volume, q_slb, 'g+', volume, q_mgd, 'b+')
    plt.xlim(0.6,1.0)
    plt.ylim(0.35,2.0)
    plt.ylabel("Grueneisen parameter")
    plt.xlabel("Relative Volume V/V0")
    plt.title("Comparing with Figure 3 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

def check_slb_fig7():
    forsterite = burnman.material() 
    forsterite.params = {       'name': 'forsterite',
                    'ref_V': 43.60e-6,
                    'ref_K': 128.0e9,
                    'K_prime': 4.2,
                    'ref_mu' : 82.0e9,
                    'mu_prime' : 1.4,
                    'n': 7.0,
                    'ref_Debye': 809.,
                    'ref_grueneisen': .99,
                    'q0': 2.1, 
                    'eta_0s' : 2.4}
    forsterite.method = slb.slb3()

    temperature = np.linspace(0., 2000., 200)
    volume = np.empty_like(temperature)
    bulk_modulus = np.empty_like(temperature)
    shear_modulus = np.empty_like(temperature)
    heat_capacity = np.empty_like(temperature)

    pressure = 1.0e5
    forsterite.set_state(pressure, 300.)
    ref_Ks = forsterite.adiabatic_bulk_modulus()

    #calculate its thermal properties
    for i in range(len(temperature)):
        forsterite.set_state(pressure, temperature[i])
        volume[i] = forsterite.molar_volume()/forsterite.params['ref_V']
        bulk_modulus[i] = forsterite.adiabatic_bulk_modulus()/ref_Ks
        shear_modulus[i] = forsterite.shear_modulus()/forsterite.params['ref_mu']
        heat_capacity[i] = forsterite.heat_capacity_p()/forsterite.params['n']
 

    #compare with figure 7
    fig1 = mpimg.imread('input_figures/slb_fig7_vol.png')
    plt.imshow(fig1, extent=[0,2200,0.99,1.08], aspect='auto')
    plt.plot(temperature, volume, 'g+')
    plt.ylim(0.99,1.08)
    plt.xlim(0,2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Volume V/V0")
    plt.title("Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

    fig1 = mpimg.imread('input_figures/slb_fig7_Cp.png')
    plt.imshow(fig1, extent=[0,2200,0.,70.], aspect='auto')
    plt.plot(temperature, heat_capacity, 'g+')
    plt.ylim(0,70)
    plt.xlim(0,2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat Capacity Cp")
    plt.title("Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()


    fig1 = mpimg.imread('input_figures/slb_fig7_K.png')
    plt.imshow(fig1, extent=[0,2200,0.6,1.02], aspect='auto')
    plt.plot(temperature, bulk_modulus, 'g+')
    plt.ylim(0.6,1.02)
    plt.xlim(0,2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Bulk Modulus K/K0")
    plt.title("Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

    fig1 = mpimg.imread('input_figures/slb_fig7_G.png')
    plt.imshow(fig1, extent=[0,2200,0.6,1.02], aspect='auto')
    plt.plot(temperature, shear_modulus, 'g+')
    plt.ylim(0.6,1.02)
    plt.xlim(0,2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Shear Modulus G/G0")
    plt.title("Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

if __name__ == "__main__":
    check_birch_murnaghan()
    check_slb_fig7()
    check_slb_fig3()
    check_mgd_shim_duffy_kenichi()
    check_mgd_fei_mao_shu_hu()
