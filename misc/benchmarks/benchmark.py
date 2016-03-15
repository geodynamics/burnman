from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))
import numpy as np
import matplotlib.pyplot as plt

import burnman

import burnman.eos.birch_murnaghan as bm
import burnman.eos.birch_murnaghan_4th as bm4
import burnman.eos.mie_grueneisen_debye as mgd
import burnman.eos.slb as slb
import burnman.eos.vinet as vinet
import matplotlib.image as mpimg


def check_birch_murnaghan():
    """
    Recreates Stixrude and Lithgow-Bertelloni (2005) Figure 1, bulk and shear modulus without thermal corrections
    """
    plt.close()

    # make a test mineral
    test_mineral = burnman.Mineral()
    test_mineral.params = {'name': 'test',
                           'V_0': 6.844e-6,
                           'K_0': 259.0e9,
                           'Kprime_0': 4.0,
                           'G_0': 175.0e9,
                           'Gprime_0': 1.7,
                           'molar_mass': .0,
                           }
    test_mineral.set_method('bm3')

    pressure = np.linspace(0., 140.e9, 100)
    volume = np.empty_like(pressure)
    bulk_modulus = np.empty_like(pressure)
    shear_modulus = np.empty_like(pressure)

    # calculate its static properties
    for i in range(len(pressure)):
        volume[i] = bm.volume(pressure[i], test_mineral.params)
        bulk_modulus[i] = bm.bulk_modulus(volume[i], test_mineral.params)
        shear_modulus[i] = bm.shear_modulus_third_order(
            volume[i], test_mineral.params)  # third order is used for the plot we are comparing against

    # compare with figure 1
    plt.plot(pressure / 1.e9, bulk_modulus /
             1.e9, pressure / 1.e9, shear_modulus / 1.e9)
    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig1.png')
    plt.imshow(fig1, extent=[0, 140, 0, 800], aspect='auto')
    plt.plot(pressure / 1.e9, bulk_modulus / 1.e9,
             'g+', pressure / 1.e9, shear_modulus / 1.e9, 'g+')
    plt.ylim(0, 800)
    plt.xlim(0, 140)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Modulus (GPa)")
    plt.title(
        "Comparing with Figure 1 of Stixrude and Lithgow-Bertelloni (2005)")

    plt.show()


def check_birch_murnaghan_4th():
    """
    Recreates the formulation of the 4th order Birch-Murnaghan EOS as in Ahmad and Alkammash, 2012; Figure 1.
    """
    plt.close()

    # make a test mineral
    test_mineral = burnman.Mineral()
    test_mineral.params = {'name': 'test',
                           'V_0': 10.e-6,
                           'K_0': 72.7e9,
                           'Kprime_0':  4.14,
                           'Kprime_prime_0': -0.0484e-9,
                           }

    test_mineral.set_method('bm4')

    pressure = np.linspace(0., 90.e9, 20)
    volume = np.empty_like(pressure)

    # calculate its static properties
    for i in range(len(pressure)):
        volume[i] = bm4.volume_fourth_order(
            pressure[i], test_mineral.params) / test_mineral.params.get('V_0')

    # compare with figure 1
    plt.plot(pressure / 1.e9, volume)
    fig1 = mpimg.imread('../../burnman/data/input_figures/Ahmad.png')
    plt.imshow(fig1, extent=[0., 90., .65, 1.], aspect='auto')
    plt.plot(pressure / 1.e9, volume, marker='o',
             color='r', linestyle='', label='BM4')
    plt.legend(loc='lower left')
    plt.xlim(0., 90.)
    plt.ylim(.65, 1.)
    plt.xlabel("Volume/V0")
    plt.ylabel("Pressure (GPa)")
    plt.title("Comparing with Figure 1 of Ahmad et al., (2012)")

    plt.show()


def check_vinet():
    """
    Recreates Dewaele et al., 2006, Figure 1, fitting a Vinet EOS to Fe data
    """
    plt.close()

    # make a test mineral
    test_mineral = burnman.Mineral()
    test_mineral.params = {'name': 'test',
                           'V_0': 6.75e-6,
                           'K_0': 163.4e9,
                           'Kprime_0': 5.38,
                           }

    test_mineral.set_method('vinet')

    pressure = np.linspace(17.7e9, 300.e9, 20)
    volume = np.empty_like(pressure)

    # calculate its static properties
    for i in range(len(pressure)):
        volume[i] = vinet.volume(pressure[i], test_mineral.params)

    # compare with figure 1
    plt.plot(pressure / 1.e9, volume / 6.02e-7)
    fig1 = mpimg.imread('../../burnman/data/input_figures/Dewaele.png')
    plt.imshow(fig1, extent=[0., 300., 6.8, 11.8], aspect='auto')
    plt.plot(pressure / 1.e9, volume / 6.02e-7, marker='o',
             color='r', linestyle='', label='Vinet Fit')
    plt.legend(loc='lower left')
    plt.xlim(0., 300.)
    plt.ylim(6.8, 11.8)
    plt.ylabel("Volume (Angstroms^3/atom")
    plt.xlabel("Pressure (GPa)")
    plt.title("Comparing with Figure 1 of Dewaele et al., (2006)")

    plt.show()


def check_mgd_shim_duffy_kenichi():
    """
    Attemmpts to recreate Shim Duffy Kenichi (2002)
    """
    plt.close()
    # Create gold material from Table 1
    gold = burnman.Mineral()
    gold.params = {'name': 'gold',
                   'V_0': 10.22e-6,
                   'K_0': 167.0e9,
                   'Kprime_0': 5.0,
                   'G_0': 0.0e9,
                   'Gprime_0': 0.0,
                   'molar_mass': .196966,
                   'n': 1.0,
                   'Debye_0': 170.,
                   'grueneisen_0': 2.97,  # this does better with gr = 2.93.  Why?
                   'q_0': 1.0}
    gold.set_method('mgd3')

    # Total pressures, pulled from Table 2
    ref_pressures = [
        np.array([0., 3.55, 7.55, 12.06, 17.16, 22.91, 29.42, 36.77, 45.11, 54.56, 65.29, 77.50, 91.42, 107.32, 125.51, 146.38, 170.38, 198.07])]
    ref_pressures.append(
        np.array([4.99, 8.53, 12.53, 17.04, 22.13, 27.88, 34.38, 41.73, 50.06, 59.50, 70.22, 82.43, 96.33, 112.22, 130.40, 151.25, 175.24, 202.90]))
    ref_pressures.append(
        np.array([12.14, 15.69, 19.68, 24.19, 29.28, 35.03, 41.53, 48.88, 57.20, 66.64, 77.37, 89.57, 103.47, 119.35, 137.53, 158.38, 182.36, 210.02]))
    ref_pressures.append(
        np.array([19.30, 22.84, 26.84, 31.35, 36.44, 42.19, 48.68, 56.03, 64.35, 73.80, 84.52, 96.72, 110.62, 126.50, 144.68, 165.53, 189.51, 217.17]))

    eos = mgd.MGD3()

    pressures = np.empty_like(ref_pressures)
    ref_dv = np.linspace(0.0, 0.34, len(pressures[0]))
    ref_volumes = (1 - ref_dv) * gold.params['V_0']
    T = np.array([300., 1000., 2000., 3000.])
    for t in range(len(pressures)):
        for i in range(len(pressures[t])):
            pressures[t][i] = eos.pressure(T[t], ref_volumes[i], gold.params)
        plt.plot(ref_dv, (pressures[t] / 1.e9 - ref_pressures[t]))
    plt.ylim(-1, 1)
    plt.ylabel("Difference in pressure (GPa)")
    plt.xlabel("1-dV/V")
    plt.title("Comparing with Shim, Duffy, and Kenichi (2002)")
    plt.show()


def check_mgd_fei_mao_shu_hu():
    """
    Benchmark agains Fei Mao Shu Hu (1991)
    """
    mgfeo = burnman.Mineral()
    mgfeo.params = {'name': 'MgFeO',
                    'V_0': 11.657e-6,
                    'K_0': 157.0e9,
                    'Kprime_0': 4.0,
                    'G_0': 0.0e9,
                    'Gprime_0': 0.0,
                    'molar_mass': .196966,
                    'n': 2.0,
                    'Debye_0': 500.,
                    'grueneisen_0': 1.50,
                    'q_0': 1.1}
    mgfeo.set_method('mgd3')

    # pulled from table 1
    temperatures = np.array(
        [300, 300, 483, 483, 483, 590, 593, 593, 593, 700, 600, 500, 650, 600,
         600, 650, 700, 737, 727, 673, 600, 543, 565, 585, 600, 628, 654, 745, 768, 747, 726, 700, 676])
    volumes = np.array(
        [77.418, 72.327, 74.427, 73.655, 72.595, 74.1, 73.834, 73.101, 70.845, 73.024, 72.630, 68.644, 72.969, 72.324, 71.857,
         72.128, 73.283, 73.337, 72.963, 71.969, 69.894, 67.430, 67.607, 67.737, 68.204, 68.518, 68.955, 70.777, 72.921, 72.476, 72.152, 71.858, 71.473])
    # change from cubic angstroms per unit cell to cubic meters per mol of
    # molecules.
    volumes = volumes / 1.e30 * 6.022141e23 / 4.0
    ref_pressures = np.array(
        [0.0, 12.23, 7.77, 9.69, 12.54, 9.21, 9.90, 11.83, 18.35, 12.68, 13.15, 25.16, 12.53, 14.01, 15.34,
         14.86, 11.99, 12.08, 13.03, 15.46, 21.44, 29.98, 29.41, 29.05, 27.36, 26.38, 24.97, 19.49, 13.39, 14.48, 15.27, 15.95, 16.94])
    ref_pressures = ref_pressures
    pressures = np.empty_like(volumes)

    eos = mgd.MGD3()

    for i in range(len(temperatures)):
        pressures[i] = eos.pressure(temperatures[i], volumes[i], mgfeo.params)

    plt.scatter(temperatures, (pressures / 1.e9 - ref_pressures))
    plt.ylim(-1, 1)
    plt.title("Comparing with Fei, Mao, Shu, and Hu (1991)")
    plt.xlabel("Temperature (K) at various volumes")
    plt.ylabel("Difference in total pressure (GPa)")
    plt.show()


def check_slb_fig3():
    """
    Benchmark grueneisen parameter against figure 3 of Stixrude and Lithgow-Bertelloni (2005b)
    """
    perovskite = burnman.Mineral()
    perovskite.params = {'name': 'perovksite',
                         'V_0': burnman.tools.molar_volume_from_unit_cell_volume(168.27, 4.),
                         'grueneisen_0': 1.63,
                         'q_0': 1.7}

    volume = np.linspace(0.6, 1.0, 100)
    grueneisen_slb = np.empty_like(volume)
    grueneisen_mgd = np.empty_like(volume)
    q_slb = np.empty_like(volume)
    q_mgd = np.empty_like(volume)

    slb_eos = slb.SLB2()
    mgd_eos = mgd.MGD2()

    # calculate its thermal properties
    for i in range(len(volume)):
        # call with dummy pressure and temperatures, they do not change it
        grueneisen_slb[i] = slb_eos.grueneisen_parameter(
            0., 0., volume[i] * perovskite.params['V_0'], perovskite.params)
        grueneisen_mgd[i] = mgd_eos.grueneisen_parameter(
            0., 0., volume[i] * perovskite.params['V_0'], perovskite.params)
        q_slb[i] = slb_eos.volume_dependent_q(
            1. / volume[i], perovskite.params)
        q_mgd[i] = perovskite.params['q_0']

    # compare with figure 7
    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig3.png')
    plt.imshow(fig1, extent=[0.6, 1.0, 0.35, 2.0], aspect='auto')
    plt.plot(volume, grueneisen_slb, 'g+', volume, grueneisen_mgd, 'b+')
    plt.plot(volume, q_slb, 'g+', volume, q_mgd, 'b+')
    plt.xlim(0.6, 1.0)
    plt.ylim(0.35, 2.0)
    plt.ylabel("Grueneisen parameter")
    plt.xlabel("Relative Volume V/V0")
    plt.title(
        "Comparing with Figure 3 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()


def check_slb_fig7_txt():
    """
    Calculates all values for forsterite and benchmarks with values from Stixrude and Lithgow-Bertelloni (personal communication)
    """
    forsterite = burnman.Mineral()
    forsterite.params = {'name': 'forsterite',
                         'V_0': 43.603e-6,
                         'K_0': 127.955e9,
                         'Kprime_0': 4.232,
                         'G_0': 81.6e9,
                         'Gprime_0': 1.4,
                         'molar_mass': .140695,
                         'n': 7.0,
                         'Debye_0': 809.183,
                         'grueneisen_0': .993,
                         'q_0': 2.093,
                         'F_0': -1.1406e5,
                         'eta_s_0': 2.364}
    forsterite.set_method('slb3')

    data = np.loadtxt(
        "../../burnman/data/input_minphys/slb_fig7.txt", skiprows=2)

    temperature = np.array(data[:, 2])
    pressure = np.array(data[:, 0])
    rho = np.array(data[:, 3])
    frho = np.empty_like(rho)
    rho_comp = np.empty_like(rho)
    Kt = np.array(data[:, 4])
    fKt = np.empty_like(Kt)
    Kt_comp = np.empty_like(Kt)
    Ks = np.array(data[:, 5])
    fKs = np.empty_like(Ks)
    Ks_comp = np.empty_like(Ks)
    G = np.array(data[:, 6])
    fG = np.empty_like(G)
    G_comp = np.empty_like(G)
    VB = np.array(data[:, 7])
    fVB = np.empty_like(VB)
    VB_comp = np.empty_like(VB)
    VS = np.array(data[:, 8])
    fVS = np.empty_like(VS)
    VS_comp = np.empty_like(VS)
    VP = np.array(data[:, 9])
    fVP = np.empty_like(VP)
    VP_comp = np.empty_like(VP)
    vol = np.array(data[:, 10])
    fvol = np.empty_like(vol)
    vol_comp = np.empty_like(vol)
    alpha = np.array(data[:, 11])
    falpha = np.empty_like(alpha)
    alpha_comp = np.empty_like(alpha)
    Cp = np.array(data[:, 12])
    fCp = np.empty_like(Cp)
    Cp_comp = np.empty_like(Cp)
    gr = np.array(data[:, 13])
    gr_comp = np.empty_like(gr)
    gibbs = np.array(data[:, 14])
    gibbs_comp = np.empty_like(gibbs)
    entropy = np.array(data[:, 15])
    entropy_comp = np.empty_like(gibbs)
    enthalpy = np.array(data[:, 16])
    enthalpy_comp = np.empty_like(gibbs)

    for i in range(len(temperature)):
        forsterite.set_state(pressure[i], temperature[i])
        rho_comp[i] = 100. * (forsterite.density / 1000. - rho[i]) / rho[i]
        Kt_comp[i] = 100. * (
            forsterite.isothermal_bulk_modulus / 1.e9 - Kt[i]) / Kt[i]
        Ks_comp[i] = 100. * (
            forsterite.adiabatic_bulk_modulus / 1.e9 - Ks[i]) / Ks[i]
        G_comp[i] = 100. * (forsterite.shear_modulus / 1.e9 - G[i]) / G[i]
        VB_comp[i] = 100. * (forsterite.v_phi / 1000. - VB[i]) / VB[i]
        VS_comp[i] = 100. * (forsterite.v_s / 1000. - VS[i]) / VS[i]
        VP_comp[i] = 100. * (forsterite.v_p / 1000. - VP[i]) / VP[i]
        vol_comp[i] = 100. * (forsterite.molar_volume * 1.e6 - vol[i]) / vol[i]
        alpha_comp[i] = 100. * (
            forsterite.thermal_expansivity / 1.e-5 - alpha[i]) / (alpha[-1])
        Cp_comp[i] = 100. * (forsterite.heat_capacity_p /
                             forsterite.params['molar_mass'] / 1000. - Cp[i]) / (Cp[-1])
        gr_comp[i] = (forsterite.grueneisen_parameter - gr[i]) / gr[i]
        gibbs_comp[i] = 100. * (
            forsterite.molar_gibbs / 1.e6 - gibbs[i]) / gibbs[i]
        entropy_comp[i] = 100. * (
            forsterite.molar_entropy - entropy[i]) / (entropy[i] if entropy[i] != 0. else 1.)
        enthalpy_comp[i] = 100. * (
            forsterite.molar_enthalpy / 1.e6 - enthalpy[i]) / (enthalpy[i] if enthalpy[i] != 0. else 1.)

    plt.plot(temperature, rho_comp, label=r'$\rho$')
    plt.plot(temperature, Kt_comp, label=r'$K_S$')
    plt.plot(temperature, Ks_comp, label=r'$K_T$')
    plt.plot(temperature, G_comp, label=r'$G$')
    plt.plot(temperature, VS_comp, label=r'$V_S$')
    plt.plot(temperature, VP_comp, label=r'$V_P$')
    plt.plot(temperature, VB_comp, label=r'$V_\phi$')
    plt.plot(temperature, vol_comp, label=r'$V$')
    plt.plot(temperature, alpha_comp, label=r'$\alpha$')
    plt.plot(temperature, Cp_comp, label=r'$c_P$')
    plt.plot(temperature, gr_comp, label=r'$\gamma$')
    plt.plot(temperature, gibbs_comp, label=r'Gibbs')
    plt.plot(temperature, enthalpy_comp, label=r'Enthalpy')
    plt.plot(temperature, entropy_comp, label=r'Entropy')

    plt.xlim([0, 2750])
    plt.ylim([-0.001, 0.001])
    plt.xticks([0, 800, 1600, 2200])
    plt.xlabel("Temperature (K)")
    plt.ylabel("Percent Difference from HeFESTo")
    plt.legend(loc="center right")
#    plt.savefig("output_figures/benchmark1.pdf")
    plt.show()


def check_slb_fig7():
    """
    Calculates all values for forsterite and benchmarks with figure 7 from Stixrude and Lithgow-Bertelloni (2005)
    """
    forsterite = burnman.Mineral()
    forsterite.params = {'name': 'forsterite',
                         'V_0': 43.60e-6,
                         'K_0': 128.0e9,
                         'Kprime_0': 4.2,
                         'G_0': 82.0e9,
                         'Gprime_0': 1.4,
                         'n': 7.0,
                         'molar_mass': .140695,
                         'Debye_0': 809.,
                         'grueneisen_0': .99,
                         'q_0': 2.1,
                         'eta_s_0': 2.4}
    forsterite.set_method('slb3')

    temperature = np.linspace(0., 2000., 200)
    volume = np.empty_like(temperature)
    bulk_modulus = np.empty_like(temperature)
    shear_modulus = np.empty_like(temperature)
    heat_capacity = np.empty_like(temperature)

    pressure = 1.0e5
    forsterite.set_state(pressure, 300.)
    Ks_0 = forsterite.adiabatic_bulk_modulus

    # calculate its thermal properties
    for i in range(len(temperature)):
        forsterite.set_state(pressure, temperature[i])
        volume[i] = forsterite.molar_volume / forsterite.params['V_0']
        bulk_modulus[i] = forsterite.adiabatic_bulk_modulus / Ks_0
        shear_modulus[i] = forsterite.shear_modulus / forsterite.params['G_0']
        heat_capacity[i] = forsterite.heat_capacity_p / forsterite.params['n']

    # compare with figure 7
    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig7_vol.png')
    plt.imshow(fig1, extent=[0, 2200, 0.99, 1.08], aspect='auto')
    plt.plot(temperature, volume, 'g+')
    plt.ylim(0.99, 1.08)
    plt.xlim(0, 2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Volume V/V0")
    plt.title(
        "Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig7_Cp.png')
    plt.imshow(fig1, extent=[0, 2200, 0., 70.], aspect='auto')
    plt.plot(temperature, heat_capacity, 'g+')
    plt.ylim(0, 70)
    plt.xlim(0, 2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat Capacity Cp")
    plt.title(
        "Comparing with adiabatic_bulk_modulus7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig7_K.png')
    plt.imshow(fig1, extent=[0, 2200, 0.6, 1.02], aspect='auto')
    plt.plot(temperature, bulk_modulus, 'g+')
    plt.ylim(0.6, 1.02)
    plt.xlim(0, 2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Bulk Modulus K/K0")
    plt.title(
        "Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()

    fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig7_G.png')
    plt.imshow(fig1, extent=[0, 2200, 0.6, 1.02], aspect='auto')
    plt.plot(temperature, shear_modulus, 'g+')
    plt.ylim(0.6, 1.02)
    plt.xlim(0, 2200)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Shear Modulus G/G0")
    plt.title(
        "Comparing with Figure 7 of Stixrude and Lithgow-Bertelloni (2005)")
    plt.show()


def check_averaging():
    """
    Reproduce Figure 1a from Watt et. al. 1976 to check the Voigt, Reuss,
    Voigt-Reuss-Hill, and Hashin-Shtrikman bounds for an elastic composite
    """
    voigt = burnman.averaging_schemes.Voigt()
    reuss = burnman.averaging_schemes.Reuss()
    voigt_reuss_hill = burnman.averaging_schemes.VoigtReussHill()
    hashin_shtrikman_upper = burnman.averaging_schemes.HashinShtrikmanUpper()
    hashin_shtrikman_lower = burnman.averaging_schemes.HashinShtrikmanLower()

    # create arrays for sampling in volume fraction
    volumes = np.linspace(0.0, 1.0, 100)
    v_bulk_modulus = np.empty_like(volumes)
    v_shear_modulus = np.empty_like(volumes)
    r_bulk_modulus = np.empty_like(volumes)
    r_shear_modulus = np.empty_like(volumes)
    vrh_bulk_modulus = np.empty_like(volumes)
    vrh_shear_modulus = np.empty_like(volumes)
    hsu_bulk_modulus = np.empty_like(volumes)
    hsu_shear_modulus = np.empty_like(volumes)
    hsl_bulk_modulus = np.empty_like(volumes)
    hsl_shear_modulus = np.empty_like(volumes)

    # MgO bulk and shear moduli taken from Landolt-Boernstein
    # - Group III Condensed Matter Volume 41B, 1999, pp 1-3
    K2 = 152.  # Bulk modulus, GPa
    G2 = 155.  # Shear modulus, GPa

    # AgCl bulk and shear moduli (estimated from plot)
    G1 = G2 * 0.07
    K1 = K2 * 0.27

    for i in range(len(volumes)):
        v_bulk_modulus[i] = voigt.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        v_shear_modulus[i] = voigt.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        r_bulk_modulus[i] = reuss.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        r_shear_modulus[i] = reuss.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        vrh_bulk_modulus[i] = voigt_reuss_hill.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        vrh_shear_modulus[i] = voigt_reuss_hill.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        hsu_bulk_modulus[i] = hashin_shtrikman_upper.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        hsu_shear_modulus[i] = hashin_shtrikman_upper.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        hsl_bulk_modulus[i] = hashin_shtrikman_lower.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        hsl_shear_modulus[i] = hashin_shtrikman_lower.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

    fig = mpimg.imread('../../burnman/data/input_figures/watt_1976_a1.png')
    plt.imshow(fig, extent=[0, 1.0, 0.25, 1.0], aspect='auto')
    plt.plot(volumes, v_bulk_modulus / K2, 'g-')
    plt.plot(volumes, r_bulk_modulus / K2, 'g-')
    plt.plot(volumes, vrh_bulk_modulus / K2, 'g-')
    plt.plot(volumes, hsu_bulk_modulus / K2, 'g-')
    plt.plot(volumes, hsl_bulk_modulus / K2, 'g-')
    plt.ylim(0.25, 1.00)
    plt.xlim(0, 1.0)
    plt.xlabel("Volume fraction")
    plt.ylabel("Averaged bulk modulus")
    plt.title("Comparing with Figure 1 of Watt et al 1976")
    plt.show()

    fig = mpimg.imread('../../burnman/data/input_figures/watt_1976_a2.png')
    plt.imshow(fig, extent=[0, 1.0, 0.0, 1.0], aspect='auto')
    plt.plot(volumes, v_shear_modulus / G2, 'g-')
    plt.plot(volumes, r_shear_modulus / G2, 'g-')
    plt.plot(volumes, vrh_shear_modulus / G2, 'g-')
    plt.plot(volumes, hsu_shear_modulus / G2, 'g-')
    plt.plot(volumes, hsl_shear_modulus / G2, 'g-')
    plt.ylim(0.0, 1.00)
    plt.xlim(0, 1.0)
    plt.xlabel("Volume fraction")
    plt.ylabel("Averaged shear modulus")
    plt.title("Comparing with Figure 1 of Watt et al 1976")
    plt.show()

    # also check against some numerical values given in Berryman (1995) for
    # porous glass
    K = 46.3
    G = 30.5
    # the value for porosity=0.46 in the table appears to be a typo.  Remove
    # it here
    porosity = np.array(
        [0.0, 0.05, 0.11, 0.13, 0.25, 0.33, 0.36, 0.39, 0.44, 0.50, 0.70])
    berryman_bulk_modulus = np.array(
        [46.3, 41.6, 36.6, 35.1, 27.0, 22.5, 21.0, 19.6, 17.3, 14.8, 7.7])  # 15.5 probably a typo?
    hsu_bulk_modulus_vals = np.empty_like(porosity)
    for i in range(len(porosity)):
        hsu_bulk_modulus_vals[i] = hashin_shtrikman_upper.average_bulk_moduli(
            [porosity[i], 1.0 - porosity[i]], [0.0, K], [0.0, G])
    for i in range(len(volumes)):
        hsu_bulk_modulus[i] = hashin_shtrikman_upper.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [0.0, K], [0.0, G])
    fig = mpimg.imread('../../burnman/data/input_figures/berryman_fig4.png')
    plt.imshow(fig, extent=[0, 1.0, 0.0, 50.0], aspect='auto')
    plt.plot(volumes, hsu_bulk_modulus, 'g-')
    plt.scatter(porosity, hsu_bulk_modulus_vals, c='r')
    plt.scatter(porosity, berryman_bulk_modulus, c='y')
    plt.ylim(0.0, 50.0)
    plt.xlim(0, 1.0)
    plt.xlabel("Porosity")
    plt.ylabel("Averaged bulk modulus")
    plt.title("Comparing with Figure 4 of Berryman (1995)")
    plt.show()


def check_averaging_2():
    """
    Reproduce Figure 1 from Hashin and Shtrikman (1963) to check the
    Hashin-Shtrikman bounds for an elastic composite
    """

    hashin_shtrikman_upper = burnman.averaging_schemes.HashinShtrikmanUpper()
    hashin_shtrikman_lower = burnman.averaging_schemes.HashinShtrikmanLower()

    # create arrays for sampling in volume fraction
    volumes = np.linspace(0.0, 1.0, 100)
    hsu_bulk_modulus = np.empty_like(volumes)
    hsu_shear_modulus = np.empty_like(volumes)
    hsl_bulk_modulus = np.empty_like(volumes)
    hsl_shear_modulus = np.empty_like(volumes)

    # These values are from Hashin and Shtrikman (1963)
    K1 = 25.0
    K2 = 60.7
    G1 = 11.5
    G2 = 41.8

    for i in range(len(volumes)):
        hsu_bulk_modulus[i] = hashin_shtrikman_upper.average_bulk_moduli(
            [1.0 - volumes[i], volumes[i]], [K1, K2], [G1, G2])
        hsu_shear_modulus[i] = hashin_shtrikman_upper.average_shear_moduli(
            [1.0 - volumes[i], volumes[i]], [K1, K2], [G1, G2])

        hsl_bulk_modulus[i] = hashin_shtrikman_lower.average_bulk_moduli(
            [1.0 - volumes[i], volumes[i]], [K1, K2], [G1, G2])
        hsl_shear_modulus[i] = hashin_shtrikman_lower.average_shear_moduli(
            [1.0 - volumes[i], volumes[i]], [K1, K2], [G1, G2])

    fig = mpimg.imread(
        '../../burnman/data/input_figures/Hashin_Shtrikman_1963_fig1_K.png')
    plt.imshow(fig, extent=[0, 1.0, 1.1, K2 + 0.3], aspect='auto')
    plt.plot(volumes, hsu_bulk_modulus, 'g-')
    plt.plot(volumes, hsl_bulk_modulus, 'g-')
    plt.ylim(K1, K2)
    plt.xlim(0, 1.0)
    plt.xlabel("Volume fraction")
    plt.ylabel("Averaged bulk modulus")
    plt.title("Comparing with Figure 1 of Hashin and Shtrikman (1963)")
    plt.show()

    fig = mpimg.imread(
        '../../burnman/data/input_figures/Hashin_Shtrikman_1963_fig2_G.png')
    plt.imshow(fig, extent=[0, 1.0, 0.3, G2], aspect='auto')
    plt.plot(volumes, hsu_shear_modulus, 'g-')
    plt.plot(volumes, hsl_shear_modulus, 'g-')
    plt.ylim(G1, G2)
    plt.xlim(0, 1.0)
    plt.xlabel("Volume fraction")
    plt.ylabel("Averaged shear modulus")
    plt.title("Comparing with Figure 2 of Hashin and Shtrikman (1963)")
    plt.show()


def check_averaging_3():
    """
    Reproduce Figure 3 from Avseth et al. (2010) to check the Voigt, Reuss,
    Voigt-Reuss-Hill, and Hashin-Shtrikman bounds for an elastic composite
    """
    voigt = burnman.averaging_schemes.Voigt()
    reuss = burnman.averaging_schemes.Reuss()
    voigt_reuss_hill = burnman.averaging_schemes.VoigtReussHill()
    hashin_shtrikman_upper = burnman.averaging_schemes.HashinShtrikmanUpper()
    hashin_shtrikman_lower = burnman.averaging_schemes.HashinShtrikmanLower()

    # create arrays for sampling in volume fraction
    volumes = np.linspace(0.0, 1.0, 100)
    v_bulk_modulus = np.empty_like(volumes)
    v_shear_modulus = np.empty_like(volumes)
    r_bulk_modulus = np.empty_like(volumes)
    r_shear_modulus = np.empty_like(volumes)
    vrh_bulk_modulus = np.empty_like(volumes)
    vrh_shear_modulus = np.empty_like(volumes)
    hsu_bulk_modulus = np.empty_like(volumes)
    hsu_shear_modulus = np.empty_like(volumes)
    hsl_bulk_modulus = np.empty_like(volumes)
    hsl_shear_modulus = np.empty_like(volumes)
    hs_av_bulk_modulus = np.empty_like(volumes)
    hs_av_shear_modulus = np.empty_like(volumes)

    # Quartz bulk and shear moduli
    K2 = 37.
    G2 = 45.

    # Fluid bulk and shear moduli
    G1 = 0.00001
    K1 = 2.35

    for i in range(len(volumes)):
        v_bulk_modulus[i] = voigt.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        v_shear_modulus[i] = voigt.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        r_bulk_modulus[i] = reuss.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        r_shear_modulus[i] = reuss.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        vrh_bulk_modulus[i] = voigt_reuss_hill.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        vrh_shear_modulus[i] = voigt_reuss_hill.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        hsu_bulk_modulus[i] = hashin_shtrikman_upper.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        hsu_shear_modulus[i] = hashin_shtrikman_upper.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        hsl_bulk_modulus[i] = hashin_shtrikman_lower.average_bulk_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])
        hsl_shear_modulus[i] = hashin_shtrikman_lower.average_shear_moduli(
            [volumes[i], 1.0 - volumes[i]], [K1, K2], [G1, G2])

        hs_av_bulk_modulus[i] = 0.5 * hsl_bulk_modulus[
            i] + 0.5 * hsu_bulk_modulus[i]
        hs_av_shear_modulus[i] = 0.5 * hsl_shear_modulus[
            i] + 0.5 * hsu_shear_modulus[i]

    fig = mpimg.imread(
        '../../burnman/data/input_figures/Avseth_et_al_2010_fig3_K.png')
    plt.imshow(fig, extent=[0, 1.0, 0., 40.0], aspect='auto')
    plt.plot(volumes, v_bulk_modulus, 'g-')
    plt.plot(volumes, r_bulk_modulus, 'g-')
    plt.plot(volumes, vrh_bulk_modulus, 'g-')
    plt.plot(volumes, hsu_bulk_modulus, 'g-')
    plt.plot(volumes, hsl_bulk_modulus, 'g-')
    plt.plot(volumes, hs_av_bulk_modulus, 'g-')
    plt.ylim(0., 40.00)
    plt.xlim(0., 1.0)
    plt.xlabel("Volume fraction")
    plt.ylabel("Averaged bulk modulus")
    plt.title("Comparing with Figure 3 of Avseth et al., 2010")
    plt.show()


if __name__ == "__main__":
    check_averaging()
    check_averaging_2()
    check_averaging_3()
    check_birch_murnaghan()
    check_birch_murnaghan_4th()
    check_vinet()
    check_slb_fig7()
    check_slb_fig3()
    check_mgd_shim_duffy_kenichi()
    check_mgd_fei_mao_shu_hu()
    check_slb_fig7_txt()
