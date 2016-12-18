# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
_________________________________________________________________________
This example script produces an input file to compute normal modes in
Mineos for a lower mantle with pyrolitic or chondritic composition.
Thanks for input from Jessica Irving.
_________________________________________________________________________


"""
# Import supporting libraries
# Imports to be compatible with Python2 and Python3
from __future__ import absolute_import
from __future__ import print_function
import os
import sys  # Library used to interact with your operating system
import numpy as np  # Library used for general array
import matplotlib.pyplot as plt  # Library used for plotting
# Import BurnMan
sys.path.insert(1, os.path.abspath('../'))  # add path to burnman
import burnman
from burnman import minerals  # import mineral library seperately


def write_mineos_input(rock, min_depth=670.e3, max_depth=2890.e3, name='burnmantestrock',
                       mineos_ref='prem_noocean', plotting=False):
    """
    _________________________________________________________________________
    Writing velocities and densities to Mineos input file
    _________________________________________________________________________

    """
    # Load reference input
    lines = open('../burnman/data/input_seismic/mineos_' + mineos_ref + '.txt').readlines()
    table = []
    for line in lines[3:]:
        numbers = np.fromstring(line, sep=' ')
        table.append(numbers)
    table = np.array(table)
    ref_radius = table[:, 0]
    ref_depth = 6371.e3 - ref_radius
    ref_density = table[:, 1]
    ref_vpv = table[:, 2]
    ref_vsv = table[:, 3]
    ref_Qk = table[:, 4]
    ref_Qmu = table[:, 5]
    ref_vph = table[:, 6]
    ref_vsh = table[:, 7]
    ref_eta = table[:, 8]

    # Cutting out range to input in Mineos (currently the lower mantle)
    indrange = [x for x in range(len(ref_depth)) if ref_depth[
        x] > min_depth and ref_depth[x] < max_depth]
    # pad both ends to include up to discontinuity, bit of a hack...
    indrange.insert(0, indrange[0] - 1)
    indrange.append(indrange[-1] + 1)
    # Invert depthrange so adiabatic computations work!
    depthrange = ref_depth[indrange][::-1]

    # convert depths to pressures
    pressures = burnman.seismic.PREM().pressure(depthrange)

    # Computing adiabatic temperatures. T0 is a choice!
    T0 = 1900  # K
    temperatures = burnman.geotherm.adiabatic(pressures, T0, rock)
    # An alternative is the Brown+Shankland (1981)
    # geotherm for mapping pressure to temperature.
    # To use this include the line below.
    #temperature = burnman.geotherm.brown_shankland(pressure)
    
    print("Calculations are done for:")
    rock.debug_print()

    rock_vp, rock_vs, rock_rho = rock.evaluate(
        ['v_p', 'v_s', 'density'], pressures, temperatures)

    # WRITE OUT FILE
    f = open('mineos_' + name + '.txt', 'w')
    f.write(lines[0][:-2] + ' +  ' + name + '\n')
    for line in lines[1:3]:
        f.write(line[:-2] + '\n')
    for i in range(indrange[0]):
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
            (ref_radius[i],
             ref_density[i],
             ref_vpv[i],
             ref_vsv[i],
             ref_Qk[i],
                ref_Qmu[i],
                ref_vph[i],
                ref_vsh[i],
                ref_eta[i]))
    for i in range(indrange[0], indrange[-1]):
        ind2 = -1 - i + indrange[0]
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
            (ref_radius[i],
             rock_rho[ind2],
             rock_vp[ind2],
             rock_vs[ind2],
             ref_Qk[i],
                ref_Qmu[i],
                rock_vp[ind2],
                rock_vs[ind2],
                ref_eta[i]))
    for i in range(indrange[-1], len(ref_radius)):
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
            (ref_radius[i],
             ref_density[i],
             ref_vpv[i],
             ref_vsv[i],
             ref_Qk[i],
                ref_Qmu[i],
                ref_vph[i],
                ref_vsh[i],
                ref_eta[i]))
    f.close()

    if plotting:
        # plot vp
        plt.plot(ref_depth / 1.e3, ref_vph / 1.e3, color='g', linestyle='-', label='vp')
        plt.plot(depthrange / 1.e3, rock_vp / 1.e3, color='g', linestyle='-',
                 marker='o', markerfacecolor='g', markersize=1)

        # plot Vs
        plt.plot(ref_depth / 1.e3, ref_vsh / 1.e3, color='b', linestyle='-', label='vs')
        plt.plot(depthrange / 1.e3, rock_vs / 1.e3, color='b', linestyle='-',
                 marker='o', markerfacecolor='b', markersize=1)

        # plot density
        plt.plot(ref_depth / 1.e3, ref_density / 1.e3, color='r', linestyle='-', label='density')
        plt.plot(depthrange / 1.e3, rock_rho / 1.e3, color='r', linestyle='-',
                 marker='o', markerfacecolor='r', markersize=1)

        plt.title(mineos_ref + ' replaced with ' + name + ' between ' +
                  str(min_depth / 1.e3) + ' and ' + str(max_depth / 1.e3) + ' km')
        plt.legend(loc='lower right')
        plt.show()


if __name__ == "__main__":
    # We'll compute the velocities for different compositions at 20 points
    # within the lower mantle. Here we define the array of those depths.
    depths = np.linspace(750e3, 2700e3, 20)
    # We use PREM to convert these depths to pressure values.
    [pressures] = burnman.seismic.PREM().evaluate(['pressure'], depths)

    #-Defining the rock-
    #The object burnman.Composite expects two lists, one with the minerals
    #themselves and one with the molar fractions of the different minerals
    #making up the rock
    #Here we test the models for a Pyrolitic and Chondritic lower mantle.


    # Perovksite solide solution
    frac_mg = 0.94
    frac_fe = 0.06
    frac_al = 0.00
    mg_fe_perovskite = minerals.SLB_2011.mg_fe_perovskite()
    mg_fe_perovskite.set_composition([frac_mg, frac_fe, frac_al])

    # ferropericlase solid solution
    frac_mg = 0.8
    frac_fe = 0.2
    mg_fe_periclase = minerals.SLB_2011.ferropericlase()
    mg_fe_periclase.set_composition([frac_mg, frac_fe])

    # Ca Perovskite
    ca_perovskite = minerals.SLB_2011.ca_perovskite()

    # Pyrolitic composition
    pyr_pv = 0.75
    pyr_fp = 0.18
    pyr_capv = 0.07
    pyrolitic_mantle = burnman.Composite(
        [mg_fe_perovskite, mg_fe_periclase, ca_perovskite], [pyr_pv, pyr_fp, pyr_capv])

    # Chondritic composition
    chon_pv = 0.88
    chon_fp = 0.05
    chon_capv = 0.07
    chondritic_mantle = burnman.Composite(
        [mg_fe_perovskite, mg_fe_periclase, ca_perovskite], [chon_pv, chon_fp, chon_capv])




    print("Calculations are done for:")
    pyrolitic_mantle.debug_print()

    # writing mineos input!
    write_mineos_input(pyrolitic_mantle, name='pyrolite', plotting=True)
    write_mineos_input(chondritic_mantle, name='chondrite', plotting=True)
