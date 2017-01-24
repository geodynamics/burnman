# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import warnings
import scipy.integrate
import matplotlib.pyplot as plt
import pkgutil

from . import tools
from . import constants
from . import seismic
from . import geotherm

def write_axisem_input(rock, min_depth=670.e3, max_depth=2890.e3, T0= 1900, filename='axisem_burnmantestrock.txt',
                       axisem_ref='axisem_prem_ani_noocean.txt', plotting=False):
    """
    Writing velocities and densities to AXISEM (www.axisem.info) input file
    Default is set to replacing the lower mantle with the BurnMan rock
    Note:
        - This implementation uses PREM to convert from depths to pressures to compute at
        - This implementation assumes an adiabatic temperature profile, only T0 at min_depth can be set
        - Currently, it only honors the discontinuities already in the synthetic input file, so it is best 
            to only replace certain layers with burnman values (this should be improved in the future).
            
            
    Parameters
    ----------
    rock  :  burnman.Composite()
        Composition to implement in the model
    min_depth : float
        minimum depth to replace model (m) (default = 670 km)
    max_depth : float 
        minimum depth to replace model (m) (default = 2890 km)
    T0 : float
        Anchor temperature at min_depth for adiabatic profile (K) (default=1900)
    filename: string
        Output filename (default ='axisem_burnmantestrock.txt')
    axisem_ref: string
        Input filename (in burnman/data/input_seismic/) (default = 'axisem_prem_ani_noocean.txt')
    plotting: Boolean
        True means plot of the old model and replaced model will be shown (default = False)

    """


    
    # Load reference input
    datastream = pkgutil.get_data('burnman', 'data/input_seismic/' + axisem_ref)
    lines = [line.strip()
                 for line in datastream.decode('ascii').split('\n') if line.strip()]
    table = []
    for line in lines[18:]:
        numbers = np.fromstring(line, sep=' ')
        if len(numbers)>0:
            if line[0] != "#" and line[0] != "%":
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

    # Cutting out range to input in Axisem reference file (currently the lower mantle)
    indrange = [x for x in range(len(ref_depth)) if ref_depth[
        x] > min_depth and ref_depth[x] < max_depth]
    # pad both ends to include up to discontinuity, bit of a hack...
    indrange.insert(0, indrange[0] - 1)
    indrange.append(indrange[-1] + 1)
    # Invert depthrange so adiabatic computations work!
    depthrange = ref_depth[indrange]

    # convert depths to pressures
    pressures = seismic.PREM().pressure(depthrange)

    # Computing adiabatic temperatures. T0 is an input parameter!
    T0 = T0  # K
    temperatures = geotherm.adiabatic(pressures, T0, rock)


    
    print("Calculations are done for:")
    rock.debug_print()

    rock_vp, rock_vs, rock_rho = rock.evaluate(
        ['v_p', 'v_s', 'density'], pressures, temperatures)
    discontinuity =0
    # WRITE OUT FILE
    f = open(filename, 'w')
    print('Writing ' + filename + ' ...')
    f.write('# Input file '+ filename +' for AXISEM created using BurnMan, replacing ' + axisem_ref+ ' between ' +str(np.round(min_depth/1.e3)) + ' and ' + str(np.round(max_depth /1.e3)) +' km \n')
    f.write('NAME ' + filename + '\n')
    for line in lines[2:18]:
        f.write(line[:-1] + '\n')
    for i in range(indrange[0]):
        if i>0 and ref_radius[i] ==ref_radius[i-1]:
                discontinuity = discontinuity + 1
                f.write('#          Discontinuity   ' +str(discontinuity) + ', depth:    '+ str(np.round(ref_depth[i]/1.e3,decimals=2)) +' km \n')
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
            (ref_radius[i], ref_density[i], ref_vpv[i], ref_vsv[i], ref_Qk[i],
                ref_Qmu[i], ref_vph[i], ref_vsh[i], ref_eta[i]))


    for i in range(indrange[0], indrange[-1]):
        ind2 = -1 + i - indrange[0]
        if ref_radius[i] ==ref_radius[i-1]:
            discontinuity = discontinuity + 1
            f.write('#          Discontinuity   '+ str(discontinuity) + ', depth:    '+ str(np.round(ref_depth[i]/1.e3,decimals=2))+' km \n')
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
            (ref_radius[i], rock_rho[ind2], rock_vp[ind2], rock_vs[ind2], ref_Qk[i],
                ref_Qmu[i], rock_vp[ind2], rock_vs[ind2], ref_eta[i]))

    for i in range(indrange[-1], len(ref_radius)):
        if ref_radius[i] ==ref_radius[i-1]:
            discontinuity = discontinuity + 1
            f.write('#          Discontinuity   ' +str(discontinuity) + ', depth:    '+ str(np.round(ref_depth[i]/1.e3,decimals=2))+' km \n')
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
                (ref_radius[i], ref_density[i], ref_vpv[i], ref_vsv[i], ref_Qk[i],
                 ref_Qmu[i], ref_vph[i], ref_vsh[i], ref_eta[i]))

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

        plt.title(filename + ' = ' + axisem_ref + ' replaced  between ' +
                  str(min_depth / 1.e3) + ' and ' + str(max_depth / 1.e3) + ' km')
        plt.legend(loc='lower right')
        plt.show()


def write_mineos_input(rock, min_depth=670.e3, max_depth=2890.e3, T0 = 1900, filename='mineos_burnmantestrock.txt',
                       mineos_ref='mineos_prem_noocean.txt', plotting=False):
    """
    Writing velocities and densities to Mineos (https://geodynamics.org/cig/software/mineos/) input file
    Default is set to replacing the lower mantle with the BurnMan rock
    Note:
        - This implementation uses PREM to convert from depths to pressures to compute at
        - This implementation assumes an adiabatic temperature profile, only T0 at min_depth can be set
        - Currently, it only honors the discontinuities already in the synthetic input file, so it is best 
            to only replace certain layers with burnman values (this should be improved in the future).


    Parameters
    ----------
    rock  :  burnman.Composite()
    Composition to implement in the model
    min_depth : float
    minimum depth to replace model (m) (default = 670 km)
    max_depth : float
    minimum depth to replace model (m) (default = 2890 km)
    T0 : float
    Anchor temperature at min_depth for adiabatic profile (K) (default=1900)
    filename: string
    Output filename (default ='mineos_burnmantestrock.txt')
    axisem_ref: string
    Input filename (in burnman/data/input_seismic/) (default = 'mineos_prem_noocean.txt')
    plotting: Boolean
    True means plot of the old model and replaced model will be shown (default = False)
    
    """
    
    
    # Load reference input
    datastream = pkgutil.get_data('burnman', 'data/input_seismic/' + mineos_ref)
    lines = [line.strip()
             for line in datastream.decode('ascii').split('\n') if line.strip()]
    table=[]
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
    pressures = seismic.PREM().pressure(depthrange)

    # Computing adiabatic temperatures. T0 is a choice!
    T0 = T0  # K
    temperatures = geotherm.adiabatic(pressures, T0, rock)

    
    print("Calculations are done for:")
    rock.debug_print()

    rock_vp, rock_vs, rock_rho = rock.evaluate(
        ['v_p', 'v_s', 'density'], pressures, temperatures)

    # WRITE OUT FILE
    f = open(filename , 'w')
    print('Writing ' + filename + ' ...')
    f.write(lines[0][:-2] + ' +  ' + filename + '\n')
    for line in lines[1:3]:
        f.write(line[:-2] + '\n')
    for i in range(indrange[0]):
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
                (ref_radius[i], ref_density[i], ref_vpv[i], ref_vsv[i], ref_Qk[i],
                 ref_Qmu[i], ref_vph[i], ref_vsh[i], ref_eta[i]))

    for i in range(indrange[0], indrange[-1]):
        ind2 = -1 - i + indrange[0]
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
                (ref_radius[i], rock_rho[ind2], rock_vp[ind2], rock_vs[ind2], ref_Qk[i],
                 ref_Qmu[i], rock_vp[ind2], rock_vs[ind2], ref_eta[i]))


    for i in range(indrange[-1], len(ref_radius)):
        f.write(
            '%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' %
                (ref_radius[i], ref_density[i], ref_vpv[i], ref_vsv[i], ref_Qk[i],
                 ref_Qmu[i], ref_vph[i], ref_vsh[i], ref_eta[i]))

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

        plt.title(filename + ' = ' + mineos_ref + ' replaced  between ' +
                  str(min_depth / 1.e3) + ' and ' + str(max_depth / 1.e3) + ' km')
        plt.legend(loc='lower right')
        plt.show()
