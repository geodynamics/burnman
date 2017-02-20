from __future__ import absolute_import
from __future__ import print_function

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


import burnman
from burnman.minerals import SLB_2011
from scipy.optimize import fsolve

from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline


# Define fitting function to find the temperature along the isentrope
def isentrope(rock, pressures, entropy):
    def temperature_on_isentrope(args, S, P, rock):
        T = args[0]
        rock.set_state(P, T)
        return rock.S - S

    sol = [1600.]
    temperatures = np.empty_like(pressures)
    for i, P in enumerate(pressures):
        sol = fsolve(temperature_on_isentrope, sol,  args=(entropy, P, rock))
        temperatures[i] = sol[0]

    return temperatures

def smooth_isentrope(interp, pressures, entropy):
    def temperature_on_isentrope(args, S, P):
        T = args[0]
        return interp(P, T)[0] - S
    
    sol = [1600.]
    temperatures = np.empty_like(pressures)
    for i, P in enumerate(pressures):
        sol = fsolve(temperature_on_isentrope, sol,  args=(entropy, P))
        temperatures[i] = sol[0]

    return temperatures

def compute_pressure_gradient(pressures, densities):
    g0 = 9.81
    gravity = pressures * 0. + 10. # starting guess
    n_gravity_iterations = 5
    for i in xrange(n_gravity_iterations):    
        # Integrate the hydrostatic equation
        # Make a spline fit of densities as a function of pressures
        rhofunc = UnivariateSpline(pressures, densities)
        # Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline(pressures, gravity)
        
        # integrate the hydrostatic equation
        depths = np.ravel(odeint((lambda p, x: 1./(gfunc(x) * rhofunc(x))), 0.0, pressures))
        
        radii = 6371.e3 - depths
        
        rhofunc = UnivariateSpline(radii[::-1], densities[::-1])
        poisson = lambda p, x: 4.0 * np.pi * burnman.constants.G * rhofunc(x) * x * x
        gravity = np.ravel(odeint(poisson, g0*radii[0]*radii[0], radii))
        gravity = gravity / radii / radii
    return depths, gravity


# Declare the rock we want to use 
rock=burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')

# First we calculate the isentrope at a given potential temperature
potential_temperature = 1550. # K
n_points = 251
rock.set_state(1.e5, potential_temperature)
entropy = rock.S
pressures = np.linspace(1.e5, 25.e9, n_points)
temperatures = isentrope(rock, pressures, entropy)


volumes, densities, C_p, alphas, compressibilities = rock.evaluate(['V', 'rho',
                                                                    'heat_capacity_p',
                                                                    'thermal_expansivity',
                                                                    'isothermal_compressibility'],
                                                                   pressures,
                                                                   temperatures)

specific_heats = C_p / rock.params['molar_mass']
depths, gravity = compute_pressure_gradient(pressures, densities)

for i in xrange(1,7):
    plt.subplot(2, 4, i)
    plt.xlabel('Pressures (GPa)')

x = pressures/1.e9

plt.subplot(2, 4, 1)
plt.plot(x, temperatures)
plt.ylabel('Temperature (K)')
plt.subplot(2, 4, 2)
plt.plot(x, depths/1.e3)
plt.ylabel('Depths (km)')
plt.subplot(2, 4, 3)
plt.plot(x, gravity)
plt.ylabel('Gravity (m/s^2)')
plt.subplot(2, 4, 4)
plt.plot(x, densities)
plt.ylabel('Density (kg/m^3)')
plt.subplot(2, 4, 5)
plt.plot(x, specific_heats)
plt.ylabel('Cp (J/K/kg)')
plt.subplot(2, 4, 6)
plt.plot(x, alphas)
plt.ylabel('alpha (/K)')
plt.subplot(2, 4, 7)
plt.plot(x, compressibilities)
plt.ylabel('compressibilities (/Pa)')

grid_pressures = np.linspace(1.e5, 25.e9, 501)
grid_temperatures = np.linspace(1400., 2000., 101)
temperature_stdev = 10.


for pressure_stdev in [0., 5.e8]:

    S_interps = burnman.tools.interp_smoothed_property(rock, 'S',
                                                       grid_pressures,
                                                       grid_temperatures,
                                                       pressure_stdev,
                                                       temperature_stdev)
    interp_smoothed_S, interp_smoothed_dSdP, interp_smoothed_dSdT = S_interps
    
    V_interps = burnman.tools.interp_smoothed_property(rock, 'V',
                                                       grid_pressures,
                                                       grid_temperatures,
                                                       pressure_stdev,
                                                       temperature_stdev)
    
    interp_smoothed_V, interp_smoothed_dVdP, interp_smoothed_dVdT = V_interps

    temperatures = smooth_isentrope(interp_smoothed_S, pressures, entropy)
    densities = rock.evaluate(['rho'], pressures, temperatures)[0]
    depths, gravity = compute_pressure_gradient(pressures, densities)
    
    volumes = np.array([interp_smoothed_V(p, T)[0] for (p, T) in zip(*[pressures, temperatures])])
    dSdT = np.array([interp_smoothed_dSdT(p, T)[0] for (p, T) in zip(*[pressures, temperatures])])
    dVdT = np.array([interp_smoothed_dVdT(p, T)[0] for (p, T) in zip(*[pressures, temperatures])])
    dVdP = np.array([interp_smoothed_dVdP(p, T)[0] for (p, T) in zip(*[pressures, temperatures])])
    
    specific_heats_relaxed = temperatures * dSdT / rock.params['molar_mass']
    alphas_relaxed = dVdT / volumes
    compressibilities_relaxed = -dVdP / volumes


    plt.subplot(2, 4, 1)
    plt.plot(x, temperatures)
    plt.subplot(2, 4, 2)
    plt.plot(x, depths/1.e3)
    plt.subplot(2, 4, 3)
    plt.plot(x, gravity)
    plt.subplot(2, 4, 4)
    plt.plot(x, densities)
    plt.subplot(2, 4, 5)
    plt.plot(x, specific_heats_relaxed)
    plt.subplot(2, 4, 6)
    plt.plot(x, alphas_relaxed)
    plt.subplot(2, 4, 7)
    plt.plot(x, compressibilities_relaxed)
    
plt.show()




# depth, pressure, temperature, density, gravity, Cp (per kilo), thermal expansivity
np.savetxt('isentrope_properties.txt', X=np.array([depths, pressures, temperatures, densities, gravity, alphas, specific_heats, compressibilities]).T,
           header='POINTS: '+str(n_points)+' \ndepth (m), pressure (Pa), temperature (K), density (kg/m^3), gravity (m/s^2), thermal expansivity (/K), Cp (J/K/kg), beta (/Pa)',
           fmt='%.10e', delimiter='\t')
