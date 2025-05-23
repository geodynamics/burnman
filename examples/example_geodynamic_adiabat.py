# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_geodynamic_adiabat
----------------

This example script demonstrates how BurnMan can be used to
self-consistently calculate properties along 1D adiabatic profiles
for use in geodynamics simulations.

We use interrogate a PerplexMaterial for material properties
as a function of pressure and temperature, and calculate
both unrelaxed (short timescale) and relaxed (long timescale)
properties. The latter are particularly important for
convection studies, where reactions are fast compared with
timescales of convection.

Finally, we show how burnman can be used to smooth entropy and
volume in order to create smoothly varying relaxed properties.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.PerplexMaterial`
* :func:`burnman.Material.evaluate`
* :func:`burnman.utils.math.interp_smoothed_array_and_derivatives`


*Demonstrates:*

* creation of a PerpleX material
* calculation of relaxed and unrelaxed thermodynamic properties
* smoothing of thermodynamic properties
* self consistent 1D geophysical profile generation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline

import burnman
from burnman.utils.math import interp_smoothed_array_and_derivatives


if __name__ == "__main__":
    #  Define fitting function to find the temperature along the isentrope
    def isentrope(rock, pressures, entropy, T_guess):
        def _deltaS(T, S, P, rock):
            rock.set_state(P, T)
            return rock.S - S

        sol = T_guess
        temperatures = np.empty_like(pressures)
        for i, P in enumerate(pressures):
            sol = brentq(
                _deltaS, rock.bounds[1][0], rock.bounds[1][1], args=(entropy, P, rock)
            )
            temperatures[i] = sol

        return temperatures

    # Define function to find an isentrope given a
    # 2D entropy interpolation function
    def interp_isentrope(interp, pressures, entropy, T_guess):
        def _deltaS(args, S, P):
            T = args[0]
            return interp((P, T)) - S

        sol = [T_guess]
        temperatures = np.empty_like(pressures)
        for i, P in enumerate(pressures):
            sol = fsolve(_deltaS, sol, args=(entropy, P))
            temperatures[i] = sol[0]

        return temperatures

    # Define function to self consistently calculate depth and gravity profiles
    # from pressure and density profiles.
    def compute_depth_gravity_profiles(
        pressures, densities, surface_gravity, outer_radius
    ):
        gravity = [surface_gravity] * len(pressures)  # starting guess
        n_gravity_iterations = 5
        for i in range(n_gravity_iterations):
            # Integrate the hydrostatic equation
            # Make a spline fit of densities as a function of pressures
            rhofunc = UnivariateSpline(pressures, densities)
            # Make a spline fit of gravity as a function of depth
            gfunc = UnivariateSpline(pressures, gravity)

            # integrate the hydrostatic equation
            depths = np.ravel(
                odeint((lambda p, x: 1.0 / (gfunc(x) * rhofunc(x))), 0.0, pressures)
            )

            radii = outer_radius - depths

            rhofunc = UnivariateSpline(radii[::-1], densities[::-1])
            poisson = (
                lambda p, x: 4.0 * np.pi * burnman.constants.G * rhofunc(x) * x * x
            )
            gravity = np.ravel(
                odeint(poisson, surface_gravity * radii[0] * radii[0], radii)
            )
            gravity = gravity / radii / radii
        return depths, gravity

    # BEGIN USER INPUTS

    # Declare the rock we want to use
    rock = burnman.PerplexMaterial("../burnman/data/input_perplex/in23_1.tab")

    # First we calculate the isentrope at a given potential temperature
    potential_temperature = 1550.0  # K
    max_pressure = 25.0e9
    outer_radius = 6371.0e3
    surface_gravity = 9.81

    n_P_gridpoints = 501
    n_T_gridpoints = 101

    n_points = 301

    pressure_std_dev = 5.0e8
    temperature_smoothing_factor = 0.5
    truncate = 4.0

    # END USER INPUTS

    min_grid_pressure = rock.bounds[0][0]
    max_grid_pressure = rock.bounds[0][1]
    min_grid_temperature = rock.bounds[1][0]
    max_grid_temperature = rock.bounds[1][1]

    rock.set_state(1.0e5, potential_temperature)
    entropy = rock.S
    pressures = np.linspace(1.0e5, max_pressure, n_points)
    temperatures = isentrope(rock, pressures, entropy, potential_temperature)

    isentrope = UnivariateSpline(pressures, temperatures)

    # Properties can then be calculated along the isentrope
    properties = rock.evaluate(
        [
            "V",
            "rho",
            "molar_heat_capacity_p",
            "thermal_expansivity",
            "isothermal_compressibility_reuss",
            "p_wave_velocity",
            "shear_wave_velocity",
        ],
        pressures,
        temperatures,
    )
    (
        volumes,
        densities,
        C_p,
        alphas,
        compressibilities,
        p_wave_velocities,
        s_wave_velocities,
    ) = properties
    specific_heats = C_p / rock.params["molar_mass"]
    depths, gravity = compute_depth_gravity_profiles(
        pressures, densities, surface_gravity, outer_radius
    )

    x = pressures / 1.0e9

    plt.rcParams["figure.figsize"] = 16, 8  # inches
    fig = plt.figure()
    ax_T = fig.add_subplot(2, 4, 1)
    ax_T.plot(x, temperatures, label="unrelaxed")
    ax_T.set_ylabel("Temperature (K)")
    ax_T.set_xlabel("Pressures (GPa)")

    ax_z = fig.add_subplot(2, 4, 2)
    ax_z.plot(x, depths / 1.0e3)
    ax_z.set_ylabel("Depths (km)")
    ax_z.set_xlabel("Pressures (GPa)")

    ax_g = fig.add_subplot(2, 4, 3)
    ax_g.plot(x, gravity)
    ax_g.set_ylabel("Gravity (m/s^2)")
    ax_g.set_xlabel("Pressures (GPa)")

    ax_rho = fig.add_subplot(2, 4, 4)
    ax_rho.plot(x, densities)
    ax_rho.set_ylabel("Density (kg/m^3)")
    ax_rho.set_xlabel("Pressures (GPa)")

    ax_cp = fig.add_subplot(2, 4, 5)
    ax_cp.plot(x, specific_heats)
    ax_cp.set_ylabel("Cp (J/K/kg)")
    ax_cp.set_xlabel("Pressures (GPa)")

    ax_alpha = fig.add_subplot(2, 4, 6)
    ax_alpha.plot(x, alphas)
    ax_alpha.set_ylabel("alpha (/K)")
    ax_alpha.set_xlabel("Pressures (GPa)")

    ax_beta = fig.add_subplot(2, 4, 7)
    ax_beta.plot(x, compressibilities)
    ax_beta.set_ylabel("compressibilities (/Pa)")
    ax_beta.set_xlabel("Pressures (GPa)")

    ax_vs = fig.add_subplot(2, 4, 8)
    ax_vs.plot(x, p_wave_velocities, label="P")
    ax_vs.plot(x, s_wave_velocities, label="S")
    ax_vs.legend(loc="upper left")
    ax_vs.set_ylabel("Velocities (km/s)")
    ax_vs.set_xlabel("Pressures (GPa)")

    # Now let's calculate some relaxed material properties.
    # These are commonly used in geodynamic simulations, because
    # advective velocities due to convection are slow compared
    # with the velocity of reaction fronts at mantle temperatures

    # Before computing relaxed properties, we can optionally choose
    # to smooth the entropy and volume. In this way, we avoid peaks
    # in the relaxed thermal expansivity and compressibility which
    # tend to cause numerical problems for geodynamics software.

    # Let's first define a grid to calculate properties.
    # This grid needs to be reasonably densely sampled to
    # capture all the gradients. Obviously, the density of
    # sampling can be lower if we intend to smooth over long wavelengths
    grid_pressures = np.linspace(min_grid_pressure, max_grid_pressure, n_P_gridpoints)
    grid_temperatures = np.linspace(
        min_grid_temperature, max_grid_temperature, n_T_gridpoints
    )

    # Here we choose to smooth by convolving the entropy and volume
    # with a 2D Gaussian function with RMS widths of 0 or 0.5 GPa in pressure, and
    # a temperature smoothing of f*P_stdev*max(dT/dP), where f is a factor > 0.25.
    # The smoothing is truncated at the
    # 4 sigma level.

    for pressure_stdev in [0.0, pressure_std_dev]:
        unsmoothed_isentrope_temperatures = isentrope(grid_pressures)
        temperature_stdev = (
            temperature_smoothing_factor
            * pressure_stdev
            * np.max(np.abs(np.gradient(unsmoothed_isentrope_temperatures)))
            / (grid_pressures[1] - grid_pressures[0])
        )

        pp, TT = np.meshgrid(grid_pressures, grid_temperatures)
        mesh_shape = pp.shape
        pp = np.ndarray.flatten(pp)
        TT = np.ndarray.flatten(TT)

        # We could compute properties over the whole grid:
        # grid_entropies, grid_volumes = rock.evaluate(['S', 'V'], pp, TT)
        # However, we can save some time by computing only when temperature is close enough
        # to the unsmoothed isentrope to affect the smoothing.
        # The maximum temperature jump for this rock is about 50 K, so a reasonable Tmax is
        # ~50 + 4.*temperature_stdev. We pad a bit more (an extra 30 K) just to be sure.
        Tdiff_max = 50 + 30 + truncate * temperature_stdev
        grid_entropies = np.zeros_like(pp)
        grid_volumes = np.zeros_like(pp)
        Tdiff = np.abs(isentrope(pp) - TT)
        mask = [idx for idx, Td in enumerate(Tdiff) if Td < Tdiff_max]
        grid_entropies[mask], grid_volumes[mask] = rock.evaluate(
            ["S", "V"], pp[mask], TT[mask]
        )

        grid_entropies = grid_entropies.reshape(mesh_shape)
        grid_volumes = grid_volumes.reshape(mesh_shape)

        # Having defined the grid and calculated unsmoothed properties,
        # we now calculate the smoothed entropy and volume and derivatives with
        # respect to pressure and temperature.
        S_interps = interp_smoothed_array_and_derivatives(
            array=grid_entropies,
            x_values=grid_pressures,
            y_values=grid_temperatures,
            x_stdev=pressure_stdev,
            y_stdev=temperature_stdev,
            truncate=truncate,
        )
        interp_smoothed_S, interp_smoothed_dSdP, interp_smoothed_dSdT = S_interps

        V_interps = interp_smoothed_array_and_derivatives(
            array=grid_volumes,
            x_values=grid_pressures,
            y_values=grid_temperatures,
            x_stdev=pressure_stdev,
            y_stdev=temperature_stdev,
            truncate=truncate,
        )

        interp_smoothed_V, interp_smoothed_dVdP, interp_smoothed_dVdT = V_interps

        # Now we can calculate and plot the relaxed and smoothed properties along the isentrope
        smoothed_temperatures = interp_isentrope(
            interp_smoothed_S, pressures, entropy, potential_temperature
        )
        densities = rock.evaluate(["rho"], pressures, smoothed_temperatures)[0]
        depths, gravity = compute_depth_gravity_profiles(
            pressures, densities, surface_gravity, outer_radius
        )

        dT = 0.1
        Vpsub, Vssub = rock.evaluate(
            ["p_wave_velocity", "shear_wave_velocity"],
            pressures,
            smoothed_temperatures - dT / 2.0,
        )
        Vpadd, Vsadd = rock.evaluate(
            ["p_wave_velocity", "shear_wave_velocity"],
            pressures,
            smoothed_temperatures + dT / 2.0,
        )

        Vps = (Vpadd + Vpsub) / 2.0
        Vss = (Vsadd + Vssub) / 2.0
        dVpdT = (Vpadd - Vpsub) / dT
        dVsdT = (Vsadd - Vssub) / dT

        p_smoothT = np.array([pressures, smoothed_temperatures]).T
        volumes = interp_smoothed_V(p_smoothT)
        dSdT = interp_smoothed_dSdT(p_smoothT)
        dVdT = interp_smoothed_dVdT(p_smoothT)
        dVdP = interp_smoothed_dVdP(p_smoothT)

        specific_heats_relaxed = (
            smoothed_temperatures * dSdT / rock.params["molar_mass"]
        )
        alphas_relaxed = dVdT / volumes
        compressibilities_relaxed = -dVdP / volumes

        print(
            "Min and max relaxed property when pressure smoothing standard deviation is {0:.2f} GPa".format(
                pressure_stdev / 1.0e9
            )
        )
        print(
            "Specific heat: {0:.2e}, {1:.2e}".format(
                np.min(specific_heats_relaxed), np.max(specific_heats_relaxed)
            )
        )
        print(
            "Thermal expansivity: {0:.2e}, {1:.2e}".format(
                np.min(alphas_relaxed), np.max(alphas_relaxed)
            )
        )
        print(
            "Compressibilities: {0:.2e}, {1:.2e}\n".format(
                np.min(compressibilities_relaxed), np.max(compressibilities_relaxed)
            )
        )

        ax_T.plot(
            x,
            smoothed_temperatures,
            label="relaxed, smoothed (P_sd: {0:.1f} GPa)".format(
                pressure_stdev / 1.0e9
            ),
        )
        ax_z.plot(x, depths / 1.0e3)
        ax_g.plot(x, gravity)
        ax_rho.plot(x, densities)
        ax_cp.plot(x, specific_heats_relaxed)
        ax_alpha.plot(x, alphas_relaxed)
        ax_beta.plot(x, compressibilities_relaxed)

    ax_T.legend(loc="upper left")
    fig.set_layout_engine("tight")

    plt.show()
