from __future__ import print_function

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.optimize import fsolve
from burnman import constants
from burnman.utils import geotherm
import warnings

from .material import Material, material_property


class Layer(object):
    """
    The base class for a planetary layer.
    The user needs to set the following before properties can be computed:

    - set_material(), which sets the material of the layer,
      e.g. a mineral, solid_solution, or composite
    - set_temperature_mode(), either predefine, or set to an adiabatic profile
    - set_pressure_mode(), to set the self-consistent pressure
      (with user-defined option the pressures can be overwritten).
      To set the self-consistent pressure the pressure at the top and the
      gravity at the bottom of the layer need to be set.
    - make(), computes the self-consistent part of the layer and starts the
      settings to compute properties within the layer

    Note that the entire planet this layer sits in is not necessarily
    self-consistent, as the pressure at the top of the layer is a
    function of the density within the layer (through the gravity).
    Entire planets can be computed self-consistently with the planet class.
    Properties will be returned at the pre-defined radius array,
    although the evaluate() function can take a newly defined depthlist
    and values are interpolated between these (sufficient sampling of the layer
    is needed for this to be accurate).
    """

    def __init__(self, name=None, radii=None, verbose=False):
        self.name = name
        assert np.all(np.diff(radii) > 0)
        self.radii = radii
        self.outer_radius = max(self.radii)
        self.inner_radius = min(self.radii)
        self.thickness = self.outer_radius - self.inner_radius
        self.n_slices = len(self.radii)
        self.verbose = verbose
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None
        self.material = None
        self.pressure_mode = "self-consistent"
        self.temperature_mode = None

    def __str__(self):
        """
        Prints details of the layer
        """
        writing = (
            f"The {self.name} is made of {self.material.name}"
            f" with {self.temperature_mode} temperatures and "
            f"{self.pressure_mode} pressures\n"
        )
        return writing

    def reset(self):
        """
        Resets all cached material properties.
        It is typically not required for the user to call this function.
        """
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None

    def set_material(self, material):
        """
        Set the material of a Layer with a Material
        """
        assert isinstance(material, Material)
        self.material = material
        self.reset()

    def set_temperature_mode(
        self, temperature_mode="adiabatic", temperatures=None, temperature_top=None
    ):
        """
        Sets temperatures within the layer as user-defined values or as
        a (potentially perturbed) adiabat.

        :param temperature_mode: This can be set to 'user-defined', 'adiabatic',
            or 'perturbed-adiabatic'. 'user-defined' fixes the temperature
            with the profile input by the user. 'adiabatic' self-consistently
            computes the adiabat when setting the state of the layer.
            'perturbed-adiabatic' adds the user input array to the adiabat.
            This allows the user to apply boundary layers (for example).
        :type temperature_mode: string
        :param temperatures: The desired fixed temperatures in [K].
            Should have same length as defined radii in layer.
        :type temperatures: array of float
        :param temperature_top: Temperature at the top of the layer.
            Used if the temperature mode is chosen to be 'adiabatic' or
            'perturbed-adiabatic'. If 'perturbed-adiabatic' is chosen as the
            temperature mode, temperature_top corresponds to the true temperature
            at the top of the layer, and the reference isentrope at this radius
            is defined to lie at a temperature of
            temperature_top - temperatures[-1].
        :type temperature_top: float
        """
        self.reset()
        assert (
            temperature_mode == "user-defined"
            or temperature_mode == "adiabatic"
            or temperature_mode == "perturbed-adiabatic"
        )

        self.temperature_mode = temperature_mode

        if (
            temperature_mode == "user-defined"
            or temperature_mode == "perturbed-adiabatic"
        ):
            assert len(temperatures) == len(self.radii)
            self.usertemperatures = temperatures
        else:
            self.usertemperatures = np.zeros_like(self.radii)

        if temperature_mode == "adiabatic" or temperature_mode == "perturbed-adiabatic":
            self.temperature_top = temperature_top
        else:
            self.temperature_top = None

    def set_pressure_mode(
        self,
        pressure_mode="self-consistent",
        pressures=None,
        gravity_bottom=None,
        pressure_top=None,
        n_max_iterations=50,
        max_delta=1.0e-5,
    ):
        """
        Sets the pressure mode of the layer,
        which can either be 'user-defined', or 'self-consistent'.

        :param pressure_mode: This can be set to 'user-defined' or 'self-consistent'.
            'user-defined' fixes the pressures with the profile input
            by the user in the 'pressures' argument.
            'self-consistent' forces Layer to calculate pressures
            self-consistently. If this is selected, the user will need
            to supply values for the gravity_bottom [m/s^2]
            and pressure_top [Pa] arguments.
        :type pressure_mode: string
        :param pressures: Pressures [Pa] to set layer to
            (if the 'user-defined' pressure_mode has been selected).
            The array should be the same length as
            the layers user-defined radii array.
        :type pressures: array of floats
        :param pressure_top: Pressure [Pa] at the top of the layer.
        :type pressure_top: float
        :param gravity_bottom: Gravity [m/s^2] at the bottom of the layer.
        :type gravity_bottom: float
        :param n_max_iterations: Maximum number of iterations to reach
            self-consistent pressures.
        :type n_max_iterations: integer
        :param max_delta: Relative update to the highest pressure in the layer between
            iterations to stop iterations.
        :type max_delta: float
        """
        self.reset()
        assert pressure_mode == "user-defined" or pressure_mode == "self-consistent"
        self.pressure_mode = pressure_mode

        assert gravity_bottom is not None
        self.gravity_bottom = gravity_bottom

        if pressure_mode == "user-defined":
            assert pressures is not None
            assert len(pressures) == len(self.radii)
            self.pressures = pressures
            warnings.warn(
                "By setting the pressures in Layer they "
                "are unlikely to be self-consistent"
            )
        elif pressure_mode == "self-consistent":
            self.pressure_top = pressure_top
            self.n_max_iterations = n_max_iterations
            self.max_delta = max_delta
        else:
            raise NotImplementedError(
                f"pressure mode {pressure_mode} " "not recognised"
            )

    def make(self):
        """
        This routine needs to be called before evaluating any properties.
        If pressures and temperatures are not user-defined, they
        are computed here. This method also initializes an array of copied
        materials from which properties can be computed.
        """
        self.reset()
        if not hasattr(self, "material"):
            raise AttributeError(
                "You must set_material() for the layer " "before running make()."
            )
        if not hasattr(self, "temperature_mode"):
            raise AttributeError(
                "You must set_temperature_mode() for the "
                "layer before running make()."
            )
        if not hasattr(self, "pressure_mode"):
            raise AttributeError(
                "You must set_pressure_mode() for the layer " "before running make()."
            )

        if self.pressure_mode == "user-defined":
            self.temperatures = self._evaluate_temperature(
                self.pressures, self.temperature_top
            )
        elif self.pressure_mode == "self-consistent":
            new_press = (
                self.pressure_top + (-self.radii + max(self.radii)) * 1.0e3
            )  # initial pressure curve guess
            temperatures = self._evaluate_temperature(new_press, self.temperature_top)
            # Make it self-consistent!!!
            i = 0

            while i < self.n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(
                    new_press, temperatures, self.gravity_bottom, self.pressure_top
                )
                temperatures = self._evaluate_temperature(
                    new_press, self.temperature_top
                )
                rel_err = abs((max(ref_press) - max(new_press)) / max(new_press))
                if self.verbose:
                    print(
                        f"Iteration {i:0d} maximum relative pressure error: "
                        f"{rel_err:.1f}"
                    )

                if rel_err < self.max_delta:
                    break

            self.pressures = new_press
            self.temperatures = temperatures
        else:
            raise NotImplementedError("pressure mode not recognised")

        self.sublayers = []
        for r in range(len(self.radii)):
            self.sublayers.append(self.material.copy())
            self.sublayers[r].set_state(self.pressures[r], self.temperatures[r])

    def evaluate(self, properties, radlist=None, radius_planet=None):
        """
        Function that is used to evaluate properties
        across the layer. If radlist is not defined, values are
        returned at the internal radlist.
        If asking for different radii than the internal radlist,
        pressure and temperature values are interpolated and the
        layer material evaluated at those pressures and
        temperatures.

        :param properties: List of properties to evaluate.
        :type properties: list of strings
        :param radlist: Radii to evaluate properties at. If left empty,
            internal radii list is used.
        :type radlist: array of floats
        :param planet_radius: Planet outer radius. Used only to calculate depth.
        :type planet_radius: float

        :returns: 1D or 2D array of requested properties
            (1D if only one property was requested)
        :rtype: numpy.array
        """

        if radlist is None:
            values = np.empty([len(properties), len(self.radii)])
            for i, prop in enumerate(properties):
                if prop == "depth":
                    values[i] = radius_planet - self.radii
                else:
                    try:
                        values[i] = getattr(self, prop)
                    except:
                        values[i] = np.array(
                            [
                                getattr(self.sublayers[i], prop)
                                for i in range(len(self.sublayers))
                            ]
                        )
        else:
            func_p = interp1d(self.radii, self.pressures)
            pressures = func_p(radlist)
            func_t = interp1d(self.radii, self.temperatures)
            temperatures = func_t(radlist)
            values = np.empty([len(properties), len(radlist)])
            for i, prop in enumerate(properties):
                if prop == "depth":
                    values[i] = radius_planet - radlist
                else:
                    try:
                        values[i] = self.material.evaluate(
                            [prop], pressures, temperatures
                        )
                    except:
                        func_prop = interp1d(self.radii, getattr(self, prop))
                        values[i] = func_prop(radlist)

        if values.shape[0] == 1:
            values = values[0]
        return values

    def _evaluate_temperature(self, pressures=None, temperature_top=None):
        """
        Returns the temperatures of the layer for given pressures.
        Used by make()
        """
        if (
            self.temperature_mode == "adiabatic"
            or self.temperature_mode == "perturbed-adiabatic"
        ):
            adiabat = geotherm.adiabatic(
                pressures[::-1],
                temperature_top - self.usertemperatures[-1],
                self.material,
            )[::-1]
        else:
            adiabat = np.zeros_like(self.radii)
        return adiabat + self.usertemperatures

    def _evaluate_eos(self, pressures, temperatures, gravity_bottom, pressure_top):
        """
        Returns updated gravity and pressure
        make() loops over this until consistency is achieved.
        """
        [density] = self.material.evaluate(["density"], pressures, temperatures)
        grav = self._compute_gravity(density, gravity_bottom)
        press = self._compute_pressure(density, grav, pressure_top)
        return grav, press

    # Functions needed to compute self-consistent radii-pressures
    def _compute_gravity(self, density, gravity_bottom):
        """
        Computes the gravity of a layer
        Used by _evaluate_eos()
        """
        # Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(self.radii, density)
        # Numerically integrate Poisson's equation

        def poisson(p, x):
            return 4.0 * np.pi * constants.G * rhofunc(x) * x * x

        grav = np.ravel(
            odeint(poisson, gravity_bottom * self.radii[0] * self.radii[0], self.radii)
        )

        if self.radii[0] == 0:
            grav[0] = 0
            grav[1:] = grav[1:] / self.radii[1:] / self.radii[1:]
        else:
            grav[:] = grav[:] / self.radii[:] / self.radii[:]
        return grav

    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.
        This integrates the equation for hydrostatic equilibrium P = rho g z.
        Used by _evaluate_eos()
        """
        # flip radius, density and gravity to increasing pressure
        depthfromtop = -self.radii[::-1] + max(self.radii)
        density = density[::-1]
        gravity = gravity[::-1]
        # Make a spline fit of density as a function of depth
        rhofunc = UnivariateSpline(depthfromtop, density)
        # Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline(depthfromtop, gravity)

        # integrate the hydrostatic equation
        pressure = np.ravel(
            odeint((lambda p, x: gfunc(x) * rhofunc(x)), pressure_top, depthfromtop)
        )

        return pressure[::-1]

    @property
    def mass(self):
        """
        Calculates the mass of the layer [kg]
        """
        mass = 0.0
        radii = self.radii
        density = self.evaluate(["density"])
        rhofunc = UnivariateSpline(radii, density)
        mass = np.abs(
            quad(lambda r: 4 * np.pi * rhofunc(r) * r * r, radii[0], radii[-1])[0]
        )
        return mass

    @property
    def moment_of_inertia(self):
        """
        Returns the moment of inertia of the layer [kg m^2]
        """
        moment = 0.0
        radii = self.radii
        density = self.evaluate(["density"])
        rhofunc = UnivariateSpline(radii, density)
        moment = np.abs(
            quad(
                lambda r: 8.0 / 3.0 * np.pi * rhofunc(r) * r * r * r * r,
                radii[0],
                radii[-1],
            )[0]
        )
        return moment

    @property
    def gravity(self):
        """
        Returns gravity profile of the layer [m s^(-2)]
        """
        return self._compute_gravity(self.density, self.gravity_bottom)

    @property
    def bullen(self):
        """
        Returns the Bullen parameter across the layer.
        The Bullen parameter assess if compression as a function of pressure is
        like homogeneous, adiabatic compression.
        Bullen parameter =1  , homogeneous, adiabatic compression
        Bullen parameter > 1 , more compressed with pressure,
        e.g. across phase transitions Bullen parameter < 1,
        less compressed with pressure, e.g. across a boundary layer.
        """
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        phi = self.bulk_sound_velocity * self.bulk_sound_velocity
        try:
            dkappadP = np.gradient(kappa, edge_order=2) / np.gradient(
                self.pressures, edge_order=2
            )
            dphidr = (
                np.gradient(phi, edge_order=2)
                / np.gradient(self.radii, edge_order=2)
                / self.gravity
            )
        except:
            dkappadP = np.gradient(kappa) / np.gradient(self.pressures)
            dphidr = np.gradient(phi) / np.gradient(self.radii) / self.gravity
        bullen = dkappadP + dphidr
        return bullen

    @property
    def brunt_vasala(self):
        """
        Returns the brunt-vasala (or buoyancy) frequency, N, across the layer.
        This frequency assess the stabilty of the layer:
        N < 0, fluid will convect
        N= 0, fluid is neutral
        N > 0, fluid is stabily stratified.
        """
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        brunt_vasala = (
            self.density * self.gravity * self.gravity * (self.bullen - 1.0) / kappa
        )
        return brunt_vasala

    @property
    def pressure(self):
        """
        Returns current pressures across the layer that was set
        with :func:`~burnman.Material.set_state`.

        Aliased with :func:`~burnman.Material.P`.

        :returns: Pressures in [Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return self.pressures

    @property
    def temperature(self):
        """
        Returns current temperature  across the layer that was set with
        :func:`~burnman.Material.set_state`.

        - Aliased with :func:`~burnman.Material.T`.

        :returns: Temperatures in [K] at the predefined radii.
        :rtype: numpy.array
        """
        return self.temperatures

    @material_property
    def molar_internal_energy(self):
        """
        Returns the molar internal energies across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.Material.energy`.

        :returns: The internal energies in [J/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].molar_internal_energy
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def molar_gibbs(self):
        """
        Returns the molar Gibbs free energies across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.gibbs`.

        :returns: Gibbs energies in [J/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_gibbs for i in range(len(self.sublayers))]
        )

    @material_property
    def molar_helmholtz(self):
        """
        Returns the molar Helmholtz free energies across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.helmholtz`.

        :returns: Helmholtz energies in [J/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_helmholtz for i in range(len(self.sublayers))]
        )

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the layer.

        Needs to be implemented in derived classes.

        :returns: Molar mass in [kg/mol].
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_mass for i in range(len(self.sublayers))]
        )

    @material_property
    def molar_volume(self):
        """
        Returns molar volumes across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.V`.

        :returns: Molar volumes in [m^3/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_volume for i in range(len(self.sublayers))]
        )

    @material_property
    def density(self):
        """
        Returns the densities across this layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.rho`.

        :returns: The densities of this material in [kg/m^3] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array([self.sublayers[i].density for i in range(len(self.sublayers))])

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropies acroos the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.S`.

        :returns: Entropies in [J/K/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_entropy for i in range(len(self.sublayers))]
        )

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpies across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.H`.

        :returns: Enthalpies in [J/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].molar_enthalpy for i in range(len(self.sublayers))]
        )

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        Returns isothermal bulk moduli across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.Material.K_T`.

        :returns: Bulk moduli in [Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].isothermal_bulk_modulus_reuss
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Returns the adiabatic bulk moduli across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.K_S`.

        :returns: Adiabatic bulk modulus in [Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].isentropic_bulk_modulus_reuss
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns isothermal compressibilities across the layer
        (or inverse isothermal bulk moduli).

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_T`.

        :returns: Isothermal compressibilities in [1/Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].isothermal_compressibility_reuss
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        Returns adiabatic compressibilities across the layer
        (or inverse adiabatic bulk moduli).

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_S`.

        :returns: Adiabatic compressibilities in [1/Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].isentropic_compressibility_reuss
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def shear_modulus(self):
        """
        Returns shear moduli across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_G`.

        :returns: Shear moduli in [Pa] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].shear_modulus for i in range(len(self.sublayers))]
        )

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speeds across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_p`.

        :returns: P wave speeds in [m/s] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].p_wave_velocity for i in range(len(self.sublayers))]
        )

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speeds across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_phi`.

        :returns: Bulk sound velocities in [m/s] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].bulk_sound_velocity for i in range(len(self.sublayers))]
        )

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speeds across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_s`.

        :returns: Shear wave speeds in [m/s] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].shear_wave_velocity for i in range(len(self.sublayers))]
        )

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameters across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.gr`.

        :returns: Grueneisen parameters [unitless] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].grueneisen_parameter for i in range(len(self.sublayers))]
        )

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficients across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.alpha`.

        :returns: Thermal expansivities in [1/K] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [self.sublayers[i].thermal_expansivity for i in range(len(self.sublayers))]
        )

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volumes across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.C_v`.

        :returns: Heat capacities in [J/K/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].molar_heat_capacity_v
                for i in range(len(self.sublayers))
            ]
        )

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar_heat capacity at constant pressures across the layer.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.C_p`.

        :returns: Heat capacities in [J/K/mol] at the predefined radii.
        :rtype: numpy.array
        """
        return np.array(
            [
                self.sublayers[i].molar_heat_capacity_p
                for i in range(len(self.sublayers))
            ]
        )

    # Aliased properties
    @property
    def P(self):
        """Alias for :func:`~burnman.Layer.pressure`"""
        return self.pressure

    @property
    def T(self):
        """Alias for :func:`~burnman.Layer.temperature`"""
        return self.temperature

    @property
    def energy(self):
        """Alias for :func:`~burnman.Layer.molar_internal_energy`"""
        return self.molar_internal_energy

    @property
    def helmholtz(self):
        """Alias for :func:`~burnman.Layer.molar_helmholtz`"""
        return self.molar_helmholtz

    @property
    def gibbs(self):
        """Alias for :func:`~burnman.Layer.molar_gibbs`"""
        return self.molar_gibbs

    @property
    def V(self):
        """Alias for :func:`~burnman.Layer.molar_volume`"""
        return self.molar_volume

    @property
    def rho(self):
        """Alias for :func:`~burnman.Layer.density`"""
        return self.density

    @property
    def S(self):
        """Alias for :func:`~burnman.Layer.molar_entropy`"""
        return self.molar_entropy

    @property
    def H(self):
        """Alias for :func:`~burnman.Layer.molar_enthalpy`"""
        return self.molar_enthalpy

    @property
    def K_T(self):
        """Alias for :func:`~burnman.Layer.isothermal_bulk_modulus_reuss`"""
        return self.isothermal_bulk_modulus_reuss

    @property
    def K_S(self):
        """Alias for :func:`~burnman.Layer.isentropic_bulk_modulus_reuss`"""
        return self.isentropic_bulk_modulus_reuss

    @property
    def beta_T(self):
        """Alias for :func:`~burnman.Layer.isothermal_compressibility_reuss`"""
        return self.isothermal_compressibility_reuss

    @property
    def beta_S(self):
        """Alias for :func:`~burnman.Layer.isentropic_compressibility_reuss`"""
        return self.isentropic_compressibility_reuss

    @property
    def G(self):
        """Alias for :func:`~burnman.Layer.shear_modulus`"""
        return self.shear_modulus

    @property
    def v_p(self):
        """Alias for :func:`~burnman.Layer.p_wave_velocity`"""
        return self.p_wave_velocity

    @property
    def v_phi(self):
        """Alias for :func:`~burnman.Layer.bulk_sound_velocity`"""
        return self.bulk_sound_velocity

    @property
    def v_s(self):
        """Alias for :func:`~burnman.Layer.shear_wave_velocity`"""
        return self.shear_wave_velocity

    @property
    def gr(self):
        """Alias for :func:`~burnman.Layer.grueneisen_parameter`"""
        return self.grueneisen_parameter

    @property
    def alpha(self):
        """Alias for :func:`~burnman.Layer.thermal_expansivity`"""
        return self.thermal_expansivity

    @property
    def C_v(self):
        """Alias for :func:`~burnman.Material.molar_heat_capacity_v`"""
        return self.molar_heat_capacity_v

    @property
    def C_p(self):
        """Alias for :func:`~burnman.Material.molar_heat_capacity_p`"""
        return self.molar_heat_capacity_p


class BoundaryLayerPerturbation(object):
    """
    A class that implements a temperature perturbation model corresponding to a
    simple thermal boundary layer.
    The model takes the following form:
    T = a*exp((r - r1)/(r0 - r1)*c) + b*exp((r - r0)/(r1 - r0)*c)
    The relationships between the input parameters and a, b and c are
    given below.

    This model is a simpler version of that proposed
    by :cite:`Richter1981`.

    :param radius_bottom: The radius at the bottom of the layer (r0) [m].
    :type radius_bottom: float

    :param radius_top: The radius at the top of the layer (r1) [m].
    :type radius_top: float

    :param rayleigh_number: The Rayleigh number of convection within the layer. The
        exponential scale factor is this number to the power of 1/4
        (Ra = c^4).
    :type rayleigh_number: float

    :param temperature_change: The total difference in potential
        temperature across the layer [K]. temperature_change = (a + b)*exp(c).
    :type temperature_change: float

    :param boundary_layer_ratio: The ratio of the linear scale factors (a/b)
        corresponding to the thermal boundary layers at the top and bottom of
        the layer. A number greater than 1 implies a larger change in
        temperature across the top boundary than the bottom boundary.
    :type boundary_layer_ratio: float
    """

    def __init__(
        self,
        radius_bottom,
        radius_top,
        rayleigh_number,
        temperature_change,
        boundary_layer_ratio,
    ):
        self.r0 = radius_bottom
        self.r1 = radius_top

        self.Ra = rayleigh_number
        self.c = np.power(self.Ra, 1.0 / 4.0)

        self.a = temperature_change / (np.exp(self.c) * (1.0 + boundary_layer_ratio))
        self.b = -boundary_layer_ratio * self.a

    def temperature(self, radii):
        """
        Returns the temperature at one or more radii [K].

        :param radii: The radii at which to evaluate the temperature.
        :type radii: float or numpy.array

        :returns: The temperatures at the requested radii.
        :rtype: float or numpy.array
        """
        return self.a * np.exp(
            (radii - self.r1) / (self.r0 - self.r1) * self.c
        ) + self.b * np.exp((radii - self.r0) / (self.r1 - self.r0) * self.c)

    def dTdr(self, radii):
        """
        Returns the thermal gradient at one or more radii [K/m].

        :param radii: The radii at which to evaluate the thermal gradients.
        :type radii: float or numpy.array

        :returns: The thermal gradient at the requested radii.
        :rtype: float or numpy.array
        """
        return (
            self.c
            / (self.r0 - self.r1)
            * (
                self.a * np.exp((radii - self.r1) / (self.r0 - self.r1) * self.c)
                - self.b * np.exp((radii - self.r0) / (self.r1 - self.r0) * self.c)
            )
        )

    def set_model_thermal_gradients(self, dTdr_bottom, dTdr_top):
        """
        Reparameterizes the model based on the thermal gradients
        at the bottom and top of the model.

        :param dTdr_bottom: The thermal gradient at the bottom of the model [K/m].
            Typically negative for a cooling planet.
        :type dTdr_bottom: float

        :param dTdr_top: The thermal gradient at the top of the model [K/m].
            Typically negative for a cooling planet.
        :type dTdr_top: float
        """

        def delta_dTdrs(args, dTdr_bottom, dTdr_top):
            a, b = args
            self.a = a
            self.b = b

            return [dTdr_bottom - self.dTdr(self.r0), dTdr_top - self.dTdr(self.r1)]

        fsolve(delta_dTdrs, [self.a, self.b], args=(dTdr_bottom, dTdr_top))
