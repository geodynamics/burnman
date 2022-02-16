# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

from subprocess import Popen, PIPE, STDOUT
from os import rename

import numpy as np
from scipy.interpolate import interp2d, griddata

from .material import Material, material_property

from ..utils.misc import copy_documentation


def create_perplex_table(werami_path, project_name, outfile, n_pressures, n_temperatures, pressure_range=None, temperature_range=None):
    '''
    This function uses PerpleX's werami software to output a table file containing the following material properties.
    2 - Density (kg/m3)
    4 - Expansivity (1/K, for volume)
    5 - Compressibility (1/bar, for volume)
    10 - Adiabatic bulk modulus (bar)
    11 - Adiabatic shear modulus (bar)
    12 - Sound velocity (km/s)
    13 - P-wave velocity (Vp, km/s)
    14 - S-wave velocity (Vs, km/s)
    17 - Entropy (J/K/kg)
    18 - Enthalpy (J/kg)
    19 - Heat Capacity (J/K/kg)
    22 - Molar Volume (J/bar)
    '''

    print('Working on creating {0}x{1} P-T table file using werami. Please wait.\n'.format(
        n_pressures, n_temperatures))

    try:
        str2 = 'y\n{0} {1}\n{2} {3}\n'.format(pressure_range[0]/1.e5, pressure_range[1]/1.e5,
                                              temperature_range[0], temperature_range[1])
    except:
        print('Keeping P-T range the same as the original project range.\n')
        str2 = 'n\n'

    stdin = '{0:s}\n2\n' \
        '2\nn\n' \
        '4\nn\n' \
        '5\nn\n' \
        '10\nn\n' \
        '11\nn\n' \
        '12\nn\n' \
        '13\nn\n' \
        '14\nn\n' \
        '17\nn\n' \
        '18\nn\n' \
        '19\nn\n' \
        '22\nn\n' \
        '0\n' \
        '{1:s}' \
        '{2:d} {3:d}\n' \
        '0'.format(project_name, str2, n_pressures, n_temperatures)

    p = Popen(werami_path, stdout=PIPE, stdin=PIPE,
              stderr=STDOUT, encoding='utf8')
    stdout = p.communicate(input=stdin)[0]
    print(stdout)
    out = [s for s in stdout.split(
        '\n') if "Output has been written to the" in s][0].split()[-1]
    rename(out, outfile)
    print('Output file renamed to {0:s}'.format(outfile))
    print('Processing complete')


class PerplexMaterial(Material):
    """
    This is the base class for a PerpleX material. States of the material
    can only be queried after setting the pressure and temperature
    using set_state().

    Instances of this class are initialised with
    a 2D PerpleX tab file. This file should be in the standard format
    (as output by werami), and should have columns with the following names:
    'rho,kg/m3', 'alpha,1/K', 'beta,1/bar', 'Ks,bar', 'Gs,bar', 'v0,km/s',
    'vp,km/s', 'vs,km/s', 's,J/K/kg', 'h,J/kg', 'cp,J/K/kg', 'V,J/bar/mol'.
    The order of these names is not important.

    Properties of the material are determined by linear interpolation from
    the PerpleX grid. They are all returned in SI units on a molar basis,
    even though the PerpleX tab file is not in these units.

    This class is available as ``burnman.PerplexMaterial``.
    """

    def __init__(self, tab_file, name='Perple_X material'):
        self.name = name
        self.params = {'name': name}
        self._property_interpolators, self.params['molar_mass'], self.bounds = self._read_2D_perplex_file(
            tab_file)
        Material.__init__(self)

    def _read_2D_perplex_file(self, filename):
        with open(filename, 'r') as f:
            datastream = f.read()

        lines = [line.strip().split()
                 for line in datastream.split('\n') if line.strip()]

        if lines[2][0] != '2':
            raise Exception('This is not a 2D PerpleX table')

        bounds = [(float(lines[4][0]), float(lines[5][0]), int(lines[6][0])),
                  (float(lines[8][0]), float(lines[9][0]), int(lines[10][0]))]

        if lines[3][0] == 'P(bar)' and lines[7][0] == 'T(K)':
            Pmin, Pint, nP = bounds[0]
            Tmin, Tint, nT = bounds[1]
        elif lines[3][0] == 'T(K)' and lines[7][0] == 'P(bar)':
            Pmin, Pint, nP = bounds[1]
            Tmin, Tint, nT = bounds[0]
        else:
            raise Exception('This file does not have a recognised PerpleX structure.\n'
                            'Are the independent variables P(bar) and T(K)?')

        Pmin = Pmin*1.e5  # bar to Pa
        Pint = Pint*1.e5  # bar to Pa
        Pmax = Pmin + Pint*(nP-1.)
        Tmax = Tmin + Tint*(nT-1.)
        pressures = np.linspace(Pmin, Pmax, nP)
        temperatures = np.linspace(Tmin, Tmax, nT)
        n_properties = int(lines[11][0])
        property_list = lines[12]

        # property_table[i][j][k] returns the kth property at the ith pressure and jth temperature
        table = np.array([[float(string) for string in line]
                         for line in lines[13:13+nP*nT]])

        if lines[3][0] == 'P(bar)':
            property_table = np.swapaxes(
                table.reshape(nT, nP, n_properties), 0, 1)
        else:
            property_table = table.reshape(nP, nT, n_properties)

        ordered_property_list = ['rho,kg/m3',
                                 'alpha,1/K',
                                 'beta,1/bar',
                                 'Ks,bar',
                                 'Gs,bar',
                                 'v0,km/s',
                                 'vp,km/s',
                                 'vs,km/s',
                                 's,J/K/kg',
                                 'h,J/kg',
                                 'cp,J/K/kg',
                                 'V,J/bar/mol']
        p_indices = [i for i, p in enumerate(
            property_list) for ordered_p in ordered_property_list if p == ordered_p]

        properties = {}
        for i, p_idx in enumerate(p_indices):
            # Fill in NaNs as long as they aren't in the corners of the P-T grid
            a = np.array(property_table[:, :, [p_idx]][:, :, 0])
            x, y = np.indices(a.shape)
            a[np.isnan(a)] = griddata((x[~np.isnan(a)], y[~np.isnan(a)]),  # points we know
                                      a[~np.isnan(a)],  # values we know
                                      (x[np.isnan(a)], y[np.isnan(a)]))

            # Fill any remaining NaNs with zeros
            properties[ordered_property_list[i]] = np.nan_to_num(a, 0.)

        densities = properties['rho,kg/m3']
        volumes = 1.e-5 * properties['V,J/bar/mol']
        molar_masses = densities*volumes
        molar_mass = np.mean(molar_masses)
        property_interpolators = {'rho': interp2d(pressures, temperatures, densities.T, bounds_error=True),
                                  'alpha': interp2d(pressures, temperatures, properties['alpha,1/K'].T, bounds_error=True),
                                  'K_T': interp2d(pressures, temperatures, 1.e5 / properties['beta,1/bar'].T, bounds_error=True),
                                  'K_S': interp2d(pressures, temperatures, 1.e5 * properties['Ks,bar'].T, bounds_error=True),
                                  'G_S': interp2d(pressures, temperatures, 1.e5 * properties['Gs,bar'].T, bounds_error=True),
                                  'bulk_sound_velocity': interp2d(pressures, temperatures, 1.e3*properties['v0,km/s'].T, bounds_error=True),
                                  'p_wave_velocity': interp2d(pressures, temperatures, 1.e3*properties['vp,km/s'].T, bounds_error=True),
                                  's_wave_velocity': interp2d(pressures, temperatures, 1.e3*properties['vs,km/s'].T, bounds_error=True),
                                  'S': interp2d(pressures, temperatures, properties['s,J/K/kg'].T*molar_masses.T, bounds_error=True),
                                  'H': interp2d(pressures, temperatures, properties['h,J/kg'].T*molar_masses.T, bounds_error=True),
                                  'C_p': interp2d(pressures, temperatures, properties['cp,J/K/kg'].T*molar_masses.T, bounds_error=True),
                                  'V': interp2d(pressures, temperatures, volumes.T, bounds_error=True)}

        bounds = [[Pmin, Pmax], [Tmin, Tmax]]
        return property_interpolators, molar_mass, bounds

    @copy_documentation(Material.set_state)
    def set_state(self, pressure, temperature):

        if not np.logical_and(np.all(self.bounds[0][0] <= pressure),
                              np.all(pressure <= self.bounds[0][1])):
            raise ValueError("The set_state pressure ({0:.4f}) is outside the bounds of this rock ({1:.4f}-{2:.4f} GPa)".format(pressure,
                                                                                                                                self.bounds[
                                                                                                                                    0][0]/1.e9,
                                                                                                                                self.bounds[0][1]/1.e9))
        if not np.logical_and(np.all(self.bounds[1][0] <= temperature),
                              np.all(temperature <= self.bounds[1][1])):
            raise ValueError("The set_state temperature ({0:.1f}) is outside the bounds of this rock ({1:.1f}-{2:.1f} K)".format(temperature,
                                                                                                                                 self.bounds[
                                                                                                                                     1][0],
                                                                                                                                 self.bounds[1][1]))
        Material.set_state(self, pressure, temperature)

    """
    Properties by linear interpolation of Perple_X output
    """

    @material_property
    @copy_documentation(Material.molar_volume)
    def molar_volume(self):
        return self._property_interpolators['V'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.molar_enthalpy)
    def molar_enthalpy(self):
        return self._property_interpolators['H'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.molar_entropy)
    def molar_entropy(self):
        return self._property_interpolators['S'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.isothermal_bulk_modulus)
    def isothermal_bulk_modulus(self):
        return self._property_interpolators['K_T'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.adiabatic_bulk_modulus)
    def adiabatic_bulk_modulus(self):
        return self._property_interpolators['K_S'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.molar_heat_capacity_p)
    def molar_heat_capacity_p(self):
        return self._property_interpolators['C_p'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.thermal_expansivity)
    def thermal_expansivity(self):
        return self._property_interpolators['alpha'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.shear_modulus)
    def shear_modulus(self):
        return self._property_interpolators['G_S'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.p_wave_velocity)
    def p_wave_velocity(self):
        return self._property_interpolators['p_wave_velocity'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.bulk_sound_velocity)
    def bulk_sound_velocity(self):
        return self._property_interpolators['bulk_sound_velocity'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.shear_wave_velocity)
    def shear_wave_velocity(self):
        return self._property_interpolators['s_wave_velocity'](self.pressure, self.temperature)[0]

    """
    Properties from mineral parameters,
    Legendre transformations
    or Maxwell relations
    """

    @material_property
    @copy_documentation(Material.molar_gibbs)
    def molar_gibbs(self):
        return self.molar_enthalpy - self.temperature*self.molar_entropy

    @material_property
    @copy_documentation(Material.molar_mass)
    def molar_mass(self):
        if 'molar_mass' in self.params:
            return self.params['molar_mass']
        else:
            raise ValueError(
                "No molar_mass parameter for mineral " + self.to_string + ".")

    @material_property
    @copy_documentation(Material.density)
    def density(self):
        return self._property_interpolators['rho'](self.pressure, self.temperature)[0]

    @material_property
    @copy_documentation(Material.molar_internal_energy)
    def molar_internal_energy(self):
        return self.molar_gibbs - self.pressure * self.molar_volume + self.temperature * self.molar_entropy

    @material_property
    @copy_documentation(Material.molar_helmholtz)
    def molar_helmholtz(self):
        return self.molar_gibbs - self.pressure * self.molar_volume

    @material_property
    @copy_documentation(Material.isothermal_compressibility)
    def isothermal_compressibility(self):
        return 1. / self.isothermal_bulk_modulus

    @material_property
    @copy_documentation(Material.adiabatic_compressibility)
    def adiabatic_compressibility(self):
        return 1. / self.adiabatic_bulk_modulus

    @material_property
    @copy_documentation(Material.molar_heat_capacity_v)
    def molar_heat_capacity_v(self):
        return (self.molar_heat_capacity_p
                - self.molar_volume * self.temperature
                * self.thermal_expansivity * self.thermal_expansivity
                * self.isothermal_bulk_modulus)

    @material_property
    @copy_documentation(Material.grueneisen_parameter)
    def grueneisen_parameter(self):
        return (self.thermal_expansivity
                * self.molar_volume
                * self.adiabatic_bulk_modulus
                / self.molar_heat_capacity_p)
