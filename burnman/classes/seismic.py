# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import warnings
import scipy.integrate

from ..utils.misc import read_table
from .. import constants


class Seismic1DModel(object):
    """
    Base class for all the seismological models.
    """

    def __init__(self):
        pass

    def evaluate(self, vars_list, depth_list=None):
        """
        Returns the lists of data for a Seismic1DModel for the depths provided

        :param vars_list: Available variables depend on the seismic model,
            and can be chosen from 'pressure', 'density', 'gravity',
            'v_s', 'v_p', 'v_phi', 'G', 'K', 'QG' and 'QK'.
        :type vars_list: array of str
        :param depth_list: Array of depths [m] to evaluate seismic model at.
        :type depth_list: array of floats

        :returns: Array of values shapes as (len(vars_list),len(depth_list)).
        :rtype: numpy.array
        """
        if depth_list is None:
            depth_list = self.internal_depth_list()
        values = np.empty((len(vars_list), len(depth_list)))
        for a in range(len(vars_list)):
            values[a, :] = getattr(self, vars_list[a])(depth_list)
        return values

    def internal_depth_list(
        self, mindepth=0.0, maxdepth=1.0e99, discontinuity_interval=1.0
    ):
        """
        Returns a sorted list of depths at which this seismic data is specified.
        This allows you to compare the seismic data without interpolation.
        The depths can be bounded by the mindepth and maxdepth parameters.

        :param mindepth: Minimum depth value to be returned [m].
        :type mindepth: float
        :param maxdepth: Maximum depth value to be returned [m].
        :type maxdepth: float
        :param discontinuity_interval: Shift continuities to remove
            ambigious values for depth [m].
        :type discontinuity_interval: float

        :returns: Depths [m].
        :rtype: numpy.array
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def pressure(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Pressure(s) at given depth(s) in [Pa].
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def v_p(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: P wave velocity at given depth(s) in [m/s].
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def v_s(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: S wave velocity at given depth(s) in [m/s].
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def v_phi(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Bulk sound wave velocity at given depth(s) in [m/s].
        :rtype: float or numpy.array of floats
        """
        v_s = self.v_s(depth)
        v_p = self.v_p(depth)
        return np.sqrt(v_p * v_p - 4.0 / 3.0 * v_s * v_s)

    def density(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Density at given depth(s) in [kg/m^3].
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def G(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Shear modulus at given depth(s) in [Pa].
        :rtype: float or numpy.array of floats
        """
        return np.power(self.v_s(depth), 2.0) * self.density(depth)

    def K(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Bulk modulus at given depth(s) in [Pa]
        :rtype: float or numpy.array of floats
        """
        return np.power(self.v_phi(depth), 2.0) * self.density(depth)

    def QK(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Quality factor (dimensionless) for bulk modulus at given depth(s).
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def QG(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Quality factor (dimensionless) for shear modulus at given depth(s).
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def depth(self, pressure):
        """
        :param pressure: Pressure(s) [Pa] at which to evaluate depths.
        :type pressure: float or numpy.array of floats

        :returns: Depth(s) [m] for given pressure(s).
        :type depth: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")

    def gravity(self, depth):
        """
        :param depth: Depth(s) [m] at which to evaluate seismic model.
        :type depth: float or numpy.array of floats

        :returns: Gravity at given depths in [m/s^2].
        :rtype: float or numpy.array of floats
        """
        raise NotImplementedError("abstract method to be implemented in derived class")


class SeismicTable(Seismic1DModel):
    """
    This is a base class that gets a 1D seismic model from a table indexed and
    sorted by radius. Fill the tables in the constructor after deriving
    from this class. This class uses :class:`burnman.seismic.Seismic1DModel`

    Note: all tables need to be sorted by increasing depth.
    self.table_depth needs to be defined. Alternatively, you can also overwrite
    the _lookup function if you want to access with something else.
    """

    def __init__(self):
        Seismic1DModel.__init__(self)

        self.table_depth = []
        self.table_radius = []
        self.table_pressure = []
        self.table_gravity = []
        self.table_density = []
        self.table_vp = []
        self.table_vs = []
        self.table_QG = []
        self.table_QK = []

        self.earth_radius = 6371.0e3

    def internal_depth_list(
        self, mindepth=0.0, maxdepth=1.0e10, discontinuity_interval=1.0
    ):
        depths = np.array(
            [
                self.table_depth[x]
                for x in range(len(self.table_depth))
                if self.table_depth[x] >= mindepth and self.table_depth[x] <= maxdepth
            ]
        )
        discontinuities = np.where(depths[1:] - depths[:-1] == 0)[0]
        # Shift values at discontinities by 1 m to simplify evaluating values
        # around these.
        depths[discontinuities] = depths[discontinuities] - discontinuity_interval
        depths[discontinuities + 1] = (
            depths[discontinuities + 1] + discontinuity_interval
        )
        return depths

    def pressure(self, depth):
        if len(self.table_pressure) == 0:
            warnings.warn(
                "Pressure is not given in "
                + self.__class__.__name__
                + " and is now being computed. This will only work when density is defined for the entire planet. Use at your own risk."
            )
            self._compute_pressure()
        return self._lookup(depth, self.table_pressure)

    def gravity(self, depth):
        if len(self.table_gravity) == 0:
            warnings.warn(
                "Gravity is not given in "
                + self.__class__.__name__
                + " and is now being computed. This will only work when density is defined for the entire planet. Use at your own risk."
            )
            self._compute_gravity()
        return self._lookup(depth, self.table_gravity)

    def v_p(self, depth):
        return self._lookup(depth, self.table_vp)

    def v_s(self, depth):
        return self._lookup(depth, self.table_vs)

    def QK(self, depth):
        return self._lookup(depth, self.table_QK)

    def QG(self, depth):
        return self._lookup(depth, self.table_QG)

    def density(self, depth):
        if len(self.table_density) == 0:
            raise ValueError("Density has not been defined for this seismic model")
        return self._lookup(depth, self.table_density)

    def bullen(self, depth):
        """
        Returns the Bullen parameter only for significant arrays
        """
        assert len(depth) > 3
        v_phi = self.v_phi(depth)
        density = self.density(depth)
        phi = v_phi * v_phi
        kappa = phi * density
        try:
            dkappadP = np.gradient(kappa, edge_order=2) / np.gradient(
                self.pressure(depth), edge_order=2
            )
            dphidz = (
                np.gradient(phi, edge_order=2)
                / np.gradient(depth, edge_order=2)
                / self.gravity(depth)
            )
        except:
            dkappadP = np.gradient(kappa) / np.gradient(self.pressure(depth))
            dphidz = np.gradient(phi) / np.gradient(depth) / self.gravity(depth)
        bullen = dkappadP - dphidz
        return bullen

    def depth(self, pressure):
        if max(pressure) > max(self.table_pressure) or min(pressure) < min(
            self.table_pressure
        ):
            raise ValueError("Pressure outside range of SeismicTable")

        depth = np.interp(pressure, self.table_pressure, self.table_depth)
        return depth

    def radius(self, pressure):
        radius = np.interp(
            pressure,
            self.table_pressure[::-1],
            self.earth_radius - self.table_depth[::-1],
        )
        return radius

    def _lookup(self, depth, value_table):
        return np.interp(depth, self.table_depth, value_table)

    def _compute_gravity(self):
        # Calculate the gravity of the planet, based on a density profile.
        # There is a check of the surface value is within reason for this
        # model, otherwise values for PREM are used.

        density = self.table_density[::-1]

        if len(density) > 0:
            radii = self.table_radius[::-1]
            g = scipy.integrate.cumulative_trapezoid(
                constants.G * 4.0 * np.pi * density * radii * radii, x=radii, initial=0
            )
            g[1:] = g[1:] / radii[1:] / radii[1:]

            self.table_gravity = g[::-1]
        else:
            raise ValueError(
                "Density profile is an empty list " "in this SeismicTable instance."
            )

    def _compute_pressure(self):
        # Calculate the pressure profile based on density and gravity.
        # This integrates the equation for hydrostatic equilibrium P = rho g z.
        radii = self.table_radius
        density = self.table_density
        gravity = self.gravity(self.earth_radius - radii)

        # convert radii to depths
        depth = self.earth_radius - radii
        pressure = scipy.integrate.cumulative_trapezoid(
            gravity * density, x=depth, initial=0
        )

        self.table_pressure = pressure


class PREM(SeismicTable):
    """
    Reads PREM (1s) (input_seismic/prem.txt, :cite:`dziewonski1981`).
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        SeismicTable.__init__(self)
        table = read_table("input_seismic/prem.txt")
        table = np.array(table)
        self.table_depth = table[:, 0]
        self.table_radius = table[:, 1]
        self.table_pressure = table[:, 2]
        self.table_density = table[:, 3]
        self.table_vp = table[:, 4]
        self.table_vs = table[:, 5]
        self.table_QK = table[:, 6]
        self.table_QG = table[:, 7]


class Slow(SeismicTable):
    """
    Inserts the mean profiles for slower regions in the lower mantle
    (Lekic et al. 2012). We stitch together tables
    'input_seismic/prem_lowermantle.txt',
    'input_seismic/swave_slow.txt',
    'input_seismic/pwave_slow.txt').
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        SeismicTable.__init__(self)

        # data is: depth radius pressure density V_p V_s Q_K Q_G
        table = read_table("input_seismic/prem.txt")
        table = np.array(table)
        table2 = read_table("input_seismic/swave_slow.txt")
        table2 = np.array(table2)
        table3 = read_table("input_seismic/pwave_slow.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius - max(table2[:, 0])
        max_radius = self.earth_radius - min(table2[:, 0])

        table = np.array(
            list(filter(lambda x: (x[1] >= min_radius and x[1] <= max_radius), table))
        )

        self.table_depth = table[:, 0]
        self.table_radius = table[:, 1]
        self.table_pressure = table[:, 2]
        self.table_density = table[:, 3]
        self.table_vp = np.interp(
            self.table_depth, table3[:, 0][::-1], table3[:, 1][::-1]
        )
        self.table_vs = np.interp(
            self.table_depth, table2[:, 0][::-1], table2[:, 1][::-1]
        )


class Fast(SeismicTable):
    """
    Inserts the mean profiles for faster regions in the lower mantle
    (Lekic et al. 2012). We stitch together tables
    'input_seismic/prem_lowermantle.txt',
    'input_seismic/swave_fast.txt',
    'input_seismic/pwave_fast.txt').
    See also :class:`burnman.seismic.Seismic1DModel`.
    """

    def __init__(self):
        SeismicTable.__init__(self)

        # data is: radius pressure density V_p V_s Q_K Q_G
        table = read_table("input_seismic/prem.txt")
        table = np.array(table)
        table2 = read_table("input_seismic/swave_fast.txt")
        table2 = np.array(table2)
        table3 = read_table("input_seismic/pwave_fast.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius - max(table2[:, 0])
        max_radius = self.earth_radius - min(table2[:, 0])

        table = np.array(
            list(filter(lambda x: (x[1] >= min_radius and x[1] <= max_radius), table))
        )

        self.table_depth = table[:, 0]
        self.table_radius = table[:, 1]
        self.table_pressure = table[:, 2]
        self.table_density = table[:, 3]
        self.table_vp = np.interp(
            self.table_depth, table3[:, 0][::-1], table3[:, 1][::-1]
        )
        self.table_vs = np.interp(
            self.table_depth, table2[:, 0][::-1], table2[:, 1][::-1]
        )


class STW105(SeismicTable):
    """
    Reads STW05 (a.k.a. REF) (1s) (input_seismic/STW105.txt, :cite:`kustowski2008`).
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        SeismicTable.__init__(self)
        # radius, pressure, density, v_p, v_s
        table = read_table("input_seismic/STW105.txt")
        table = np.array(table)
        self.table_radius = table[:, 0][::-1]
        self.table_density = table[:, 1][::-1]
        self.table_vpv = table[:, 2][::-1]
        self.table_vsv = table[:, 3][::-1]
        self.table_QK = table[:, 4][::-1]
        self.table_QG = table[:, 5][::-1]
        self.table_vph = table[:, 6][::-1]
        self.table_vsh = table[:, 7][::-1]

        self.table_depth = self.earth_radius - self.table_radius

        # Voigt averages for Vs and Vp
        self.table_vs = np.sqrt(
            (2.0 * self.table_vsv * self.table_vsv + self.table_vsh * self.table_vsh)
            / 3.0
        )
        self.table_vp = np.sqrt(
            (self.table_vpv * self.table_vpv + 4.0 * self.table_vph * self.table_vph)
            / 5.0
        )


class IASP91(SeismicTable):
    """
    Reads REF/STW05 (input_seismic/STW105.txt, :cite:`kustowski2008`).
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        SeismicTable.__init__(self)
        table = read_table("input_seismic/iasp91.txt")  # depth, radius, v_p, v_s
        table = np.array(table)
        self.table_depth = table[:, 0]
        self.table_radius = table[:, 1]
        self.table_vp = table[:, 2]
        self.table_vs = table[:, 3]


class AK135(SeismicTable):
    """
    Reads AK135 (input_seismic/ak135.txt, :cite:`kennett1995`).
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        SeismicTable.__init__(self)
        table = read_table(
            "input_seismic/ak135.txt"
        )  # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[:, 0]
        self.table_radius = table[:, 1]
        self.table_density = table[:, 2]
        self.table_vp = table[:, 3]
        self.table_vs = table[:, 4]
        self.table_QG = table[:, 5]
        self.table_QK = table[:, 6]


def attenuation_correction(v_p, v_s, v_phi, Qs, Qphi):
    """
    Applies the attenuation correction following Matas et al. (2007), page 4.
    This is simplified, and there is also currently no 1D Q model implemented.
    The correction, however, only slightly reduces the velocities,
    and can be ignored for our current applications.
    Arguably, it might not be as relevant when comparing computations
    to PREM for periods of 1s as is implemented here.
    Called from :func:`burnman.main.apply_attenuation_correction`

    :param v_p: P wave velocity in [m/s].
    :type v_p: float
    :param v_s: S wave velocitiy in [m/s].
    :type v_s: float
    :param v_phi: Bulk sound velocity in [m/s].
    :type v_phi: float
    :param Qs: Shear quality factor [dimensionless].
    :type Qs: float
    :param Qphi: Bulk quality factor [dimensionless].
    :type Qphi: float

    :returns: Corrected P wave, S wave and bulk sound velocities in [m/s].
    :rtype: tuple
    """
    beta = 0.3  # Matas et al. (2007) page 4
    Qp = 3.0 / 4.0 * pow((v_p / v_s), 2.0) * Qs  # Matas et al. (2007) page 4

    cot = 1.0 / np.tan(beta * np.pi / 2.0)
    v_p *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qp  # Matas et al. (2007) page 1
    v_s *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qs
    v_phi *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qphi
    return v_p, v_s, v_phi


"""
shared variable of prem, so that other routines do not need to create
prem over and over. See geotherm for example.
"""
prem_model = PREM()
