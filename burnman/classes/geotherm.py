# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.optimize import brentq

from ..utils.misc import read_table, lookup_and_interpolate
from ..utils.math import bracket
from .seismic import prem_model


class Geotherm:
    """
    Base class for all geotherms.
    Subclasses should implement a temperatures(depths)
    method or inherit the one provided here.

    :param depths_array: Array of depths (in meters).
    :type depths_array: np.array
    :param temperatures_array: Array of temperatures (in Kelvin).
    :type temperatures_array: np.array
    """

    def __init__(self, depths_array, temperatures_array):
        """
        Initialize the Geotherm with optional depth and temperature arrays.

        :param depths_array: Array of depths (in meters).
        :type depths_array: np.array
        :param temperatures_array: Array of temperatures (in Kelvin).
        :type temperatures_array: np.array
        """
        self.depths_array = np.asarray(depths_array)
        self.temperatures_array = np.asarray(temperatures_array)

        self.minimum_depth = self.depths_array.min()
        self.maximum_depth = self.depths_array.max()

    def temperatures(self, depths):
        """
        Evaluate the geotherm at given depths.

        :param depths: Depths at which to evaluate the geotherm (in meters).
        :type depths: array-like of float
        :return: Temperatures at the given depths (in Kelvin).
        :rtype: np.array of float
        :raises ValueError: If depths are outside the range of the geotherm.
        :raises ValueError: If depths_array or temperatures_array is None.
        """
        depths = np.asarray(depths)

        z_min = depths.min()
        z_max = depths.max()
        if z_min < self.minimum_depth:
            raise ValueError(
                f"Requested depth ({z_min}) m is shallower than "
                f"minimum geotherm depth ({self.minimum_depth}) m."
            )
        if z_max > self.maximum_depth:
            raise ValueError(
                f"Requested depth ({z_max}) m is deeper than "
                f"maximum geotherm depth ({self.maximum_depth}) m."
            )
        temperatures = np.empty_like(depths)

        for i, d in enumerate(depths):
            temperatures[i] = lookup_and_interpolate(
                self.depths_array, self.temperatures_array, d
            )

        return temperatures


class BrownShankland(Geotherm):
    """
    Geotherm from Brown and Shankland, 1981 (:cite:`Brown1981`)
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/brown_81.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class Anderson(Geotherm):
    """
    Geotherm from Anderson, 1982 (:cite:`anderson1982earth`)
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/anderson_82.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class StaceyContinental(Geotherm):
    """
    Continental geotherm from Stacey, 1977 (:cite:`stacey1977`)
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/Stacey_1977_continents.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class StaceyOceanic(Geotherm):
    """
    Oceanic geotherm from Stacey, 1977 (:cite:`stacey1977`)
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/Stacey_1977_oceans.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class Katsura2022(Geotherm):
    """
    Geotherm from Katsura, 2022 (:cite:`Katsura2022`)
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/Katsura_2022.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class Plesa2022Mars6cm3(Geotherm):
    """
    Martian areotherm from Plesa, 2022 (:cite:`Plesa2022`)
    assuming 6 cm³/mol activation volume for mantle viscosity.
    """

    def __init__(self):
        """
        No arguments need to be passed to the constructor.
        """
        table = read_table("input_geotherm/Plesa_2022_Mars_V_6cm3.txt")
        Geotherm.__init__(self, table[:, 0], table[:, 1])


class GeothermFromPressures(Geotherm):
    """
    Geotherm defined by pressures and temperatures.

    :param pressures_array: Array of pressures (in Pa).
    :type pressures_array: np.array
    :param temperatures_array: Array of temperatures (in Kelvin).
    :type temperatures_array: np.array
    :param depth_to_pressure_model: Object with an object.pressure(depth)
        method, such as :class:`burnman.seismic.SeismicModel`.
        If not given, the PREM model will be used.
    :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
        other object with a pressure(depth) method.
    """

    def __init__(
        self, pressures_array, temperatures_array, depth_to_pressure_model=prem_model
    ):
        """
        Initialize the GeothermFromPressures with pressure and temperature arrays, and
        a depth-to-pressure model.

        :param pressures_array: Array of pressures (in Pa).
        :type pressures_array: np.array
        :param temperatures_array: Array of temperatures (in Kelvin).
        :type temperatures_array: np.array
        :param depth_to_pressure_model: Object with an object.pressure(depth)
            method, such as :class:`burnman.seismic.SeismicModel`.
            If not given, the PREM model will be used.
        :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
            other object with a pressure(depth) method.
        """
        self.pressures_array = np.asarray(pressures_array)
        self.temperatures_array = np.asarray(temperatures_array)
        self.depth_to_pressure_model = depth_to_pressure_model

        self.minimum_pressure = self.pressures_array.min()
        self.maximum_pressure = self.pressures_array.max()

    def temperatures_from_pressures(self, pressures):
        """
        Evaluate the geotherm at given pressures.

        :param pressures: Pressures at which to evaluate the geotherm (in Pa).
        :type pressures: array-like of float
        :return: Temperatures at the given pressures (in Kelvin).
        :rtype: np.ndarray of float
        """

        pressures = np.asarray(pressures)

        P_min = pressures.min()
        P_max = pressures.max()
        if P_min < self.minimum_pressure:
            raise ValueError(
                f"Requested pressure ({P_min}) Pa is lower than "
                f"minimum geotherm pressure ({self.minimum_pressure}) Pa."
            )
        if P_max > self.maximum_pressure:
            raise ValueError(
                f"Requested pressure ({P_max}) Pa is higher than "
                f"maximum geotherm pressure ({self.maximum_pressure}) Pa."
            )

        temperatures = np.empty_like(pressures)
        for i, P in enumerate(pressures):
            temperatures[i] = lookup_and_interpolate(
                self.pressures_array, self.temperatures_array, P
            )

        return temperatures

    def temperatures(self, depths):
        """
        Evaluate the geotherm at given depths.

        :param depths: Depths at which to evaluate the geotherm (in meters).
        :type depths: array-like of float
        :return: Temperatures at the given depths (in Kelvin).
        :rtype: np.ndarray of float
        """
        depths = np.asarray(depths)
        pressures = self.depth_to_pressure_model.pressure(depths)
        return self.temperatures_from_pressures(pressures)


class Anzellini2013(GeothermFromPressures):
    """
    Geotherm from Anzellini et al., 2013 (:cite:`Anzellini2013`) (their Figure S4).
    Mantle part is from Steinberger and Holme, 2008 (:cite:`Steinberger2008`).

    :param depth_to_pressure_model: Object with an object.pressure(depth)
        method, such as :class:`burnman.seismic.SeismicModel`.
        If not given, the PREM model will be used.
    :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
        other object with a pressure(depth) method.
    """

    def __init__(self, depth_to_pressure_model=prem_model):
        """
        Initialize the Anzellini et al., 2013 (:cite:`Anzellini2013`) geotherm.

        :param depth_to_pressure_model: Object with an object.pressure(depth)
            method, such as :class:`burnman.seismic.SeismicModel`.
            If not given, the PREM model will be used.
        :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
            other object with a pressure(depth) method, optional
        """
        table = read_table("input_geotherm/Anzellini_2013_PT.txt")
        GeothermFromPressures.__init__(
            self,
            table[:, 0],
            table[:, 1],
            depth_to_pressure_model=depth_to_pressure_model,
        )


class AdiabaticGeotherm(GeothermFromPressures):
    """
    Geotherm calculated as an adiabat for a given material, anchored
    on a specified temperature and pressure.

    In addition to the temperature method of the standard Geotherm class,
    this class also provides a temperatures_from_pressures method.

    Initialisation parameters:

    :param material: :class:`burnman.Material` object.
    :type material: burnman.Material
    :param T_0: Temperature at the anchor pressure [K].
    :type T_0: float
    :param P_0: The anchor pressure [Pa].
    :type P_0: float
    :param depth_to_pressure_model: Object with an object.pressure(depth)
        method, such as :class:`burnman.seismic.SeismicModel`.
        If not given, the PREM model will be used.
    :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
        other object with a pressure(depth) method.
    """

    def __init__(self, material, T_0, P_0, depth_to_pressure_model=prem_model):
        """
        Initialize the AdiabaticGeotherm.
        :param material: :class:`burnman.Material` object.
        :type material: burnman.Material
        :param T_0: Temperature at the anchor pressure [K].
        :type T_0: float
        :param P_0: The anchor pressure [Pa].
        :type P_0: float
        :param depth_to_pressure_model: Object with an object.pressure(depth)
            method, such as :class:`burnman.seismic.SeismicModel`.
            If not given, the PREM model will be used.
        :type depth_to_pressure_model: :class:`burnman.seismic.SeismicModel` or
            other object with a pressure(depth) method.
        """
        self.material = material
        self.T_0 = T_0
        self.P_0 = P_0
        self.depth_to_pressure_model = depth_to_pressure_model

    def temperatures_from_pressures(self, pressures):
        """
        Evaluate the adiabat geotherm at given pressures.

        :param pressures: Pressures at which to evaluate the geotherm (in Pa).
        :type pressures: array-like of float
        :return: Temperatures at the given pressures (in Kelvin).
        :rtype: np.ndarray of float
        """
        return adiabatic_profile(pressures, self.material, self.T_0, self.P_0)


def adiabatic_profile(pressures, rock, T_0, P_0=None):
    """
    This function calculates an adiabat for a rock
    anchored on a specified temperature and pressure.

    A starting guess is provided by integrating:

    .. math::
        \\frac{\\partial T}{\\partial P} = \\frac{ \\gamma  T}{ K_s }

    from the anchor point. In this expression, :math:`\\gamma`
    is the Grueneisen parameter and :math:`K_s` is
    the adiabatic bulk modulus. The exact adiabatic profile is then found
    by finding the temperature profile along which the entropy is constant.

    :param pressures: The list of pressures in :math:`[Pa]` at which
        to evaluate the geotherm.
    :type pressures: list of floats

    :param rock: Composite for which we compute the adiabat.
        From this material we must compute average Grueneisen parameters
        and adiabatic bulk moduli for each pressure/temperature.
    :type rock: :class:`burnman.composite`

    :param T_0: An anchor temperature, corresponding to the temperature at the
        anchor pressure. :math:`[K]`
    :type T_0: float

    :param P_0: An anchor pressure. If not given, the first pressure in
        the pressures list will be used.
    :type P_0: float, optional

    :returns: The list of temperatures for each pressure. :math:`[K]`
    :rtype: numpy.array of floats
    """
    pressures = np.asarray(pressures)

    if P_0 is None:
        P_0 = pressures[0]

    rock.set_state(P_0, T_0)
    S0 = rock.S

    def delta_S(T, P, rock, S0):
        return S0 - rock.evaluate(["S"], [P], [T])[0][0]

    temperatures = np.empty_like(pressures)

    # estimate the first temperature
    if P_0 != pressures[0]:
        args = (pressures[0], rock, S0)
        sol = bracket(
            fn=delta_S,
            x0=T_0 + (rock.gr * T_0 / rock.K_S) * (pressures[0] - P_0),
            dx=1.0,
            args=args,
        )
        temperatures[0] = brentq(delta_S, sol[0], sol[1], args=args)
    else:
        temperatures[0] = T_0

    for i in range(1, len(pressures)):
        args = (pressures[i], rock, S0)
        sol = bracket(
            fn=delta_S,
            x0=(
                temperatures[i - 1]
                + (rock.gr * temperatures[i - 1] / rock.K_S)
                * (pressures[i] - pressures[i - 1])
            ),
            dx=1.0,
            args=args,
        )
        temperatures[i] = brentq(delta_S, sol[0], sol[1], args=args)

    return temperatures
