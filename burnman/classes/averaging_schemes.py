# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
import warnings


class AveragingScheme(object):
    """
    Base class defining an interface for determining average
    elastic properties of a rock.  Given a list of volume
    fractions for the different mineral phases in a rock,
    as well as their bulk and shear moduli, an averaging
    will give back a single scalar values for the averages.
    New averaging schemes should define the functions
    average_bulk_moduli and average_shear_moduli, as
    specified here.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli :math:`K` for a composite. This defines the interface
        for this method, and is not implemented in the base class.

        :param volumes: List of the volume of each phase in the composite
            :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli : List of bulk moduli of each phase in the composite
            :math:`[Pa]`.
        :type bulk_moduli: list of floats
        :param shear_moduli : List of shear moduli of each phase in the composite
            :math:`[Pa]`.
        :type shear_moduli: list of floats

        :returns: The average bulk modulus :math:`K`. :math:`[Pa]`
        :rtype: float
        """
        raise NotImplementedError("")

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli :math:`G` for a composite.  This defines the interface
        for this method, and is not implemented in the base class.

        :param volumes: List of the volume of each phase in the composite
            :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli of each phase in the composite
            :math:`[Pa]`.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli of each phase in the composite
            :math:`[Pa]`.
        :type shear_moduli: list of floats

        :returns: The average shear modulus :math:`G`. :math:`[Pa]`
        :rtype: float
        """
        raise NotImplementedError("")

    def average_density(self, volumes, densities):
        """
        Average the densities of a composite, given a list of volume
        fractions and densitites. This is implemented in the base class,
        as how to calculate it is not dependent on the geometry of the rock.
        The formula for density is given by

        .. math::
            \\rho = \\frac{\\Sigma_i \\rho_i V_i }{\\Sigma_i V_i}

        :param volumes: List of the volume of each phase in the composite
            :math:`[m^3]`.
        :type volumes: list of floats
        :param densities: List of densities of each phase in the composite
            :math:`[kg/m^3]`.
        :type densities: list of floats

        :returns: Density :math:`\\rho` :math:`[kg/m^3]`.
        :rtype: float
        """
        total_mass = np.sum(np.array(densities) * np.array(volumes))
        total_vol = np.sum(np.array(volumes))  # should sum to one
        density = total_mass / total_vol
        return density

    def average_thermal_expansivity(self, volumes, alphas):
        """
        Averages the thermal expansion coefficient of the mineral :math:`\\alpha`
        :math:`[1/K]`.

        :param volumes: List of volume fractions of each phase
            in the composite (should sum to 1.0).
        :type volumes: list of floats
        :param alphas: List of thermal expansivities :math:`\\alpha`
            of each phase in the composite. :math:`[1/K]`
        :type alphas: list of floats

        :returns: Thermal expansivity of the composite :math:`\\alpha`. :math:`[1/K]`
        :rtype: float
        """
        total_vol = np.sum(np.array(volumes))
        return np.sum(np.array(alphas) * np.array(volumes)) / total_vol

    def average_heat_capacity_v(self, fractions, c_v):
        # TODO: double-check that the formula we use is appropriate here.
        """
        Averages the heat capacities at constant volume :math:`C_V` by molar fractions
        as in eqn. (16) in :cite:`Ita1992`.

        :param fractions: List of molar fractions of each phase
            in the composite (should sum to 1.0).
        :type fractions: list of floats
        :param c_v: List of heat capacities at constant volume :math:`C_V`
            of each phase in the composite. :math:`[J/K/mol]`
        :type c_v: list of floats

        :returns: Heat capacity at constant volume of the composite :math:`C_V`
            :math:`[J/K/mol]`.
        :rtype: float
        """
        return np.sum(np.array(fractions) * np.array(c_v))

    def average_heat_capacity_p(self, fractions, c_p):
        # TODO: double-check that the formula we use is correct.
        """
        Averages the heat capacities at constant pressure :math:`C_P`
        by molar fractions.

        :param fractions: List of molar fractions of each phase in the composite
            (should sum to 1.0).
        :type fractions: list of floats
        :param c_p: List of heat capacities at constant pressure :math:`C_P` of each
            phase in the composite :math:`[J/K/mol]`.
        :type c_p: list of floats

        :returns: Heat capacity at constant pressure :math:`C_P` of the composite
            :math:`[J/K/mol]`.
        :rtype: float
        """
        return np.sum(np.array(fractions) * np.array(c_p))


class VoigtReussHill(AveragingScheme):
    """
    Class for computing the Voigt-Reuss-Hill average for elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli`
    and
    :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli`
    functions.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the Voigt-Reuss-Hill average,
        given by:

        .. math::
            K_{VRH} = \\frac{K_V + K_R}{2}

        This is simply a shorthand for an arithmetic average of the bounds given
        by :class:`burnman.averaging_schemes.voigt`
        and :class:`burnman.averaging_schemes.reuss`.

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
            Not used in this average.
        :type shear_moduli: list of floats

        :returns: The Voigt-Reuss-Hill average bulk modulus :math:`K_{VRH}`
            :math:`[Pa]`.
        :rtype: float
        """
        return voigt_reuss_hill_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli :math:`G` of a composite with the
        Voigt-Reuss-Hill average, given by:

        .. math::
            G_{VRH} = \\frac{G_V + G_R}{2}

        This is simply a shorthand for an arithmetic average of the bounds given
        by :class:`burnman.averaging_schemes.voigt`
        and :class:`burnman.averaging_schemes.reuss`.

        .. math::
            G_V = \\Sigma_i V_i G_i

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
            Not used in this average.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
        :type shear_moduli: list of floats

        :returns: The Voigt-Reuss-Hill average shear modulus :math:`G_{VRH}`
            :math:`[Pa]`.
        :rtype: float
        """
        return voigt_reuss_hill_function(volumes, shear_moduli)


class Voigt(AveragingScheme):
    """
    Class for computing the Voigt (iso-strain) bound for elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli` and
    :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli` functions.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite :math:`K` with the Voigt (iso-strain)
        bound, given by:

        .. math::
            K_V = \\Sigma_i V_i K_i

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
            Not used in this average.
        :type shear_moduli: list of floats

        :returns: The Voigt average bulk modulus :math:`K_R` :math:`[Pa]`.
        :rtype: float
        """
        return voigt_average_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of a composite with the Voigt (iso-strain)
        bound, given by:

        .. math::
            G_V = \\Sigma_i V_i G_i

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
            Not used in this average.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
        :type shear_moduli: list of floats

        :returns: The Voigt average shear modulus :math:`G_V` :math:`[Pa]`.
        :rtype: float
        """
        return voigt_average_function(volumes, shear_moduli)


class Reuss(AveragingScheme):
    """
    Class for computing the Reuss (iso-stress) bound for elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli` and
    :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli` functions.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the Reuss (iso-stress)
        bound, given by:

        .. math::
            K_R = \\left(\\Sigma_i \\frac{V_i}{K_i} \\right)^{-1}

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
            Not used in this average.
        :type shear_moduli: list of floats

        :returns: The Reuss average bulk modulus :math:`K_R` :math:`[Pa]`.
        :rtype: float
        """
        return reuss_average_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of a composite with the Reuss (iso-stress)
        bound, given by:

        .. math::
            G_R = \\left( \\Sigma_i \\frac{V_i}{G_i} \\right)^{-1}

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K`
            of each phase in the composite :math:`[Pa]`.
            Not used in this average.
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite :math:`[Pa]`.
        :type shear_moduli: list of floats

        :returns: The Reuss average shear modulus :math:`G_R` :math:`[Pa]`.
        :rtype: float
        """
        return reuss_average_function(volumes, shear_moduli)


class HashinShtrikmanUpper(AveragingScheme):
    """
    Class for computing the upper Hashin-Shtrikman bound for elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli`
    and :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli`
    functions. Implements formulas from :cite:`Watt1976`.
    The Hashin-Shtrikman bounds
    are tighter than the Voigt and Reuss bounds because they make the
    additional assumption that the orientation of the phases are statistically
    isotropic.  In some cases this may be a good assumption, and in others it
    may not be.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the upper Hashin-Shtrikman bound.
        Implements Formulas from :cite:`Watt1976`, which are too lengthy to reproduce
        here.

        :param volumes: List of the volumes of each phase
            in the composite. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The upper Hashin-Shtrikman average bulk modulus :math:`K`.
            :math:`[Pa]`
        :rtype: float
        """

        K_n = max(bulk_moduli)
        G_n = max(shear_moduli)

        vol_frac = volumes / sum(volumes)

        alpha_n = -3.0 / (3.0 * K_n + 4.0 * G_n)
        A_n = 0
        for i in range(len(vol_frac)):
            if bulk_moduli[i] != K_n:
                A_n += vol_frac[i] / (1.0 / (bulk_moduli[i] - K_n) - alpha_n)

        K_upper = K_n + A_n / (1.0 + alpha_n * A_n)
        return K_upper

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of a composite with the upper Hashin-Shtrikman bound.
        Implements Formulas from :cite:`Watt1976`, which are too lengthy to reproduce
        here.

        :param volumes: List of the volumes of each phase
            in the composite. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The upper Hashin-Shtrikman average shear modulus :math:`G`.
            :math:`[Pa]`
        :rtype: float
        """

        K_n = max(bulk_moduli)
        G_n = max(shear_moduli)

        vol_frac = volumes / sum(volumes)

        beta_n = -3.0 * (K_n + 2.0 * G_n) / (5.0 * G_n * (3.0 * K_n + 4.0 * G_n))
        B_n = 0
        for i in range(len(vol_frac)):
            if shear_moduli[i] != G_n:
                B_n += vol_frac[i] / (1.0 / (2.0 * (shear_moduli[i] - G_n)) - beta_n)

        G_upper = G_n + (0.5) * B_n / (1.0 + beta_n * B_n)
        return G_upper


class HashinShtrikmanLower(AveragingScheme):
    """
    Class for computing the lower Hashin-Shtrikman bound for elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli`
    and :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli`
    functions.
    Implements Formulas from :cite:`Watt1976`.  The Hashin-Shtrikman bounds
    are tighter than the Voigt and Reuss bounds because they make the
    additional assumption that the orientation of the phases are statistically
    isotropic.  In some cases this may be a good assumption, and in others it
    may not be.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the lower Hashin-Shtrikman bound.
        Implements Formulas from :cite:`Watt1976`, which are too lengthy to reproduce
        here.

        :param volumes: List of the volumes of each phase
            in the composite. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The lower Hashin-Shtrikman average bulk modulus :math:`K`.
            :math:`[Pa]`
        :rtype: float
        """

        K_1 = min(bulk_moduli)
        G_1 = min(shear_moduli)

        vol_frac = volumes / sum(volumes)

        alpha_1 = -3.0 / (3.0 * K_1 + 4.0 * G_1)
        A_1 = 0
        for i in range(len(vol_frac)):
            if bulk_moduli[i] != K_1:
                A_1 += vol_frac[i] / (1.0 / (bulk_moduli[i] - K_1) - alpha_1)

        K_lower = K_1 + A_1 / (1.0 + alpha_1 * A_1)
        return K_lower

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of a composite with the lower Hashin-Shtrikman bound.
        Implements Formulas from :cite:`Watt1976`,
        which are too lengthy to reproduce here.

        :param volumes: List of the volumes of each phase
            in the composite. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The lower Hashin-Shtrikman average shear modulus :math:`G`.
            :math:`[Pa]`
        :rtype: float
        """

        K_1 = min(bulk_moduli)
        G_1 = min(shear_moduli)

        vol_frac = volumes / sum(volumes)

        beta_1 = -3.0 * (K_1 + 2.0 * G_1) / (5.0 * G_1 * (3.0 * K_1 + 4.0 * G_1))
        B_1 = 0
        for i in range(len(vol_frac)):
            if shear_moduli[i] != G_1:
                B_1 += vol_frac[i] / (1.0 / (2.0 * (shear_moduli[i] - G_1)) - beta_1)

        G_lower = G_1 + (0.5) * B_1 / (1.0 + beta_1 * B_1)
        return G_lower


class HashinShtrikmanAverage(AveragingScheme):
    """
    Class for computing arithmetic mean of the Hashin-Shtrikman bounds on
    elastic properties.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`,
    and implements the
    :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli`
    and
    :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli`
    functions.
    """

    def __init__(self):
        self.upper = HashinShtrikmanUpper()
        self.lower = HashinShtrikmanLower()

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the arithmetic mean of the upper
        and lower Hashin-Shtrikman bounds.

        :param volumes: List of the volumes of each phase
            in the composite. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. Not used in this average. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The arithmetic mean of the Hashin-Shtrikman bounds on bulk modulus
            :math:`K` :math:`[Pa]`.
        :rtype: float
        """
        return (
            self.upper.average_bulk_moduli(volumes, bulk_moduli, shear_moduli)
            + self.lower.average_bulk_moduli(volumes, bulk_moduli, shear_moduli)
        ) / 2.0

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the arithmetic mean of the upper
        and lower Hashin-Shtrikman bounds.

        :param volumes: List of the volumes of each phase
            in the composite. [m^3].
        :type volumes: list of floats
        :param bulk_moduli: List of bulk moduli :math:`K` of each phase
            in the composite. Not used in this average. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: List of shear moduli :math:`G` of each phase
            in the composite. :math:`[Pa]`
        :type shear_moduli: list of floats

        :returns: The arithmetic mean of the Hashin-Shtrikman bounds on shear modulus
            :math:`G` :math:`[Pa]`.
        :rtype: float
        """
        return (
            self.upper.average_shear_moduli(volumes, bulk_moduli, shear_moduli)
            + self.lower.average_shear_moduli(volumes, bulk_moduli, shear_moduli)
        ) / 2.0


def voigt_average_function(phase_volume, X):
    """
    Calculates the Voigt (iso-strain) average.  Rather like
    resistors in series.  Called by voigt and
    voigt_reuss_hill classes.

    :param phase_volume: Phase volumes.
    :type phase_volume: List of floats

    :param X: Phase moduli.
    :type X: List of floats

    :returns: Voigt-averaged modulus
    :rtype: float
    """
    vol_frac = phase_volume / np.sum(phase_volume)
    X_voigt = sum(f * x for f, x in zip(vol_frac, X))
    return X_voigt


def reuss_average_function(phase_volume, X):
    """
    Calculates the Reuss (iso-stress) average.  Rather like
    resistors in parallel.  Called by reuss and
    voigt_reuss_hill classes.


    :param phase_volume: Phase volumes.
    :type phase_volume: List of floats

    :param X: Phase moduli.
    :type X: List of floats

    :returns: Reuss-averaged modulus
    :rtype: float
    """
    vol_frac = phase_volume / np.sum(phase_volume)
    for f, x in zip(vol_frac, X):
        if x <= 0 and np.abs(f) > np.finfo(float).eps:
            warnings.warn("Oops, called reuss_average with Xi<=0!")
            return 0.0
    X_reuss = 1.0 / sum(f / x for f, x in zip(vol_frac, X))
    return X_reuss


def voigt_reuss_hill_function(phase_volume, X):
    """
    Calculates the Voigt-Reuss-Hill average (arithmetic mean
    of Voigt and Reuss bounds).  Called by
    voigt_reuss_hill class.

    :param phase_volume: Phase volumes.
    :type phase_volume: List of floats

    :param X: Phase moduli.
    :type X: List of floats

    :returns: Voigt-Reuss-Hill-averaged modulus
    :rtype: float
    """
    X_vrh = (
        voigt_average_function(phase_volume, X)
        + reuss_average_function(phase_volume, X)
    ) / 2.0
    return X_vrh
