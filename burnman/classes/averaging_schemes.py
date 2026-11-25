# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import warnings
from scipy.optimize import brentq


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
        Average the bulk moduli :math:`K` for a composite. This defines the
        interface for this method, and is not implemented in the base class.

        :param volumes: List of the volume of each phase
            in the composite :math:`[m^3]`.
        :type volumes: list of floats

        :param bulk_moduli: List of bulk moduli of each phase
            in the composite :math:`[Pa]`.
        :type bulk_moduli: list of floats

        :param shear_moduli: List of shear moduli of each phase
            in the composite :math:`[Pa]`.
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


class SelfConsistent(AveragingScheme):
    """
    Self-consistent elastic moduli for an aggregate of spherical grains.

    The aggregate is isotropic, and it is assumed that no phase plays
    a special role as an interconnected medium.

    The spherical self-consistent equations of Budiansky (1965),
    https://doi.org/10.1016/0022-5096(65)90011-6, are

    .. math::

        \\sum_i f_i\frac{K_i-K}{K_i+4G/3}=0, \\qquad
        \\sum_i f_i\frac{G_i-G}{G_i+\\zeta}=0,

        \\zeta=\frac{G(9K+8G)}{6(K+2G)}, \\qquad
        f_i=\\frac{V_i}{\\sum_j V_j}.

    Average bulk and shear moduli are computed by
    :func:`self_consistent_function`.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average bulk moduli using spherical self-consistent homogenization.

        :param volumes: Volumes of each phase. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: Bulk moduli of each phase. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: Shear moduli of each phase. :math:`[Pa]`
        :type shear_moduli: list of floats
        :returns: Self-consistent bulk modulus. :math:`[Pa]`
        :rtype: float
        """
        return self_consistent_function(volumes, bulk_moduli, shear_moduli)[0]

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average shear moduli using spherical self-consistent homogenization.

        :param volumes: Volumes of each phase. :math:`[m^3]`
        :type volumes: list of floats
        :param bulk_moduli: Bulk moduli of each phase. :math:`[Pa]`
        :type bulk_moduli: list of floats
        :param shear_moduli: Shear moduli of each phase. :math:`[Pa]`
        :type shear_moduli: list of floats
        :returns: Self-consistent shear modulus. :math:`[Pa]`
        :rtype: float
        """
        return self_consistent_function(volumes, bulk_moduli, shear_moduli)[1]


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
        Average the shear moduli of a composite with the arithmetic mean of the upper
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


def self_consistent_function(phase_volume, bulk_moduli, shear_moduli):
    """
    Estimate the effective bulk and shear moduli of a multiphase rock
    using the spherical self-consistent method.

    Each phase is described by its bulk and shear moduli and the volume it
    occupies in the rock. Phase properties should be evaluated at the same
    pressure and temperature. Volumes are converted internally to volume
    fractions, so they do not need to sum to one.

    The model represents each phase as spherical grains surrounded by a
    uniform elastic medium with the unknown effective properties of the
    whole rock. It is "self-consistent" because the combined response of
    these grains must reproduce the properties assumed for the surrounding
    medium. All phases are treated on an equal footing, and the resulting
    aggregate is isotropic. The deformation of a grain depends on both its
    own stiffness and that of its surroundings, which couples the effective
    bulk and shear moduli.

    The coupled equations are given in :class:`SelfConsistent`, following
    Budiansky (1965). To solve them, the function first calculates the
    effective bulk modulus :math:`K` for a trial effective shear modulus
    :math:`G`:

    .. math::

        K(G) = \frac{\\sum_i f_i K_i/(K_i+4G/3)}
                    {\\sum_i f_i/(K_i+4G/3)},
        \\qquad f_i = \frac{V_i}{\\sum_j V_j}.

    Here :math:`K_i`, :math:`G_i`, and :math:`f_i` are the moduli and
    volume fraction of phase :math:`i`. Substituting :math:`K(G)` into
    the shear equation leaves a single unknown:

    .. math::

        \\sum_i f_i\frac{G_i-G}{G_i+\\zeta}=0,
        \\qquad \\zeta = \frac{G(9K+8G)}{6(K+2G)}.

    A bracketed root finder adjusts :math:`G` between the smallest and
    largest phase shear moduli until this equation is satisfied. The
    corresponding :math:`K(G)` completes the solution.

    Fluids can be represented by zero shear modulus, and empty pores by
    zero bulk and shear moduli. If the model predicts that the mixture
    cannot support shear, the function returns :math:`G=0` and the
    corresponding limiting bulk modulus. For an elastic solid mixed with
    empty spherical pores, this model predicts that both moduli vanish
    at a porosity of 0.5. Zero-volume phases are ignored.

    :param phase_volume: Volume occupied by each phase, or its volume
        fraction. Entries must be finite and nonnegative, with at least
        one positive entry. All volumes must use the same units.
    :type phase_volume: list of floats or numpy.ndarray
    :param bulk_moduli: Bulk moduli of the corresponding phases. Values
        for phases with positive volume must be finite and nonnegative.
        Use the same pressure unit for all bulk and shear moduli,
        normally :math:`[Pa]`.
    :type bulk_moduli: list of floats or numpy.ndarray
    :param shear_moduli: Shear moduli of the corresponding phases, in
        the same units as ``bulk_moduli``. Values for phases with positive
        volume must be finite and nonnegative.
    :type shear_moduli: list of floats or numpy.ndarray
    :returns: Effective bulk and shear moduli, :math:`(K, G)`, in the same
        units as the input moduli.
    :rtype: tuple of floats
    :raises ValueError: If the arrays have incompatible shapes, the volumes
        are invalid, or an occupied phase has invalid elastic moduli.
    """
    volumes = np.asarray(phase_volume, dtype=float)
    bulk = np.asarray(bulk_moduli, dtype=float)
    shear = np.asarray(shear_moduli, dtype=float)
    if (
        volumes.ndim != 1
        or volumes.size == 0
        or bulk.shape != volumes.shape
        or shear.shape != volumes.shape
    ):
        raise ValueError(
            "Phase volumes and moduli must be nonempty 1D arrays of equal length."
        )
    if not np.all(np.isfinite(volumes)) or np.any(volumes < 0.0):
        raise ValueError("Phase volumes must be finite and nonnegative.")
    occupied = volumes > 0.0
    if not np.any(occupied):
        raise ValueError("At least one phase must have positive volume.")
    # Ignore absent phases before inspecting their elastic properties.
    volumes, bulk, shear = volumes[occupied], bulk[occupied], shear[occupied]
    if (
        not np.all(np.isfinite(bulk))
        or not np.all(np.isfinite(shear))
        or np.any(bulk < 0.0)
        or np.any(shear < 0.0)
    ):
        raise ValueError(
            "Occupied phases must have finite, nonnegative bulk and shear moduli."
        )
    if volumes.size == 1:
        return float(bulk[0]), float(shear[0])
    # Normalization makes the root tolerances scale-independent.
    fractions = volumes / np.max(volumes)
    fractions /= np.sum(fractions)
    scale = max(np.max(bulk), np.max(shear))
    if scale == 0.0:
        return 0.0, 0.0
    bulk, shear = bulk / scale, shear / scale
    zero_bulk = np.sum(fractions[bulk == 0.0])
    zero_shear = np.sum(fractions[shear == 0.0])

    def effective_bulk(g):
        if g == 0.0 and zero_bulk > 0.0:
            return 0.0
        denominators = bulk + 4.0 * g / 3.0
        # Weighted form avoids subtracting nearly equal harmonic terms.
        weights = fractions * (np.min(denominators) / denominators)
        return np.dot(weights, bulk) / np.sum(weights)

    # As G -> 0+, zeta/G -> (3-p)/(2+p), where p is the volume
    # fraction of zero-bulk phases. Thus the limiting shear residual is
    # 1-5*q/(3-p), with q the zero-shear fraction. This handles dry pores
    # (critical porosity 1/2) and fluids with K>0 (critical fraction 3/5)
    # without an arbitrary small stiffness or a spurious positive root.
    zero_limit = 1.0 - 5.0 * zero_shear / (3.0 - zero_bulk)
    if zero_shear > 0.0 and zero_limit <= 0.0:
        return float(scale * effective_bulk(0.0)), 0.0

    def shear_residual(g):
        if g == 0.0:
            return zero_limit
        k = effective_bulk(g)
        zeta = g * ((9.0 * k + 8.0 * g) / (6.0 * (k + 2.0 * g)))
        return np.sum(fractions * ((shear - g) / (shear + zeta)))

    lower, upper = np.min(shear), np.max(shear)
    if lower == upper:
        g = lower
    else:
        g = brentq(
            shear_residual,
            lower,
            upper,
            xtol=np.finfo(float).tiny,
            rtol=1.0e-12,
            maxiter=200,
        )
    return float(scale * effective_bulk(g)), float(scale * g)


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
