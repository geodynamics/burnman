# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


class Material:
    """
    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type :class:`~burnman.mineral.Mineral` and their molar
    fractions. This class is available as ``burnman.Material``.

    The user needs to call set_method() (once in the beginning) and set_state()
    before querying the material with unroll() or density().

    Attributes
    ----------
    pressure : float
        The current pressure as set by :func:`~burnman.Material.set_state`.
    temperature : float
        The current temperature as set by :func:`~burnman.Material.set_state`.
    """

    def __init__(self):
        self.pressure = None
        self.temperature = None

    def set_method(self, method):
        """
        Set the averaging method. See :doc:`averaging` for details.

        Notes
        -----
        Needs to be implemented in derived classes.
        """
        raise NotImplementedError("need to implement set_method() in derived class!")

    def to_string(self):
        """
        Returns a human-readable name of this material. The default implementation will return the name of the class,
        which is a reasonable default.

        Returns
        -------
        name : string
            Name of this material.
        """
        return "'" + self.__class__.__name__ + "'"

    def debug_print(self):
        """
        Print a human-readable representation of this Material.
        """
        (frs,mins) = self.unroll()
        if len(mins)==1:
            print mins[0].to_string()
        else:
            print "Material %s:" % self.to_string()
            for (fr,mi) in zip(frs,mins):
                print "  %g of phase %s" % (fr, mi.to_string())

    def set_state(self, pressure, temperature):
        """
        Set the material to the given pressure and temperature.

        Parameters
        ----------
        pressure : float
            The desired pressure in Pa.
        temperature : float
            The desired temperature in K.
        """
        self.pressure = pressure
        self.temperature = temperature

    def unroll(self):
        """
        Unroll this material into a list of :class:`burnman.Mineral` and their molar fractions. All averaging schemes
        then operate on this list of minerals. Note that the return value of this function may depend on the current
        state (temperature, pressure).

        Notes
        -----
        Needs to be implemented in derived classes.

        Returns
        -------
        fractions : list of float
            List of molar fractions, should sum to 1.0.
        minerals : list of :class:`burnman.Mineral`
            List of minerals.
        """
        raise NotImplementedError("need to implement unroll() in derived class!")
        return ()

    def density(self):
        """
        Returns the density of this material. Note that the return value of this function may depend on the current
        state (temperature, pressure).

        Notes
        -----
        Needs to be implemented in derived classes.

        Returns
        -------
        density : float
            The density of this material in kg/m^3

        """
        raise NotImplementedError("need to implement density() in derived class!")
        return None
