# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


class material:
    """
    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type burnman.mineral and their molar
    fractions.

    The user needs to call set_method() once in the beginning and set_state()
    before querying the material with unroll() or density().
    """

    def set_method(self, method):
        """
        """
        raise NotImplementedError("need to implement set_method() in derived class!")

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        """
        Set the material to the given pressure and temperature
        """
        self.pressure = pressure
        self.temperature = temperature

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        raise NotImplementedError("need to implement unroll() in derived class!")
        return ()

    def density(self):
        raise NotImplementedError("need to implement density() in derived class!")
        return inf            
