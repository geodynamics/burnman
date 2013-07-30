# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
BurnMan
=======

BurnMan is a lower mantle shear velocity generator constrained by mineral physics.



"""

from main import *
from partitioning import calculate_partition_coefficient,calculate_phase_percents
import minerals
import seismic
from composite import composite
from minerals_base import material
import averaging_schemes
