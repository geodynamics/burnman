# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
BurnMan
=======

BurnMan is an open source mineral physics toolbox written in Python to
determine seismic velocities for the lower mantle. BurnMan calculates the
isotropic thermoelastic moduli by solving the equations-of-state for a
mixture of minerals defined by the user. The user may select from a list of
minerals applicable to the lower mantle included or easily define one of
their own.

Features:
  - form composites of arbitrary combination of :doc:`minerals`
  - extensive :doc:`mineral_database`
  - easy plotting and comparison of seismic profiles using matplotlib
  - many examples highlighting different features of BurnMan
  - different thermoelastic models, choice between second or third order accuracy
  - different averaging schemes
  - different geotherms
  - extensible: all parts can be replaced by user-written modules if desired



"""

from version import version as __version__

#classes for representing rocks and minerals
from mineral import Mineral
from material import Material
from composite import Composite
from mineral_helpers import *

#high level functions
from main import *

#mineral library
import minerals

#central user tools
import seismic
import averaging_schemes
import geotherm

#miscellaneous
import tools
from partitioning import calculate_partition_coefficient,calculate_phase_percents
