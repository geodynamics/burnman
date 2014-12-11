# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""BurnMan
=======

BurnMan is an open source mineral physics toolbox written in Python which
determines seismic velocities for the lower mantle. BurnMan calculates the
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

Please cite:

  - Cottaar S., Heister, T., Rose, I., and Unterborn, C., 2014, BurnMan: A
    lower mantle mineral physics toolkit, Geochemistry, Geophysics, and
    Geosystems, 15(4), 1164-1179
  
Acknowledgement and Support:

  - This project was initiated at, and follow-up research support was received
    through, Cooperative Institute of Deep Earth Research, CIDER (NSF FESD
    grant 1135452) -- see www.deep-earth.org

  - We thank all the fellow members of the Cider Mg/Si team for their input:
    Valentina Magni, Yu Huang, JiaChao Liu, Marc Hirschmann, and Barbara
    Romanowicz.

  - We thank Lars Stixrude for providing benchmarking calculations.

  - We thank CIG (www.geodynamics.org) for support and accepting our donation
    of BurnMan as an official project.

  - We also welcomed helpful discussions with Zack Geballe, Motohiko Murakami,
    Bill McDonough, Quentin Williams, Wendy Panero, and Wolfgang Bangerth.

"""

from version import version as __version__

#classes for representing rocks and minerals
from mineral import Mineral
from material import Material
from composite import Composite
from mineral_helpers import *

#high level functions
from main import *
from model import Model

#mineral library
import minerals

#central user tools
import seismic
import averaging_schemes
import geotherm

#miscellaneous
import tools
from partitioning import calculate_partition_coefficient,calculate_phase_percents
