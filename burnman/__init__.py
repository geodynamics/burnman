# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
BurnMan
=======

BurnMan is an open source mineral physics toolbox written in Python which
determines the velocities of seismic waves in mineral assemblages at
high pressure and temperature. It was designed to calculate seismic
velocities in the lower mantle, but is equally suited to any part
of the solid Earth (or indeed any of the terrestrial planets).
BurnMan calculates the isotropic thermoelastic moduli by solving the
equations-of-state for a mixture of minerals defined by the user. The user
may select from an extensive list of minerals obtained from published
databases. Alternatively, they can easily define their own minerals.

Features:

  - a range of thermoelastic models, choice between second or third order accuracy
  - a range of thermodynamic models for mineral endmembers
  - consistent, comprehensive treatment of minerals with solid solutions
  - form composites of arbitrary combination of :doc:`materials`
  - extensive :doc:`mineral_database`
  - easy plotting and comparison of seismic profiles using matplotlib
  - many examples highlighting different features of BurnMan
  - different averaging schemes for seismic velocities in composite materials
  - a catalogue of published geotherms
  - extensible: all parts can be replaced by user-written modules if desired

Please cite:

  - Cottaar S., Heister, T., Rose, I., and Unterborn, C., 2014, BurnMan: A
    lower mantle mineral physics toolkit, Geochemistry, Geophysics, and
    Geosystems, 15(4), 1164-1179 `(link) <http://dx.doi.org/10.1002/2013GC005122>`_

Acknowledgement and Support:

  - This project was initiated at, and follow-up research support was received
    through, Cooperative Institute of Deep Earth Research, CIDER (NSF FESD
    grant 1135452) -- see `www.deep-earth.org <http://www.deep-earth.org>`_

  - We thank all the fellow members of the CIDER Mg/Si team for their input:
    Valentina Magni, Yu Huang, JiaChao Liu, Marc Hirschmann, and Barbara
    Romanowicz.

  - We thank Lars Stixrude for providing benchmarking calculations.

  - We thank CIG (`www.geodynamics.org <http://www.geodynamics.org>`_) for support and accepting our donation
    of BurnMan as an official project.

  - We also welcomed helpful discussions with Zack Geballe, Motohiko Murakami,
    Bill McDonough, Quentin Williams, Wendy Panero, and Wolfgang Bangerth.

"""
from __future__ import absolute_import

from .version import version as __version__

# classes for representing rocks and minerals:
from .mineral import Mineral
from .material import Material
from .perplex import PerplexMaterial
from .composite import Composite
from .solutionmodel import SolutionModel
from .solidsolution import SolidSolution
from .combinedmineral import CombinedMineral
from .mineral_helpers import *

# high level functions
from .main import *
from .model import Model

# mineral library
from . import minerals

# central user tools
from . import seismic
from . import output_seismo
from . import averaging_schemes
from . import eos

from . import processchemistry
from . import chemicalpotentials
from . import geotherm

# miscellaneous
from . import tools
from . import nonlinear_fitting
from .partitioning import calculate_partition_coefficient, calculate_phase_percents
