# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Mineral database

  - :mod:`~burnman.minerals.SLB_2005`
  - :mod:`~burnman.minerals.SLB_2011_ZSB_2013`
  - :mod:`~burnman.minerals.SLB_2011`
  - :mod:`~burnman.minerals.Murakami_etal_2012`
  - :mod:`~burnman.minerals.Murakami_2013`
  - :mod:`~burnman.minerals.Matas_etal_2007`
  - :mod:`~burnman.minerals.HP_2011_ds62`
  - :mod:`~burnman.minerals.HP_2011_fluids`
  - :mod:`~burnman.minerals.HHPH_2013`
  - :mod:`~burnman.minerals.other`
"""
from __future__ import absolute_import

# Stixrude and Lithgow-Bertelloni
from . import SLB_2011
from . import SLB_2011_ZSB_2013
from . import SLB_2005

# Murakami and coworkers
from . import Murakami_etal_2012
from . import Murakami_2013

# Matas and coworkers
from . import Matas_etal_2007

# Holland, Powell and coworkers
from . import HP_2011_ds62
from . import HP_2011_fluids
from . import HHPH_2013

# Kurnosov et al. 2017
from . import KMFBZ_2017
# Other
from . import other
