# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.
# Ian Rose ian.r.rose@gmail.com

"""

CIDER 2014 BurnMan Tutorial --- step 2
--------------------------------------

In this second part of the tutorial we try to get a closer fit to our
1D seismic reference model.  In the simple Mg, Si, and O model that
we used in step 1 there was one free parameter, namely phase_1_fraction,
which goes between zero and one.

In this script we want to explore how good of a fit to PREM we can get
by varying this fraction. We create a simple function that calculates
a misfit between PREM and our mineral model as a function of phase_1_fraction,
and then plot this misfit function to try to find a best model.

This script may be run by typing

    python step_2.py

"""
from __future__ import absolute_import
from __future__ import print_function

# The imports here are identical to those from step 1
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
import os
import sys
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman import minerals

if __name__ == '__main__':
    # Again, we load the PREM seismic model and query it for pressure,
    # density, and elastic properties at lower mantle depths. This, too,
    # is identical to step 1
    n_depths = 20
    min_depth = 850.e3
    max_depth = 2800.e3
    depths = np.linspace(min_depth, max_depth, n_depths)

    seismic_model = burnman.seismic.PREM()
    pressure, seis_rho, seis_vphi, seis_vs = seismic_model.evaluate(
        ['pressure', 'density', 'v_phi', 'v_s'], depths)
    temperature = burnman.geotherm.brown_shankland(pressure)

    """
    This is the main workhorse of this script.  We define the function ``misfit''
    which takes a single parameter phase_1_fraction.  We create the rock we want,
    calculate its elastic properties, and then calculate an L2 misfit between
    the calculated profiles for Vs, Vphi, and density and those for PREM.

    Here again is our model with stishovite and wuestite.  Instead of that, you
    will want to copy the ``rock'' you created in step 1.  If you experimented with
    other rocks than those with Mg perovskite and periclase, you can also try those.
    """

    def misfit(phase_1_fraction):

        # Here we define the rock as before.
        phase_2_fraction = 1.0 - phase_1_fraction
        rock = burnman.Composite(
            [minerals.SLB_2011.stishovite(), minerals.SLB_2011.wuestite()], [phase_1_fraction, phase_2_fraction])

        # Just as in step 1, we want to set which equation of state we use,
        # then call burnman.velocities_from_rock, which evaluates the
        # elastic properties and seismic velocities at the predefined
        # pressures and temperatures
        rock.set_method('slb3')
        density, vphi, vs = rock.evaluate(
            ['density', 'v_phi', 'v_s'], pressure, temperature)

        # Since we will call this misfit function many times, we may be interested
        # in a status report.  These lines print some debug output so we
        # can keep track of what the script is doing.
        print("Calculations are done for:")
        rock.debug_print()

        # Here we integrate an L2 difference with depth between our calculated seismic
        # profiles and PREM.  We then return those misfits.
        [vs_err, vphi_err, rho_err] = burnman.compare_l2(
            depths, [vs, vphi, density], [seis_vs, seis_vphi, seis_rho])

        return vs_err, vphi_err, rho_err

    """

    With the misfit function now defined, we can call it many times for different
    phase_1_fraction values, and see how good of a fit we can get.

    """

    # We create the array ``fraction'', which has 101 fractions between
    # zero and one, and call our misfit function for each of those fractions.
    fraction = np.linspace(0.0, 1.0, 101)
    errs = np.array([misfit(f) for f in fraction])
    vs_misfit = errs[:, 0]
    vphi_misfit = errs[:, 1]
    rho_misfit = errs[:, 2]

    # Finally, we plot the misfits against the phase_1_fraction. You will probably
    # find that it is difficult to fit shear wave speed, bulk sound speed, and
    # density all at the same time.

    plt.plot(fraction, vs_misfit, "r-x", label=("Vs misfit"))
    plt.plot(fraction, vphi_misfit, "b-x", label=("Vphi misfit"))
    plt.plot(fraction, rho_misfit, "g-x", label=("Density misfit"))
    plt.yscale('log')
    plt.xlabel('Fraction Phase 1')
    plt.ylabel('Misfit')
    plt.legend()
    plt.show()
