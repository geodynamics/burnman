# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

example_solid_solution
----------------------
    
This example shows how to create different solid solution models and output
thermodynamic and thermoelastic quantities.

There are three main types of solid solution currently implemented in 
BurnMan:

1. Ideal solid solutions
2. Symmmetric solid solutions
3. Asymmetric solid solutions

These solid solutions can potentially deal with:

* Disordered endmembers (more than one element on a crystallographic site)
* Site vacancies
* More than one valence/spin state of the same element on a site

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.solidsolution.SolidSolution`
* :class:`burnman.solutionmodel.SolutionModel`


*Demonstrates:*

* Different ways to define a solid solution
* How to set composition and state
* How to output thermodynamic and thermoelastic properties

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":
    '''
    First, let's pick a starting pressure and temperature
    '''
    P=1.e5
    T=1000.

    """
    Now we can make a simple assemblage
    """
    forsterite=minerals.SLB_2011.forsterite()
    enstatite=minerals.SLB_2011.enstatite()

    class simple_mantle(burnman.Assemblage):
        def __init__(self):
            self.name='simple_mantle'
            phases = [forsterite, enstatite]
            burnman.Assemblage.__init__(self, phases)

    assemblage=simple_mantle()

    assemblage.set_phase_fractions([0.8, 0.2])
    assemblage.set_composition([[], []])
    assemblage.set_state(P, T, "all")
    print assemblage.G


    olivine=minerals.SLB_2011.mg_fe_olivine()
    orthopyroxene=minerals.SLB_2011.orthopyroxene()

    class mantle_with_solid_solutions(burnman.Assemblage):
        def __init__(self):
            self.name='Mantle with solid solutions'
            phases = [olivine, orthopyroxene]
            burnman.Assemblage.__init__(self, phases)

    assemblage=mantle_with_solid_solutions()

    assemblage.set_phase_fractions([0.8, 0.2])
    assemblage.set_composition([[0.9,0.1], [0.94, 0.06, 0.0, 0.0]])
    assemblage.set_state(P, T, "all")
    print assemblage.G

    """
    Here's a simpler way of doing the same thing
    """
    assemblage=burnman.Assemblage([olivine, orthopyroxene])
    assemblage.set_phase_fractions([0.8, 0.2])
    assemblage.set_composition([[0.9,0.1], [0.94, 0.06, 0.0, 0.0]])
    assemblage.set_state(P, T, "all")
    print assemblage.G
