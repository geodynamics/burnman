# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

example_composite
----------------------
    
This example shows how to create composites and output
thermodynamic and thermoelastic quantities.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`


*Demonstrates:*

* Different ways to define an instance of the Composite class
* How to set phase fractions, composition and state
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
    Now we can make a simple composite assemblage
    """
    forsterite=minerals.SLB_2011.forsterite()
    enstatite=minerals.SLB_2011.enstatite()

    class simple_mantle(burnman.Composite):
        def __init__(self):
            self.name='simple_mantle'
            phases = [forsterite, enstatite]
            burnman.Composite.__init__(self, phases)

    assemblage=simple_mantle()

    assemblage.set_fractions([0.8, 0.2])
    print assemblage.composition()
    assemblage.set_state(P, T)
    print 'Assemblage 1 (endmembers only)'
    print assemblage.phase_names
    print ''

    olivine=minerals.SLB_2011.mg_fe_olivine()
    orthopyroxene=minerals.SLB_2011.orthopyroxene()

    class mantle_with_solid_solutions(burnman.Composite):
        def __init__(self):
            self.name='Mantle with solid solutions'
            phases = [olivine, orthopyroxene]
            burnman.Composite.__init__(self, phases)

    assemblage=mantle_with_solid_solutions()

    assemblage.set_fractions([0.8, 0.2])
    assemblage.set_composition([[0.9,0.1], [0.94, 0.06, 0.0, 0.0]])
    assemblage.set_state(P, T)
    print 'Assemblage 2 (with solid solutions), method 1'
    print assemblage.phase_names
    print ''

    """
    Here's a briefer way of doing the same thing
    """
    olivine=minerals.SLB_2011.mg_fe_olivine()
    orthopyroxene=minerals.SLB_2011.orthopyroxene()
    assemblage=burnman.Composite([olivine, orthopyroxene])
    assemblage.set_composition([[0.9,0.1], [0.94, 0.06, 0.0, 0.0]])
    assemblage.set_fractions([0.8, 0.2])
    assemblage.set_state(P, T)
    print 'Assemblage 2 (with solid solutions), method 2'
    print assemblage.phase_names
    print ''
    '''
    Allow fractions to be set as molar (default), volume or mass
    Volume and mass fractions only possible *after* phase compositions have been set
    Volume fractions also require state to have been set
    The fractions stored and used by burnman are *always* molar fractions - volume and mass fractions are converted before use. Note that mass fractions will change with set_composition, and volume fractions will change with set_state...
    '''
    assemblage.set_fractions([0.8, 0.2], 'molar')
    print 'Assemblage composition when molar phase fractions are', assemblage.molar_fractions
    print assemblage.composition()
    print ''

    print 'Converting between fraction types'
    assemblage.set_fractions([0.8, 0.2], 'mass')
    print 'Molar fractions from mass fractions:', assemblage.molar_fractions

    assemblage.set_fractions([0.8, 0.2], 'volume')
    print 'Molar fractions from volume fractions:', assemblage.molar_fractions
    print ''


    """
    An assemblage can also be queried for chemical potential information
    """
    fayalite=minerals.HP_2011_ds62.fa()
    magnetite=minerals.HP_2011_ds62.mt()
    quartz=minerals.HP_2011_ds62.q()
    O2=minerals.HP_2011_fluids.O2()
    O2.set_state(1.e5, T)

    FMQ=burnman.Composite([fayalite, magnetite, quartz])
    FMQ.set_state(P, T)
    print 'Chemical potentials for the FMQ buffer at', P/1.e9, "GPa and", T, "K"
    print FMQ.phase_names
    components=['FeO', 'SiO2', 'O2']

    for i, component in enumerate(components):
        print component + ':', FMQ.chemical_potentials(components)[i]/1.e3, 'kJ/mol'

    print 'log10(fO2):', np.log10(FMQ.fugacity(O2))



    """
    Now we deal with elastic and thermodynamic properties of composites
    where the individual phases are treated as isotropic. 

    Averaging is conducted over phases (endmembers and solid solutions)
    which have unique isotropic bulk and shear moduli derived in a 
    different way to composites.

    """

    
