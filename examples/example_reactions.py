# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2023 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_reactions
-----------------

There are several reasons why we might want to calculate
a list of balanced endmember reactions between phases.
This is useful, for example, as an initial step in
designing conventional thermobarometers or in developing
a model for disequilibrium geodynamics.

This brief example demonstrates two BurnMan functions that
allow users to easily calculate reactions from a
list of endmember formulae, or from BurnMan Composite objects.

*Uses:*

* :func:`burnman.tools.chemistry.reactions_from_stoichiometric_matrix`
* :func:`burnman.tools.chemistry.reactions_from_formulae`

*Demonstrates:*

* Calculating balanced reactions from a set of endmembers.
"""

import numpy as np
from burnman import Composite
from burnman.minerals import SLB_2022
from burnman.tools.chemistry import reactions_from_stoichiometric_matrix
from burnman.tools.chemistry import reactions_from_formulae
from pprint import pprint

if __name__ == "__main__":
    # Let's first look at melting reactions between solid and liquid mantle
    # endmembers in the CaO-MgO-SiO2-H2O system. This system includes
    # olivine (forsterite), orthopyroxene (enstatite and orthodiopside),
    # clinopyroxene (clinoenstatite and diopside),
    # liquid (forsterite, diopside, quartz and water endmembers) and
    # fluid (H2O endmember).

    # We start by declaring a dictionary with names and formulae of endmembers.
    endmembers = {
        "foO": "Mg2SiO4",
        "enOpx": "Mg2Si2O6",
        "diOpx": "MgCaSi2O6",
        "enCpx": "Mg2Si2O6",
        "diCpx": "MgCaSi2O6",
        "foL": "Mg2SiO4",
        "diL": "MgCaSi2O6",
        "qL": "SiO2",
        "wL": "H2O",
        "wFluid": "H2O",
    }

    # Split this dictionary into names and formulae
    names = list(endmembers.keys())
    formulae = list(endmembers.values())

    # The next line does all the heavy lifting
    # and calculates all of the forward and reverse reactions
    # We can either return the reactions as a matrix...
    reactions_matrix = reactions_from_formulae(formulae, names, return_strings=False)

    # ...or as a list of easily readable strings
    reactions = reactions_from_formulae(formulae, names, return_strings=True)
    n_forward = int(len(reactions) / 2)

    print("Calculating upper mantle reactions in CMSH system...")
    print(
        f"There are {n_forward} forward reactions between the {len(names)} endmembers in this system."
    )
    print("... as matrix:")
    pprint(np.array(reactions_matrix, dtype=int)[:n_forward])

    print("... as strings:")
    for r in reactions[:n_forward]:
        print(r)
    print("")

    # Ok, let's look at another example.
    # This time, we create a Composite material using three
    # phases from the Stixrude and Lithgow-Bertelloni (2022) dataset
    assemblage = Composite(
        [SLB_2022.bridgmanite(), SLB_2022.ferropericlase(), SLB_2022.ca_perovskite()]
    )

    # The Composite material can be queried for the elements,
    # chemical formulae and names of all the endmembers in all the phases
    elements = assemblage.elements
    formulae = assemblage.endmember_formulae
    names = assemblage.endmember_names

    print("Calculating lower mantle reactions in NCFMAS system for assemblage:")
    pprint(f"{[phase.name for phase in assemblage.phases]}")
    print("Endmembers:")
    pprint(names)
    print("Elements:")
    pprint(elements)

    # The stoichiometric matrix can easily be accessed
    # The values in the matrix M[i,j] correspond to the number of atoms of
    # element j in endmember i.
    print("Stoichiometric matrix:")
    M = assemblage.stoichiometric_matrix
    pprint(np.array(M, dtype=int))

    # Now we use the second reaction calculating function from
    # the tools module. This takes the stoichiometric matrix as input
    # and returns the list of reactions.
    R = reactions_from_stoichiometric_matrix(M)
    n_forward = int(R.shape[0] / 2)
    print(
        f"There are {n_forward} forward reactions between the {len(names)} endmembers in this system."
    )
    print(R[:n_forward])
    print(reactions_from_formulae(formulae, names, return_strings=True)[:n_forward])
    print()

    # Our final example is very similar to the previous one,
    # but here we add two additional phases to the assemblage.
    # This example is designed to demonstrate that even a relatively small
    # number of endmembers can result in a lot of reactions.
    assemblage = Composite(
        [
            SLB_2022.bridgmanite(),
            SLB_2022.ferropericlase(),
            SLB_2022.ca_perovskite(),
            SLB_2022.calcium_ferrite_structured_phase(),
            SLB_2022.new_aluminous_phase(),
        ]
    )

    elements = assemblage.elements
    formulae = assemblage.endmember_formulae
    names = assemblage.endmember_names

    print("Calculating lower mantle reactions in NCFMAS system for assemblage:")
    pprint(f"{[phase.name for phase in assemblage.phases]}")

    # Here we calculate the reactions and store them as a list of readable strings.
    # However, we don't print them to standard output, but just print the number of
    # forward reactions.
    reactions = reactions_from_formulae(formulae, names, return_strings=True)
    n_forward = int(len(reactions) / 2)
    print(
        f"There are {n_forward} forward reactions between the {len(names)} endmembers in this system."
    )
