.. _ref-mineral-datasets:

Mineral Datasets
================

BurnMan includes a large collection of mineral datasets in
the the :mod:`burnman.minerals` module that
provide thermodynamic and thermoelastic properties
for minerals commonly found in the Earth and planetary interiors.

See :ref:`ref-api-mineral-database` for a list of built-in
minerals and solutions available in BurnMan.

The current HPx-eos igneous sets are available under their upstream names:

* ``ig51W24``: Weller et al. (2024), for anhydrous alkaline and subalkaline
  systems. Includes all 12 solutions, including ``liq_W24d``, nepheline,
  kalsilite, leucite and melilite.
* ``ig51G25``: Green et al. (2025), the corrected H18 set for hydrous and
  anhydrous subalkaline systems. Includes all 14 solutions, including
  ``liq_G25w``, ``fl_G25``, biotite, muscovite, amphibole, epidote and cordierite.

Both use THERMOCALC endmember dataset 6.36, ``HPx_ds636``, as
distributed with the W24 and G25 releases. The dataset parameters are distinct
from those published by Holland and Powell (2011); that citation describes
the equations of state used for the solids and liquids. Both sets include
``fmq`` and ``nno`` as ``CombinedMineral`` oxygen buffers. The supplied
``tc-ds636.txt`` records generation on 3 May 2023.
Dataset 6.36 includes solids, liquids and PS1994 water;
other gases and aqueous endmembers are omitted. Its ``cov()`` function
retains the full 291-endmember formation-energy covariance matrix.
Solution class names and endmember order match the upstream files.
For example::

    from burnman.minerals import ig51W24, ig51G25

    dry_melt = ig51W24.liq_W24d()
    hydrous_melt = ig51G25.liq_G25w()
    fluid = ig51G25.fl_G25()

These modules were generated from the authors' `W24 release of 28 June 2025
<https://hpxeosandthermocalc.org/wp-content/uploads/2025/06/tc-thermoinput-igw24-2025-06-28.zip>`_
and `G25 release uploaded 7 May 2025
<https://hpxeosandthermocalc.org/wp-content/uploads/2025/05/tc-thermoinput-igg25-h18corr-2024-12-21.zip>`_.
Their references are `Weller et al. (2024)
<https://doi.org/10.1093/petrology/egae098>`_ and `Green et al. (2025)
<https://doi.org/10.1093/petrology/egae079>`_.
To regenerate them, extract the archives into
``burnman/data/input_raw_endmember_datasets`` and run from the repository root::

    python burnman/data/input_raw_endmember_datasets/HP636data_to_burnman.py \
        burnman/data/input_raw_endmember_datasets/tc-thermoinput-igW24-2025-06-28/tc-ds636.txt
    python burnman/data/input_raw_endmember_datasets/hpx_eos_to_burnman.py

The solution porter also regenerates the older HPx-eos sets, whose original
source archives must be present. It reads the THERMOCALC 3.51 headers and
explicit activity denominators, preserves variable site multiplicities,
and converts temperature-dependent interactions from kJ/mol/K to J/mol/K.
Endmember objects with identical labels and different recipes receive
distinct Python identifiers, while identical definitions are shared.

The HPx-eos datasets include the metabasite melt ``mb50NCKFMASHTO.L``,
the metapelite melts ``mp50KFMASH.liq``, ``mp50MnNCKFMASHTO.liq`` and
``mp50NCKFMASHTO.liq``, and the aqueous fluid ``ig50NCKFMASHTOCr.fl``.
The porter preserves their activity expressions, including melt sites
whose multiplicities depend on composition, and their DQF energy corrections.

The fluid uses ``HP_2011_fluids.H2O``: the Pitzer-Sterner (1994) water
equation of state with the Holland-Powell (2011, Table 2a) ideal-gas
enthalpy, entropy and heat capacity reference. Water selects the stable
liquid or vapour root with the lowest Gibbs energy. Its supported range
is 500–1700 K and 0 < P <= 5 GPa; this also limits the fluid solution.
``ig51G25.fl_G25`` uses the equivalent water reference from ``HPx_ds636``
and has the same limits.
The EOS was ported from BurnMan's ``water_filter`` branch, with root
selection and the HP thermal reference following ``burnman_cpp``.
``Pitzer_Sterner_1994.H2O_Pitzer_Sterner`` retains that branch's original
Debye thermal reference, which has different absolute energies and is
not the reference used by the HPx-eos fluid.

The original H18 ``liq`` models in ``ig50NCKFMASHTOCr`` and
``ig50NCKFMASTOCr`` are excluded because their authors withdrew them
after finding an algebraic error. Use ``ig51G25.liq_G25w`` for the corrected
replacement, together with the solutions from ``ig51G25`` and dataset 6.36.
See the `authors' igneous dataset notes
<https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/hpx-eos-igneous-sets/>`_.
