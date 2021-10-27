

.. _ref-example-tutorial:

The BurnMan Tutorial
********************

.. toctree::
  :maxdepth: 1

  tutorial_01_material_classes.ipynb
  tutorial_02_composition_class.ipynb
  tutorial_03_layers_and_planets.ipynb
  tutorial_04_fitting.ipynb
  tutorial_05_equilibrium.ipynb

CIDER Tutorial 2014
*******************

The tutorial for BurnMan presented at CIDER 2014 consists of three separate units:
  - :mod:`step 1 <contrib.cider_tutorial_2014.step_1>`,
  - :mod:`step 2 <contrib.cider_tutorial_2014.step_2>`, and
  - :mod:`step 3 <contrib.cider_tutorial_2014.step_3>`.

.. _ref-example-tut1:

.. automodule:: contrib.cider_tutorial_2014.step_1

When run (without putting in a more realistic composition), the program produces the following image:

.. image:: figures/tut-step1.png

Your goal in this tutorial is to improve this awful fit...



.. automodule:: contrib.cider_tutorial_2014.step_2

Without changing any input, the program should produce the following image showing the misfit as a function of perovskite content:

.. image:: figures/tut-step2.png



.. automodule:: contrib.cider_tutorial_2014.step_3

After changing the standard deviations for :math:`K_{0}^{'}` and :math:`G_{0}^{'}` to 0.2, the following figure of velocities for 1000 realizations is produced:

.. image:: figures/tut-step3.png
