Examples
********

BurnMan comes with a small tutorial in the tutorial/ folder, and large
collection of example programs under examples/. Below you can find a summary
of the different examples. They are grouped into :ref:`ref-example-tutorial`,
:ref:`ref-example-simple`, and :ref:`ref-example-advanced`. We suggest
starting with the tutorial before moving on to the simpler examples,
especially if you are new to using BurnMan.

Finally, we also include the scripts that were used for all computations and
figures in the 2014 BurnMan paper in the misc/ folder, see
:ref:`ref-example-paper`.

.. _ref-example-tutorial:

Tutorial
========

The tutorial for BurnMan currently consists of three separate units:
  - :mod:`step 1 <tutorial.step_1>`,
  - :mod:`step 2 <tutorial.step_2>`, and
  - :mod:`step 3 <tutorial.step_3>`.

.. _ref-example-tut1:

.. #CIDER 2014 BurnMan Tutorial --- step 1
.. #--------------------------------------

.. automodule:: tutorial.step_1

When run (without putting in a more realistic composition), the program produces the following image:

.. image:: tut-step1.png
Your goal in this tutorial is to improve this awful fit....


link to source code: `tutorial/step_2.py <../../../tutorial/step_2.py>`_




.. automodule:: tutorial.step_2

Whithout changing any input, the program should produce the following image showing the misfit as a function of perovskite content:

.. image:: tut-step2.png

link to source code: `tutorial/step_2.py <../../../tutorial/step_2.py>`_



.. automodule:: tutorial.step_3

After changing the standard deviations for :math:`K\prime` and :math:`G\prime` to 0.2, the following figure of velocities for 1000 realizations is produced:

.. image:: tut-step3.png


.. _ref-example-simple:

Simple Examples
===============

The following is a list of simple examples:
  - :mod:`~examples.example_beginner`,
  - :mod:`~examples.example_geotherms`,
  - :mod:`~examples.example_seismic`,
  - :mod:`~examples.example_composition`, and
  - :mod:`~examples.example_averaging`.


.. automodule:: examples.example_beginner

*Resulting figure:*

.. image:: example_beginner.png



.. automodule:: examples.example_geotherms

*Resulting figure:*

.. image:: example_geotherm.png

.. automodule:: examples.example_seismic

*Resulting figures:*

.. image:: example_seismic.png

.. image:: example_seismic2.png

.. automodule:: examples.example_composition  

*Resulting figure:*

.. image:: example_composition.png

.. automodule:: examples.example_averaging

*Resulting figure:*

.. image:: example_averaging.png


.. _ref-example-advanced:

More Advanced Examples
======================

Advanced examples:
  - :mod:`~examples.example_spintransition`,
  - :mod:`~examples.example_user_input_material`,
  - :mod:`~examples.example_optimize_pv`, and
  - :mod:`~examples.example_compare_all_methods`.

.. automodule:: examples.example_spintransition


*Resulting figure:*

.. image:: example_spintransition.png

.. automodule:: examples.example_user_input_material

.. automodule:: examples.example_optimize_pv

*Resulting figure:*

.. image:: example_opt_pv.png


.. automodule:: examples.example_compare_all_methods  


*Resulting figure:*

.. image:: example_compare_all_methods.png


.. _ref-example-paper:

Reproducing Cottaar, Heister, Rose and Unterborn (2014)
=======================================================

.. automodule:: misc.paper_averaging

.. automodule:: misc.paper_benchmark
   :exclude-members: check_slb_fig7_txt

.. automodule:: misc.paper_fit_data

.. automodule:: misc.paper_incorrect_averaging
   :exclude-members: ferropericlase,perovskite

.. automodule:: misc.paper_opt_pv
   :exclude-members: calc_shear_velocities,error

.. automodule:: misc.paper_onefit
   :exclude-members: array_to_rock, make_rock, output_rock, realization_to_array

.. automodule:: misc.paper_uncertain
   :exclude-members: my_perovskite

Misc or work in progress
========================


.. automodule:: examples.example_inv_murakami
.. broken: .. automodule:: example_optimize_slb2011
.. automodule:: examples.example_compare_enstpyro     
.. automodule:: examples.example_partition_coef

.. broken: .. automodule:: example_premite_isothermal

.. automodule:: examples.example_fit_data             
.. automodule:: examples.example_grid                 
.. automodule:: examples.example_woutput

