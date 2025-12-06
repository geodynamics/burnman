.. _ref-relaxed_materials:

Relaxed Material Classes
========================

Overview
--------

Some materials have isochemical degrees of freedom, meaning that some macroscopic
characteristic of the material can change without changing its overall composition.
For example, a solid solution may be able to change how Fe and Mg are distributed
amongst its constituent sites, or Fe may be able to change its spin state, or the material
itself may be able to change by distortion of the crystal lattice. In a multiphase
material, isochemical degrees of freedom may arise from changes in the proportions or
compositions of the constituent phases (e.g., Fe-Mg exchange). If these isochemical
variables can change on timescales shorter than the timescales of changes in pressure
and temperature, then the material is said to be relaxed with respect to these
isochemical degrees of freedom.

BurnMan provides classes to model relaxed materials, including
:class:`burnman.RelaxedSolution`, :class:`burnman.RelaxedAnisotropicSolution`,
and :class:`burnman.RelaxedComposite`. These classes extend their unrelaxed
counterparts by adding the capability to minimize the Gibbs energy of the system
by varying one or more isochemical variables at fixed pressure and temperature.
These modified classes are not simply a combination of the unrelaxed classes
with a Gibbs minimiser; they also modify the way thermodynamic properties
are calculated to account for the relaxation of the isochemical variables.

The properties of particular interest in relaxed materials are those
that depend on the second derivatives of the Gibbs energy, such as
the bulk modulus, thermal expansivity, and heat capacity. In unrelaxed
materials, these properties are calculated directly from the second derivatives
of the Gibbs energy with respect to pressure and temperature, keeping all
molar fractions constant. By contrast, in relaxed materials, the
second derivatives are affected by simultaneous changes in the
isochemical variables. To account for this, BurnMan uses variational calculus
to derive expressions for the relaxed second derivatives of the Gibbs energy,
which are then used to compute the thermodynamic properties of the relaxed material.

Mathematical Derivation
-----------------------

Let the extent of the isochemical reactions be given by the vector :math:`\xi`.
Let us also define the partial derivative function
:math:`F(P,T,X,\xi) = \partial G(P,T,X,\xi) / \partial \xi`,
where :math:`G` is the Gibbs free energy of the system.
For simplicity, let us group the pressure, temperature, and
composition variables into a single vector :math:`Z = (P,T,X)`
(we aren't actually interested in the composition, but keep it here for completeness). 
The relaxed state is given by the condition :math:`F(Z,\xi) = 0`.
Finally, we define a function :math:`\phi(Z) = \xi`,
which gives the extent of the isochemical reactions at equilibrium. We now have:

.. math::
    F = F(Z, \phi(Z)) = \frac{\partial G}{\partial \xi} (Z, \phi(Z))

Taking the total derivative of this expression with respect to :math:`Z`, we have:

.. math::
    \frac{d F}{d Z} = \frac{\partial F}{\partial Z} + \frac{\partial F}{\partial \xi}
    \frac{\partial \phi}{\partial Z}

Now, since :math:`F(Z, \phi(Z)) = 0` for all :math:`Z`, we also have :math:`d F / d Z = 0`, 
and therefore:

.. math::
     \frac{\partial F}{\partial \xi}
    \frac{\partial \phi}{\partial Z} = - \frac{\partial F}{\partial Z}
    
Note that the matrix :math:`\partial F / \partial \xi` is the Hessian matrix of second derivatives
of the Gibbs energy with respect to the isochemical variables:

.. math::
    \frac{\partial F}{\partial \xi} = \frac{\partial^2 G}{\partial \xi^2}

It is therefore symmetric. Let us define the pseudo-inverse of this matrix as
:math:`R`, so that :math:`R \partial F / \partial \xi = I`, where :math:`I` is the identity matrix.
Then:

.. math::
    \frac{\partial \phi}{\partial Z} = - R \frac{\partial F}{\partial Z}


Finally, define a function :math:`g(Z) = G(Z, \phi(Z))`, which gives the Gibbs energy
of the system at equilibrium. Using the chain rule, we can write the total derivative of this function with respect to :math:`Z` as:

.. math::
    \frac{d g}{d Z} = \frac{\partial G}{\partial Z} + \frac{\partial G}{\partial \xi}
    \frac{\partial \phi}{\partial Z}

Since :math:`F = \partial G / \partial \xi = 0`, the second term on the right-hand side
vanishes, and we have:

.. math::
    \frac{d g}{d Z} = \frac{\partial G(Z, \phi(Z))}{\partial Z}

Taking the second derivative of :math:`g` with respect to :math:`Z`:

.. math::
    \frac{d^2 g}{d Z^2} = \frac{\partial^2 G}{\partial Z^2} +
    \frac{\partial^2 G}{\partial Z \partial \xi} \frac{\partial \phi}{\partial Z}


Now, substituting in our expression for :math:`\partial \phi / \partial Z`:

.. math::
    \frac{d^2 g}{d Z^2} = \frac{\partial^2 G}{\partial Z^2} -
    \frac{\partial^2 G}{\partial Z \partial \xi} R
    \frac{\partial F}{\partial Z}

or, equivalently:

.. math::
    \frac{d^2 g}{d Z^2} = \frac{\partial^2 G}{\partial Z^2} -
    \frac{\partial^2 G}{\partial Z \partial \xi} 
    R
    \frac{\partial^2 G}{\partial \xi \partial Z}

This expression gives the relaxed second derivatives of the Gibbs energy
with respect to pressure, temperature, and composition. These second derivatives
can then be used to calculate the relaxed thermodynamic properties of the system,
such as the bulk modulus, thermal expansivity, and heat capacity:

.. math::
    K_{P,relaxed} = -V \left(\frac{\partial^2 g}{\partial P^2}\right)^{-1}, \quad
    \alpha_{relaxed} = \frac{1}{V} \frac{\partial^2 g}{\partial P \partial T}, \quad
    C_{P,relaxed} = - T \frac{\partial^2 g}{\partial T^2}


Relaxed Material Classes
------------------------

.. toctree::
  :maxdepth: 2

  relaxed_materials_01_relaxed_solution
  relaxed_materials_02_relaxed_anisotropic_solution
  relaxed_materials_03_relaxed_composite
