.. _ref-anisotropic-materials-hyperelastic-mineral:

Hyperelastic Mineral Class
--------------------------

The ``HyperelasticMineral`` class extends the ``AnisotropicMineral`` class from hydrostatic
states to finite, non-hydrostatic elastic deformations. Its independent
variables are the cell vectors and temperature. Derivatives of the 
Helmholtz free energy with respect to those variables
determines the stress, entropy, and
other thermodynamic and thermoelastic properties.
The underlying AnisotropicMineral supplies both the reference cell
and the hydrostatic properties at each volume and temperature.

The following description of the
``HyperelasticMineral`` class uses fixed composition, tension-positive Cauchy
stress :math:`\boldsymbol{\sigma}`, and compression-positive pressure.
Helmholtz energy :math:`\mathcal{F}`, entropy :math:`S` and heat capacities
are **molar** quantities; :math:`V` is the molar volume.
Stresses and stiffnesses have units of Pa. Repeated Cartesian indices are
summed, and :math:`:` denotes contraction over two indices.
Superscripts :math:`T` and :math:`S` on elastic tensors
mean isothermal and isentropic, respectively. Fourth-order tensors use
blackboard-bold symbols, such as :math:`\mathbb{B}`, :math:`\mathbb{C}` and :math:`\mathbb{S}`.

Cell geometry and the hydrostatic reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\mathbf{M}` have the three molar cell vectors as its **columns**.
BurnMan stores these as **rows** in the ``cell_vectors`` attribute, so that
``M = mineral.cell_vectors.T``. For a positively oriented cell,

.. math::

    V = \det\mathbf{M}

where :math:`\mathbf{M}` is the cell tensor.
At constant volume, only :math:`\det\mathbf{M}` is fixed; the individual cell lengths and angles are able to vary. When the state
of the ``HyperelasticMineral`` is set with ``set_state_with_molar_cell_vectors``, the
cell vectors are stored in ``cell_vectors`` and the molar volume is computed. The 
underlying ``AnisotropicMineral`` is set to the same volume and temperature.
The unrotated cell tensor of the ``AnisotropicMineral`` 
(:math:`\mathbf{M}^{\mathrm h}(V,T)`) and its scalar pressure
(:math:`P^{\mathrm h}(V,T)`) are computed. The deformation of the ``HyperelasticMineral`` away from this hydrostatic cell can then be described by the **isochoric deformation gradient** :math:`\mathbf{F}_{\mathrm{iso}}` and the **right Hencky strain** :math:`\mathbf{H}`:

.. math::

    \mathbf{F}_{\mathrm{iso}} = \mathbf{M}(\mathbf{M}^{\mathrm h})^{-1}
      = \mathbf{R}\mathbf{U}, \qquad
    \mathbf{H} = \ln\mathbf{U}
      = \tfrac12\ln(\mathbf{F}_{\mathrm{iso}}^{\mathsf T}
                         \mathbf{F}_{\mathrm{iso}}).

Here :math:`\mathbf{R}` is a proper rotation, :math:`\mathbf{U}` is the
symmetric positive-definite right stretch, and the logarithm is a matrix
logarithm. Because both cells have the same volume,

.. math::

    \det\mathbf{F}_{\mathrm{iso}}=1, \qquad
    \operatorname{tr}\mathbf{H}=0, \qquad
    \mathbf{M}=\mathbf{R}\exp(\mathbf{H})\mathbf{M}^{\mathrm h}.

The **right, isochoric Hencky strain** can be returned by the ``hencky_strain`` property. The total deformation gradient :math:`\mathbf{F}=\mathbf{M}\mathbf{M}_0^{-1}` is returned by the ``deformation_gradient_tensor`` property.
Finally, the rotation matrix :math:`\mathbf{R}` is returned by the ``rotation_matrix`` property if the ``set_output_frame("current")`` method has been called.

Helmholtz equation of state
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The constitutive equation for the non-hydrostatic part of the Helmholtz free energy is chosen to be quadratic in the isochoric Hencky strain at fixed volume and temperature:

.. math::
    :label: hyperelastic-helmholtz

    \mathcal{F}(\mathbf{M},T) = \mathcal{F}^{\mathrm h}(V,T)
                              + \mathcal{F}^{\mathrm{el}}(\mathbf{M},T),
    \qquad
    \mathcal{F}^{\mathrm{el}} = \frac{V}{2}
      H_{ij} \mathbb{B}^{T,\mathrm h}_{ijkl}(V,T) H_{kl}.

:math:`\mathbb{B}^{T,\mathrm h}` is the hydrostatic mineral's **spatial
isothermal stress--strain stiffness**, expressed in the same unrotated axes
as :math:`\mathbf{H}`. The hydrostatic cell and stiffness come from the
anisotropic equation of state :cite:`Myhill2022`. At
:math:`\mathbf{H}=\mathbf{0}`, all thermo-physical properties reduce to the
hydrostatic values, with :math:`\boldsymbol{\sigma}=-P^{\mathrm h}\mathbf{I}`. The molar volume :math:`V` converts the energy density into an energy per mole.

Why the hydrostatic spatial stiffness appears
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The appearance of the spatial stiffness :math:`\mathbb{B}^{T,\mathrm h}` in the hyperelastic model
may appear surprising; normally, one would
expand the Helmholtz energy in terms of Green--Lagrange strain :math:`\mathbf{E}`
(:math:`\mathbf{E}=(\mathbf{F}^{\mathsf T}\mathbf{F}-\mathbf{I})/2`), and therefore use the material stiffness :math:`\mathbb{C}^T`.
This form actually arises directly from
an expansion in terms of :math:`\mathbf{E}`, as we now show. The following derivation 
can be deduced from the description in :cite:`Wallace1967`.

Choose one hydrostatic state as a reference, with pressure :math:`P^{T,\mathrm h}`
and molar volume :math:`V_0`. Write :math:`\mathbb{C}^{T,\mathrm h}` for the Hessian
of :math:`\Psi=\mathcal{F}/V_0` with respect to Green--Lagrange strain :math:`\mathbf{E}`.
At fixed temperature, the linear coefficient is the reference second
Piola--Kirchhoff stress, equal to :math:`-P^{T,\mathrm h}\mathbf{I}`. Therefore,

.. math::

    \Psi=\Psi^{T,\mathrm h}-P^{T,\mathrm h}\operatorname{tr}\mathbf{E}
       +\tfrac12\mathbf{E}:\mathbb{C}^{T,\mathrm h}:\mathbf{E}
       +O(\|\mathbf{E}\|^3).

The Green--Lagrange strain can itself be expressed in terms of the
isochoric Hencky strain :math:`\mathbf{H}`:

.. math::

    \mathbf{E}=\tfrac12(\exp(2\mathbf{H})-\mathbf{I})
       =\mathbf{H}+\mathbf{H}^2+O(\|\mathbf{H}\|^3).

Even when :math:`\mathbf{H}` is traceless, :math:`\mathbf{E}` generally is not.
Substituting this expansion and retaining terms through second order gives

.. math::

    \Psi=\Psi^{T,\mathrm h}-P^{T,\mathrm h}\operatorname{tr}\mathbf{H}
       -P^{T,\mathrm h}\operatorname{tr}(\mathbf{H}^2)
       +\tfrac12\mathbf{H}:\mathbb{C}^{T,\mathrm h}:\mathbf{H}
       +O(\|\mathbf{H}\|^3).

The hydrostatic stress correction is :cite:`Wallace1967` (Eq. 2.36):

.. math::

    \mathbb{B}^{T,\mathrm h}_{ijkl}
      =\mathbb{C}^{T,\mathrm h}_{ijkl}
       -P^{T,\mathrm h}(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}
                        -\delta_{ij}\delta_{kl}).

Substitution, using symmetric :math:`\mathbf{H}`, gives

.. math::

    \Psi=\Psi^{T,\mathrm h}-P^{T,\mathrm h}\left[\operatorname{tr}\mathbf{H}
                     +\tfrac12(\operatorname{tr}\mathbf{H})^2\right]
          +\tfrac12\mathbf{H}:\mathbb{B}^{T,\mathrm h}:\mathbf{H}
          +O(\|\mathbf{H}\|^3).

For an isochoric deformation the trace terms vanish.
Multiplying by :math:`V` then yields the quadratic elastic term in
:eq:`hyperelastic-helmholtz`. BurnMan adopts this quadratic dependence at
each hydrostatic :math:`V,T`. Neglecting higher powers of the isochoric
Hencky strain is a **constitutive assumption**. It does not imply that
the Green--Lagrange energy Hessian remains constant as the 
isochoric strain increases, either in the model, or in the real system.

Stress, entropy and numerical derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a symmetric
increment :math:`\boldsymbol{\varepsilon}`, define

.. math::

    \mathbf{M}(\boldsymbol{\varepsilon})
       =\exp(\boldsymbol{\varepsilon})\mathbf{M}, \qquad
    \widehat{\mathcal{F}}(\boldsymbol{\varepsilon},T)
       =\mathcal{F}(\mathbf{M}(\boldsymbol{\varepsilon}),T).

Derivatives are evaluated for the deformed cell at :math:`\boldsymbol{\varepsilon}=\mathbf{0}`. They are not derivatives with respect to :math:`\mathbf{H}`. 
Mechanical work and the entropy definition give

.. math::

    d\mathcal{F}=-S\,dT+V\boldsymbol{\sigma}:d\boldsymbol{\varepsilon},
    \qquad
    \sigma_{ij}=\frac{1}{V}\widehat{\mathcal{F}}_{,\varepsilon_{ij}},
    \qquad
    S=-\left(\frac{\partial\mathcal{F}}{\partial T}\right)_{\mathbf{M}}.

Holding the physical cell fixed is essential. A temperature change at fixed
cell generally changes :math:`\mathbf{M}^{\mathrm h}`, and hence
:math:`\mathbf{H}`. A strain perturbation can change :math:`V`,
:math:`\mathbf{M}^{\mathrm h}` and :math:`\mathbb{B}^{T,\mathrm h}`.
For every perturbed state BurnMan
recomputes the hydrostatic reference and the polar decomposition; it does
not add strain increments directly to :math:`\mathbf{H}`, since
strains generally do not commute.

For example, for a symmetric direction :math:`\mathbf{D}`, let
:math:`\mathbf{M}_{\pm}=\exp(\pm h\mathbf{D})\mathbf{M}`. The central
differences used for the elastic contribution have the form

.. math::

    \boldsymbol{\sigma}:\mathbf{D}
      \simeq -P^{\mathrm h}\operatorname{tr}\mathbf{D}
       +\frac{\mathcal{F}^{\mathrm{el}}(\mathbf{M}_+,T)
                    -\mathcal{F}^{\mathrm{el}}(\mathbf{M}_-,T)}{2hV},

    S\simeq S^{\mathrm h}
       -\frac{\mathcal{F}^{\mathrm{el}}(\mathbf{M},T+\Delta T)
                    -\mathcal{F}^{\mathrm{el}}(\mathbf{M},T-\Delta T)}{2\Delta T}.

The volume in the stress denominator belongs to the unperturbed state;
the volumes inside :math:`\mathcal{F}^{\mathrm{el}}` belong to their respective
perturbed states. The hydrostatic contributions are obtained from the
underlying equation of state rather than from differences of the full
Helmholtz energy.

The reported ``pressure`` is defined as the negative one-third of the
trace of the Cauchy stress,
:math:`P=-\operatorname{tr}\boldsymbol{\sigma}/3`. It need not equal
:math:`P^{\mathrm h}(V,T)` away from hydrostatic conditions.

Elastic stiffnesses
~~~~~~~~~~~~~~~~~~~

BurnMan's ``isothermal_stiffness_tensor`` is Wallace's
:math:`\mathbb{B}^T`, the local response of the Cauchy stress to symmetric spatial
strain at fixed temperature and zero incremental rotation
:cite:`Wallace1967` (Eqs. 2.35--2.39):

.. math::

    d\sigma_{ij}=\mathbb{B}^T_{ijkl}\,d\varepsilon_{kl}, \qquad
    \mathbb{B}^T_{ijkl}=\frac{1}{V}
         \widehat{\mathcal{F}}_{,\varepsilon_{ij}\varepsilon_{kl}}
                          -\sigma_{ij}\delta_{kl}.

The first term on the RHS of the second equation is the Hessian of the **local exponential strain
parameterization** defined in the previous section. The last term accounts for the changing volume
in the definition of Cauchy stress. It follows that

.. math::

    \mathbb{B}^T_{ijkl}-\mathbb{B}^T_{klij}
       =\sigma_{kl}\delta_{ij}-\sigma_{ij}\delta_{kl}.

Thus :math:`\mathbb{B}^T` has both minor symmetries but generally lacks
major symmetry under non-hydrostatic stress. Symmetrizing its Voigt matrix
would change the physical stress response,
and so this is not done in BurnMan.

Relation to material elasticity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a fixed original reference, let :math:`\mathbf{F}` be the total
deformation gradient, :math:`J=\det\mathbf{F}`, and
:math:`\Psi=\mathcal{F}/V_0`. The work-conjugate stresses are

.. math::

    \mathbf{P}=\frac{\partial\Psi}{\partial\mathbf{F}}
       =J\boldsymbol{\sigma}\mathbf{F}^{-\mathsf T}, \qquad
    \mathbf{S}_{\mathrm{PK2}}=\frac{\partial\Psi}{\partial\mathbf{E}}
       =J\mathbf{F}^{-1}\boldsymbol{\sigma}\mathbf{F}^{-\mathsf T}.

The material Hessian
:math:`\mathbb{C}^T_{IJKL}=\partial^2\Psi/\partial E_{IJ}\partial E_{KL}`
must first be pushed forward into the current configuration:

.. math::

    \mathbb{C}^{E,T}_{ijkl}=\frac{1}{J}
       F_{iI}F_{jJ}F_{kK}F_{lL}\mathbb{C}^T_{IJKL}.

Wallace's stress correction then gives

.. math::
    :label: hyperelastic-wallace-b

    \mathbb{B}^T_{ijkl}=\mathbb{C}^{E,T}_{ijkl}
       +\tfrac12\left(
          \sigma_{ik}\delta_{jl}+\sigma_{il}\delta_{jk}
         +\sigma_{jk}\delta_{il}+\sigma_{jl}\delta_{ik}
         -2\sigma_{ij}\delta_{kl}\right).

When the current configuration is chosen as the reference,
:math:`\mathbf{F}=\mathbf{I}` and :math:`\mathbb{C}^{E,T}=\mathbb{C}^T` at that
state. This choice does not require the total prior deformation to be small.
For a general incremental displacement gradient
:math:`\mathbf{L}=\mathbf{D}+\mathbf{W}`, with symmetric :math:`\mathbf{D}`
and antisymmetric :math:`\mathbf{W}`, the stress increment also includes
rotation:

.. math::

    d\boldsymbol{\sigma}
      =\mathbb{B}^T:\mathbf{D}
        +\mathbf{W}\boldsymbol{\sigma}-\boldsymbol{\sigma}\mathbf{W}
        +\boldsymbol{\pi}\,dT.

Seismic wave propagation
^^^^^^^^^^^^^^^^^^^^^^^^

For small, locally plane, adiabatic waves about a homogeneous prestressed
state, Wallace's equation-of-motion coefficients are
:cite:`Wallace1967` (Eqs. 2.24 and 2.28)

.. math::

    \mathbb{A}^S_{ijkl}=\mathbb{C}^{E,S}_{ijkl}+\delta_{ik}\sigma_{jl}, \qquad
    T_{ik}=\mathbb{A}^S_{ijkl}n_jn_l, \qquad
    T_{ik}u_k=\rho v^2u_i.

Wallace denotes :math:`\mathbb{A}^S` by :math:`\mathbb{S}_{ijkl}`
(the symbol S is overused in thermodynamics; this tensor should
not be confused with compliance tensors, entropy, or the second
Piola--Kirchhoff stress). In this expression,
:math:`\mathbf{n}` is a unit propagation direction in the current frame,
:math:`\rho` is the current mass density, and :math:`\mathbf{u}` is the
polarization. :math:`\mathbb{C}^{E,S}` is the pushed-forward Green-strain Hessian
of internal energy at fixed entropy.

Eliminating :math:`\mathbb{C}^{E,S}` with :eq:`hyperelastic-wallace-b`, using its
isentropic counterpart, gives the form used by the ``christoffel_tensor`` method:

.. math::

    T^{B}_{ik}=\mathbb{B}^S_{ijkl}n_jn_l, \qquad
    T_{ik}=\tfrac12(T^{B}_{ik}+T^{B}_{ki})
        +\tfrac12\left[(\mathbf{n}\cdot\boldsymbol{\sigma}\mathbf{n})
                                  \delta_{ik}-\sigma_{ik}\right].

The Christoffel tensor :math:`\mathbf{T}` can then be used to compute wave speeds
and polarizations with the ``wave_velocities`` method:

.. math::

    \mathbf{T}\mathbf{u}^{(a)}=\lambda_a\mathbf{u}^{(a)},
    \qquad
    v_a=\sqrt{\frac{\lambda_a}{\rho}},
    \qquad
    \|\mathbf{u}^{(a)}\|=1,
    \qquad a=1,2,3.

Each eigenvalue :math:`\lambda_a` gives a wave speed :math:`v_a`, and its
normalized eigenvector :math:`\mathbf{u}^{(a)}` gives the polarization.

Under hydrostatic stress the additional term vanishes and :math:`\mathbf{T}^B` is symmetric, and so either :math:`\mathbf{T}^B` or :math:`\mathbf{T}` can be used. Under non-hydrostatic stress, :math:`\mathbf{T}` must be used, as done in the ``wave_velocities`` method.

Maxwell relations and thermal response
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fixed-cell thermal stress and molar isometric heat capacity are

.. math::

    \pi_{ij}=\left(\frac{\partial\sigma_{ij}}{\partial T}\right)_{\mathbf{M}}
       =\frac{1}{V}\widehat{\mathcal{F}}_{,\varepsilon_{ij}T}
       =-\frac{1}{V}\left(\frac{\partial S}
                                      {\partial\varepsilon_{ij}}\right)_T,
    \qquad
    C_\varepsilon=T\left(\frac{\partial S}{\partial T}\right)_{\mathbf{M}}.

The Maxwell relation leads to the entropy--strain derivative:

.. math::

    dS=\frac{C_\varepsilon}{T}\,dT
                      -V\boldsymbol{\pi}:d\boldsymbol{\varepsilon}.

At constant entropy, substituting the resulting temperature increment into
the stress increment gives the isentropic stiffness,

.. math::

    \mathbb{B}^S_{ijkl}=\mathbb{B}^T_{ijkl}+\frac{VT}{C_\varepsilon}\pi_{ij}\pi_{kl}.

For zero incremental rotation, the thermal expansivity tensor at fixed
Cauchy stress and the corresponding molar heat capacity are

.. math::

    \mathbb{S}^T=(\mathbb{B}^T)^{-1}, \qquad
    \boldsymbol{\alpha}=-\mathbb{S}^T:\boldsymbol{\pi}, \qquad
    C_\sigma=C_\varepsilon-VT\boldsymbol{\pi}:\boldsymbol{\alpha}.

Here the inverse acts on symmetric tensors. Because :math:`\mathbb{S}^T`
need not have major symmetry, define separately
:math:`\widetilde\alpha_{kl}=-\pi_{ij}\mathbb{S}^T_{ijkl}`. Then

.. math::

    \mathbb{S}^S_{ijkl}=\mathbb{S}^T_{ijkl}
            -\frac{VT}{C_\sigma}\alpha_{ij}\widetilde\alpha_{kl}.

Replacing :math:`\widetilde{\boldsymbol{\alpha}}` by
:math:`\boldsymbol{\alpha}` is valid only when the two contractions agree,
as they do for a major-symmetric compliance. The volumetric expansivity is
:math:`\alpha_V=\operatorname{tr}\boldsymbol{\alpha}`; the local Reuss
compressibility is :math:`\beta_R^T=\delta_{ij}\mathbb{S}^T_{ijkl}\delta_{kl}` and
:math:`K_R^T=1/\beta_R^T`.

The relevant properties are ``thermal_stress``, ``thermal_expansivity_tensor``,
``molar_isometric_heat_capacity``, ``molar_isostress_heat_capacity`` and the
isothermal/isentropic stiffness and compliance tensors. Isometric means
fixed cell shape **and** volume. It is a stronger constraint than fixing
volume while allowing the hydrostatic reference cell to change shape.


Rotations and output frames
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rigid rotation changes tensor components but leaves the energy, volume,
entropy and wave speeds for correspondingly rotated directions unchanged.
For a rotation :math:`\mathbf{Q}`,

.. math::

    \mathbf{M}'=\mathbf{Q}\mathbf{M}, \qquad
    \boldsymbol{\sigma}'=\mathbf{Q}\boldsymbol{\sigma}\mathbf{Q}^{\mathsf T},
    \qquad
    \mathbb{B}'_{ijkl}=Q_{ia}Q_{jb}Q_{kc}Q_{ld}\mathbb{B}_{abcd}.

``set_output_frame("current")`` reports spatial tensors in the frame of the
supplied cell vectors; this is the default. ``set_output_frame("crystal")``
reports them in the crystallographic frame defined by the cell and the
underlying mineral's frame convention. This changes the representation,
not the physical state. Input to ``set_state_with_molar_cell_vectors`` is
always in the current frame, regardless of the selected output frame.
The right ``hencky_strain`` stays in the hydrostatic reference axes and
is not a spatial tensor to be rotated with the output frame.

Using the model
~~~~~~~~~~~~~~~

Given an existing ``AnisotropicMineral`` named ``hydrostatic_mineral``,
the following applies a volume-preserving simple shear at 1000 K:

.. code-block:: python

    import numpy as np
    from burnman import HyperelasticMineral

    hydrostatic_mineral.set_state(5.0e9, 1000.0)
    M_h = hydrostatic_mineral.unrotated_cell_vectors.T.copy()
    mineral = HyperelasticMineral(hydrostatic_mineral)

    F = np.eye(3)
    F[0, 1] = 0.05
    mineral.set_state_with_molar_cell_vectors((F @ M_h).T, 1000.0)

    print("Molar volume (m^3/mol):", mineral.V)
    print("Mean pressure (GPa):", mineral.pressure / 1.0e9)
    print("Cauchy stress (GPa):\n", mineral.cauchy_stress / 1.0e9)
    velocities, polarizations = mineral.wave_velocities([1.0, 0.0, 0.0])
    print("Wave speeds (km/s):", velocities / 1.0e3)

The supplied hydrostatic object is retained and its state is updated during
hyperelastic calculations. Use separate hydrostatic instances for independent
hyperelastic models. The repository example
``examples/example_hyperelastic_mineral.py`` constructs a complete model,
solves shear at fixed normal traction, and plots stresses, relative volume
change and cell shape.

``burnman.tools.eos.check_hyperelastic_eos_consistency`` checks energy
derivatives, Maxwell relations, thermal conversions and rotation behaviour.
Numerical derivatives have finite step-size and cancellation errors; exact
thermodynamic identities should therefore be tested with numerical tolerances.

Class reference
~~~~~~~~~~~~~~~

.. autoclass:: burnman.classes.hyperelastic.HyperelasticMineral
    :no-members:
    :no-undoc-members:
    :no-inherited-members:
    :no-index:
