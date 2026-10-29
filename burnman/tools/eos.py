# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import warnings
from scipy.linalg import logm, expm
from ..utils.anisotropy import (
    contract_compliances,
    contract_stiffnesses,
    contract_strains,
    contract_stresses,
    expand_strains,
)


def check_eos_consistency(
    m,
    P=1.0e9,
    T=300.0,
    tol=1.0e-4,
    verbose=False,
    including_shear_properties=True,
    equilibration_function=None,
):
    """
    Checks that numerical derivatives of the Gibbs energy of a mineral
    or composite under given conditions are equal to those provided
    analytically by the equation of state.

    :param m: The material for which the equation of state
        is to be checked for consistency.
    :type m: :class:`burnman.Mineral` or :class:`burnman.Composite`

    :param P: The pressure at which to check consistency.
    :type P: float

    :param T: The temperature at which to check consistency.
    :type T: float

    :param tol: The fractional tolerance for each of the checks.
    :type tol: float

    :param verbose: Decide whether to print information about each check.
    :type verbose: bool

    :param including_shear_properties: Decide whether to check shear information,
        which is pointless for liquids and equations of state
        without shear modulus parameterizations.
    :type including_shear_properties: bool

    :param equilibration_function: Function to internally equilibrate object.
        Called after every set_state.
        Takes the mineral object as its only argument.
    :type equilibration_function: function

    :returns: Boolean stating whether all checks have passed.
    :rtype: bool
    """
    if equilibration_function is None:

        def equilibration_function(mineral):
            pass

    dT = 1.0
    dP = 1000.0

    m.set_state(P, T)
    equilibration_function(m)
    G0 = m.gibbs
    S0 = m.S
    V0 = m.V

    expr = ["G = F + PV", "G = H - TS", "G = E - TS + PV"]
    eq = [
        [m.gibbs, (m.helmholtz + P * m.V)],
        [m.gibbs, (m.H - T * m.S)],
        [m.gibbs, (m.internal_energy - T * m.S + P * m.V)],
    ]

    m.set_state(P, T + dT)
    equilibration_function(m)
    G1 = m.gibbs
    S1 = m.S
    V1 = m.V

    m.set_state(P + dP, T)
    equilibration_function(m)
    G2 = m.gibbs
    V2 = m.V

    # T derivatives
    m.set_state(P, T + 0.5 * dT)
    equilibration_function(m)
    expr.extend(["S = -dG/dT", "alpha = 1/V dV/dT", "C_p = T dS/dT"])
    eq.extend(
        [
            [m.S, -(G1 - G0) / dT],
            [m.alpha, (V1 - V0) / dT / (V1 + V0) * 2.0],
            [
                m.heat_capacity_p,
                (T + 0.5 * dT) * (S1 - S0) / dT,
            ],
        ]
    )

    # P derivatives
    m.set_state(P + 0.5 * dP, T)
    equilibration_function(m)
    expr.extend(["V = dG/dP", "K_T = -V dP/dV"])
    eq.extend(
        [
            [m.V, (G2 - G0) / dP],
            [m.isothermal_bulk_modulus_reuss, -0.5 * (V2 + V0) * dP / (V2 - V0)],
        ]
    )

    expr.extend(
        ["C_v = Cp - alpha^2*K_T*V*T", "K_S = K_T*Cp/Cv", "gr = alpha*K_T*V/Cv"]
    )
    eq.extend(
        [
            [
                m.heat_capacity_v,
                m.heat_capacity_p
                - m.alpha * m.alpha * m.isothermal_bulk_modulus_reuss * m.V * T,
            ],
            [
                m.isentropic_bulk_modulus_reuss,
                m.isothermal_bulk_modulus_reuss * m.heat_capacity_p / m.heat_capacity_v,
            ],
            [
                m.gr,
                m.alpha * m.isothermal_bulk_modulus_reuss * m.V / m.heat_capacity_v,
            ],
        ]
    )

    if including_shear_properties:
        expr.append("Vphi = np.sqrt(K_S/rho)")
        eq.append(
            [m.bulk_sound_velocity, np.sqrt(m.isentropic_bulk_modulus_reuss / m.rho)]
        )

        expr.extend(["Vp = np.sqrt((K_S + 4G/3)/rho)", "Vs = np.sqrt(G_S/rho)"])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            eq.extend(
                [
                    [
                        m.p_wave_velocity,
                        np.sqrt(
                            (m.isentropic_bulk_modulus_reuss + 4.0 * m.G / 3.0) / m.rho
                        ),
                    ],
                    [m.shear_wave_velocity, np.sqrt(m.G / m.rho)],
                ]
            )
            if len(w) == 1:
                print(w[0].message)
                print(
                    "\nYou can suppress this message by setting the "
                    "parameter\nincluding_shear_properties to False "
                    "when calling check_eos_consistency.\n"
                )
        note = ""
    else:
        note = " (not including shear properties)"

    consistencies = [
        np.abs(e[0] - e[1]) < np.abs(tol * e[1]) + np.finfo("float").eps for e in eq
    ]
    eos_is_consistent = np.all(consistencies)

    if verbose:
        print("Checking EoS consistency for {0:s}{1}".format(m.to_string(), note))
        print("Expressions within tolerance of {0:2f}".format(tol))
        for i, c in enumerate(consistencies):
            print("{0:10s} : {1:5s}".format(expr[i], str(c)))
            if not c:
                print(eq[i])
        if eos_is_consistent:
            print(
                "All EoS consistency constraints satisfied for {0:s}".format(
                    m.to_string()
                )
            )
        else:
            print(
                "Not satisfied all EoS consistency constraints for {0:s}".format(
                    m.to_string()
                )
            )

    return eos_is_consistent


def check_anisotropic_eos_consistency(
    m, P=1.0e9, T=300.0, tol=1.0e-4, verbose=False, equilibration_function=None
):
    """
    Checks that numerical derivatives of the Gibbs energy of an anisotropic mineral
    under given conditions are equal to those provided
    analytically by the equation of state.

    :param m: The anisotropic mineral for which the equation of state
        is to be checked for consistency.
    :type m: :class:`burnman.AnisotropicMineral`

    :param P: The pressure at which to check consistency.
    :type P: float

    :param T: The temperature at which to check consistency.
    :type T: float

    :param tol: The fractional tolerance for each of the checks.
    :type tol: float

    :param verbose: Decide whether to print information about each check.
    :type verbose: bool

    :param equilibration_function: Function to internally equilibrate object.
        Called after every set_state.
        Takes the mineral object as its only argument.
    :type equilibration_function: function

    :returns: Boolean stating whether all checks have passed.
    :rtype: bool
    """
    if equilibration_function is None:

        def equilibration_function(mineral):
            pass

    dT = 1.0
    dP = 1000.0

    m.set_state(P, T)
    equilibration_function(m)
    G0 = m.gibbs
    S0 = m.S
    V0 = m.V

    expr = ["G = F + PV", "G = H - TS", "G = E - TS + PV"]
    eq = [
        [m.gibbs, (m.helmholtz + P * m.V)],
        [m.gibbs, (m.H - T * m.S)],
        [m.gibbs, (m.internal_energy - T * m.S + P * m.V)],
    ]

    m.set_state(P, T + dT)
    equilibration_function(m)
    G1 = m.gibbs
    S1 = m.S
    V1 = m.V

    m.set_state(P + dP, T)
    equilibration_function(m)
    G2 = m.gibbs
    V2 = m.V

    # T derivatives
    m.set_state(P, T + 0.5 * dT)
    equilibration_function(m)
    expr.extend(["S = -dG/dT", "alpha = 1/V dV/dT", "C_p = T dS/dT"])
    eq.extend(
        [
            [m.S, -(G1 - G0) / dT],
            [m.alpha, (V1 - V0) / dT / m.V],
            [m.heat_capacity_p, (T + 0.5 * dT) * (S1 - S0) / dT],
        ]
    )

    # P derivatives
    m.set_state(P + 0.5 * dP, T)
    equilibration_function(m)
    expr.extend(["V = dG/dP", "K_T = -V dP/dV"])
    eq.extend(
        [
            [m.V, (G2 - G0) / dP],
            [m.isothermal_bulk_modulus_reuss, -0.5 * (V2 + V0) * dP / (V2 - V0)],
        ]
    )

    expr.extend(["C_v = Cp - alpha^2*K_T*V*T", "K_S = K_T*Cp/Cv"])
    eq.extend(
        [
            [
                m.heat_capacity_v,
                m.heat_capacity_p
                - m.alpha * m.alpha * m.isothermal_bulk_modulus_reuss * m.V * T,
            ],
            [
                m.isentropic_bulk_modulus_reuss,
                m.isothermal_bulk_modulus_reuss * m.heat_capacity_p / m.heat_capacity_v,
            ],
        ]
    )

    # Third derivative
    m.set_state(P + 0.5 * dP, T)
    equilibration_function(m)
    b0 = m.isothermal_compressibility_tensor
    F0 = m.deformation_gradient_tensor

    m.set_state(P + 0.5 * dP, T + dT)
    equilibration_function(m)
    b1 = m.isothermal_compressibility_tensor
    F1 = m.deformation_gradient_tensor

    m.set_state(P, T + 0.5 * dT)
    equilibration_function(m)
    a0 = m.thermal_expansivity_tensor
    F2 = m.deformation_gradient_tensor

    m.set_state(P + dP, T + 0.5 * dT)
    equilibration_function(m)
    a1 = m.thermal_expansivity_tensor
    F3 = m.deformation_gradient_tensor

    m.set_state(P + 0.5 * dP, T + 0.5 * dT)
    equilibration_function(m)

    # Compare with derivatives of F in the unrotated frame.
    Q = m.rotation_matrix.T
    beta1 = m.isothermal_compressibility_tensor
    alpha1 = m.thermal_expansivity_tensor

    beta1 = np.einsum("mi, nj, ij->mn", Q, Q, beta1)
    alpha1 = np.einsum("mi, nj, ij->mn", Q, Q, alpha1)

    if m.orthotropic:
        beta0 = -(logm(F3) - logm(F2)) / dP
        alpha0 = (logm(F1) - logm(F0)) / dT

        expr.extend(
            [f"SI = -d(lnm(F))/dP ({i}{j})" for i in range(3) for j in range(i, 3)]
        )
        expr.extend(
            [f"alpha = d(lnm(F))/dT ({i}{j})" for i in range(3) for j in range(i, 3)]
        )
    else:
        expr.extend(
            [
                f"SI = -0.5(dF/dP*F^-1 + (dF/dP*F^-1)^T) ({i}{j})"
                for i in range(3)
                for j in range(i, 3)
            ]
        )
        expr.extend(
            [
                f"alpha = 0.5(dF/dT*F^-1 + (dF/dT*F^-1)^T) ({i}{j})"
                for i in range(3)
                for j in range(i, 3)
            ]
        )
        invF = np.linalg.inv(m.deformation_gradient_tensor)
        dFdP = (F3 - F2) / dP
        LP = np.einsum("ij,kj->ik", dFdP, invF)
        beta0 = -0.5 * (LP + LP.T)
        dFdT = (F1 - F0) / dT
        LT = np.einsum("ij,kj->ik", dFdT, invF)
        alpha0 = 0.5 * (LT + LT.T)

    eq.extend([[beta0[i, j], beta1[i, j]] for i in range(3) for j in range(i, 3)])
    eq.extend([[alpha0[i, j], alpha1[i, j]] for i in range(3) for j in range(i, 3)])

    expr.extend(
        [f"d(alpha)/dP = -d(beta_T)/dT ({i}{j})" for i in range(3) for j in range(i, 3)]
    )
    eq.extend(
        [
            [(a1[i, j] - a0[i, j]) / dP, -(b1[i, j] - b0[i, j]) / dT]
            for i in range(3)
            for j in range(i, 3)
        ]
    )

    # Consistent Phi
    expr.extend(["dPsidf_Voigt[:3,:3] == 1"])
    eq.extend(
        [
            [
                np.sum(m.isothermal_compliance_tensor[:3, :3]),
                m.isothermal_compressibility_reuss,
            ]
        ]
    )

    # Consistent inverses
    expr.extend([f"S_T = inv(C_T) ({i}{j})" for i in range(6) for j in range(i, 6)])
    S_T = m.isothermal_compliance_tensor
    S_T2 = np.linalg.inv(m.isothermal_stiffness_tensor)
    eq.extend([[S_T[i, j], S_T2[i, j]] for i in range(6) for j in range(i, 6)])

    expr.extend([f"S_N = inv(C_N) ({i}{j})" for i in range(6) for j in range(i, 6)])
    S_N = m.isentropic_compliance_tensor
    S_N2 = np.linalg.inv(m.isentropic_stiffness_tensor)
    eq.extend([[S_N[i, j], S_N2[i, j]] for i in range(6) for j in range(i, 6)])

    # Consistent isotropic and anisotropic properties
    expr.extend(
        [
            "V = det(M)",
            "alpha_v = tr(alpha)",
            "beta_RT = sum(I S_T I)",
            "beta_RS = sum(I S_S I)",
        ]
    )
    eq.extend(
        [
            [m.V, np.linalg.det(m.cell_vectors)],
            [m.alpha, np.trace(m.thermal_expansivity_tensor)],
            [
                m.isothermal_compressibility_reuss,
                np.sum(m.isothermal_compliance_tensor[:3, :3]),
            ],
            [
                m.isentropic_compressibility_reuss,
                np.sum(m.isentropic_compliance_tensor[:3, :3]),
            ],
        ]
    )

    expr.append("Vphi = np.sqrt(K_S/rho)")
    eq.append([m.bulk_sound_velocity, np.sqrt(m.isentropic_bulk_modulus_reuss / m.rho)])

    consistencies = [
        np.abs(e[0] - e[1]) < np.abs(tol * e[1]) + np.finfo("float").eps for e in eq
    ]
    eos_is_consistent = np.all(consistencies)

    if verbose:
        print("Checking EoS consistency for {0:s}".format(m.to_string()))
        print("Expressions within tolerance of {0:2f}".format(tol))
        for i, c in enumerate(consistencies):
            print("{0:10s} : {1:5s}".format(expr[i], str(c)))
            if not c:
                print(eq[i])
        if eos_is_consistent:
            print(
                "All EoS consistency constraints satisfied for {0:s}".format(
                    m.to_string()
                )
            )
        else:
            print(
                "Not satisfied all EoS consistency constraints for {0:s}".format(
                    m.to_string()
                )
            )

    return eos_is_consistent


_test_deformation_gradient = np.array(
    [[0.9, 0.1, 0.05], [0.1, 0.91, 0.04], [0.02, 0.03, 0.89]], dtype=float
)


def check_hyperelastic_eos_consistency(
    m,
    deformation_gradient=_test_deformation_gradient,
    T=300.0,
    tol=1.0e-4,
    verbose=False,
    *,
    strain_step=2.0e-4,
    temperature_step=0.5,
):
    """
    Check Helmholtz derivatives, Maxwell relations and frame indifference.

    The checked cell matrix is ``M = deformation_gradient @ M_reference``,
    where cell vectors are rows of ``M.T``.
    All stresses and tangents are compared in the same spatial frame.
    Crystal-frame output is transformed back into the spatial frame at each
    perturbed state; the selected output frame is left unchanged.
    Matrix comparisons use Voigt notation with engineering shear strains
    (2 epsilon_23, 2 epsilon_13, 2 epsilon_12) and unscaled shear stresses.

    The molar Helmholtz differential is
    :math:`dF = -S dT + V \\sigma : d\\epsilon`. Consequently,
    :math:`\\pi = (\\partial\\sigma/\\partial T)_M
    = -(1/V)(\\partial S/\\partial\\epsilon)_T` and
    :math:`C_S = C_T + VT \\pi\\otimes\\pi/C_\\epsilon`.
    The spatial stiffness and compliance tensors need not have major symmetry.

    An existing cell/temperature state is restored, even if exceptions are raised.
    An initially unset mineral is left at the checked state.

    :param m: Hyperelastic mineral to check.
    :param deformation_gradient: Finite (3, 3) matrix with positive determinant,
        declared relative to ``m.hydrostatic.cell_vectors_0``.
    :param T: Temperature in Kelvin.
    :param tol: Fractional tolerance, relative to the largest absolute entry
        in each compared scalar or tensor. Entropy and energy comparisons
        use C_eps and max(T C_eps, V max(abs(C_T))) as characteristic scales,
        so zero reference values do not require zero numerical error.
    :param verbose: Print each result and its maximum absolute error.
    :param strain_step: Positive dimensionless central-difference half-step.
    :param temperature_step: Positive temperature half-step in Kelvin, less
        than T. The perturbed states must be within the model's valid range.
    :returns: Whether every check passes.
    :rtype: bool
    """
    deformation_gradient = np.asarray(deformation_gradient, dtype=float)
    if (
        deformation_gradient.shape != (3, 3)
        or not np.all(np.isfinite(deformation_gradient))
        or np.linalg.det(deformation_gradient) <= 0.0
    ):
        raise ValueError(
            "deformation_gradient must be finite, (3, 3), and have positive determinant"
        )
    if not np.isfinite(T) or T <= 0.0:
        raise ValueError("T must be finite and positive")
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive")
    if not np.isfinite(strain_step) or strain_step <= 0.0:
        raise ValueError("strain_step must be finite and positive")
    if not np.isfinite(temperature_step) or not 0.0 < temperature_step < T:
        raise ValueError("temperature_step must be finite, positive, and less than T")

    cell_vectors = (deformation_gradient @ m.hydrostatic.cell_vectors_0.T).T
    try:
        # Restore physical input vectors, not a rotated display of them.
        original_cell_vectors = getattr(m, "_cell_vectors_current", None)
        if original_cell_vectors is None:
            original_cell_vectors = m.cell_vectors
        original_cell_vectors, original_T = original_cell_vectors.copy(), m.T
    except AttributeError:
        original_cell_vectors, original_T = cell_vectors.copy(), T

    # Unit engineering strains: shear entries in the full tensor are 1/2.
    strains = np.array([expand_strains(e) for e in np.eye(6)])

    def set_state(cell_vectors, temperature):
        m.set_state_with_molar_cell_vectors(cell_vectors, temperature)

    def spatial_tensor(property_name):
        tensor = np.array(getattr(m, property_name), copy=True)
        if getattr(m, "output_frame") == "crystal":
            # Q maps the supplied spatial cell into the crystal frame and
            # varies with strain. Undo it before forming a derivative: the
            # derivative of moving-frame components includes extra frame terms.
            Q = m.rotation_to_crystallographic_frame.T
            if tensor.ndim == 2:
                return Q @ tensor @ Q.T
            return np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, tensor)
        return tensor

    results = []

    def compare(name, actual, expected, scale=0.0, roundoff_scale=0.0):
        actual, expected = np.asarray(actual), np.asarray(expected)
        scale = max(scale, np.max(np.abs(actual)), np.max(np.abs(expected)))
        error = np.max(np.abs(actual - expected))
        passed = bool(
            np.all(np.isfinite(actual))
            and np.all(np.isfinite(expected))
            and error <= tol * scale + 64.0 * np.finfo(float).eps * roundoff_scale
        )
        results.append((name, passed, error))

    de, dT = strain_step, temperature_step
    try:
        set_state(cell_vectors, T)
        V, F, S = m.V, m.helmholtz, m.entropy
        sigma = spatial_tensor("cauchy_stress")
        C_full = spatial_tensor("full_isothermal_stiffness_tensor")
        CS_full = spatial_tensor("full_isentropic_stiffness_tensor")
        compliance = spatial_tensor("full_isothermal_compliance_tensor")
        compliance_S = spatial_tensor("full_isentropic_compliance_tensor")
        C, CS = contract_stiffnesses(C_full), contract_stiffnesses(CS_full)
        Ceps = m.molar_isometric_heat_capacity
        Csigma = m.molar_isostress_heat_capacity
        pi = spatial_tensor("thermal_stress")
        alpha = spatial_tensor("thermal_expansivity_tensor")
        sigma_v, pi_v = contract_stresses(sigma), contract_stresses(pi)
        alpha_v = contract_strains(alpha)
        entropy_scale = abs(Ceps)
        energy_scale = max(T * entropy_scale, V * np.max(np.abs(C)))

        compare("V = det(M)", V, np.linalg.det(cell_vectors))
        compare("U = F + TS", m.internal_energy, F + T * S, scale=energy_scale)
        compare("P = -tr(sigma)/3", m.pressure, -np.trace(sigma) / 3.0)
        compare("symmetric Cauchy stress", sigma, sigma.T)
        compare("symmetric thermal stress", pi, pi.T)
        compare("symmetric thermal expansivity", alpha, alpha.T)
        for name, tensor in [
            ("C_T", C_full),
            ("C_S", CS_full),
            ("S_T", compliance),
            ("S_S", compliance_S),
        ]:
            compare(name + " first minor symmetry", tensor, tensor.swapaxes(0, 1))
            compare(name + " second minor symmetry", tensor, tensor.swapaxes(2, 3))
        compare("S_T C_T = I", contract_compliances(compliance) @ C, np.eye(6))
        compare("S_S C_S = I", contract_compliances(compliance_S) @ CS, np.eye(6))

        # The volume factor in dF/de contributes sigma (x) I. Testing
        # C_T == C_T.T alone would incorrectly reject nonhydrostatic states.
        curvature = C + np.outer(sigma_v, contract_stresses(np.eye(3)))
        compare("Maxwell: C_T + sigma (x) I is symmetric", curvature, curvature.T)

        set_state(cell_vectors, T + dT)
        F_plus, S_plus = m.helmholtz, m.entropy
        sigma_plus_T = spatial_tensor("cauchy_stress")
        set_state(cell_vectors, T - dT)
        F_minus, S_minus = m.helmholtz, m.entropy
        sigma_minus_T = spatial_tensor("cauchy_stress")
        dSdT = (S_plus - S_minus) / (2.0 * dT)
        dsigma_dT = (sigma_plus_T - sigma_minus_T) / (2.0 * dT)
        compare(
            "S = -dF/dT at fixed M",
            S,
            -(F_plus - F_minus) / (2.0 * dT),
            scale=entropy_scale,
        )
        compare("C_eps = T dS/dT at fixed M", Ceps, T * dSdT)
        compare("pi = d(sigma)/dT at fixed M", pi, dsigma_dT)

        stress_from_energy = np.empty(6)
        dSde = np.empty(6)
        stress_derivatives = np.empty((6, 6))
        for k, E in enumerate(strains):
            set_state((expm(de * E) @ cell_vectors.T).T, T)
            F_plus, S_plus = m.helmholtz, m.entropy
            sigma_plus = spatial_tensor("cauchy_stress")
            set_state((expm(-de * E) @ cell_vectors.T).T, T)
            F_minus, S_minus = m.helmholtz, m.entropy
            sigma_minus = spatial_tensor("cauchy_stress")
            stress_from_energy[k] = (F_plus - F_minus) / (2.0 * de * V)
            dSde[k] = (S_plus - S_minus) / (2.0 * de)
            stress_derivatives[:, k] = contract_stresses(sigma_plus - sigma_minus) / (
                2.0 * de
            )

        compare("sigma = (1/V) dF/dstrain", sigma_v, stress_from_energy)
        compare("C_T = d(sigma)/dstrain", C, stress_derivatives)
        compare("Maxwell: pi = -(1/V) dS/dstrain", pi_v, -dSde / V)
        compare(
            "Maxwell: d(sigma)/dT = -(1/V) dS/dstrain",
            contract_stresses(dsigma_dT),
            -dSde / V,
        )
        compare("C_T : alpha = -pi", C @ alpha_v, -pi_v)
        # Compare thermal increments themselves, allowing roundoff from
        # subtraction: a large elastic modulus must not hide a wrong update.
        compare(
            "C_sigma = C_eps - VT pi:alpha",
            Csigma - Ceps,
            -V * T * pi_v @ alpha_v,
            roundoff_scale=max(abs(Csigma), abs(Ceps)),
        )
        compare(
            "C_S = C_T + VT pi (x) pi / C_eps",
            CS - C,
            V * T / Ceps * np.outer(pi_v, pi_v),
            roundoff_scale=max(np.max(np.abs(CS)), np.max(np.abs(C))),
        )
        # Independent chain rule along dS = 0, using measured derivatives.
        compare(
            "C_S from constant-entropy derivative",
            CS,
            stress_derivatives - np.outer(contract_stresses(dsigma_dT), dSde) / dSdT,
        )

        # Follow the first-order isostress path: dstrain = alpha dT.
        # Maxwell gives dS = C_eps/T dT - V pi:dstrain.
        set_state((expm(dT * alpha) @ cell_vectors.T).T, T + dT)
        S_plus = m.entropy
        set_state((expm(-dT * alpha) @ cell_vectors.T).T, T - dT)
        S_minus = m.entropy
        compare(
            "C_sigma = T dS/dT along isostress direction",
            Csigma,
            T * (S_plus - S_minus) / (2.0 * dT),
        )

        # Active rotation of the actual cell matrix; all public spatial
        # tensors must transform with this same Q, including at finite strain.
        Q = expm(np.array([[0.0, -0.37, -0.23], [0.37, 0.0, -0.19], [0.23, 0.19, 0.0]]))
        set_state((Q @ cell_vectors.T).T, T)
        for name, actual, expected, scale in [
            ("volume", m.V, V, 0.0),
            ("Helmholtz energy", m.helmholtz, F, energy_scale),
            ("entropy", m.entropy, S, entropy_scale),
            ("isometric heat capacity", m.molar_isometric_heat_capacity, Ceps, 0.0),
            ("isostress heat capacity", m.molar_isostress_heat_capacity, Csigma, 0.0),
        ]:
            compare("objective " + name, actual, expected, scale=scale)
        for name, actual, expected in [
            ("Cauchy stress", spatial_tensor("cauchy_stress"), sigma),
            ("thermal stress", spatial_tensor("thermal_stress"), pi),
            (
                "thermal expansivity",
                spatial_tensor("thermal_expansivity_tensor"),
                alpha,
            ),
        ]:
            compare("objective " + name, actual, Q @ expected @ Q.T)
        for name, actual, expected in [
            (
                "isothermal stiffness",
                spatial_tensor("full_isothermal_stiffness_tensor"),
                C_full,
            ),
            (
                "isentropic stiffness",
                spatial_tensor("full_isentropic_stiffness_tensor"),
                CS_full,
            ),
            (
                "isothermal compliance",
                spatial_tensor("full_isothermal_compliance_tensor"),
                compliance,
            ),
            (
                "isentropic compliance",
                spatial_tensor("full_isentropic_compliance_tensor"),
                compliance_S,
            ),
        ]:
            compare(
                "objective " + name,
                actual,
                np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, expected),
            )
    finally:
        set_state(original_cell_vectors, original_T)

    if verbose:
        print(f"Checking hyperelastic EoS consistency for {m.to_string()}")
        for name, passed, error in results:
            print(f"{name}: {passed} (max absolute error: {error:.6g})")

    return all(passed for _, passed, _ in results)
