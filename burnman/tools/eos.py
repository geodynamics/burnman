# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import warnings
from scipy.linalg import logm


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
    under given conditions are equal to those provided
    analytically by the equation of state.

    :param m: The mineral for which the equation of state
        is to be checked for consistency.
    :type m: :class:`burnman.Mineral`

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
        [m.gibbs, (m.molar_internal_energy - T * m.S + P * m.V)],
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
            [m.molar_heat_capacity_p, (T + 0.5 * dT) * (S1 - S0) / dT],
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
                m.molar_heat_capacity_v,
                m.molar_heat_capacity_p
                - m.alpha * m.alpha * m.isothermal_bulk_modulus_reuss * m.V * T,
            ],
            [
                m.isentropic_bulk_modulus_reuss,
                m.isothermal_bulk_modulus_reuss
                * m.molar_heat_capacity_p
                / m.molar_heat_capacity_v,
            ],
            [
                m.gr,
                m.alpha
                * m.isothermal_bulk_modulus_reuss
                * m.V
                / m.molar_heat_capacity_v,
            ],
        ]
    )

    expr.append("Vphi = np.sqrt(K_S/rho)")
    eq.append([m.bulk_sound_velocity, np.sqrt(m.isentropic_bulk_modulus_reuss / m.rho)])

    if including_shear_properties:
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
        [m.gibbs, (m.molar_internal_energy - T * m.S + P * m.V)],
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
            [m.molar_heat_capacity_p, (T + 0.5 * dT) * (S1 - S0) / dT],
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
                m.molar_heat_capacity_v,
                m.molar_heat_capacity_p
                - m.alpha * m.alpha * m.isothermal_bulk_modulus_reuss * m.V * T,
            ],
            [
                m.isentropic_bulk_modulus_reuss,
                m.isothermal_bulk_modulus_reuss
                * m.molar_heat_capacity_p
                / m.molar_heat_capacity_v,
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

    beta0 = -(logm(F3) - logm(F2)) / dP
    alpha0 = (logm(F1) - logm(F0)) / dT

    Q = m.deformed_coordinate_frame
    beta1 = m.isothermal_compressibility_tensor
    alpha1 = m.thermal_expansivity_tensor

    beta1 = np.einsum("mi, nj, ij->mn", Q, Q, beta1)
    alpha1 = np.einsum("mi, nj, ij->mn", Q, Q, alpha1)

    if m.orthotropic:
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
