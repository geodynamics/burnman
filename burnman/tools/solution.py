# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
from ..classes.combinedmineral import CombinedMineral
from ..classes.solution import Solution
from ..classes.solutionmodel import (
    IdealSolution,
    SymmetricRegularSolution,
    AsymmetricRegularSolution,
    SubregularSolution,
    logish,
)
from ..utils.chemistry import site_occupancies_to_strings
from ..utils.math import complete_basis
from ..constants import gas_constant as R


def _decompose_3D_matrix(W):
    """
    Decomposes a 3D matrix W_ijk where E = W_ijk p_i p_j p_k
    into a subregular form where
    E = G_i p_i + WB_ij (1 - p_j + p_i) / 2 + WT_ijk p_i p_j p_k,
    and i < j < k.

    :param W: 3D interaction matrix.
    :type W: numpy.array

    :returns: The 1D array G_i, the 2D upper triangular array WB_ij and
        the ternary terms in a list where each item is in the form [i, j, k, WT_ijk]
    :rtype: tuple
    """

    n_mbrs = len(W)
    # New endmember components
    # W_iii needs to be copied, otherwise just a view onto W
    new_endmember_excesses = np.copy(np.einsum("iii->i", W))

    # Removal of endmember components from 3D representation
    W -= (
        np.einsum(
            "i, j, k->ijk", new_endmember_excesses, np.ones(n_mbrs), np.ones(n_mbrs)
        )
        + np.einsum(
            "i, j, k->ijk", np.ones(n_mbrs), new_endmember_excesses, np.ones(n_mbrs)
        )
        + np.einsum(
            "i, j, k->ijk", np.ones(n_mbrs), np.ones(n_mbrs), new_endmember_excesses
        )
    ) / 3.0

    # Transformed 2D components
    # (i=j, i=k, j=k)
    new_binary_matrix = (
        np.einsum("jki, jk -> ij", W, np.identity(n_mbrs))
        + np.einsum("jik, jk -> ij", W, np.identity(n_mbrs))
        + np.einsum("ijk, jk -> ij", W, np.identity(n_mbrs))
    ).round(decimals=12)

    # Wb is the 3D matrix corresponding to the terms in the binary matrix,
    # such that the two following print statements produce the same answer
    # for a given array of endmember proportions
    Wb = (
        np.einsum("ijk, ij->ijk", W, np.identity(n_mbrs))
        + np.einsum("ijk, jk->ijk", W, np.identity(n_mbrs))
        + np.einsum("ijk, ik->ijk", W, np.identity(n_mbrs))
    )

    # Remove binary component from 3D representation
    # The extra terms are needed because the binary term in the formulation
    # of a subregular solution model given by
    # Helffrich and Wood includes ternary components (the sum_k X_k part)..
    W -= (
        Wb
        + (
            np.einsum("ij, k", new_binary_matrix, np.ones(n_mbrs))
            - np.einsum("ij, ik->ijk", new_binary_matrix, np.identity(n_mbrs))
            - np.einsum("ij, jk->ijk", new_binary_matrix, np.identity(n_mbrs))
        )
        / 2.0
    )

    # Find the 3D components Wijk by adding the elements at
    # the six equivalent positions in the matrix
    new_ternary_terms = []
    for i in range(n_mbrs):
        for j in range(i + 1, n_mbrs):
            for k in range(j + 1, n_mbrs):
                val = (
                    W[i, j, k]
                    + W[j, k, i]
                    + W[k, i, j]
                    + W[k, j, i]
                    + W[j, i, k]
                    + W[i, k, j]
                ).round(decimals=12)
                if np.abs(val) > 1.0e-12:
                    new_ternary_terms.append([i, j, k, val])

    return (new_endmember_excesses, new_binary_matrix, new_ternary_terms)


def _subregular_matrix_conversion(new_basis, Wijk, n_mbrs):
    """
    Converts the arrays reguired to describe a subregular solution
    from one endmember basis to another.

    The excess nonconfigurational energies of the subregular solution model
    are described as follows:
    E = G_i p_i + WB_ij (1 - p_j + p_i) / 2 + WT_ijk p_i p_j p_k,
    and i < j < k.

    :param new_basis: The new endmember basis, given as amounts of the old endmembers.
    :type new_basis: 2D numpy array

    :param binary_matrix: The upper triangular matrix WB_ij.
    :type binary_matrix: 2D numpy array

    :param ternary_terms: The ternary terms in a list where each
        item is in the form [i, j, k, WT_ijk]
    :type ternary_terms: list of lists of length 4

    :param endmember_excesses: The array G_i
    :type endmember_excesses: 1D numpy array

    :returns: The 1D array G_i, the 2D upper triangular array WB_ij and
        the ternary terms in a list where each item is in the form [i, j, k, WT_ijk]
    :rtype: tuple
    """
    # Transformation to new 3D representation
    A = new_basis.T
    Wn = np.einsum("il, jm, kn, ijk -> lmn", A, A, A, Wijk)
    Wn = Wn[:n_mbrs, :n_mbrs, :n_mbrs]
    new_endmember_excesses, new_binary_terms, new_ternary_terms = _decompose_3D_matrix(
        Wn
    )

    return (new_endmember_excesses, new_binary_terms, new_ternary_terms)


def transform_solution_to_new_basis(
    solution,
    new_basis,
    n_mbrs=None,
    solution_name=None,
    endmember_names=None,
    molar_fractions=None,
):
    """
    Transforms a solution model from one endmember basis to another.
    Returns a new Solution object.

    :param solution: The original solution object.
    :type solution: :class:`burnman.Solution` object

    :param new_basis: The new endmember basis, given as amounts of the old endmembers.
    :type new_basis: 2D numpy array

    :param n_mbrs: The number of endmembers in the new solution
        (defaults to the length of new_basis).
    :type n_mbrs: float, optional

    :param solution_name: A name corresponding to the new solution.
    :type solution_name: str, optional

    :param endmember_names: A list corresponding to the names of the new endmembers.
    :type endmember_names: list of str, optional

    :param molar_fractions: Fractions of the new endmembers in the new solution.
    :type molar_fractions: numpy.array, optional

    :returns: The transformed solution.
    :rtype: :class:`burnman.Solution` object
    """
    new_basis = np.array(new_basis)
    if n_mbrs is None:
        n_mbrs, n_all_mbrs = new_basis.shape
    else:
        _, n_all_mbrs = new_basis.shape

    if solution_name is None:
        name = "child solution"
    else:
        name = solution_name

    solution_model = solution.solution_model

    # Use type here to avoid inheritance problems
    solution_type = type(solution_model)
    if solution_type == IdealSolution:
        ESV_modifiers = [[0.0, 0.0, 0.0] for _ in new_basis]

    elif (
        solution_type == AsymmetricRegularSolution
        or solution_type == SymmetricRegularSolution
    ):
        A = complete_basis(new_basis).T

        old_alphas = solution.solution_model.alphas
        alphas = np.einsum("i, ij", solution.solution_model.alphas, A)
        inv_diag_alphas = np.diag(1.0 / alphas)
        B = np.einsum("ij, jk, kl->il", np.diag(old_alphas), A, inv_diag_alphas)
        alphas = list(alphas[0:n_mbrs])
        Qe = np.einsum("ik, ij, kl->jl", solution.solution_model.We, B, B)
        Qs = np.einsum("ik, ij, kl->jl", solution.solution_model.Ws, B, B)
        Qv = np.einsum("ik, ij, kl->jl", solution.solution_model.Wv, B, B)

        def new_interactions(Q, n_mbrs):
            return [
                [
                    float(
                        (Q[i, j] + Q[j, i] - Q[i, i] - Q[j, j])
                        * (alphas[i] + alphas[j])
                        / 2.0
                    )
                    for j in range(i + 1, n_mbrs)
                ]
                for i in range(n_mbrs - 1)
            ]

        energy_interaction = new_interactions(Qe, n_mbrs)
        entropy_interaction = new_interactions(Qs, n_mbrs)
        volume_interaction = new_interactions(Qv, n_mbrs)

        ESV_modifiers = [
            [Qe[i, i] * alphas[i], Qs[i, i] * alphas[i], Qv[i, i] * alphas[i]]
            for i in range(n_mbrs)
        ]

    elif solution_type == SubregularSolution:
        full_basis = complete_basis(new_basis)

        def new_interactions(W, n_mbrs):
            return [
                [[W[i, j], W[j, i]] for j in range(i + 1, n_mbrs)]
                for i in range(n_mbrs - 1)
            ]

        # N.B. initial endmember_excesses are zero
        Emod, We, ternary_e = _subregular_matrix_conversion(
            full_basis, solution.solution_model.Wijke, n_mbrs
        )
        Smod, Ws, ternary_s = _subregular_matrix_conversion(
            full_basis, solution.solution_model.Wijks, n_mbrs
        )
        Vmod, Wv, ternary_v = _subregular_matrix_conversion(
            full_basis, solution.solution_model.Wijkv, n_mbrs
        )
        energy_interaction = new_interactions(We, n_mbrs)
        entropy_interaction = new_interactions(Ws, n_mbrs)
        volume_interaction = new_interactions(Wv, n_mbrs)

        ESV_modifiers = [[Emod[i], Smod[i], Vmod[i]] for i in range(n_mbrs)]

    else:
        raise Exception(
            "The function to change basis for the "
            "{0} solution type has not yet been "
            "implemented.".format(solution_type)
        )

    # Create site formulae
    new_occupancies = np.array(new_basis).dot(
        solution.solution_model.endmember_occupancies
    )
    new_noccs = np.array(new_basis).dot(solution.solution_model.endmember_noccupancies)
    new_multiplicities = np.array(new_basis).dot(
        solution.solution_model.site_multiplicities
    )
    site_formulae = site_occupancies_to_strings(
        solution.solution_model.sites, new_multiplicities, new_occupancies
    )

    S_conf = -R * (new_noccs * (logish(new_noccs) - logish(new_multiplicities))).sum(-1)

    # Create endmembers
    endmembers = []
    for i, vector in enumerate(new_basis[:n_mbrs]):

        # Add ideal configurational entropy to endmembers
        ESV_modifiers[i][1] += S_conf[i]

        nonzero_indices = np.nonzero(vector)[0]
        if len(nonzero_indices) == 1:
            endmembers.append(
                [solution_model.endmembers[nonzero_indices[0]][0], site_formulae[i]]
            )
        else:
            mbr = CombinedMineral(
                [solution_model.endmembers[idx][0] for idx in nonzero_indices],
                [vector[idx] for idx in nonzero_indices],
                ESV_modifiers[i],
            )
            mbr.params["formula"] = {
                key: value
                for (key, value) in mbr.params["formula"].items()
                if value > 1.0e-12
            }
            endmembers.append([mbr, site_formulae[i]])

    if endmember_names is not None:
        for i in range(n_mbrs):
            endmembers[i][0].params["name"] = endmember_names[i]
            endmembers[i][0].name = endmember_names[i]

    if n_mbrs == 1:
        endmembers[0][0].name = name
        endmembers[0][0].parent = solution
        endmembers[0][0].basis = new_basis
        return endmembers[0][0]
    else:
        if solution_type == IdealSolution:
            new_solution_model = IdealSolution(endmembers=endmembers)
        elif solution_type == SymmetricRegularSolution:
            new_solution_model = type(solution_model)(
                endmembers=endmembers,
                energy_interaction=energy_interaction,
                volume_interaction=volume_interaction,
                entropy_interaction=entropy_interaction,
            )
        elif solution_type == SubregularSolution:
            new_solution_model = type(solution_model)(
                endmembers=endmembers,
                energy_interaction=energy_interaction,
                volume_interaction=volume_interaction,
                entropy_interaction=entropy_interaction,
                energy_ternary_terms=ternary_e,
                volume_ternary_terms=ternary_v,
                entropy_ternary_terms=ternary_s,
            )
        else:
            new_solution_model = type(solution_model)(
                endmembers=endmembers,
                energy_interaction=energy_interaction,
                volume_interaction=volume_interaction,
                entropy_interaction=entropy_interaction,
                alphas=alphas,
            )
        new_solution = Solution(
            name=name,
            solution_model=new_solution_model,
            molar_fractions=molar_fractions,
        )
        new_solution.parent = solution
        new_solution.basis = new_basis
        return new_solution
