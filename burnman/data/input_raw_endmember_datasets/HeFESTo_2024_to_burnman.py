# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This is a script that converts a version of the
# HeFESTo data format into the standard burnman
# format (printed to stdout)
# This script is used to generate the 2024 dataset,
# taking parameters from https://github.com/stixrude/HeFESTo_parameters_010123
# via git clone https://github.com/stixrude/HeFESTo_parameters_010123.git
# The SLB_2024 BurnMan file is then created via the commands
# python HeFESTo_2024_to_burnman.py > ../../minerals/SLB_2024.py; black ../../minerals/SLB_2024.py
# and from . import SLB_2024 added to the minerals __init__.py


import os
import pathlib
import pandas as pd
import numpy as np
from burnman.utils.chemistry import dictionarize_formula
from burnman.constants import gas_constant
from itertools import groupby
import re


def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


def rfloat(x, m=1.0):
    return round(float(x) * m, ndigits=10)


# Hard-coded into HeFESTo
"""
Ferropericlase:
Crystallography of wustites shows that the locations of the tetrahedral site
and the octahedral vacancies are highly correlated, with the tetrahedral
site being surrounded by a tetrahedral array of octahedral vacancies
(Gavarri et al. 1979 ).Because of this highly structured arrangement,
we ignore the contribution to the ideal entropy of mixing from the
tetrahedral site and the vacant octahedral sites.

...the formula of mag : [Fe]Fe2O4, with the square brackets indicating
that the first site does not participate in mixing. Thus, the ideal entropy
of mixing of an equimolar solution of wu and mag is 2Rln2 per 4 oxygens,
where R is the gas constant.

Spinel:
Pure magnetite is known to have an inverse arrangement of cations at
ambient conditions, rather than normal, as we have assumed.

Bridgmanite:
hepv and hlpv where only the cation on the
smaller B site undergoes the high spin to low spin transition...
"""

hardcoded_site_occupancies = {
    "wu": "[Fe]2[Fe]2",
    "wuls": "[Fels]2[Fels]2",
    "mag": "[Fef]0[Fef]2",  # first site not counted
    "smag": "[Fe][Fef][Fef]",
    "hepv": "[Fef][Fef]",
    "hlpv": "[Fef][Fefls]",
    "fapv": "[Fef][Al]",
    "hppv": "[Fef][Fef]",
    "hmag": "[Fe][Fef][Fef]",
    "acm": "[Na][Fef][Si]2",
    "andr": "[Ca]3[Fef][Fef]",
    "hem": "[Fef][Fef]",
    "lppv": "[Fefls][Fefls]",
}

magnetic_structural_parameters = {
    "fea": {"p": 0.4},
    "mag": {"p": 0.4},
    "smag": {"p": 0.4},
    "hmag": {"p": 0.4},
}

hefesto_path = "HeFESTo_parameters_010123"
mbrdir = pathlib.Path(hefesto_path)
phasedir = pathlib.Path(f"{hefesto_path}/phase")

param_dict = {}
with open("hefesto_parameter_names.txt") as file:
    lines = [line.rstrip() for line in file]
    for line in lines:
        param_dict[" ".join(line.split()[1:])] = line.split()[0]

ignored_mbrs = []

solution_aliases = {
    "c2c": "c2c_pyroxene",
    "cf": "calcium_ferrite_structured_phase",
    "cpx": "clinopyroxene",
    "gt": "garnet",
    "il": "ilmenite",
    "mw": "ferropericlase",
    "nal": "new_aluminous_phase",
    "ol": "olivine",
    "opx": "orthopyroxene",
    "plg": "plagioclase",
    "ppv": "post_perovskite",
    "pv": "bridgmanite",
    "ri": "ringwoodite",
    "sp": "mg_fe_aluminous_spinel",
    "wa": "wadsleyite",
}

# solutions
sol_params = {}
for f in phasedir.iterdir():
    # Do not process liquid!!
    if os.stat(f).st_size > 0 and str(f).split("/")[-1] != "liq":
        name = str(f).split("/")[-1]
        d = pd.read_csv(f, sep="\\s+")

        endmembers = list(d.columns)

        n_mbrs = len(endmembers)

        data = d.values.tolist()

        n_data = len(data)
        v0 = n_data - n_mbrs

        # Preprocess interaction matrices
        # to get rid of ignored endmembers
        use_indices = np.array(
            [i for i, mbr in enumerate(endmembers) if mbr not in ignored_mbrs]
        )
        use_endmembers = np.array(
            [mbr for mbr in endmembers if mbr not in ignored_mbrs]
        )
        n_use_mbrs = len(use_indices)

        energy_matrix = np.array(
            np.array(data)[np.ix_(use_indices, use_indices)], dtype=float
        )
        volume_matrix = np.array(
            np.array(data)[np.ix_(v0 + use_indices, use_indices)], dtype=float
        )

        if n_data == 2 * n_mbrs + 1:
            energy_interactions = []
            energy_interactions_2 = []
            volume_interactions = []
            for i in range(n_use_mbrs - 1):
                j = i + 1
                energy_interactions.append(
                    list(np.around(energy_matrix[i][j:] * 1000.0, decimals=6))
                )
                energy_interactions_2.append(
                    list(np.around(energy_matrix[i][:j] * 1000.0, decimals=6))
                )
                volume_interactions.append(
                    list(np.around(volume_matrix[i][j:] * 1.0e-6, decimals=20))
                )

            max_abs_e = np.max(
                np.abs([item for sublist in energy_interactions for item in sublist])
            )
            max_abs_v = np.max(
                np.abs([item for sublist in volume_interactions for item in sublist])
            )
            max_abs_lower_triangular_e2 = np.max(
                np.abs([item for sublist in energy_interactions_2 for item in sublist])
            )
            assert max_abs_lower_triangular_e2 < 1.0e-3

            # print(endmembers, n_mbrs)
            # print(energy_interactions)
            if max_abs_v > 1.0e-20:
                # print(volume_interactions)
                pass

            alphas = []
            mbr_sites = []
            for mbr in use_endmembers:
                with open(f"{hefesto_path}/{mbr}") as file:
                    lines = [line.rstrip() for line in file]
                    i_laar = [i for i, line in enumerate(lines) if "Laar" in line][0]
                    alphas.append(float(lines[i_laar].split()[0]))

                    formula_raw = lines[0].split()[0]

                    # Catch for MgTs mixing
                    if mbr == "mgts":
                        formula_raw = "Mg_1Al_1Si_2O_6"

                    # Catch for plag mixing
                    if mbr == "ab":
                        formula_raw = "Na_1Al_2Si_2O_8"

                    # Catch for mt mixing
                    if mbr == "mag":
                        formula_raw = "(Vac_1Fe_1)Fe_2O_4"

                    starts = [i for i, char in enumerate(formula_raw) if char.isupper()]
                    mixstarts = [i for i, char in enumerate(formula_raw) if char == "("]
                    mixends = [i for i, char in enumerate(formula_raw) if char == ")"]

                    for i in range(len(mixstarts)):
                        ms = mixstarts[i]
                        me = mixends[i]
                        starts = [i for i in starts if i < ms or i > me]

                    starts.extend(mixstarts)
                    starts.sort()
                    sites = [
                        formula_raw[i:j] for i, j in zip(starts, starts[1:] + [None])
                    ]
                    mbr_sites.append(sites)

            n_sites = [len(sites) for sites in mbr_sites]

            assert all_equal(n_sites)

            active_sites = []
            for i in range(n_sites[0]):
                occupancies = [sites[i] for sites in mbr_sites]
                if not all_equal(occupancies):
                    active_sites.append(i)

            site_formula_strings = []
            totals = []
            for mbr_site in mbr_sites:
                totals.append([])
                site_formula_strings.append("")
                for i in active_sites:
                    # Get number of atoms on sites
                    n_atoms = [int(s) for s in re.findall(r"\d+", mbr_site[i])]

                    species_starts = [
                        i for i, char in enumerate(mbr_site[i]) if char.isupper()
                    ]
                    species_ends = [
                        i for i, char in enumerate(mbr_site[i]) if char == "_"
                    ]
                    total_atoms = np.sum(n_atoms)
                    totals[-1].append(total_atoms)

                    site_string = "["
                    if len(n_atoms) > 1:
                        species = ""
                        for j in range(len(n_atoms)):
                            species += mbr_site[i][species_starts[j] : species_ends[j]]
                            species += f"{n_atoms[j]}/{total_atoms}"
                        site_string += species
                    else:
                        site_string += mbr_site[i][species_starts[0] : species_ends[0]]

                    site_string += "]"
                    if total_atoms != 1:
                        site_string += f"{total_atoms}"

                    site_formula_strings[-1] += site_string

            for i in range(len(totals[0])):
                t = [total[i] for total in totals]
                assert all_equal(t)

            solution_is_symmetric = np.max(np.abs(np.array(alphas) - 1.0)) < 1.0e-5

            sol_params[name] = {
                "name": name,
                "n_mbrs": n_use_mbrs,
                "mbr_names": use_endmembers,
                "mbr_site_formulae": site_formula_strings,
                "alphas": alphas,
            }

            if max_abs_e < 1.0:
                sol_params[name]["solution_type"] = "IdealSolution"
            elif solution_is_symmetric:
                sol_params[name]["solution_type"] = "SymmetricRegularSolution"
            else:
                sol_params[name]["solution_type"] = "AsymmetricRegularSolution"

            if max_abs_e > 1.0:
                sol_params[name]["energy_interaction"] = energy_interactions
            if max_abs_v > 1.0e-20:
                sol_params[name]["volume_interaction"] = volume_interactions

        else:
            pass

# endmembers
mbr_params = {}
mbr_mods = {}
mbr_names = {}
for f in mbrdir.iterdir():
    if os.stat(f).st_size > 0:
        name = str(f).split("/")[-1]

        process = False
        try:
            with open(f) as file:
                lines = [line.rstrip().split() for line in file]
                if lines[-1][-1] == "prime":
                    process = True
                if lines[0][-1] == "Liquid":
                    process = False
        except IsADirectoryError:
            pass

    if process:
        form = lines[0][0].replace("_", "").replace("(", "").replace(")", "")
        formula = dictionarize_formula(form)
        full_name = "_".join(lines[0][1:])
        # formula = formula_to_string(formula)
        # print(formula)
        # print(name)
        idict = {}
        for line in lines[1:]:
            key = param_dict[" ".join(line[1:])]
            idict[key] = line[0]

        if "T_crit" not in idict:
            idict["T_crit"] = "0.00000"

        assert idict["debye"] == "1.00000"

    if idict["zpp"] != "0.00000":
        process = False

    if process:
        mbr_params[name] = {
            "name": full_name,
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": rfloat(idict["F_0"], 1.0e3),
            "V_0": rfloat(idict["V_0"], 1.0e-6),
            "K_0": rfloat(idict["K_0"], 1.0e9),
            "Kprime_0": rfloat(idict["Kprime_0"]),
            "Debye_0": rfloat(idict["Theta_0"]),
            "grueneisen_0": rfloat(idict["gamma_0"]),
            "q_0": rfloat(idict["q_0"]),
            "G_0": rfloat(idict["G_0"], 1.0e9),
            "Gprime_0": rfloat(idict["Gprime_0"]),
            "eta_s_0": rfloat(idict["eta_s_0"]),
            "n": rfloat(idict["n"]),
            "Z": rfloat(idict["Z"]),
            "molar_mass": rfloat(idict["molar_mass"], 1.0e-3),
        }

        # electronic terms
        if np.abs(rfloat(idict["bel_0"])) > 1.0e-10:
            mbr_params[name]["equation_of_state"] = "slb3-conductive"
            mbr_params[name]["bel_0"] = rfloat(idict["bel_0"])
            mbr_params[name]["gel"] = rfloat(idict["gel"])

        # The positivity constraint for T_crit below
        # is hardcoded into HeFESTo.
        if idict["T_crit"] != "0.00000" and float(idict["T_crit"]) > 0.0:
            if name in magnetic_structural_parameters:
                # The magnetic model in SLB2024 is identical
                # to the chs model already coded in BurnMan
                # but with the zero Gibbs energy and entropy
                # assigned to the fully ordered form,
                # not the disordered form,.
                p = magnetic_structural_parameters[name]["p"]
                Scrit = float(idict["S_crit"])
                magnetic_moment = np.exp(Scrit / gas_constant) - 1.0
                A = (518.0 / 1125.0) + (11692.0 / 15975.0) * ((1.0 / p) - 1.0)
                fT0 = (1.0 / A) * float(idict["T_crit"]) * 79.0 / (140.0 * p)
                E0 = Scrit * fT0

                mbr_mods[name] = [
                    [
                        "magnetic_chs",
                        {
                            "structural_parameter": p,
                            "curie_temperature": [float(idict["T_crit"]), 0.0],
                            "magnetic_moment": [magnetic_moment, 0.0],
                        },
                    ],
                    [
                        "linear",
                        {
                            "delta_E": E0,
                            "delta_S": Scrit,
                            "delta_V": 0.0,
                        },
                    ],
                ]
            else:
                mbr_mods[name] = [
                    [
                        "landau_slb_2022",
                        {
                            "Tc_0": float(idict["T_crit"]),
                            "S_D": float(idict["S_crit"]),
                            "V_D": float(idict["V_crit"]) / 1.0e6,
                        },
                    ]
                ]
        mbr_names[name] = (
            full_name.lower()
            .replace("-", "_")
            .replace(" ", "_")
            .replace("'", "")
            .replace('\\"', "")
        )

# WRITE FILE

print(
    "# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit\n"
    "# for the Earth and Planetary Sciences\n"
    "# Copyright (C) 2012 - 2024 by the BurnMan team,\n"
    "# released under the GNU GPL v2 or later.\n"
    "\n\n"
    '"""\n'
    "SLB_2024\n"
    "Minerals from Stixrude & Lithgow-Bertelloni 2024 and references therein\n"
    f"File autogenerated using HeFESTO_2024_to_burnman.py and {hefesto_path}\n"
    '"""\n\n'
    "from __future__ import absolute_import\n"
    "\n\n"
    "from ..eos.slb import SLB3\n"
    "from ..eos import debye\n"
    "from ..eos import birch_murnaghan as bm\n"
    "from ..classes.mineral import Mineral\n"
    "from ..classes.solution import Solution, RelaxedSolution\n"
    "from ..classes.solutionmodel import (\n"
    "    IdealSolution,\n"
    "    SymmetricRegularSolution,\n"
    "    AsymmetricRegularSolution,\n"
    ")\n\n"
)


st_class = (
    "class SLB3Stishovite(SLB3):\n"
    "    def shear_modulus(self, pressure, temperature, volume, params):\n"
    '        """\n'
    "        Returns the shear modulus. :math:`[Pa]`\n"
    "        SLB2024 first calculate the shear modulus at T_0 using\n"
    "        the third order Birch-Murnaghan equation of state,\n"
    "        then apply the softening relevant to the temperature\n"
    "        and pressure of interest, and then apply the usual thermal term.\n"
    '        """\n'
    "\n"
    "        G_bare = bm.shear_modulus_third_order(volume, params)\n"
    "        Pc = 51.6e9 + 11.1 * (temperature - 300.0) * 1.0e6  # Eq B4\n"
    "        PcQ = (\n"
    "            Pc + 50.7e9\n"
    "        )  # From the HeFESTo code in stishtran.f. PcQ is T-independent in paper.\n"
    "        Cs0 = 128.0e9 * PcQ / Pc\n"
    "        Delta = 0.0  # 100.e9 in the HeFESTo code in stishtran.f, 0. in paper.\n"
    "        PminusPc = pressure - Pc\n"
    "        if PminusPc < 0.0:\n"
    "            Cs = Cs0 * PminusPc / (pressure - PcQ)  # Eq B2\n"
    "        else:\n"
    "            Ast = (Cs0 - Delta) * (PcQ - Pc)\n"
    "            Cs = Cs0 - Ast / (PcQ - pressure + 3.0 * PminusPc)  # Eq B3\n"
    "\n"
    "        # VRH-like average of the bare shear modulus with softened (C11-C12)/2.\n"
    "        G = 0.5 * (13.0 / 15.0 * G_bare + 2.0 / 15.0 * Cs) + 0.5 / (\n"
    "            13.0 / 15.0 / G_bare + 2.0 / 15.0 / Cs\n"
    "        )  # Eq B5\n"
    "\n"
    "        debye_T = self._debye_temperature(params['V_0'] / volume, params)\n"
    "        eta_s = self._isotropic_eta_s(params['V_0'] / volume, params)\n"
    "\n"
    "        E_th = debye.thermal_energy(temperature, debye_T, params['n'])\n"
    "        E_th_ref = debye.thermal_energy(params['T_0'], debye_T, params['n'])\n"
    "\n"
    "        return G - eta_s * (E_th - E_th_ref) / volume\n"
)

fper_relaxed_class = (
    "class ferropericlase_relaxed(RelaxedSolution):\n"
    "    def __init__(self):\n"
    '        """RelaxedSolution model for ferropericlase (mw).\n'
    "        Endmembers (and site species distributions) are given in the order:\n"
    "        - pe ([Mg]2[Mg]2)\n"
    "        - wu and wuls ([Fe,Fels]2[Fe,Fels]2)\n"
    "        - anao ([Na]2[Al]2)\n"
    "        - mag ([Fe]0[Fef]2)\n"
    "        The entropy from the first site in magnetite is not counted.\n"
    "        The equilibrium spin state is calculated automatically.\n"
    '        """\n'
    "        solution = ferropericlase()\n"
    "        vrel = [[0.0, 0.5, -0.5, 0.0, 0.0]]\n"
    "        vunrel = [\n"
    "            [1.0, 0.0, 0.0, 0.0, 0.0],\n"
    "            [0.0, 0.5, 0.5, 0.0, 0.0],\n"
    "            [0.0, 0.0, 0.0, 1.0, 0.0],\n"
    "            [0.0, 0.0, 0.0, 0.0, 1.0],\n"
    "        ]\n"
    "        RelaxedSolution.__init__(self, solution, vrel, vunrel)\n"
)

bdg_relaxed_class = (
    "class bridgmanite_relaxed(RelaxedSolution):\n"
    "    def __init__(self):\n"
    '        """RelaxedSolution model for bridgmanite (pv).\n'
    "        Only the spin transition is relaxed.\n"
    "        Endmembers (and site species distributions) are given in the order:\n"
    "        - mgpv ([Mg][Si])\n"
    "        - fepv ([Fe][Si])\n"
    "        - alpv ([Al][Al])\n"
    "        - hepv + hlpv ([Fef][Fef,Fefls])\n"
    "        - fapv ([Fef][Al])\n"
    "        - crpv ([Cr][Cr])\n"
    "        The equilibrium spin state is calculated automatically.\n"
    '        """\n'
    "        solution = bridgmanite()\n"
    "        vrel = [[0.0, 0., 0., 0.5, -0.5, 0.0, 0.0]]\n"
    "        vunrel = [\n"
    "            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\n"
    "            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],\n"
    "            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],\n"
    "            [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0],\n"
    "            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],\n"
    "            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],\n"
    "        ]\n"
    "        RelaxedSolution.__init__(self, solution, vrel, vunrel)\n"
)

print('"""\n' "ENDMEMBERS\n" '"""\n')
for key, mbr_prm in sorted(mbr_params.items()):
    if key == "st":
        print(st_class)
        mbr_prm = mbr_prm.__repr__().replace("'slb3'", "SLB3Stishovite()")
        mbr_prm = mbr_prm.replace("slb3", "fault")
    print(
        f"\nclass {key}(Mineral):\n"
        "    def __init__(self):\n"
        f"        self.params = {mbr_prm}\n"
    )
    if key in mbr_mods:
        print(
            f"        self.property_modifiers = {mbr_mods[key]}",
        )
        print("")
    print("        Mineral.__init__(self)")
    print("")

print('"""\n' "SOLUTIONS\n" '"""\n')

for key, prm in sorted(sol_params.items()):
    for i in range(prm["n_mbrs"]):
        if prm["mbr_names"][i] in hardcoded_site_occupancies:
            prm["mbr_site_formulae"][i] = hardcoded_site_occupancies[
                prm["mbr_names"][i]
            ]

    docstring = '"""'
    docstring += f'{prm["solution_type"]} model for {solution_aliases[key]} ({key}).\n'
    docstring += "Endmembers (and site species distributions) are given in the order:\n"
    for i in range(prm["n_mbrs"]):
        docstring += f'- {prm["mbr_names"][i]} ({prm["mbr_site_formulae"][i]})\n'
    if key == "mw":
        docstring += "The entropy from the first site in magnetite is not counted.\n"
    docstring += '"""'

    print(
        f"\nclass {solution_aliases[key]}(Solution):\n"
        "    def __init__(self, molar_fractions=None):\n"
        f"        {docstring}\n"
        f'        self.name = "{solution_aliases[key]}"\n'
        f"        self.solution_model = {prm['solution_type']}(\n"
        "            endmembers=["
    )
    for i in range(prm["n_mbrs"]):
        print(
            f'                [{prm["mbr_names"][i]}(), "{prm["mbr_site_formulae"][i]}"],'
        )
    print("            ],")

    if prm["solution_type"] == "AsymmetricRegularSolution":
        print(f"            alphas={prm['alphas']},")

    if prm["solution_type"] != "IdealSolution":
        print(f'            energy_interaction={prm["energy_interaction"]},')
        if "volume_interaction" in prm:
            print(f'            volume_interaction={prm["volume_interaction"]},')

    print("        )\n")
    print("        Solution.__init__(self, molar_fractions=molar_fractions)\n")

    if key == "mw":
        print(fper_relaxed_class)
    if key == "pv":
        print(bdg_relaxed_class)


print('\n"""\n' "ENDMEMBER ALIASES\n" '"""\n')

for key, value in sorted(mbr_names.items()):
    print(f"{value.replace('(', '').replace(')', '')} = {key}")

print("")


print('\n"""\n' "SOLUTION ALIASES\n" '"""\n')

for key, value in sorted(solution_aliases.items()):
    if key != "sp":
        print(f"{key} = {value}")
