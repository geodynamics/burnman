# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This is a script that converts a HPx-eos data format into the standard burnman
# format (printed to stdout). For example, the dataset file
# tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASHTOCr.txt
# can be processed  via the commands
# python hpx_eos_to_burnman.py > ../../minerals/ig50NCKFMASHTOCr.py; black ../../minerals/ig50NCKFMASHTOCr.py
# the text "from . import ig50NCKFMASHTOCr must then be added to the minerals __init__.py

import numpy as np
import logging
import sys
from sympy import Symbol, prod, sympify
from sympy.parsing.sympy_parser import parse_expr
from burnman.constants import gas_constant
from burnman.minerals import HGP_2018_ds633, HP_2011_ds62
from datetime import date

ds = [
    ["tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASHTOCr.txt", "HGP_2018_ds633"],
    ["tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASTOCr.txt", "HGP_2018_ds633"],
    ["tc-thermoinput-metabasite-2022-01-30/tc-mb50NCKFMASHTO.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50KFMASH.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50MnNCKFMASHTO.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50NCKFMASHTO.txt", "HP_2011_ds62"],
]

ignore_solutions = ["liq", "L", "fl"]  # currently unreadable Temkin models

# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
current_year = str(date.today())[:4]


def reverse_polish(lines):
    # get symbols
    syms = []
    for line in lines:
        n_terms = int(line[0])
        i = 2
        for term in range(n_terms):
            n_add = int(line[i])
            syms.extend(line[i + 2 : i + 2 * n_add + 2][::2])
            i += 2 * n_add + 2
    unique_syms = list(np.unique(sorted(syms)))
    sympy_syms = {u: Symbol(u) for u in unique_syms}

    # get expressions
    expr = 0.0
    for line in lines:
        n_terms = int(line[0])
        i = 2
        terms = []
        for term in range(n_terms):
            n_add = int(line[i])
            en = parse_expr(line[i - 1])  # this is the constant
            for j in range(n_add):
                en += parse_expr(line[i + 1 + 2 * j]) * sympy_syms[line[i + 2 + 2 * j]]
            terms.append(en)
            # syms.extend(l[i+2:i+2*n_add + 2][::2])
            # print("a", n_add)
            i += 2 * n_add + 2
        product = prod(terms)
        expr += product
    return expr


def ordering_modifier(property_modifier, ordered):
    mtype, params = property_modifier
    if mtype == "bragg_williams":
        if ordered:
            return np.array([0.0, 0.0, 0.0])
        else:
            n = params["n"]
            if params["factor"] > 0.0:
                f = [params["factor"], params["factor"]]
            else:
                f = [1.0, -params["factor"]]

            S_disord = (
                -gas_constant
                * (
                    f[0] * (np.log(1.0 / (n + 1.0)) + n * np.log(n / (n + 1.0)))
                    + f[1]
                    * (n * np.log(1.0 / (n + 1.0)) + n * n * np.log(n / (n + 1.0)))
                )
                / (n + 1.0)
            )
            return np.array([params["deltaH"], S_disord, params["deltaV"]])
    else:
        if params["T_0"] < params["Tc_0"]:
            Q_0 = np.power((params["Tc_0"] - params["T_0"]) / params["Tc_0"], 0.25)
        else:
            Q_0 = 0.0

        if ordered:
            E_ord = (
                params["Tc_0"] * params["S_D"] * (Q_0 * Q_0 - np.power(Q_0, 6.0) / 3.0)
                - (
                    params["S_D"] * params["Tc_0"] * 2.0 / 3.0
                    - params["V_D"] * params["P_0"]
                )
                - params["P_0"] * params["V_D"] * Q_0 * Q_0
            )
            S_ord = params["S_D"] * (Q_0 * Q_0 - 1.0)
            V_ord = params["V_D"] * (Q_0 * Q_0 - 1.0)
            return np.array([E_ord, S_ord, V_ord])

        else:
            E_disord = (
                params["Tc_0"] * params["S_D"] * (Q_0 * Q_0 - np.power(Q_0, 6.0) / 3.0)
                - params["P_0"] * params["V_D"] * Q_0 * Q_0
            )
            S_disord = params["S_D"] * Q_0 * Q_0
            V_disord = params["V_D"] * Q_0 * Q_0
            return np.array([E_disord, S_disord, V_disord])


for solution_file, mbr_dataset in ds:
    logging.debug(f"{solution_file}, {mbr_dataset}")
    dataset = eval(mbr_dataset)
    out_ss = ""
    out_make = ""

    i = 0
    data = []
    with open(solution_file, "r") as file:
        # Read each line in the file one by one
        process = True
        for line in file:
            # Process the line (for example, print it)
            ln = line.strip().split()
            if len(ln) > 0:
                if ln[0] == "verbatim" or ln[0] == "header":
                    process = not process
                if (
                    ln[0][0] != "%"
                    and ln[0][0] != "*"
                    and ln[0] != "header"
                    and ln[0] != "verbatim"
                    and process
                ):
                    data.append(ln)

        start_indices = []
        for i_d, d in enumerate(data):
            if len(d) == 3:
                if d[0] != "check" and d[0] != "%" and d[0] not in ignore_solutions:
                    if int(d[1]) > 1.1:
                        start_indices.append(i_d)

    n_data = len(data)

    noods = []
    for ind in start_indices:
        i = ind
        name = data[i][0]
        n_mbrs = int(data[i][1])
        logging.debug(f"Solution: {name}")
        logging.debug(f"Number of endmembers: {n_mbrs}")

        vars = data[i + 1 : i + n_mbrs]
        logging.debug(vars)
        i += n_mbrs

        mbr_proportions = []
        mbr_names = []
        for j in range(n_mbrs):
            n_lines = int(data[i][1])
            mbr_names.append(data[i][0].split("(")[1][:-1])
            mbr_proportions.append(data[i : i + n_lines])
            i += n_lines
        logging.debug(mbr_proportions)

        formulation = data[i][0]
        logging.debug(f"formulation: {formulation}")
        i += 1
        n_ints = int((n_mbrs * (n_mbrs - 1.0)) / 2.0)
        interactions = data[i : i + n_ints]
        logging.debug(interactions)

        nonideal_entropies = False
        nonideal_volumes = False

        We = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        Ws = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        Wv = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        for ints in interactions:
            m = ints[0].split("(")[1].replace(")", "").split(",")
            m0 = mbr_names.index(m[0])
            m1 = mbr_names.index(m[1])
            We[m0][m1 - m0 - 1] = float(ints[1]) * 1.0e3
            Ws[m0][m1 - m0 - 1] = -float(ints[2])
            Wv[m0][m1 - m0 - 1] = float(ints[3]) * 1.0e-5

            if np.abs(Ws[m0][m1 - m0 - 1]) > 1.0e-10:
                nonideal_entropies = True
            if np.abs(Wv[m0][m1 - m0 - 1]) > 1.0e-10:
                nonideal_volumes = True

        i += n_ints

        if formulation == "asf":
            text_alphas = data[i : i + n_mbrs]
            logging.debug(text_alphas)
            i += n_mbrs

            Ae = [0.0 for i in range(n_mbrs)]
            As = [0.0 for i in range(n_mbrs)]
            Av = [0.0 for i in range(n_mbrs)]
            for alpha in text_alphas:
                try:
                    m = alpha[0].split("(")[1].replace(")", "").split(",")[0]
                except IndexError:
                    m = alpha[0]
                m0 = mbr_names.index(m)
                Ae[m0] = float(alpha[1])
                As[m0] = -float(alpha[2])
                Av[m0] = float(alpha[3]) * 1.0e6

        n_site_species = int(data[i][0])

        logging.debug(f"Number of site species: {n_site_species}")
        i += 1

        site_fractions = []
        for j in range(n_site_species):
            n_lines = int(data[i][1])
            site_fractions.append(data[i : i + n_lines])
            i += n_lines

        summed = 0.0
        sites = [[]]
        for site_species_and_expression in site_fractions:
            site_species = site_species_and_expression[0][0]
            expression = site_species_and_expression
            expression[0] = expression[0][2:]
            expression = reverse_polish(expression)
            summed += expression
            summed = sympify(summed)
            sites[-1].append(site_species)
            if summed == 1:
                summed = 0
                sites.append([])
        sites = sites[:-1]
        site_species = {}
        for n_site, site in enumerate(sites):
            for species in site:
                site_species[species] = n_site
        n_sites = len(sites)

        endmember_site_fractions_and_checks = []
        for j in range(n_mbrs):
            n_lines = 1
            while True:
                if (
                    data[i + n_lines][0] == "make"
                    or data[i + n_lines][0] == "delG(tran)"
                    or data[i + n_lines][0] == "check"
                    or data[i + n_lines][0] == "delG(od)"
                    or data[i + n_lines][0] == "delG(make)"
                    or data[i + n_lines][0] == "delG(mod)"
                    or data[i + n_lines][0] == "delG(rcal)"
                    or data[i + n_lines][0] == "DQF"
                ):
                    n_lines += 1
                else:
                    break
            endmember_site_fractions_and_checks.append(data[i : i + n_lines])
            i += n_lines

        mbr_initializations = []
        for f in endmember_site_fractions_and_checks:
            mbr = f[0][0]
            if len(f) == 1:
                mbr = f"{mbr_dataset}.{mbr}()"
            elif len(f) == 2:
                mbr = f"{mbr_dataset}.{mbr}()"
                assert f[1][0] == "check"
            else:
                make_idx = [g[0] for g in f].index("make")
                make = f[make_idx][1:]
                delG = [np.array(g[1:4], dtype=float) for g in f if "delG" in g[0]]
                delG = np.array([0.0, 0.0, 0.0]) if len(delG) == 0 else delG[0]
                delG *= [1.0e3, -1.0e3, 1.0e-5]
                n_make = int(make[0])
                el = 1
                combined_mbrs = "["
                combined_amounts = "["
                for j in range(n_make):
                    if make[el] == "ordered":
                        el += 1
                        m_string = make[el]
                        m_amount = float(parse_expr(make[el + 1]))
                        mineral = getattr(dataset, m_string)()
                        property_modifier = mineral.property_modifiers[0]
                        delG += ordering_modifier(property_modifier, True) * m_amount
                        combined_mbrs += f"{make[el]}_nood, "
                        noods.append(make[el])
                    elif make[el] == "disordered":
                        el += 1
                        m_string = make[el]
                        m_amount = float(parse_expr(make[el + 1]))
                        mineral = getattr(dataset, m_string)()
                        property_modifier = mineral.property_modifiers[0]
                        delG += ordering_modifier(property_modifier, False) * m_amount
                        combined_mbrs += f"{make[el]}_nood, "
                        noods.append(make[el])
                    elif make[el] == "equilibrium":
                        el += 1
                        combined_mbrs += f"{mbr_dataset}.{make[el]}(), "
                    else:
                        combined_mbrs += f"{mbr_dataset}.{make[el]}(), "
                    combined_amounts += f"{float(parse_expr(make[el+1]))}, "
                    el += 2
                combined_mbrs = combined_mbrs[:-2]
                combined_amounts = combined_amounts[:-2]
                combined_mbrs += "]"
                combined_amounts += "]"
                out_make += f'{mbr} = CombinedMineral({combined_mbrs}, {combined_amounts}, {list(delG)}, "{mbr}")\n'
            n_occs = int(f[0][2])
            site_occ = [[] for ln in range(n_sites)]
            site_occ_list = f[0][3:]
            ss = [[] for idx in range(n_sites)]
            site_amounts = [[] for idx in range(n_sites)]
            for j in range(n_occs):
                site_sp = site_occ_list[2 * j]
                i_site = site_species[site_sp]
                # rename for insertion into string
                if site_sp[0] == "x":
                    site_sp = site_sp[1:]
                    site_sp = "".join(filter(str.isalpha, site_sp)).title()

                ss[i_site].append(site_sp)
                site_amounts[i_site].append(parse_expr(site_occ_list[2 * j + 1]))
            total_site_amounts = [sum(amounts) for amounts in site_amounts]
            site_fractions = [
                [amount / total_site_amounts[j] for amount in amounts]
                for j, amounts in enumerate(site_amounts)
            ]

            site_formula = ""
            for m, site in enumerate(ss):
                site_formula += "["
                for n, species in enumerate(site):
                    site_formula += species
                    if site_fractions[m][n] != 1:
                        site_formula += str(site_fractions[m][n])
                site_formula += "]"
                if str(total_site_amounts[m]) != "1":
                    site_formula += str(total_site_amounts[m])

            mbr_initializations.append(f'[{mbr}, "{site_formula}"]')

        out_ss += f"class {name}(Solution):\n"
        out_ss += "    def __init__(self, molar_fractions=None):\n"
        out_ss += f'        self.name = "{name}"\n'
        if formulation == "asf":
            out_ss += "        self.solution_model = AsymmetricRegularSolution(\n"
        else:
            out_ss += "        self.solution_model = SymmetricRegularSolution(\n"
        out_ss += "            endmembers=[\n"
        for mbr_initialization in mbr_initializations:
            out_ss += f"                {mbr_initialization},\n"
        out_ss += "            ],\n"
        if formulation == "asf":
            out_ss += f"            alphas={Ae},\n"
        out_ss += f"            energy_interaction={We},\n"
        if nonideal_entropies:
            out_ss += f"            entropy_interaction={Ws},\n"
        if nonideal_volumes:
            out_ss += f"            volume_interaction={Wv},\n"

        out_ss += "        )\n"
        out_ss += "        Solution.__init__(self, molar_fractions=molar_fractions)\n\n"

        logging.debug(site_fractions)
        logging.debug(endmember_site_fractions_and_checks)
        logging.debug("")

    out_noods = ""
    noods = np.unique(sorted(noods))
    for nood in noods:
        mineral = getattr(dataset, nood)()
        params = mineral.params
        out_noods += f"{nood}_nood = Mineral({params})\n\n"

    # Preamble
    solution_dataset = solution_file.split("-")[-1].split(".")[0]
    underline = "^" * len(solution_dataset)
    preamble = (
        "# This file is part of BurnMan - a thermoelastic\n"
        "# and thermodynamic toolkit for the Earth and Planetary Sciences\n"
        f"# Copyright (C) 2012 - {current_year} by the BurnMan team, released under the GNU\n"
        "# GPL v2 or later.\n"
        "\n"
        '"""\n'
        f"{solution_dataset}\n"
        f"{underline}\n"
        "\n"
        "HPx-eos solutions using endmembers from\n"
        f"dataset {mbr_dataset}.\n"
        "The values in this document are all in S.I. units,\n"
        "unlike those in the original THERMOCALC file.\n"
        "This file is autogenerated using process_HPX_eos.py\n"
        '"""\n\n'
        f"from numpy import array, nan\n"
        f"from . import {mbr_dataset}\n"
        "from ..classes.mineral import Mineral\n"
        "from ..classes.solution import Solution\n"
        "from ..classes.solutionmodel import SymmetricRegularSolution\n"
        "from ..classes.solutionmodel import AsymmetricRegularSolution\n"
        "from ..classes.combinedmineral import CombinedMineral\n\n"
    )

    output_string = f"{preamble}\n{out_noods}\n{out_make}\n{out_ss}"

    output_file = f"../../minerals/{solution_dataset}.py"

    with open(output_file, "w") as text_file:
        print(f"from . import {solution_dataset}")
        text_file.write(output_string)

print("Copy the above to ../../minerals/__init__.py:")
