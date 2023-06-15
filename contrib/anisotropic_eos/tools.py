import numpy as np
from tabulate import tabulate


def print_table_for_mineral_constants(mineral, indices):
    constants = []
    for i, j in indices:
        constants.append(mineral.anisotropic_params["c"][i - 1, j - 1, :, :])

    constants = np.array(constants)

    mn_pairs = []
    for n in range(constants.shape[2]):
        for m in range(constants.shape[1]):
            if not np.all(constants[:, m, n] == 0):
                mn_pairs.append((m, n))

    rows = [[f"$c_{{pq{m}{n}}}$" for (m, n) in mn_pairs]]
    for ci, (i, j) in enumerate(indices):
        row = [f"$c_{{{i}{j}}}$"]
        row.extend([f"{constants[ci, m, n]:.4e}" for (m, n) in mn_pairs])
        rows.append(row)

    print(tabulate(rows, headers="firstrow", tablefmt="latex_raw"))


def print_table_for_mineral_constants_2(mineral, param_list, indices):
    constants = []
    for i, j in indices:
        cs = []
        for param in param_list:
            cs.append(mineral.anisotropic_params[param][i - 1, j - 1])
        constants.append(cs)
    constants = np.array(constants)

    param_list = [p + "_" for p in param_list]
    rows = [["$p$", "$q$"]]
    rows[0].extend(
        [f'${param.split("_")[0]}_{{{param.split("_")[1]}pq}}$' for param in param_list]
    )
    for ci, (i, j) in enumerate(indices):
        row = [f"{i}", f"{j}"]
        row.extend(
            [
                f"{constants[ci, i]:.4e}"
                if constants[ci, i] != 0 and constants[ci, i] != 1
                else "-"
                for i in range(len(param_list))
            ]
        )
        rows.append(row)

    print(tabulate(rows, headers="firstrow", tablefmt="latex_raw"))
