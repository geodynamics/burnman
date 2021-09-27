import numpy as np
from tabulate import tabulate

def print_table_for_mineral_constants(mineral, indices):

    constants = []
    for (i, j) in indices:
        constants.append(mineral.c[i-1, j-1, :, :])

    constants = np.array(constants)

    mn_pairs = []
    for n in range(constants.shape[2]):
        for m in range(constants.shape[1]):
            if not np.all(constants[:, m, n] == 0):
                mn_pairs.append((m, n))

    rows = [[f'$c_{{ij{m}{n}}}$' for (m, n) in mn_pairs]]
    for ci, (i, j) in enumerate(indices):
        row = [f'$c_{{{i}{j}}}$']
        row.extend([f'{constants[ci, m, n]:.4e}' for (m, n) in mn_pairs])
        rows.append(row)

    print(tabulate(rows, headers='firstrow', tablefmt='latex_raw'))
