import pandas as pd
import numpy as np
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume


def get_data():
    d = pd.read_csv(
        "data/plagioclase_stiffness_tensor_Brown_2016.dat",
        comment="#",
        header=0,
        index_col=0,
        sep="\\s+",
    )

    rows = list(d.index[1:])
    d = np.array(d.to_numpy()[1:], dtype=float)
    p_an = d[0, ::2] / 100.0  # fractions

    CN = np.empty((len(p_an), 6, 6))
    CN_err = np.empty((len(p_an), 6, 6))
    for irow in range(1, 22):
        name = rows[irow]
        i, j = int(name[1]) - 1, int(name[2]) - 1
        CN[:, i, j] = d[irow, ::2]
        CN[:, j, i] = d[irow, ::2]
        CN_err[:, i, j] = d[irow, 1::2]
        CN_err[:, j, i] = d[irow, 1::2]

    SN = np.linalg.inv(CN)
    beta_RN = np.sum(SN[:, :3, :3], axis=(1, 2))
    KRN = 1.0 / beta_RN
    psiN = np.einsum("ijk, i->ijk", SN, KRN)
    KRN_err = np.sqrt(np.sum(np.power(CN_err[:, :3, :3], 2.0), axis=(1, 2))) / 9.0
    data = {
        "CN": {
            "P": p_an * 0.0 + 1.0e5,
            "T": p_an * 0.0 + 298.15,
            "p_an": p_an,
            "psiN": psiN,
            "SN": SN / 1.0e9,
            "CN": CN * 1.0e9,
            "CN_err": CN_err * 1.0e9,
            "KRN": KRN * 1.0e9,
            "KRN_err": KRN_err,
        }
    }

    d = pd.read_csv(
        "data/plagioclase_cell_tensor_Brown_2016.dat",
        comment="#",
        header=0,
        index_col=0,
        sep="\\s+",
    ).to_dict("list")
    data["cell"] = {k: np.array(v) for k, v in d.items()}

    # convert units
    f_V = molar_volume_from_unit_cell_volume(1.0, 4.0)
    f_ln = np.cbrt(f_V)
    for prp in ["a", "a_err", "b", "b_err", "c", "c_err"]:
        data["cell"][prp] = data["cell"][prp] * f_ln
    for prp in ["V", "V_err"]:
        data["cell"][prp] = data["cell"][prp] * f_V
    for prp in ["rho", "rho_err"]:
        data["cell"][prp] = data["cell"][prp] * 1.0e6

    mask = data["cell"]["Z"] > 6.0
    data["cell"]["c"][mask] = data["cell"]["c"][mask] / 2.0
    data["cell"]["c_err"][mask] = data["cell"]["c_err"][mask] / 2.0
    data["cell"]["V"][mask] = data["cell"]["V"][mask] / 2.0
    data["cell"]["V_err"][mask] = data["cell"]["V_err"][mask] / 2.0

    data["cell"]["P"] = data["cell"]["p_an"] * 0.0 + 1.0e5
    data["cell"]["T"] = data["cell"]["p_an"] * 0.0 + 298.15

    d = np.loadtxt("data/plagioclase_compressibilities_Brown_2016.dat")[:, 1:]

    b = np.array([d[:, 2 * i + 1] for i in range(6)]) / 1.0e9
    b_err = np.array([d[:, 2 * i + 2] for i in range(6)]) / 1.0e9
    b = np.moveaxis(b, 1, 0)
    b_err = np.moveaxis(b_err, 1, 0)

    # Brown reports the sum in Voigt form (Section 3.1)
    b[:, 3:] = b[:, 3:] / 2
    b_err[:, 3:] = b_err[:, 3:] / 2

    data["beta"] = {
        "P": d[:, 0] * 0.0 + 1.0e5,
        "T": d[:, 0] * 0.0 + 298.15,
        "p_an": d[:, 0] / 100.0,
        "b": b,
        "b_err": b_err,
        "bTR": np.sum(b[:, :3], axis=1),
        "bTR_err": np.sqrt(np.sum(np.power(b_err[:, :3], 2.0), axis=1)),
    }

    return data


def get_C1_data():
    data = get_data()
    for dataset, d in data.items():
        if dataset == "beta":
            mask = d["p_an"] < 0.47
        else:
            mask = d["p_an"] < 0.50
        data[dataset] = {k: v[mask] for k, v in d.items()}
    return data


if __name__ == "__main__":
    data = get_C1_data()
    print(data)
