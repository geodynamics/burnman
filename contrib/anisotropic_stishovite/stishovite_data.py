import numpy as np
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume


def add_to_data(data, publication, compilation):
    if publication not in data:
        data[publication] = {}
    for name, header, dataset in compilation:
        data[publication][name] = {
            header[i]: dataset[:, i] for i in range(len(dataset[0]))
        }


def get_data():
    data = {}

    # Zhang data
    CT_data = np.loadtxt("data/Zhang_2021_stishovite_elastic_tensor.dat")
    unit_cell_data = np.loadtxt("data/Zhang_2021_stishovite_unit_cell.dat")
    CT_header = [
        "P",
        "P_err",
        "C11",
        "C11_err",
        "C12",
        "C12_err",
        "C13",
        "C13_err",
        "C22",
        "C22_err",
        "C23",
        "C23_err",
        "C33",
        "C33_err",
        "C44",
        "C44_err",
        "C55",
        "C55_err",
        "C66",
        "C66_err",
        "rho",
        "rho_err",
    ]
    unit_cell_header = [
        "P",
        "P_err",
        "a",
        "a_err",
        "b",
        "b_err",
        "c",
        "c_err",
        "V",
        "V_err",
        "rho",
        "rho_err",
    ]

    isstv = unit_cell_data[:, 2] == unit_cell_data[:, 4]

    compilation = [
        ["CT", CT_header, CT_data],
        ["stv", unit_cell_header, unit_cell_data[isstv]],
        ["poststv", unit_cell_header, unit_cell_data[~isstv]],
    ]
    add_to_data(data, "Zhang_2021", compilation)
    for d in [
        data["Zhang_2021"]["CT"],
        data["Zhang_2021"]["stv"],
        data["Zhang_2021"]["poststv"],
    ]:
        d["T"] = 0.0 * d["P"] + 298.15
        d["T_err"] = 0.0 * d["P"] + 5.0

    # Andrault data
    stv_data = np.loadtxt("data/Andrault_2003_stishovite_unit_cell.dat")
    poststv_data = np.loadtxt("data/Andrault_2003_post_stishovite_unit_cell.dat")
    stv_header = ["P", "P_err", "a", "a_err", "c", "c_err", "V", "V_err"]
    poststv_header = [
        "P",
        "P_err",
        "a",
        "a_err",
        "b",
        "b_err",
        "c",
        "c_err",
        "V",
        "V_err",
    ]
    compilation = [
        ["stv", stv_header, stv_data],
        ["poststv", poststv_header, poststv_data],
    ]
    add_to_data(data, "Andrault_2003", compilation)
    d = data["Andrault_2003"]["stv"]
    d["b"] = d["a"]
    d["b_err"] = d["a_err"]
    for d in [data["Andrault_2003"]["stv"], data["Andrault_2003"]["poststv"]]:
        d["T"] = 0.0 * d["P"] + 298.15
        d["T_err"] = 0.0 * d["P"] + 5.0

    # Fischer data (swap a and b)
    F_data = np.genfromtxt("data/Fischer_2018_stishovite.dat", dtype=str)
    file, phase, beamline = F_data.T[:3]
    mask = phase == "stishovite"

    F_data = np.genfromtxt("data/Fischer_2018_stishovite.dat")[:, 3:]
    header = [
        "P",
        "P_err",
        "P_Pt",
        "P_Pt_err",
        "T",
        "T_err",
        "T_KBr",
        "T_KBr_err",
        "a_KBr",
        "a_KBr_err",
        "a_Pt",
        "a_Pt_err",
        "b",
        "b_err",
        "a",
        "a_err",
        "c",
        "c_err",
    ]
    compilation = [["stv", header, F_data[mask]], ["poststv", header, F_data[~mask]]]
    add_to_data(data, "Fischer_2018", compilation)

    for d in [data["Fischer_2018"]["stv"], data["Fischer_2018"]["poststv"]]:
        d["V"] = d["a"] * d["b"] * d["c"]
        d["V_err"] = d["V"] * np.sqrt(
            (
                np.power(d["a_err"] / d["a"], 2.0)
                + np.power(d["b_err"] / d["b"], 2.0)
                + np.power(d["c_err"] / d["c"], 2.0)
            )
        )

    # Ito data (only stishovite)
    I_data = np.loadtxt("data/Ito_1974_stishovite_room_pressure.dat")
    header = ["T", "a", "a_err", "b", "b_err", "c", "c_err", "V", "V_err"]
    compilation = [["stv", header, I_data]]
    add_to_data(data, "Ito_1974", compilation)

    d = data["Ito_1974"]["stv"]
    d["T_err"] = d["T"] * 0.0 + 5.0
    d["P"] = d["T"] * 0.0 + 0.0001
    d["P_err"] = d["T"] * 0.0 + 0.000001

    # Wang data (only stishovite)
    W_data = np.loadtxt("data/Wang_et_al_2012_stv.dat")
    header = [
        "id",
        "T",
        "P",
        "P_err",
        "a",
        "a_err",
        "c",
        "c_err",
        "V",
        "V_err",
        "VAu",
        "VAu_err",
    ]
    compilation = [["stv", header, W_data]]
    add_to_data(data, "Wang_2012", compilation)

    d = data["Wang_2012"]["stv"]
    d["T_err"] = d["T"] / 100.0
    d["b"] = d["a"]
    d["b_err"] = d["a_err"]

    # Nishihara data (only stishovite)
    N_data = np.loadtxt("data/Nishihara_et_al_2005_stv.dat")
    header = ["P", "T", "a", "a_err", "c", "c_err", "V", "V_err"]
    compilation = [["stv", header, N_data]]
    add_to_data(data, "Nishihara_2005", compilation)

    d = data["Nishihara_2005"]["stv"]
    d["T_err"] = d["T"] / 100.0
    d["P_err"] = d["P"] * 0.05 + 0.2  # conservative estimate of the pressure error
    d["b"] = d["a"]
    d["b_err"] = d["a_err"]

    # Convert to SI units
    d = data["Zhang_2021"]["CT"]
    d["P"] = d["P"] * 1.0e9
    d["P_err"] = d["P_err"] * 1.0e9

    f_ln = np.cbrt(molar_volume_from_unit_cell_volume(1.0, 2.0))

    for d in [
        data["Andrault_2003"]["stv"],
        data["Andrault_2003"]["poststv"],
        data["Ito_1974"]["stv"],
        data["Fischer_2018"]["stv"],
        data["Fischer_2018"]["poststv"],
        data["Zhang_2021"]["stv"],
        data["Zhang_2021"]["poststv"],
        data["Wang_2012"]["stv"],
        data["Nishihara_2005"]["stv"],
    ]:
        d["V"] = molar_volume_from_unit_cell_volume(d["V"], 2.0)
        d["V_err"] = molar_volume_from_unit_cell_volume(d["V_err"], 2.0)
        d["P"] = d["P"] * 1.0e9
        d["P_err"] = d["P_err"] * 1.0e9

        for ln in ["a", "b", "c"]:
            d[f"{ln}"] = d[f"{ln}"] * f_ln
            d[f"{ln}_err"] = d[f"{ln}_err"] * f_ln

    # Molar mass (M) is apparently rounded by Zhang to get rho.
    # As rho is a derived value, we correct rho here and
    # return the derived volumes
    M_Zhang = 0.0601
    M = 0.06008
    f = M / M_Zhang
    data["Zhang_2021"]["CT"]["rho"] = 1.0e3 * f * data["Zhang_2021"]["CT"]["rho"]
    data["Zhang_2021"]["CT"]["rho_err"] = (
        1.0e3 * f * data["Zhang_2021"]["CT"]["rho_err"]
    )
    data["Zhang_2021"]["CT"]["V"] = M / data["Zhang_2021"]["CT"]["rho"]
    data["Zhang_2021"]["CT"]["V_err"] = M * (
        data["Zhang_2021"]["CT"]["rho_err"]
        / np.power(data["Zhang_2021"]["CT"]["rho"], 2.0)
    )
    for phase in ["stv", "poststv"]:
        data["Zhang_2021"][phase]["rho"] = 1.0e3 * f * data["Zhang_2021"][phase]["rho"]
        data["Zhang_2021"][phase]["rho_err"] = (
            1.0e3 * f * data["Zhang_2021"][phase]["rho_err"]
        )
        data["Zhang_2021"][phase]["V2"] = M / data["Zhang_2021"][phase]["rho"]
        data["Zhang_2021"][phase]["V2_err"] = M * (
            data["Zhang_2021"][phase]["rho_err"]
            / np.power(data["Zhang_2021"][phase]["rho"], 2.0)
        )

    return data


def common_data():
    # First, let's read in the PVT equation of state data
    data = get_data()

    d = {}
    d["PTV"] = {"stv": {}, "poststv": {}}
    d["PTV_err"] = {"stv": {}, "poststv": {}}
    d["abc"] = {"stv": {}, "poststv": {}}
    d["abc_err"] = {"stv": {}, "poststv": {}}

    for pub in [
        "Ito_1974",
        "Andrault_2003",
        "Nishihara_2005",
        "Wang_2012",
        "Fischer_2018",
        "Zhang_2021",
    ]:
        for phase in ["stv", "poststv"]:
            try:
                Z_data = data[pub][phase]
                d["PTV"][phase][pub] = np.array(
                    [Z_data["P"], Z_data["T"], Z_data["V"]]
                ).T
                d["PTV_err"][phase][pub] = np.array(
                    [Z_data["P_err"], Z_data["T_err"], Z_data["V_err"]]
                ).T
            except:
                pass

            try:
                Z_data = data[pub][phase]
                d["abc"][phase][pub] = np.array(
                    [Z_data["a"], Z_data["b"], Z_data["c"]]
                ).T
                d["abc_err"][phase][pub] = np.array(
                    [Z_data["a_err"], Z_data["b_err"], Z_data["c_err"]]
                ).T
            except:
                pass
    return d


def get_stiffnesses(CN):
    """ """
    n = len(CN["C11"])
    Cs = np.zeros((n, 6, 6))
    Cs_err = np.zeros((n, 6, 6))
    for ij in ["11", "12", "13", "22", "23", "33", "44", "55", "66"]:
        i, j = (int(ij[0]) - 1, int(ij[1]) - 1)
        Cs[:, i, j] = CN[f"C{ij}"]
        Cs[:, j, i] = CN[f"C{ij}"]
        Cs_err[:, i, j] = CN[f"C{ij}_err"]
        Cs_err[:, j, i] = CN[f"C{ij}_err"]
    K_NR = 1.0 / np.sum(np.linalg.inv(Cs)[:, :3, :3], axis=(1, 2))
    return Cs, Cs_err, K_NR


CN_data = get_data()["Zhang_2021"]["CT"]
CN_GPa, CN_err_GPa, KNR_GPa = get_stiffnesses(CN_data)
P_for_CN = CN_data["P"]
T_for_CN = CN_data["T"]
V_for_CN = CN_data["V"]
SN_invGPa = np.linalg.inv(CN_GPa)

# approximate error in K_NR
KNR_err_GPa = np.sqrt(np.sum(np.power(CN_err_GPa[:, :3, :3], 2.0), axis=(1, 2))) / 9.0

beta_NR_invGPa = np.sum(SN_invGPa[:, :3, :3], axis=(1, 2))
SNoverbetaNR_obs = np.einsum("ijk, i->ijk", SN_invGPa, 1.0 / beta_NR_invGPa)
betaNRoverSN_err_obs = np.ones_like(SNoverbetaNR_obs) * 0.2
lnV_for_CN = np.log(V_for_CN)
