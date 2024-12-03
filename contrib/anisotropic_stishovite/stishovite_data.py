import numpy as np
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume
from copy import copy
from burnman.calibrants.tools import pressure_to_pressure
from burnman import calibrants
from scipy.integrate import cumulative_simpson


def add_to_data(data, publication, compilation):
    if publication not in data:
        data[publication] = {}
    for name, header, dataset in compilation:
        data[publication][name] = {
            header[i]: dataset[:, i] for i in range(len(dataset[0]))
        }


Au_Fei = calibrants.Fei_2007.Au()
Pt_Fei = calibrants.Fei_2007.Pt()
Pt_Holmes = calibrants.Holmes_1989.Pt()
Au_Anderson = calibrants.Anderson_1989.Au()
Au_Tsuchiya = calibrants.Tsuchiya_2003.Au()
Pt_Dorogokupets = calibrants.Dorogokupets_2007.Pt()
Pt_Zha = calibrants.Zha_2008.Pt()
Au_Matsui = calibrants.Matsui_2010.Au()
Pt_Matsui = calibrants.Matsui_2009.Pt()
Au_Dorogokupets = calibrants.Dorogokupets_2007.Au()

calibrant_conversions_none = {}

calibrant_conversions_Dorogokupets = {
    "Andrault_2003": [Pt_Holmes, Pt_Dorogokupets],
    "Murakami_2003": [Pt_Holmes, Pt_Dorogokupets],
    "Nishihara_2005": [Au_Anderson, Au_Dorogokupets],
    "Wang_2012": [Au_Tsuchiya, Au_Dorogokupets],
    "Grocholski_2013": [Au_Tsuchiya, Au_Dorogokupets],
    "Fischer_2018": [Pt_Dorogokupets, Pt_Dorogokupets],
    "Zhang_2021": [Au_Fei, Au_Dorogokupets],
    "Zhang_2023": [Au_Fei, Au_Dorogokupets],
}

calibrant_conversions_Holmes_Tsuchiya = {
    "Andrault_2003": [Pt_Holmes, Pt_Holmes],
    "Murakami_2003": [Pt_Holmes, Pt_Holmes],
    "Nishihara_2005": [Au_Anderson, Au_Tsuchiya],
    "Wang_2012": [Au_Tsuchiya, Au_Tsuchiya],
    "Grocholski_2013": [Au_Tsuchiya, Au_Tsuchiya],
    "Fischer_2018": [Pt_Dorogokupets, Pt_Holmes],
    "Zhang_2021": [Au_Fei, Au_Tsuchiya],
    "Zhang_2023": [Au_Fei, Au_Tsuchiya],
}

calibrant_conversions_Matsui = {
    "Andrault_2003": [Pt_Holmes, Pt_Matsui],
    "Murakami_2003": [Pt_Holmes, Pt_Matsui],
    "Nishihara_2005": [Au_Anderson, Au_Matsui],
    "Wang_2012": [Au_Tsuchiya, Au_Matsui],
    "Grocholski_2013": [Au_Tsuchiya, Au_Matsui],
    "Fischer_2018": [Pt_Dorogokupets, Pt_Matsui],
    "Zhang_2021": [Au_Fei, Au_Matsui],
    "Zhang_2023": [Au_Fei, Au_Matsui],
}

calibrant_conversions_Fei = {
    "Andrault_2003": [Pt_Holmes, Pt_Fei],
    "Murakami_2003": [Pt_Holmes, Pt_Fei],
    "Nishihara_2005": [Au_Anderson, Au_Fei],
    "Wang_2012": [Au_Tsuchiya, Au_Fei],
    "Grocholski_2013": [Au_Tsuchiya, Au_Fei],
    "Fischer_2018": [Pt_Dorogokupets, Pt_Fei],
    "Zhang_2021": [Au_Fei, Au_Fei],
    "Zhang_2023": [Au_Fei, Au_Fei],
}

calibrant_conversions = calibrant_conversions_Fei


def get_data(unify_calibrants=True):
    data = {}

    # Zhang et al., 2021 elastic data
    CT_data = np.loadtxt("data/Zhang_2021_stishovite_elastic_tensor.dat")
    idx = np.argsort(CT_data[:, 0])
    CT_data = CT_data[idx]
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

    # Zhang et al., 2021 XRD data
    unit_cell_data = np.loadtxt("data/Zhang_2021_stishovite_unit_cell.dat")
    idx = np.argsort(unit_cell_data[:, 0])
    unit_cell_data = unit_cell_data[idx]
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
    idx = np.argsort(stv_data[:, 0])
    stv_data = stv_data[idx]
    poststv_data = np.loadtxt("data/Andrault_2003_post_stishovite_unit_cell.dat")
    idx = np.argsort(poststv_data[:, 0])
    poststv_data = poststv_data[idx]
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

    # Murakami data
    poststv_data = np.loadtxt("data/Murakami_2003_post_stishovite_unit_cell.dat")
    idx = np.argsort(poststv_data[:, 0])
    poststv_data = poststv_data[idx]
    poststv_header = [
        "runid",
        "P",
        "a",
        "a_err",
        "b",
        "b_err",
        "c",
        "c_err",
        "V",
        "V_err",
    ]
    compilation = [["poststv", poststv_header, poststv_data]]
    add_to_data(data, "Murakami_2003", compilation)
    d = data["Murakami_2003"]["poststv"]
    d["T"] = 0.0 * d["P"] + 298.15
    d["T_err"] = 0.0 * d["P"] + 5.0

    d["V"] = d["V"] * 1.0e-6 / molar_volume_from_unit_cell_volume(1.0, 2.0)
    d["V_err"] = d["V_err"] * 1.0e-6 / molar_volume_from_unit_cell_volume(1.0, 2.0)
    d["P"] = d["P"]
    d["P_err"] = d["P"] * 0.01

    f_ln = np.cbrt(molar_volume_from_unit_cell_volume(1.0, 2.0))
    for ln in ["a", "b", "c"]:
        d[f"{ln}"] = d[f"{ln}"]  # * f_ln
        d[f"{ln}_err"] = d[f"{ln}_err"]  # * f_ln

    # Grocholski data
    all_data = np.loadtxt(
        "data/Grocholski_2013_stishovite_post_stishovite_unit_cell.dat"
    )
    idx = np.argsort(all_data[:, 0])
    all_data = all_data[idx]
    stv_header = [
        "P",
        "P_err",
        "a_p",
        "a_p_err",
        "a",
        "a_err",
        "b",
        "b_err",
        "c",
        "c_err",
        "V",
        "V_err",
    ]
    stv_mask = all_data[:, 0] < 50.0
    compilation = [
        ["stv", stv_header, copy(all_data[stv_mask])],
        ["poststv", stv_header, copy(all_data[~stv_mask])],
    ]
    add_to_data(data, "Grocholski_2013", compilation)
    for d in [data["Grocholski_2013"]["stv"], data["Grocholski_2013"]["poststv"]]:
        d["T"] = 0.0 * d["P"] + 298.15
        d["T_err"] = 0.0 * d["P"] + 5.0

        # Volumes are reported as molar volumes, here we change to unit cell volumes
        # We correct back later, with all the other datasets
        d["V"] = d["V"] * 1.0e-6 / molar_volume_from_unit_cell_volume(1.0, 2.0)
        d["V_err"] = d["V_err"] * 1.0e-6 / molar_volume_from_unit_cell_volume(1.0, 2.0)

    # Zhang data
    all_data = np.loadtxt("data/Zhang_2023_stishovite_unit_cell.dat")
    idx = np.argsort(all_data[:, 0])
    all_data = all_data[idx]
    stv_header = [
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
        "unique_refl",
        "Rint",
        "R1",
    ]
    stv_mask = all_data[:, 0] < 50.0
    compilation = [
        ["stv", stv_header, copy(all_data[stv_mask])],
        ["poststv", stv_header, copy(all_data[~stv_mask])],
    ]
    add_to_data(data, "Zhang_2023", compilation)
    for d in [data["Zhang_2023"]["stv"], data["Zhang_2023"]["poststv"]]:
        d["T"] = 0.0 * d["P"] + 298.15
        d["T_err"] = 0.0 * d["P"] + 5.0

    # Sun et al., 2019 data
    S_data = np.loadtxt("data/Sun_et_al_2019_post_stishovite_volumes.dat")
    idx = np.argsort(S_data[:, 0])
    S_data = S_data[idx]
    poststv_header = ["P", "P_err", "T", "T_err", "V", "V_err", "V_NaCl", "V_NaCl_err"]
    compilation = [
        ["poststv", poststv_header, S_data],
    ]
    add_to_data(data, "Sun_2019", compilation)

    # Fischer et al., 2018 data (swap a and b)
    F_data = np.genfromtxt("data/Fischer_2018_stishovite.dat", dtype=str)
    file, phase, beamline = F_data.T[:3]
    mask = phase == "stishovite"

    F_data = np.genfromtxt("data/Fischer_2018_stishovite.dat")[:, 3:]
    idx = np.argsort(F_data[:, 0])
    F_data = F_data[idx]

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

    # Wang et al., 2012 data (only stishovite)
    W_data = np.loadtxt("data/Wang_et_al_2012_stv.dat")
    idx = np.argsort(W_data[:, 2])
    W_data = W_data[idx]
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

    # Nishihara et al., 2005 data (only stishovite)
    N_data = np.loadtxt("data/Nishihara_et_al_2005_stv.dat")
    idx = np.argsort(N_data[:, 0])
    N_data = N_data[idx]

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

    # Check that P and V are in the correct units
    # Unify calibrants if desired
    for study in [
        "Andrault_2003",
        "Murakami_2003",
        "Nishihara_2005",
        "Wang_2012",
        "Grocholski_2013",
        "Fischer_2018",
        "Sun_2019",
        "Zhang_2021",
        "Zhang_2023",
    ]:

        for phase in ["stv", "poststv"]:
            try:
                P = data[study][phase]["P"]
                T = data[study][phase]["T"]
                V = data[study][phase]["V"]
                if max(P) > 1.0e6:
                    print("P too small", study)
                    exit()
                if min(V) < 1.0:
                    print("V too big", study)
                    exit()
                if unify_calibrants and study in calibrant_conversions:
                    print(f"Changing calibrant for {study} data ({phase})")
                    P = [
                        pressure_to_pressure(
                            calibrant_conversions[study][0],
                            calibrant_conversions[study][1],
                            P[i] * 1.0e9,
                            T[i],
                        )
                        for i in range(len(P))
                    ]
                    data[study][phase]["P"] = np.array(P) / 1.0e9

            except KeyError:
                pass

    for study, phase in [
        ("Ito_1974", "stv"),
        ("Andrault_2003", "stv"),
        ("Andrault_2003", "poststv"),
        ("Murakami_2003", "poststv"),
        ("Nishihara_2005", "stv"),
        ("Wang_2012", "stv"),
        ("Grocholski_2013", "stv"),
        ("Grocholski_2013", "poststv"),
        ("Fischer_2018", "stv"),
        ("Fischer_2018", "poststv"),
        ("Sun_2019", "poststv"),
        ("Zhang_2021", "stv"),
        ("Zhang_2021", "poststv"),
        ("Zhang_2023", "stv"),
        ("Zhang_2023", "poststv"),
    ]:
        d = data[study][phase]

        d["V"] = molar_volume_from_unit_cell_volume(d["V"], 2.0)
        d["V_err"] = molar_volume_from_unit_cell_volume(d["V_err"], 2.0)
        d["P"] = d["P"] * 1.0e9
        d["P_err"] = d["P_err"] * 1.0e9

        try:
            for ln in ["a", "b", "c"]:
                d[f"{ln}"] = d[f"{ln}"] * f_ln
                d[f"{ln}_err"] = d[f"{ln}_err"] * f_ln
        except KeyError as e:
            if study == "Sun_2019":
                pass
            else:
                print(study, phase)
                print(e)

    # Molar mass (M) is apparently rounded by Zhang_2021 to get rho.
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

    #### TODO Add Zhang 2023 here?
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
            except KeyError:
                pass

            try:
                Z_data = data[pub][phase]
                d["abc"][phase][pub] = np.array(
                    [Z_data["a"], Z_data["b"], Z_data["c"]]
                ).T
                d["abc_err"][phase][pub] = np.array(
                    [Z_data["a_err"], Z_data["b_err"], Z_data["c_err"]]
                ).T
            except KeyError:
                pass
    return d


def other_data():
    # First, let's read in the PVT equation of state data
    data = get_data()

    d = {}
    d["PTV"] = {"stv": {}, "poststv": {}}
    d["PTV_err"] = {"stv": {}, "poststv": {}}
    d["abc"] = {"stv": {}, "poststv": {}}
    d["abc_err"] = {"stv": {}, "poststv": {}}

    for pub in [
        "Murakami_2003",
        "Grocholski_2013",
        "Zhang_2023",
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
            except KeyError:
                pass

            try:
                Z_data = data[pub][phase]
                d["abc"][phase][pub] = np.array(
                    [Z_data["a"], Z_data["b"], Z_data["c"]]
                ).T
                d["abc_err"][phase][pub] = np.array(
                    [Z_data["a_err"], Z_data["b_err"], Z_data["c_err"]]
                ).T
            except KeyError:
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
CN_err_GPa[CN_err_GPa < 1.0] = 1.0

KTR_GPa = KNR_GPa
P_for_CN = CN_data["P"]
T_for_CN = CN_data["T"]
V_for_CN = CN_data["V"]
SN_invGPa = np.linalg.inv(CN_GPa)
P_for_CN_new = cumulative_simpson(KTR_GPa * 1.0e9, x=-np.log(V_for_CN), initial=0)

# approximate error in K_NR
KNR_err_GPa = np.sqrt(np.sum(np.power(CN_err_GPa[:, :3, :3], 2.0), axis=(1, 2))) / 9.0

beta_NR_invGPa = np.sum(SN_invGPa[:, :3, :3], axis=(1, 2))
SNoverbetaNR_obs = np.einsum("ijk, i->ijk", SN_invGPa, 1.0 / beta_NR_invGPa)
betaNRoverSN_err_obs = np.ones_like(SNoverbetaNR_obs) * 0.2
lnV_for_CN = np.log(V_for_CN)

# Jiang data
J2009_data = np.loadtxt("data/Jiang_2009_elastic.dat", skiprows=1)
J2009_stiffness_matrices = np.zeros((10, 6, 6))
J2009_pressures = J2009_data[:, 0] * 1.0e9
J2009_temperatures = J2009_pressures * 0.0 + 298.15
for i in range(10):
    # C11, C33, C13, C12, C44, C66
    J2009_stiffness_matrices[i] = [
        [J2009_data[i, 1], J2009_data[i, 4], J2009_data[i, 3], 0, 0, 0],
        [J2009_data[i, 4], J2009_data[i, 1], J2009_data[i, 3], 0, 0, 0],
        [J2009_data[i, 3], J2009_data[i, 3], J2009_data[i, 2], 0, 0, 0],
        [0, 0, 0, J2009_data[i, 5], 0, 0],
        [0, 0, 0, 0, J2009_data[i, 5], 0],
        [0, 0, 0, 0, 0, J2009_data[i, 6]],
    ]
J2009_stiffness_matrices = J2009_stiffness_matrices * 1.0e9
J2009_compliance_matrices = np.linalg.inv(J2009_stiffness_matrices)
J2009_KRT = 1.0 / np.sum(J2009_compliance_matrices[:, :3, :3], axis=(1, 2))
J2009_dPsidfs = np.einsum("i, ijk->ijk", J2009_KRT, J2009_compliance_matrices)

C2024 = np.loadtxt("data/Chen_2024_stishovite_velocities.dat")
C2024_pressures = C2024[:, 0] * 1.0e9
C2024_pressures_unc = C2024[:, 1] * 1.0e9
C2024_temperatures = C2024[:, 2] + 273.15
C2024_temperatures_unc = C2024[:, 3] + 273.15
C2024_Vp = C2024[:, 4] * 1.0e3
C2024_Vp_unc = C2024[:, 5] * 1.0e3
C2024_Vs = C2024[:, 6] * 1.0e3
C2024_Vs_unc = C2024[:, 7] * 1.0e3
C2024_Vphi = (
    np.sqrt(np.power(C2024[:, 4], 2.0) - 4.0 / 3.0 * np.power(C2024[:, 6], 2.0)) * 1.0e3
)
C2024_Vphi_unc = (
    (1 / C2024_Vphi)
    * np.sqrt(
        (C2024[:, 4] * C2024[:, 5]) ** 2
        + ((4.0 / 3.0) * C2024[:, 6] * C2024[:, 7]) ** 2
    )
    * 1.0e6
)
C2024_KS = C2024[:, 10] * 1.0e9
C2024_KS_unc = C2024[:, 11] * 1.0e9
C2024_G = C2024[:, 12] * 1.0e9
C2024_G_unc = C2024[:, 13] * 1.0e9
