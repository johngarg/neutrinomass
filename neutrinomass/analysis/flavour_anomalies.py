from neutrinomass.analysis.utils import PMNS, NEUTRINO_MASSES, VEV, GF, CKM, ALPHA_EM
import numpy as np

import flavio
from wilson import Wilson

_, MV2, MV3 = NEUTRINO_MASSES
_, U2, U3 = np.conjugate(PMNS.cols())


def x(xi, m0):
    prefactor = xi / np.sqrt(2 * m0)
    return (
        prefactor * (U2[0, 0] * np.sqrt(MV2) + 1j * np.sqrt(MV3) * U3[0, 0]),
        prefactor * (U2[0, 1] * np.sqrt(MV2) + 1j * np.sqrt(MV3) * U3[0, 1]),
        prefactor * (U2[0, 2] * np.sqrt(MV2) + 1j * np.sqrt(MV3) * U3[0, 2]),
    )


def Z(xi, m0):
    prefactor = 1 / (np.sqrt(2 * m0) * xi)
    return (
        prefactor * (U2[0, 0] * np.sqrt(MV2) - 1j * np.sqrt(MV3) * U3[0, 0]),
        prefactor * (U2[0, 1] * np.sqrt(MV2) - 1j * np.sqrt(MV3) * U3[0, 1]),
        prefactor * (U2[0, 2] * np.sqrt(MV2) - 1j * np.sqrt(MV3) * U3[0, 2]),
    )


def run_point():
    arg_xi = 1.5
    xi = 1 * (np.cos(arg_xi) + 1j * np.sin(arg_xi))
    m_eta = 5e4
    lambda_kappa = 0.05
    m0 = lambda_kappa * VEV ** 2 / ((16 * np.pi ** 2) ** 2 * m_eta ** 2)
    z23 = 1
    w22 = 1
    w32 = -1

    arg_y33 = 2
    y33 = 2 * (np.cos(arg_y33) + 1j * np.sin(arg_y33))

    m_Pi7 = 1000
    x32 = x(xi, m0)[2]
    clequ1_3332 = np.conjugate(x32) * y33 / (2 * m_Pi7 ** 2)

    m_zeta = 15000
    z22_minus_z23 = Z(xi, m0)[1]
    z22 = (z22_minus_z23 - z23 * w32) / w22
    clq3_2232 = np.conjugate(z22) * z23 / (4 * m_zeta ** 2)

    m_Phi = 2000
    cqu1_2322 = -2 * w22 * w32 / (9 * m_Phi ** 2)

    wsmeft = Wilson(
        {
            "lequ1_3332": clequ1_3332,
            "lequ3_3332": 0.25 * clequ1_3332,
            "lq3_2223": clq3_2232,
            "lq1_2223": 3 * clq3_2232,
            "qu1_2322": cqu1_2322,
            "qu8_2322": -3 * cqu1_2322 / 2,
        },
        scale=(m_Pi7 + m_Phi) / 2,
        eft="SMEFT",
        basis="Warsaw",
    )

    wwet = wsmeft.match_run(scale=4.2, eft="WET", basis="flavio")
    rdstar = flavio.np_prediction("Rtaul(B->D*lnu)", wsmeft)
    rd = flavio.np_prediction("Rtaul(B->Dlnu)", wsmeft)

    print("---------------")
    print("Values")
    print("---------------")

    print(f"lequ1_3332: {clequ1_3332 * 1e6}")
    print(f"C9mumu: {wwet['C9_bsmumu']}")
    print(f"C10mumu: {wwet['C10_bsmumu']}")
    print(f"C9ee: {wwet['C9_bsee']}")
    print(f"C9tautau: {wwet['C9_bstautau']}")
    print(f"RD: {rd}")
    print(f"RD*: {rdstar}")
    print(f"CSL: {wwet['CSL_bctaunutau']}")

    print("---------------")
    print("With parameters")
    print("---------------")
    print(f"xi: {xi}")
    print(f"m_eta: {m_eta/1e3} TeV")
    print(f"m_zeta: {m_zeta/1e3} TeV")
    print(f"m_Pi7: {m_Pi7/1e3} TeV")

    print(f"z22: {z22}")
    print(f"z13: {-Z(xi, m0)[0]}")
    print(f"z23: {z23}")
    print(f"z33: {-Z(xi, m0)[2]}")

    print(f"x12: {x(xi, m0)[0]}")
    print(f"x22: {x(xi, m0)[1]}")
    print(f"x32: {x32}")
    print(f"y33: {y33}")
