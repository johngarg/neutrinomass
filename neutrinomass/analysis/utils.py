#!/usr/bin/env python3

import sympy as sym
import numpy as np
import functools as ft

# CONSTANTS
pert_bound = np.sqrt(4.0 * np.pi)
GF = 1.16638e-5  # Inverse GeV squared
ALPHA_EM = 1.0 / 127.0  # EM alpha at M_z
mb = 4.18  # Mass of the b-quark in GeV


## Measured neutrino oscillation params from NuFit (http://www.nu-fit.org/?q=node/166)
## arXiv:1611.01514,  NuFIT 3.2 (2018), www.nu-fit.org
THETA12 = np.deg2rad(33.82)
THETA23 = np.deg2rad(49.6)
THETA13 = np.deg2rad(8.61)
DELTA_CP = np.deg2rad(215.0)

# classes
class UnitaryMatrix(object):
    """
    A unitary matrix, by default the PMNS matrix.
    """

    def __init__(
        self,
        c12=np.cos(THETA12),
        c13=np.cos(THETA13),
        c23=np.cos(THETA23),
        s12=np.sin(THETA12),
        s13=np.sin(THETA13),
        s23=np.sin(THETA23),
        δ_CP=DELTA_CP,
        α_1=0.0,
        α_2=0.0,
        symb=False,
    ):
        """
        Creates a unitary matrix in the parametrisation of eq. 1.1 in 1611.01514.
        Conventions for Majorana phases from from eq. 8 of 1710.00715.
        """
        self.symb = symb

        if not symb:
            # numpy
            dtype = np.complex128

            matrix_1 = np.matrix(
                [[1.0, 0.0, 0.0], [0.0, c23, s23], [0.0, -s23, c23]], dtype=dtype
            )
            matrix_2 = np.matrix(
                [
                    [c13, 0.0, s13 * np.exp(-1j * δ_CP)],
                    [0.0, 1.0, 0.0],
                    [-s13 * np.exp(1j * δ_CP), 0.0, c13],
                ],
                dtype=dtype,
            )
            matrix_3 = np.matrix(
                [[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]], dtype=dtype
            )
            # P matrix contains the Majorana phases
            matrix_p = np.matrix(
                [
                    [np.exp(1j * α_1), 0.0, 0.0],
                    [0.0, np.exp(1j * α_2), 0.0],
                    [0.0, 0.0, 1.0],
                ],
                dtype=dtype,
            )
        if symb:
            # sympy matrices
            matrix_1 = sym.Matrix([[1.0, 0.0, 0.0], [0.0, c23, s23], [0.0, -s23, c23]])
            matrix_2 = sym.Matrix(
                [
                    [c13, 0.0, s13 * sym.exp(-1j * δ_CP)],
                    [0.0, 1.0, 0.0],
                    [-s13 * sym.exp(1j * δ_CP), 0.0, c13],
                ]
            )
            matrix_3 = sym.Matrix([[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]])
            # P matrix contains the Majorana phases
            matrix_p = sym.Matrix(
                [
                    [sym.exp(1j * α_1), 0.0, 0.0],
                    [0.0, sym.exp(1j * α_2), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            )

        if not symb:
            self.val = ft.reduce(np.dot, [matrix_1, matrix_2, matrix_3, matrix_p])
        else:
            self.val = ft.reduce(
                lambda x, y: x * y, [matrix_1, matrix_2, matrix_3, matrix_p]
            )

    def cols(self):
        return self.val.T

    def get_dagger(self):
        return self.val.getH()


# Measured squared neutrino mass differences
# atm is \sqrt{Δm_{21}^2} and sol is \sqrt{Δm_{3l}^2} (in GeV)
measured_mv_diffs = {"sol": np.sqrt(7.40e-5) * 1e-9, "atm": np.sqrt(2.494e-3) * 1e-9}

# NO: v1                 v2                               v3
#          <- Δm_sol ->      <------   Δm_atm   ------>

# IO: v3                v1          v2
#        <------   Δm_atm   ------>
#                         <- Δm_sol ->
NEUTRINO_MASSES = [0.0, measured_mv_diffs["sol"], measured_mv_diffs["atm"]]

PMNS = UnitaryMatrix()

# Measured CKM matrix parameters
THETA12_CKM = np.deg2rad(13.04)
THETA13_CKM = np.deg2rad(0.201)
THETA23_CKM = np.deg2rad(2.38)
DELTA_CP_CKM = np.deg2rad(1.20)

CKM = UnitaryMatrix(
    c12=np.cos(THETA12_CKM),
    c13=np.cos(THETA13_CKM),
    c23=np.cos(THETA23_CKM),
    s12=np.sin(THETA12_CKM),
    s13=np.sin(THETA13_CKM),
    s23=np.sin(THETA23_CKM),
    δ_CP=DELTA_CP_CKM,
    α_1=0.0,
    α_2=0.0,
    symb=False,
).val

VEV = 246.0 / np.sqrt(2)
