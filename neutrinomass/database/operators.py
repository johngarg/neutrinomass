#!/usr/bin/env python3

from neutrinomass.completions import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.database.utils import table_data
from neutrinomass.database import MVDF

import pickle
from collections import Counter
import os


def print_paper_table():
    operator_latex = {}
    operator_latex[
        "1"
    ] = r"$L^{i} L^{j} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "2"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "3a"
    ] = r"$L^{i} L^{j} Q^{k} \bar{d} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "3b"
    ] = r"$L^{i} L^{j} Q^{k} \bar{d} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "4a"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \bar{u}^{\dagger} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "4b"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \bar{u}^{\dagger} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "5a"
    ] = r"$L^{i} L^{j} Q^{k} \bar{d} H^{l} H^{m} \tilde{H}^{n}  \cdot \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "5b"
    ] = r"$\mathcal{O}_1 \cdot Q^{i} \bar{d} \tilde{H}^{j}  \cdot \epsilon_{i j}$"
    operator_latex[
        "5c"
    ] = r"$\mathcal{O}_{3a} \cdot H^{i} \tilde{H}^{j}  \cdot \epsilon_{i j}$"
    operator_latex[
        "5d"
    ] = r"$\mathcal{O}_{3b} \cdot H^{i} \tilde{H}^{j}  \cdot \epsilon_{i j}$"
    operator_latex[
        "6a"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \bar{u}^{\dagger} H^{l} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "6b"
    ] = r"$\mathcal{O}_1 \cdot \tilde{Q}^{i} \bar{u}^{\dagger} \tilde{H}^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "6c"
    ] = r"$\mathcal{O}_{4a} \cdot H^{i} \tilde{H}^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "6d"
    ] = r"$\mathcal{O}_{4b} \cdot H^{i} \tilde{H}^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "7"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \tilde{Q}^{k} H^{l} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "8"
    ] = r"$L^{i} \bar{e}^{\dagger} \bar{u}^{\dagger} \bar{d} H^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "9"
    ] = r"$L^{i} L^{j} L^{k} L^{l} \bar{e} \bar{e}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "10"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} Q^{l} \bar{d}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "11a"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \bar{d} \bar{d}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "11b"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \bar{d} \bar{d}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "12a"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{u}^{\dagger}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "12b"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{u}^{\dagger}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "13"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} \tilde{Q}^{l} \bar{u}^{\dagger}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "14a"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "14b"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "15"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u}^{\dagger} \bar{d}  \cdot \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "16"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} \bar{u}^{\dagger} \bar{d}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "17"
    ] = r"$L^{i} L^{j} \bar{u}^{\dagger} \bar{d} \bar{d} \bar{d}^{\dagger}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "18"
    ] = r"$L^{i} L^{j} \bar{u} \bar{u}^{\dagger} \bar{u}^{\dagger} \bar{d}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "19"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \bar{u}^{\dagger} \bar{d} \bar{d}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "20"
    ] = r"$L^{i} \bar{e}^{\dagger} \tilde{Q}^{j} \bar{u}^{\dagger} \bar{u}^{\dagger} \bar{d}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "21a"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} Q^{l} \bar{u} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "21b"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} Q^{l} \bar{u} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "22a"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{e} \bar{e}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "22b"
    ] = r"$\mathcal{O}_2 \cdot \tilde{L}^i \bar{e}^\dagger H^j \epsilon_{ij}$"
    operator_latex[
        "23a"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} \tilde{Q}^{l} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "23b"
    ] = r"$\mathcal{O}_2 \cdot \tilde{Q}^i \bar{d}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "24a"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \bar{d} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "24b"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \bar{d} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "24c"
    ] = r"$\mathcal{O}_{3a} \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "24d"
    ] = r"$\mathcal{O}_{3b} \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "24e"
    ] = r"$\mathcal{O}_{11a} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "24f"
    ] = r"$\mathcal{O}_{11b} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "25a"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \bar{u} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "25b"
    ] = r"$\mathcal{O}_{3a} \cdot Q^i \bar{u} H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "25c"
    ] = r"$\mathcal{O}_{3b} \cdot Q^i \bar{u} H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "26a"
    ] = r"$L^{i} L^{j} \tilde{L}^{k} \bar{e}^{\dagger} Q^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "26b"
    ] = r"$L^{i} L^{j} \tilde{L}^{k} \bar{e}^{\dagger} Q^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "26c"
    ] = r"$\mathcal{O}_{3a} \cdot \tilde{L}^i \bar{e}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "26d"
    ] = r"$\mathcal{O}_{3b} \cdot \tilde{L}^i \bar{e}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "27a"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{d} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "27b"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{d} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "27c"
    ] = r"$\mathcal{O}_{3a} \cdot \tilde{Q}^i \bar{d}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "27d"
    ] = r"$\mathcal{O}_{3b} \cdot \tilde{Q}^i \bar{d}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28a"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "28b"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "28c"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j n} \epsilon_{l m}$"
    operator_latex[
        "28d"
    ] = r"$\mathcal{O}_{3a} \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28e"
    ] = r"$\mathcal{O}_{3b} \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28f"
    ] = r"$\mathcal{O}_{4a} \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28g"
    ] = r"$\mathcal{O}_{4b} \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28h"
    ] = r"$\mathcal{O}_{14a} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "28i"
    ] = r"$\mathcal{O}_{14b} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "29a"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "29b"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} \bar{u} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "29c"
    ] = r"$\mathcal{O}_{4a} \cdot Q^i \bar{u} H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "29d"
    ] = r"$\mathcal{O}_{4b} \cdot Q^i \bar{u} H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "30a"
    ] = r"$L^{i} L^{j} \tilde{L}^{k} \bar{e}^{\dagger} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "30b"
    ] = r"$L^{i} L^{j} \tilde{L}^{k} \bar{e}^{\dagger} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "30c"
    ] = r"$\mathcal{O}_{4a} \cdot \tilde{L}^i \bar{e}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "30d"
    ] = r"$\mathcal{O}_{4b} \cdot \tilde{L}^i \bar{e}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "31a"
    ] = r"$\mathcal{O}_{4a} \cdot \tilde{Q}^i \bar{d}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "31b"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "31c"
    ] = r"$\mathcal{O}_{4b} \cdot \tilde{Q}^i \bar{d}^\dagger H^j \cdot \epsilon_{ij}$"
    operator_latex[
        "32a"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{u}^{\dagger} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "32b"
    ] = r"$L^{i} L^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} \bar{u}^{\dagger} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "32c"
    ] = r"$\mathcal{O}_{4a} \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "32d"
    ] = r"$\mathcal{O}_{4b} \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex["32e"] = r"$\mathcal{O}_{12a} \cdot H^i \tilde{H}^j$"
    operator_latex["32f"] = r"$\mathcal{O}_{12b} \cdot H^i \tilde{H}^j$"
    operator_latex[
        "33"
    ] = r"$\mathcal{O}_1 \cdot \bar{e} \bar{e} \bar{e}^{\dagger} \bar{e}^{\dagger}$"
    operator_latex[
        "34"
    ] = r"$L^{i} \bar{e} \bar{e}^{\dagger} \bar{e}^{\dagger} Q^{j} \bar{d} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "35"
    ] = r"$L^{i} \bar{e} \bar{e}^{\dagger} \bar{e}^{\dagger} \tilde{Q}^{j} \bar{u}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "36"
    ] = r"$\bar{e}^{\dagger} \bar{e}^{\dagger} Q^{i} Q^{j} \bar{d} \bar{d} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "37"
    ] = r"$\bar{e}^{\dagger} \bar{e}^{\dagger} Q^{i} \tilde{Q}^{j} \bar{u}^{\dagger} \bar{d} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "38"
    ] = r"$\bar{e}^{\dagger} \bar{e}^{\dagger} \tilde{Q}^{i} \tilde{Q}^{j} \bar{u}^{\dagger} \bar{u}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "39a"
    ] = r"$\mathcal{O}_1 \cdot L^{i} L^{j} \tilde{L}^{k} \tilde{L}^{l} \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "39b"
    ] = r"$L^{i} L^{j} L^{k} L^{l} \tilde{L}^{m} \tilde{L}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m p} \epsilon_{n q}$"
    operator_latex[
        "39c"
    ] = r"$L^{i} L^{j} L^{k} L^{l} \tilde{L}^{m} \tilde{L}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k p} \epsilon_{m q}$"
    operator_latex[
        "39d"
    ] = r"$\mathcal{O}_1 \cdot L^{i} L^{j} \tilde{L}^{k} \tilde{L}^{l} \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "40a"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k p} \epsilon_{m q}$"
    operator_latex[
        "40b"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j p} \epsilon_{k q} \epsilon_{m n}$"
    operator_latex[
        "40c"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i n} \epsilon_{j p} \epsilon_{k q} \epsilon_{l m}$"
    operator_latex[
        "40d"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m p} \epsilon_{n q}$"
    operator_latex[
        "40e"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k p} \epsilon_{n q}$"
    operator_latex[
        "40f"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j n} \epsilon_{l p} \epsilon_{m q}$"
    operator_latex[
        "40g"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k p} \epsilon_{l q}$"
    operator_latex[
        "40h"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l p} \epsilon_{n q}$"
    operator_latex[
        "40i"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i m} \epsilon_{j p} \epsilon_{k q} \epsilon_{l n}$"
    operator_latex[
        "40j"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j p} \epsilon_{l n} \epsilon_{m q}$"
    operator_latex[
        "40k"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j p} \epsilon_{l m} \epsilon_{n q}$"
    operator_latex[
        "40l"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} Q^{m} \tilde{Q}^{n} H^{p} H^{q} \cdot  \epsilon_{i k} \epsilon_{j p} \epsilon_{l q} \epsilon_{m n}$"
    operator_latex[
        "41a"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{d} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "41b"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{d} \bar{d}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "42a"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "42b"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "43a"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j n} \epsilon_{l m}$"
    operator_latex[
        "43b"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "43c"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "43d"
    ] = r"$L^{i} L^{j} L^{k} \tilde{L}^{l} \bar{u}^{\dagger} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m n}$"
    operator_latex[
        "44a"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} Q^{k} \tilde{Q}^{l} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "44b"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} Q^{k} \tilde{Q}^{l} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "44c"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} Q^{k} \tilde{Q}^{l} H^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "44d"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} Q^{k} \tilde{Q}^{l} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "45"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} \bar{d} \bar{d}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "46"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} \bar{u} \bar{u}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "47a"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k p} \epsilon_{l q}$"
    operator_latex[
        "47b"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i m} \epsilon_{j p} \epsilon_{k n} \epsilon_{l q}$"
    operator_latex[
        "47c"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i p} \epsilon_{j q} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "47d"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k p} \epsilon_{m q}$"
    operator_latex[
        "47e"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i n} \epsilon_{j p} \epsilon_{k l} \epsilon_{m q}$"
    operator_latex[
        "47f"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i j} \epsilon_{k n} \epsilon_{l p} \epsilon_{m q}$"
    operator_latex[
        "47g"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j p} \epsilon_{k n} \epsilon_{m q}$"
    operator_latex[
        "47h"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i j} \epsilon_{k p} \epsilon_{l q} \epsilon_{m n}$"
    operator_latex[
        "47i"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i l} \epsilon_{j p} \epsilon_{k q} \epsilon_{m n}$"
    operator_latex[
        "47j"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i p} \epsilon_{j q} \epsilon_{k l} \epsilon_{m n}$"
    operator_latex[
        "47k"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m p} \epsilon_{n q}$"
    operator_latex[
        "47l"
    ] = r"$L^{i} L^{j} Q^{k} Q^{l} \tilde{Q}^{m} \tilde{Q}^{n} H^{p} H^{q}  \cdot  \epsilon_{i j} \epsilon_{k l} \epsilon_{m p} \epsilon_{n q}$"
    operator_latex[
        "48"
    ] = r"$L^{i} L^{j} \bar{d} \bar{d} \bar{d}^{\dagger} \bar{d}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "49"
    ] = r"$L^{i} L^{j} \bar{u} \bar{u}^{\dagger} \bar{d} \bar{d}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "50a"
    ] = r"$L^{i} L^{j} \bar{u}^{\dagger} \bar{d} \bar{d} \bar{d}^{\dagger} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "50b"
    ] = r"$\mathcal{O}_{17} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "51"
    ] = r"$L^{i} L^{j} \bar{u} \bar{u} \bar{u}^{\dagger} \bar{u}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "52a"
    ] = r"$L^{i} L^{j} \bar{u} \bar{u}^{\dagger} \bar{u}^{\dagger} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "52b"
    ] = r"$\mathcal{O}_{18} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "53"
    ] = r"$L^{i} L^{j} \bar{u}^{\dagger} \bar{u}^{\dagger} \bar{d} \bar{d} \tilde{H}^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "54a"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} Q^{k} \tilde{Q}^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "54b"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} Q^{k} \tilde{Q}^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j l} \epsilon_{k n}$"
    operator_latex[
        "54c"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} Q^{k} \tilde{Q}^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j k} \epsilon_{l n}$"
    operator_latex[
        "54d"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} Q^{k} \tilde{Q}^{l} \bar{d} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "55a"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "55b"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j l} \epsilon_{k n}$"
    operator_latex[
        "55c"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "55d"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \tilde{Q}^{k} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "56"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \bar{d} \bar{d} \bar{d}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "57"
    ] = r"$L^{i} \bar{e}^{\dagger} \tilde{Q}^{j} \bar{u}^{\dagger} \bar{d} \bar{d}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "58"
    ] = r"$L^{i} \bar{e}^{\dagger} \tilde{Q}^{j} \bar{u} \bar{u}^{\dagger} \bar{u}^{\dagger} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "59a"
    ] = r"$L^{i} \bar{e}^{\dagger} Q^{j} \bar{u}^{\dagger} \bar{d} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "59b"
    ] = r"$\mathcal{O}_{19} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "59c"
    ] = r"$\mathcal{O}_{8} \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "60a"
    ] = r"$L^{i} \bar{e}^{\dagger} \tilde{Q}^{j} \bar{u}^{\dagger} \bar{u}^{\dagger} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "60b"
    ] = r"$\mathcal{O}_8 \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "60c"
    ] = r"$\mathcal{O}_{20} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "61a"
    ] = r"$\mathcal{O}_1 \cdot L^i \bar{e} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex["61b"] = r"$\mathcal{O}_2 \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "62a"
    ] = r"$L^{i} L^{j} L^{k} L^{l} \bar{e} \bar{e} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex["62b"] = r"$\mathcal{O}_9 \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "63a"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} Q^{l} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j n} \epsilon_{l m}$"
    operator_latex[
        "63b"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} Q^{l} \bar{d} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "63c"
    ] = r"$\mathcal{O}_2 \cdot Q^i \bar{d} \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "63d"
    ] = r"$\mathcal{O}_{10} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "64a"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "64b"
    ] = r"$L^{i} L^{j} L^{k} \bar{e} \tilde{Q}^{l} \bar{u}^{\dagger} H^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j n} \epsilon_{l m}$"
    operator_latex[
        "64c"
    ] = r"$\mathcal{O}_2 \cdot \tilde{Q}^i \bar{u}^\dagger \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "64d"
    ] = r"$\mathcal{O}_{13} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "65a"
    ] = r"$L^{i} L^{j} \bar{e} \bar{e}^{\dagger} \bar{u}^{\dagger} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "65b"
    ] = r"$\mathcal{O}_{16} \cdot H^i \tilde{H}^j \cdot \epsilon_{ij}$"
    operator_latex[
        "71"
    ] = r"$\mathcal{O}_1 \cdot Q^{i} \bar{u} H^{j} \cdot  \epsilon_{i j}$"
    operator_latex[
        "75"
    ] = r"$\mathcal{O}_8 \cdot Q^{i} \bar{u} H^j \cdot \epsilon_{i j}$"
    operator_latex[
        "76"
    ] = r"$\bar{e}^\dagger \bar{e}^\dagger \bar{u}^\dagger \bar{u}^\dagger \bar{d} \bar{d}$"
    operator_latex[
        "77"
    ] = r"$\mathcal{O}_{1} \cdot \tilde{L}^{i} \bar{e}^{\dagger} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "78"
    ] = r"$\mathcal{O}_{1} \cdot \tilde{Q}^{i} \bar{d}^{\dagger} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "1p"
    ] = r"$\mathcal{O}_{1} \cdot \tilde{H}^{i} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "8p"
    ] = r"$\mathcal{O}_{8} \cdot \tilde{H}^{i} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "1pp"
    ] = r"$\mathcal{O}_{1} \cdot\tilde{H}^{i} H^{j} \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "1ppp"
    ] = r"$\mathcal{O}_{1} \cdot \tilde{H}^{i} H^{j} \tilde{H}^{k} H^{l} \tilde{H}^{m} H^{n} \cdot \epsilon_{ij} \epsilon_{kl} \epsilon_{mn}$"
    operator_latex[
        "7p"
    ] = r"$\mathcal{O}_{7} \cdot \tilde{H}^{i} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "8pp"
    ] = r"$\mathcal{O}_{8} \cdot \tilde{H}^{i} H^{j} \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "71p"
    ] = r"$\mathcal{O}_{71} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "76p"
    ] = r"$\mathcal{O}_{76} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "77p"
    ] = r"$\mathcal{O}_{77} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "78p"
    ] = r"$\mathcal{O}_{78} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "79a"
    ] = r"$\mathcal{O}_{61a} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "79b"
    ] = r"$\mathcal{O}_{2} \cdot \tilde{H}^{i} H^{j}  \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "80a"
    ] = r"$\mathcal{O}_{5a} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "80b"
    ] = r"$\mathcal{O}_{5b} \cdot \tilde{H}^{i} H^{j}  \cdot \epsilon_{ij}$"
    operator_latex[
        "80c"
    ] = r"$\mathcal{O}_{3a} \cdot \tilde{H}^{i} H^{j}  \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "80d"
    ] = r"$\mathcal{O}_{3b} \cdot \tilde{H}^{i} H^{j}  \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "81a"
    ] = r"$\mathcal{O}_{6a} \cdot \tilde{H}^{i} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "81b"
    ] = r"$\mathcal{O}_{6b} \cdot \tilde{H}^{i} H^{j} \cdot \epsilon_{ij}$"
    operator_latex[
        "81c"
    ] = r"$\mathcal{O}_{4a} \cdot \tilde{H}^{i} H^{j}  \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "81d"
    ] = r"$\mathcal{O}_{4b} \cdot \tilde{H}^{i} H^{j}  \tilde{H}^{k} H^{l} \cdot \epsilon_{ij} \epsilon_{kl}$"
    operator_latex[
        "82"
    ] = r"$L^{i} \tilde{L}^{j} \bar{e}^{\dagger} \bar{e}^{\dagger} \bar{u}^{\dagger} \bar{d} H^{k} H^{l} \cdot \epsilon_{ik} \epsilon_{jl}$"
    operator_latex[
        "D1"
    ] = r"$(DL)^{i} L^{j} {\bar{u}^{\dagger}} \bar{d}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "D2a"
    ] = r"$(DL)^{i} L^{j} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D2b"
    ] = r"$(DL)^{i} L^{j} (DH)^{k} H^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D2c"
    ] = r"$(DL)^{i} L^{j} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D3"
    ] = r"$L^{i} {\bar{e}^{\dagger}} H^{j} H^{k} (DH)^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D4a"
    ] = r"$L^{i} L^{j} (DL)^{k} (D\bar{e}) H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D4b"
    ] = r"$L^{i} L^{j} (DL)^{k} (D\bar{e}) H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D5a"
    ] = r"$L^{i} L^{j} (DL)^{k} \tilde{L}^{l} H^{m} H^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "D5b"
    ] = r"$L^{i} L^{j} (DL)^{k} \tilde{L}^{l} H^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "D5c"
    ] = r"$L^{i} L^{j} (DL)^{k} \tilde{L}^{l} H^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "D5d"
    ] = r"$L^{i} L^{j} (DL)^{k} \tilde{L}^{l} H^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "D6a"
    ] = r"$L^{i} L^{j} \bar{e} {\bar{e}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D6b"
    ] = r"$L^{i} L^{j} \bar{e} {\bar{e}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D7a"
    ] = r"$(DL)^{i} L^{j} Q^{k} (D\bar{d}) H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D7b"
    ] = r"$(DL)^{i} L^{j} Q^{k} (D\bar{d}) H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D7c"
    ] = r"$(DL)^{i} L^{j} Q^{k} (D\bar{d}) H^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D8a"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i n} \epsilon_{j k} \epsilon_{l m}$"
    operator_latex[
        "D8b"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i n} \epsilon_{j l} \epsilon_{k m}$"
    operator_latex[
        "D8c"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m n}$"
    operator_latex[
        "D8d"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j k} \epsilon_{l n}$"
    operator_latex[
        "D8e"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j l} \epsilon_{k n}$"
    operator_latex[
        "D8f"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i m} \epsilon_{j n} \epsilon_{k l}$"
    operator_latex[
        "D8g"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "D8h"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k n} \epsilon_{l m}$"
    operator_latex[
        "D8i"
    ] = r"$L^{i} L^{j} Q^{k} \tilde{Q}^{l} (DH)^{m} H^{n}  \cdot  \epsilon_{i j} \epsilon_{k l} \epsilon_{m n}$"
    operator_latex[
        "D9a"
    ] = r"$L^{i} L^{j} \bar{d} {\bar{d}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D9b"
    ] = r"$L^{i} L^{j} \bar{d} {\bar{d}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D10a"
    ] = r"$(DL)^{i} L^{j} {\bar{u}^{\dagger}} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D10b"
    ] = r"$(DL)^{i} L^{j} {\bar{u}^{\dagger}} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D10c"
    ] = r"$(DL)^{i} L^{j} {\bar{u}^{\dagger}} \bar{d} H^{k} \tilde{H}^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D11"
    ] = r"$(DL)^{i} L^{j} (D\bar{u}^{\dagger}) (D\bar{d})  \cdot  \epsilon_{i j}$"
    operator_latex[
        "D12a"
    ] = r"$L^{i} L^{j} \bar{u} {\bar{u}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D12b"
    ] = r"$L^{i} L^{j} \bar{u} {\bar{u}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D13a"
    ] = r"$(DL)^{i} L^{j} \tilde{Q}^{k} (D{\bar{u}^{\dagger}}) H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D13b"
    ] = r"$(DL)^{i} L^{j} \tilde{Q}^{k} (D{\bar{u}^{\dagger}}) H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D14a"
    ] = r"$L^{i} {\bar{e}^{\dagger}} Q^{j} \bar{d} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D14b"
    ] = r"$L^{i} {\bar{e}^{\dagger}} Q^{j} \bar{d} (DH)^{k} H^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D14c"
    ] = r"$L^{i} {\bar{e}^{\dagger}} Q^{j} \bar{d} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D15"
    ] = r"$(DL)^{i} {\bar{e}^{\dagger}} (D \bar{u}^{\dagger}) \bar{d} H^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "D16a"
    ] = r"$L^{i} {\bar{e}^{\dagger}} \tilde{Q}^{j} {\bar{u}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D16b"
    ] = r"$L^{i} {\bar{e}^{\dagger}} \tilde{Q}^{j} {\bar{u}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D16c"
    ] = r"$L^{i} {\bar{e}^{\dagger}} \tilde{Q}^{j} {\bar{u}^{\dagger}} (DH)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D17"
    ] = r"${\bar{e}^{\dagger}} {\bar{e}^{\dagger}} {\bar{u}^{\dagger}} \bar{d} (DH)^{i} H^{j}  \cdot  \epsilon_{i j}$"
    operator_latex[
        "D18a"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j m} \epsilon_{l n}$"
    operator_latex[
        "D18b"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i k} \epsilon_{j l} \epsilon_{m n}$"
    operator_latex[
        "D18c"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i m} \epsilon_{j l} \epsilon_{k n}$"
    operator_latex[
        "D18d"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i j} \epsilon_{k m} \epsilon_{l n}$"
    operator_latex[
        "D18e"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i n} \epsilon_{j l} \epsilon_{k m}$"
    operator_latex[
        "D18f"
    ] = r"$(DL)^{i} L^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j n} \epsilon_{k m}$"
    operator_latex[
        "D19a"
    ] = r"$(DL)^{i} L^{j} (D^{3} H)^{k} H^{l}  \cdot  \epsilon_{i j} \epsilon_{k l}$"
    operator_latex[
        "D19b"
    ] = r"$(DL)^{i} L^{j} (D^{3} H)^{k} H^{l}  \cdot  \epsilon_{i l} \epsilon_{j k}$"
    operator_latex[
        "D19c"
    ] = r"$(DL)^{i} L^{j} (D^{3} H)^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D20"
    ] = r"$L^{i} {\bar{e}^{\dagger}} H^{j} H^{k} H^{l} (DH)^{m} \tilde{H}^{n}  \cdot  \epsilon_{i l} \epsilon_{j m} \epsilon_{k n}$"
    operator_latex[
        "D21"
    ] = r"$(DL)^{i} (D\bar{e}^{\dagger}) H^{j} H^{k} (DH)^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"
    operator_latex[
        "D22"
    ] = r"${\bar{e}^{\dagger}} {\bar{e}^{\dagger}} (DH)^{i} (DH)^{j} H^{k} H^{l}  \cdot  \epsilon_{i k} \epsilon_{j l}$"

    models_dict = pickle.load(
        open(os.path.join(os.path.dirname(__file__), "models.p"), "rb")
    )
    op_labels = set(models_dict.keys())
    filtered_dict = Counter(list(MVDF["op"]))

    for k, v in models_dict.items():
        if k not in filtered_dict:
            filtered_dict[k] = 0

    models_dict = {k: str(v) for k, v in models_dict.items()}
    filtered_dict = {k: str(v) for k, v in filtered_dict.items()}

    for label, latex in operator_latex.items():

        if label in EFF_OPERATORS:
            op = EFF_OPERATORS[label]
        elif label in DERIV_EFF_OPERATORS:
            op = DERIV_EFF_OPERATORS[label]

        loops, scale = table_data(op)

        if label not in op_labels:
            models, filtered = "", ""
        else:
            models = models_dict[label]
            filtered = filtered_dict[label]

        row = fr"${format_primes(label)}$ & {latex} & {models} & {filtered} & {loops} & {scale} \\"
        print(row)


def format_primes(string):
    if string.endswith("ppp"):
        return string[:-3] + r"^{\prime\prime\prime}"
    if string.endswith("pp"):
        return string[:-2] + r"^{\prime\prime}"
    if string.endswith("p"):
        return string[:-1] + r"^{\prime}"
    return string
