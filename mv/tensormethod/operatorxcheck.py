#!/usr/bin/env python3

"""Script to check that no operators are missing from the numbered list."""

from mv.tensormethod.parse_hs import parse
from mv.tensormethod.parse_hs import H7_LNV_NF3
from mv.tensormethod.parse_hs import H9_LNV_NF3
from mv.tensormethod.parse_hs import H11_LNV_NF3
from mv.tensormethod.contract import invariants
from mv.tensormethod.sm import H, L


def list_invariants(hs):
    """Returns a list of non-vanishing invariants, grouped by field content."""
    out = []
    for fields in parse(hs):
        # don't worry about O1pp
        if fields == [H, H, H, H, H, H.conj, H.conj, H.conj, L, L]:
            print("Skipped!")
            continue

        try:
            invs = invariants(*fields)
        except:
            breakpoint()

        # flatten list
        for inv in invs:

            if inv.contains_derivative:
                continue

            bl = inv.BL_classification
            if bl == -1:
                out.append(inv)

    return out
