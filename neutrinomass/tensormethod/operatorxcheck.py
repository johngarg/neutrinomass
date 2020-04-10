#!/usr/bin/env python3

"""Script to check that no operators are missing from the numbered list."""

from neutrinomass.tensormethod.parse_hs import parse
from neutrinomass.tensormethod.parse_hs import H6_NF3
from neutrinomass.tensormethod.parse_hs import H7_LNV_NF3
from neutrinomass.tensormethod.parse_hs import H9_LNV_NF3
from neutrinomass.tensormethod.parse_hs import H11_LNV_NF3
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.sm import H, L


def list_invariants(hs):
    """Returns a list of non-vanishing invariants, grouped by field content."""
    out = []
    for fields in parse(hs):
        # don't worry about O1pp
        if fields == [H, H, H, H, H, H.conj, H.conj, H.conj, L, L]:
            print("Skipped!")
            continue

        invs = invariants(*fields)

        # flatten list
        for inv in invs:

            # if inv.contains_derivative:
            #     continue

            bl = inv.BL_classification
            if bl == -1:
                out.append(inv)

    return out
