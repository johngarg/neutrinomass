#!/usr/bin/env python3

import sympy
import neutrinomass.tensormethod as tm

# you can hard code in the vanishing Higgs index contractions


class Eps:
    def __init__(self, i, j):
        indices = [i, j]
        if sorted(indices) == indices:
            self.fst = i
            self.scd = j
            self.is_neg = False
        else:
            self.fst = j
            self.scd = i
            self.is_neg = True

        self.indices = (self.fst, self.scd)

    def __lt__(self, other):
        return self.indices < other.indices

    def __repr__(self):
        maybe_neg = "-" if self.is_neg else ""
        return f"{maybe_neg}eps{self.fst, self.scd}"

    def str_no_neg(self):
        return f"eps{self.fst, self.scd}"


class Sum:
    def __init__(self, *args):
        self.args = sorted(args)

    def __repr__(self):
        return " + ".join(str(a) for a in self.args)


class Prod:
    def __init__(self, *args):
        self.args = sorted(args)

    def __repr__(self):
        string = "*".join(a.str_no_neg() for a in self.args)
        if self.is_neg:
            string = "-" + string
        return string

    @property
    def is_neg(self):
        neg = 0
        for arg in self.args:
            if arg.is_neg:
                neg += 1
        return not (neg % 2 == 0)
