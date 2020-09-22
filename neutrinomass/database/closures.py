#!/usr/bin/env python3

"""The file contains functions that calculcate operator-closure estimates for
the new-physics scale associated with ΔL = 2 operators under the assumption that
they give rise to radiative neutrino masses.

"""

from typing import Union, List

import math
import sympy
from functools import reduce
from matchpy import Operation, Symbol, Arity, match, Pattern, Wildcard, substitute
from neutrinomass.completions import EffectiveOperator, EFF_OPERATORS
from neutrinomass.tensormethod import H, L, Q, Field, IndexedField
from neutrinomass.tensormethod.core import ISOSPIN, GENERATION, Operator

Op = Operation.new("Op", Arity.variadic, commutative=True, associative=True)
Const = Operation.new("Const", Arity.unary)
c = Operation.new("c", Arity.unary)  # conjugate

# fields
e = Symbol("e")
nu = Symbol("nu")
u = Symbol("u")
d = Symbol("d")
h0 = Symbol("h0")
hp = Symbol("hp")
db = Symbol("db")
ub = Symbol("ub")
eb = Symbol("eb")
W = Symbol("W")

# masses
ye = Symbol("ye")
yu = Symbol("yu")
yd = Symbol("yd")
v = Symbol("v")

# couplings and symbols
nunu = Symbol("nunu")
g2 = Symbol("g2")
loop = Symbol("loop")
loopv2 = Symbol("loopv2")

# replacement rules for closures
RST = Wildcard.dot("rst")
RULES = {
    # neutrinos
    Op(nu, nu, RST): Op(Const(nunu), RST),
    # free loops
    Op(c(h0), h0, RST): Op(Const(loopv2), RST),
    Op(c(hp), hp, RST): Op(Const(loop), RST),
    Op(c(eb), eb, RST): Op(Const(loop), RST),
    Op(c(db), db, RST): Op(Const(loop), RST),
    Op(c(ub), ub, RST): Op(Const(loop), RST),
    Op(c(d), d, RST): Op(Const(loop), RST),
    Op(c(u), u, RST): Op(Const(loop), RST),
    Op(c(e), e, RST): Op(Const(loop), RST),
    Op(c(nu), nu, RST): Op(Const(loop), RST),
    # masses
    Op(e, eb, RST): Op(Const(v), Const(ye), Const(loop), RST),
    Op(u, ub, RST): Op(Const(v), Const(yu), Const(loop), RST),
    Op(d, db, RST): Op(Const(v), Const(yd), Const(loop), RST),
    Op(c(e), c(eb), RST): Op(Const(v), Const(ye), Const(loop), RST),
    Op(c(u), c(ub), RST): Op(Const(v), Const(yu), Const(loop), RST),
    Op(c(d), c(db), RST): Op(Const(v), Const(yd), Const(loop), RST),
    # make W
    Op(e, nu, RST): Op(Const(g2), W, nu, nu, RST),
    Op(c(eb), db, c(ub), RST): Op(
        Const(v),
        Const(g2),
        Const(loop),
        Const(loop),
        Const(loopv2),
        Const(yu),
        Const(yd),
        Const(ye),
        nu,
        RST,
    ),
    Op(c(eb), c(d), u, RST): Op(
        Const(v),
        Const(g2),
        Const(loop),
        Const(loop),
        Const(loopv2),
        Const(yu),
        Const(yd),
        Const(ye),
        nu,
        RST,
    ),
    Op(db, c(ub), c(hp), RST): Op(
        Const(v),
        Const(g2),
        Const(loop),
        Const(loop),
        Const(loopv2),
        Const(yu),
        Const(yd),
        RST,
    ),
    # remove W
    Op(W, hp, RST): Op(Const(v), Const(loop), RST),
    Op(W, u, RST): Op(d, Const(loop), RST),
    Op(W, c(d), RST): Op(c(u), Const(loop), RST),
    Op(W, c(ub), db, RST): Op(
        Const(v), Const(v), Const(yu), Const(yd), Const(loop), Const(loop), RST
    ),
    # remove hp
    Op(hp, c(eb), RST): Op(Const(loop), Const(ye), nu, RST),
    Op(hp, c(u), RST): Op(Const(yd), db, Const(loop), RST),
    Op(c(hp), db, RST): Op(Const(yd), c(u), Const(loop), RST),
    Op(hp, d, RST): Op(Const(yu), c(ub), Const(loop), RST),
    Op(hp, ub, RST): Op(Const(yu), c(d), Const(loop), RST),
    Op(c(hp), c(ub), RST): Op(Const(yu), d, Const(loop), RST),
    Op(hp, u, c(d), RST): Op(Const(g2), Const(loop), Const(loop), RST),
    # make hp
    Op(c(eb), nu, RST): Op(Const(ye), hp, nu, nu, RST),
    # vev
    Op(h0, RST): Op(Const(v), RST),
    Op(c(h0), RST): Op(Const(v), RST),
}


def apply_rules(rules, subject):
    for k, v in rules.items():
        for substitution in match(subject, Pattern(k)):
            subject = substitute(Pattern(v), substitution)
            return subject

    return subject


def fixed_point(start, rules=RULES, max_iterations=10, verbose=False):
    old, new = None, start
    counter = 1
    if verbose:
        print(start)
    while new != old:
        # Check if max iterations reached
        if counter > max_iterations:
            print("Maximum iterations reached on fixed_point")
            return new

        old = new
        new = apply_rules(rules, old)
        if verbose:
            print(new)
        counter += 1

    return new


def clean(op):
    lst = list(op)
    out = []
    while lst:
        item = lst.pop(0)
        if item.head.name == "Const":
            s = str(item.operands[0])
            out.append(s)
        else:
            # raise Exception(f"Non Const symbol {item} encountered in reduced operator.")
            return op

    return out


h = [hp, h0]
l = [nu, e]
q = [u, d]
hc = [h0, hp]
lc = [e, nu]
qc = [d, u]
eps = [[0, -1], [1, 0]]


def parse_operator(eff_op: Union[EffectiveOperator, Operator]):
    """Parse the operator `eff_op` into matchpy symbols with the SU(2) structure
    expanded.

    This is an unfortunately messy function and needs to be rewritten in a
    cleaner way.

    """
    if isinstance(eff_op, EffectiveOperator):
        operator = eff_op.operator
    else:
        operator = eff_op

    fields, epsilons = [], []
    n_indices = 0
    for expr in operator.tensors:
        if isinstance(expr, Field):
            if expr.derivs:
                expr = expr.strip_derivs_with_indices()

            i = expr.indices_by_type["Isospin"]
            if expr.is_conj:
                label = expr.label.lower()[:-1]
                label = label + ("c" if i else "")
            else:
                label = expr.label.lower()

            if i:
                label += f"[{i[0]}]"

            if expr.is_conj:
                label = f"c({label})"

            fields.append(label)

        else:
            epsilons.append(expr.indices)
            n_indices += 2

    field_string = "*".join(fields)
    eval_str = ""
    indices = [("_i", "_j"), ("_k", "_l"), ("_m", "_n"), ("_p", "_q")]
    for (i, j), (a, b) in zip(epsilons, indices):
        field_string = field_string.replace(str(-i), a)
        field_string = field_string.replace(str(-j), b)
        field_string += f"*eps[{a}][{b}]"

    loop_ranges = [1, 1, 1, 1, 1, 1, 1, 1]
    for i in range(n_indices):
        loop_ranges[i] += 1

    _a, _b, _c, _d, _e, _f, _g, _h = loop_ranges
    res = []
    for s1, _i in zip(["_i", "_i"], range(_a)):
        for s2, _j in zip(["_j", "_j"], range(_b)):
            for s3, _k in zip(["_k", "_k"], range(_c)):
                for s4, _l in zip(["_l", "_l"], range(_d)):
                    for s5, _m in zip(["_m", "_m"], range(_e)):
                        for s6, _n in zip(["_n", "_n"], range(_f)):
                            for s7, _p in zip(["_p", "_p"], range(_g)):
                                for s8, _q in zip(["_q", "_q"], range(_h)):
                                    res.append(
                                        field_string.replace(s1, str(_i))
                                        .replace(s2, str(_j))
                                        .replace(s3, str(_k))
                                        .replace(s4, str(_l))
                                        .replace(s5, str(_m))
                                        .replace(s6, str(_n))
                                        .replace(s7, str(_p))
                                        .replace(s8, str(_q))
                                    )

    out = []
    for elem in res:
        new_elem = elem.replace("*eps[0][1]", "").replace("*eps[1][0]", "")
        if not "eps" in new_elem:
            out.append(new_elem)

    return [eval(f'Op({elem.replace("*", ",")})') for elem in out]


def neutrino_mass_estimate(eff_op: Union[EffectiveOperator, List[Op]], verbose=False):
    if isinstance(eff_op, EffectiveOperator):
        clean_lst = [clean(fixed_point(op)) for op in parse_operator(eff_op)]
    else:
        clean_lst = [clean(fixed_point(op)) for op in eff_op]

    out = []
    for lst in clean_lst:
        if not isinstance(lst, list):
            if verbose:
                print(f"Skipping structure {lst} since it should be negligible")
            continue

        # make sure only two neutrinos are in the diagram
        assert lst.count("nunu") == 1
        lst.remove("nunu")

        n_vevs = lst.count("v")
        assert n_vevs % 2 == 0
        n_loopv2 = (n_vevs - 2) // 2

        for _ in range(n_vevs):
            lst.remove("v")

        # need to account for seesaw case
        prod = reduce(lambda x, y: x * y, [sympy.Symbol(i) for i in lst]) if lst else 1
        prod *= sympy.Symbol("v") * sympy.Symbol("v") / sympy.Symbol("Λ")
        for _ in range(n_loopv2):
            prod *= sympy.Symbol("loopv2")
        out.append(prod)

    return out


def numerical_np_scale_estimate(expr):
    """Returns log10 of estimate of Λ in TeV."""
    vev = 174
    mv = 5e-11
    loop = 1.0 / (16 * math.pi ** 2)
    subs_list = [
        (sympy.Symbol("v"), vev),
        (sympy.Symbol("loop"), loop),
        (sympy.Symbol("loopv2"), (loop + vev ** 2 / sympy.Symbol("Λ") ** 2),),
        (sympy.Symbol("g2"), 0.6295 ** 2),
        (sympy.Symbol("yu"), 172.76 / vev),
        (sympy.Symbol("yd"), 4.18 / vev),
        (sympy.Symbol("ye"), 1.78 / vev),
    ]
    m = expr.subs(subs_list)
    sol = sympy.solve(m - mv, sympy.Symbol("Λ"))[0]
    scale = abs(sol) * 1e-3
    return scale


def numerical_mv(expr):
    """Returns estimate of neutrino-mass scale in GeV assuming a NP scale of 1
    TeV.

    """
    vev = 174
    mv = 5e-11
    loop = 1.0 / (16 * math.pi ** 2)
    subs_list = [
        (sympy.Symbol("v"), vev),
        (sympy.Symbol("loop"), loop),
        (sympy.Symbol("loopv2"), (loop + vev ** 2 / sympy.Symbol("Λ") ** 2),),
        (sympy.Symbol("g2"), 0.6295 ** 2),
        (sympy.Symbol("yu"), 172.76 / vev),
        (sympy.Symbol("yd"), 4.18 / vev),
        (sympy.Symbol("ye"), 1.78 / vev),
        (sympy.Symbol("Λ"), 1000),
    ]
    return expr.subs(subs_list)
