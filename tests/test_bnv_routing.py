"""Dimension-eight benchmark from Heeck and Sokhashvili, arXiv:2611.08621."""
from functools import reduce
from operator import mul

import pytest
from sympy import Rational

from neutrinomass.completions.amplitudes import derivative_placement_polynomial
from neutrinomass.completions.completions import (
    deriv_operator_completion_stream, exact_completions, partition_leaves,
)
from neutrinomass.completions.core import EffectiveOperator, VectorLikeDiracFermion
from neutrinomass.completions.routing import (
    InvalidMomentumGraph,
    classify_fermion_pair,
    species_states,
    validate_routes,
)
from neutrinomass.tensormethod.contract import su2_singlets_type
from neutrinomass.tensormethod.core import D, Field
from neutrinomass.tensormethod.sm import H, L, db, ub


@pytest.mark.parametrize(
    "use_baryon_labels,seed,generate",
    [
        (True, 1, exact_completions), (False, 1, exact_completions),
        (True, 0, deriv_operator_completion_stream), (True, 0, exact_completions),
    ],
    ids=["standard-quarks", "no-3b", "higgs-stream", "higgs-exact"],
)
def test_dhdduL_exact_mediators_and_fermion_momentum(use_baryon_labels, seed, generate):
    fields = [H, db.conj, db.conj, ub.conj, L]
    assert sum(f.charges["3b"] for f in fields) == 3
    if not use_baryon_labels:
        # With a common neutral spectator assignment, local and routed
        # candidates enter the same exact-equivalence bucket.
        fields = [
            Field(
                f.label, f.dynkin, charges={**f.charges, "3b": 0},
                is_conj=f.is_conj, comm=f.comm, nf=f.nf, latex=f.latex,
            )
            for f in fields
        ]
    fields[seed] = D(fields[seed], "11" if seed == 0 else "10")
    tensor = reduce(mul, (f.fresh_indices() for f in fields))
    weak_singlet, = su2_singlets_type(tensor, "i")
    operator = EffectiveOperator("DHdduL", weak_singlet)

    # Each field is written as a colour triplet or singlet, with its signed Y.
    s2 = ((1, 0), 1, -Rational(1, 3))
    s3 = ((1, 0), 1, Rational(2, 3))
    f1 = ((0, 0), 1, 0)
    f8 = ((1, 0), 2, Rational(1, 6))
    f13 = ((1, 0), 2, Rational(7, 6))
    expected = {(s2, f1), (s2, f8), (s3, f1), (s3, f13)}
    found = set()
    verified_cuts = set()
    higgs_placements = set()
    for completion in generate(operator):
        assert completion.canonical_topology == "1s4f_2"
        scalar, fermion = sorted(completion.exotics, key=lambda f: f.is_fermion)
        assert scalar.is_scalar and fermion.is_fermion
        pair = []
        for field in (scalar, fermion):
            if field.colour_irrep == (0, 1):
                field = field.conj
            pair.append((field.colour_irrep, field.isospin_irrep[0] + 1, field.y))
        pair = tuple(pair)
        found.add(pair)
        audit = validate_routes(completion, derivative_order=1)
        contribution, = completion.momentum_contributions
        assert contribution.numerator_kind == "momentum"
        assert audit.momentum_edges == (tuple(sorted(contribution.edge)),)
        if pair != (s3, f13):
            continue
        leaves = {
            leaf.node: (leaf.field.strip_derivs_with_indices() if leaf.field.derivs else leaf.field)
            for leaf in partition_leaves(completion.partition)
        }
        side = set(contribution.cut_side)
        assert side and side < leaves.keys()
        assert derivative_placement_polynomial(completion) == {
            ((node, 1),): 1 for node in side
        }
        cut = frozenset(
            tuple(sorted(leaves[node].label for node in group))
            for group in (side, leaves.keys() - side)
        )
        verified_cuts.add(cut)
        if any(leaf.field.label == "DH" for leaf in partition_leaves(completion.partition)):
            higgs_placements.add(cut)

    assert found == expected
    assert verified_cuts == {frozenset({("H", "ub†"), ("L", "db†", "db†")})}
    if seed == 0 and generate is deriv_operator_completion_stream:
        assert higgs_placements == verified_cuts


def test_endpoint_charge_policy_keeps_other_quantum_numbers():
    def fermion(y, x, three_b):
        return VectorLikeDiracFermion(
            "F", "u0 c0", charges={"y": y, "x": x, "3b": three_b}
        )

    left = fermion(1, 2, 3)
    different_b = fermion(1, 2, 0).conj
    with pytest.raises(InvalidMomentumGraph):
        classify_fermion_pair(left, different_b, species_states({left}))

    ignored = ("3b",)
    states = species_states({left}, ignored)
    assert classify_fermion_pair(left, different_b, states, ignored) == "momentum"
    for wrong in (fermion(2, 2, 0).conj, fermion(1, 3, 0).conj):
        with pytest.raises(InvalidMomentumGraph):
            classify_fermion_pair(left, wrong, states, ignored)

    neutral = VectorLikeDiracFermion("N", "u0", charges={"y": 0, "3b": 3})
    with pytest.raises(InvalidMomentumGraph, match="Wrong vectorlike Weyl partner"):
        classify_fermion_pair(
            neutral, neutral.dirac_partner().conj, species_states({neutral}, ignored), ignored
        )
