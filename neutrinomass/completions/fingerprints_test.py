#!/usr/bin/env python3

from collections import Counter

from neutrinomass.completions import EFF_OPERATORS
from neutrinomass.completions.completions import collect_completions, operator_completions
from neutrinomass.completions.fingerprints import (
    completion_digest,
    interaction_fingerprint,
)
from neutrinomass.tensormethod import L, eps


def test_interaction_fingerprint_normalises_free_indices():
    first = (
        L("u0 i0 g8_")
        * L("u1 i1 g9_")
        * eps("-u0 -u1")
        * eps("-i0 -i1")
    )
    second = (
        L("u4 i4 g20_")
        * L("u5 i5 g21_")
        * eps("-u4 -u5")
        * eps("-i4 -i5")
    )

    assert interaction_fingerprint(first) == interaction_fingerprint(second)


def test_operator_one_completion_baseline():
    completions = list(operator_completions(EFF_OPERATORS["1"]))
    collected = collect_completions(completions)

    assert len(completions) == 8
    assert Counter(
        (completion.topology, completion.canonical_topology)
        for completion in completions
    ) == Counter({("2s2f_1", "2s2f_1"): 4, ("2s2f_2", "2s2f_2"): 4})
    assert sorted(map(len, collected.values())) == [2, 2, 4]
    assert completion_digest(completions) == (
        "7870cabf4e444ff946702b676463d7d497eb91b8da2a08aa91e07129c839d365"
    )
