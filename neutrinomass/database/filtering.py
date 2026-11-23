"""Versioned filtering contracts and explicit exclusion witnesses.

The published contract intentionally preserves the old loop-table prefix
lookup, conjugation identity and repeated factors. It is a comparison
baseline, not a proof of the completeness of all physical mass mechanisms.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
import math

from neutrinomass.utils.functions import conjugate_field

ONE_LOOP_SCALE_TEV = 605520000000.0 / (16 * math.pi**2)


def conjugacy_key(field):
    return min(field, conjugate_field(field))


def published_model_key(fields):
    from neutrinomass.database.database import EXOTICS
    return tuple(sorted(
        f"published:{int(EXOTICS[field])}" if field in EXOTICS
        else "new:" + conjugacy_key(field)
        for field in fields
    ))


@lru_cache(maxsize=1)
def published_one_loop_templates():
    from neutrinomass.database.database import EXOTICS
    from neutrinomass.database.heavyloops import generate_models
    lookup = defaultdict(list)
    for field in EXOTICS:
        lookup[field.rsplit(",", 1)[0]].append(field)
    result = {}
    for model in generate_models():
        # Exact comma-delimited prefix; these templates contain four components.
        choices = [lookup[field] for field in model]
        if any(not group for group in choices):
            continue
        if any(len(group) != 1 for group in choices):
            raise ValueError("Ambiguous published one-loop baryon assignment")
        fields = tuple(group[0] for group in choices)
        result.setdefault(published_model_key(fields), fields)
    return tuple(sorted(result.items()))


def _subsets(fields):
    return (key for n in range(1, len(fields)+1) for key in combinations(fields, n))


def filter_published_registry(registry, scales):
    if set(registry) != set(scales):
        raise ValueError("Registry and scale inventory differ")
    survivors = {operator: {} for operator in registry}
    witnesses = []
    removed_mass, removed_loop = Counter(), Counter()
    groups = defaultdict(list)
    active = {}
    for operator, scale in scales.items():
        groups[scale].append(operator)
    loops = dict(published_one_loop_templates())
    for scale in sorted(groups, reverse=True):
        additions = {}
        for operator in sorted(groups[scale]):
            for fields, data in sorted(registry[operator].items()):
                key = published_model_key(fields)
                excluded = next((s for s in _subsets(key) if s in active), None)
                if excluded is not None:
                    removed_mass[operator] += 1
                    witnesses.append(dict(operator=operator, fields=fields,
                        code="higher_scale_published_model", mechanism=active[excluded],
                        matched_key=excluded))
                    continue
                # The historical mass stage precedes the one-loop stage.
                additions.setdefault(key, dict(operator=operator, fields=fields,
                                                scale_tev=scale))
                loop = next((s for s in _subsets(key) if s in loops), None)
                if scale < ONE_LOOP_SCALE_TEV and loop is not None:
                    removed_loop[operator] += 1
                    witnesses.append(dict(operator=operator, fields=fields,
                        code="published_one_loop_weinberg", mechanism=dict(
                            fields=loops[loop], scale_tev=ONE_LOOP_SCALE_TEV), matched_key=loop))
                else:
                    survivors[operator][fields] = data
        active.update(additions)
    return dict(mode="published", schema_version=2, scale_unit="TeV",
        survivors=survivors, witnesses=witnesses,
        removed_by_mass=dict(removed_mass), removed_by_one_loop_weinberg=dict(removed_loop),
        one_loop_scale_tev=ONE_LOOP_SCALE_TEV,
        distinct_input_models=len({published_model_key(f) for m in registry.values() for f in m}),
        distinct_surviving_models=len({published_model_key(f) for m in survivors.values() for f in m}))


def multiset_contains(available, required):
    """Repeated occurrences of a species are not additional admitted copies."""
    counts = Counter(available)
    return all(counts[key] >= count for key, count in Counter(required).items())


def filter_species_registry(registry, scales, *, mechanism_search, loop_search):
    """Filter species models with complete candidate-specific search callbacks.

    A search returns a checked witness or None after exhausting its specified
    finite search. It must raise for unsupported/incomplete work. In
    particular, absence from the published representation registry is not a
    permitted substitute for a complete search.
    """
    if set(registry) != set(scales):
        raise ValueError("Registry and scale inventory differ")
    survivors, witnesses = {}, []
    for operator in sorted(registry):
        survivors[operator] = {}
        higher = tuple(sorted(op for op in scales if scales[op] > scales[operator]))
        for species, data in sorted(registry[operator].items()):
            witness = mechanism_search(operator, species, higher)
            if witness is None and scales[operator] < ONE_LOOP_SCALE_TEV:
                witness = loop_search(species)
            if witness is None:
                survivors[operator][species] = data
                continue
            if not multiset_contains(species, witness["required_species"]):
                raise ValueError("Mechanism requires species absent from the candidate")
            if witness.get("verified") is not True:
                raise ValueError("Unverified exclusion witness")
            witnesses.append(dict(operator=operator, species=species, mechanism=witness))
    return dict(mode="species-aware", schema_version=2, scale_unit="TeV",
                survivors=survivors, witnesses=witnesses)
