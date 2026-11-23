"""Stable ordered generation ranges; no scientific class is pruned here."""

from math import factorial

from neutrinomass.completions.completions import (
    derivative_generation_context, get_topology_data, colour_singlets,
    operator_strip_derivs, construct_operator,
)
from neutrinomass.completions.core import EffectiveOperator
from neutrinomass.database.rebuild import operator_registry, operator_kind, DERIVATIVE_KIND


def generation_plan(operator_name,assignment_size=256):
    if assignment_size < 1:
        raise ValueError("Assignment chunk size must be positive")
    operator=operator_registry()[operator_name]
    producers=[]
    if operator_kind(operator_name) == DERIVATIVE_KIND:
        _,combos,_,_,_=derivative_generation_context(operator)
        producers.extend(("local",i,combo) for i,combo in enumerate(combos)
                         if combo.operator.simplify() != 0)
        fields,epsilons,_=operator_strip_derivs(operator.operator).values()
        producers.append(("routed",None,EffectiveOperator(operator_name,construct_operator(fields,epsilons))))
    else:
        producers.append(("regular",None,operator))
    units=[]
    for producer,placement,combo in producers:
        kind=combo.topology_type
        assignments=factorial(kind["n_scalars"])*factorial(kind["n_fermions"])
        colours=len(colour_singlets([combo.operator],overcomplete=True))
        for topology in get_topology_data(**kind):
            for colour in range(colours):
                for start in range(0,assignments,assignment_size):
                    units.append(dict(producer=producer,placement=placement,
                        topology=topology["canonical_topology"],colour=colour,
                        assignment_start=start,assignment_stop=min(start+assignment_size,assignments)))
    return units
