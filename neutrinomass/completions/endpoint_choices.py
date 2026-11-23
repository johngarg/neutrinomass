"""Choose mass or kinetic transport on every internal Weyl line explicitly."""

from contextlib import contextmanager
from contextvars import ContextVar
from itertools import combinations
from copy import deepcopy

_kinetic_edges=ContextVar("physical_kinetic_edges",default=frozenset())


@contextmanager
def kinetic_edges(edges):
    token=_kinetic_edges.set(frozenset(tuple(sorted(e)) for e in edges))
    try:
        yield
    finally:
        _kinetic_edges.reset(token)


def transported_endpoint(field,edge):
    if tuple(sorted(edge)) not in _kinetic_edges.get():
        return field
    from neutrinomass.completions.core import VectorLikeDiracFermion,MajoranaFermion
    if isinstance(field,VectorLikeDiracFermion):
        return field.dirac_partner().conj
    if isinstance(field,MajoranaFermion):
        return field.majorana_partner().conj
    raise ValueError("Kinetic numerator requested on a nonfermion line")


def physical_edge_completion_stream(operator,*,second_derivative_basis=None,unique_projection=None):
    from neutrinomass.completions.completions import (
        operator_strip_derivs,construct_operator,partition_stream,partition_leaves,
        canonical_propagator_cut,partition_completion,
    )
    from neutrinomass.completions.core import EffectiveOperator,FailedCompletion
    from neutrinomass.completions.physical_routes import local_physical_expansions
    fields,epsilons,order=operator_strip_derivs(operator.operator).values()
    stripped=EffectiveOperator(operator.name,construct_operator(fields,epsilons))
    for partition in partition_stream(stripped):
        leaves={leaf.node:leaf.field for leaf in partition_leaves(partition["partition"])}
        fermions=[]
        for edge in partition["graph"].edges:
            if set(edge)&leaves.keys():
                continue
            cut=canonical_propagator_cut(partition["partition"],partition["graph"],edge)
            if sum(leaves[node].is_fermion for node in cut)%2:
                fermions.append(tuple(sorted(edge)))
        fermions.sort()
        for count in range(order%2,min(order,len(fermions))+1,2):
            for selected in combinations(fermions,count):
                branch=dict(partition,graph=deepcopy(partition["graph"]),completion_operator=operator)
                with kinetic_edges(selected):
                    completion=partition_completion(branch)
                if isinstance(completion,FailedCompletion):
                    continue
                for candidate in local_physical_expansions(completion,order,
                        second_derivative_basis=second_derivative_basis):
                    actual={tuple(sorted(p.edge)) for p in candidate.momentum_contributions
                            if p.numerator_kind == "momentum"}
                    if actual != set(selected):
                        raise ValueError("Constructed kinetic endpoints differ from the selected physical branch")
                    if unique_projection is not None:
                        candidate.lorentz_projection=unique_projection
                    yield candidate
