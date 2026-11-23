"""Candidate-specific heavy T1/T3 Weinberg loops.

These are the finite box/triangle templates of arXiv:1204.5862. Their roles
are occurrences, not automatically distinct species. Every vertex is rebuilt
with the admitted copies before its Bose/Fermi and gauge checks. The external
legs are LLHH; no extra symmetry or zero baryon assignment is imposed.
"""

from collections import defaultdict
from itertools import combinations,product,permutations

import networkx as nx

from neutrinomass.completions.core import Completion,VectorLikeDiracFermion,MajoranaFermion
from neutrinomass.completions.species import CandidateAssignment
from neutrinomass.completions.topologies import Leaf
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.completions.equivalence import field_label_parts
from neutrinomass.completions.routing import state_key
from neutrinomass.completions.amplitudes import induced_amplitude_operator
from neutrinomass.database.serialization import tensor_to_data,tensor_from_data,dumps_completion
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.filtering import ONE_LOOP_SCALE_TEV
from neutrinomass.tensormethod import L,H
from neutrinomass.tensormethod.core import IndexedField,Operator
from neutrinomass.tensormethod.contract import unsimplified_invariants


def paired_endpoints(field):
    if not field.is_fermion:
        return ((field.conj,"scalar"),)
    if isinstance(field,VectorLikeDiracFermion):
        mass=field.dirac_partner()
    elif isinstance(field,MajoranaFermion):
        mass=field.majorana_partner()
    else:
        raise ValueError("UNSUPPORTED_ONE_LOOP_PARTICLE_NATURE")
    return ((mass,"mass"),(field.conj,"momentum"))


def vertex_allowed(fields):
    if sum(field.mass_dim for field in fields)>4:
        return False
    if any(sum(field.charges.get(charge,0) for field in fields) != 0
           for charge in {key for field in fields for key in field.charges}):
        return False
    fermions=[field for field in fields if field.is_fermion]
    return not fermions or (len(fermions)==2 and fermions[0].lorentz_irrep == fermions[1].lorentz_irrep)


def concrete_vertices(fields):
    if not vertex_allowed(fields):
        return ()
    concrete={field.label:field for field in fields if isinstance(field,IndexedField)}
    result=[]
    for term in unsimplified_invariants(*(field.field if isinstance(field,IndexedField) else field
                                         for field in fields),ignore=[]):
        tensors=[]
        for tensor in term.tensors:
            if isinstance(tensor,IndexedField) and tensor.label in concrete:
                data=tensor_to_data(concrete[tensor.label])
                data["indices"]=[str(index) for index in tensor.indices]
                tensor=tensor_from_data(data)
            tensors.append(tensor)
        term=Operator(*tensors)
        if term.safe_simplify() != 0:
            # fun_eval reconstructs generic IndexedFields without their charge
            # and conjugation kwargs in SymPy 1.2. The amplitude contraction
            # freshens its own temporary tensors after mapping these ports.
            result.append(term)
    return tuple(result)


def cycle_assignments(states,external):
    """Close quantum-number-compatible cycles before constructing invariants."""
    size=len(external)
    choices=tuple((left,right,kind) for left in states for right,kind in paired_endpoints(left))
    def extend(edges):
        if len(edges)==size:
            if (sum(edge[2]=="momentum" for edge in edges)%2 == 0
                    and vertex_allowed((edges[-1][1],edges[0][0],*external[0]))):
                yield tuple(edges)
            return
        for edge in choices:
            if vertex_allowed((edges[-1][1],edge[0],*external[len(edges)])):
                yield from extend(edges+[edge])
    for edge in choices:
        yield from extend([edge])


def loop_records(fields):
    states=CandidateAssignment(fields,()).states
    shapes=[]
    for leptons in combinations(range(4),2):
        shapes.append(tuple((L,) if i in leptons else (H,) for i in range(4)))
    shapes.append(((L,),(L,),(H,H)))
    vertex_cache={}
    for external in shapes:
        size=len(external)
        for edges in cycle_assignments(states,external):
            vertices=[]
            for i in range(size):
                inputs=(edges[i-1][1],edges[i][0],*external[i])
                key=tuple((f.label,f.dynkin,tuple(sorted(f.charges.items())),
                           type(f).__name__,getattr(f,"comm",None)) for f in inputs)
                if key not in vertex_cache:
                    # Bound the cache independently of the candidate's size.
                    if len(vertex_cache)>=256:
                        vertex_cache.clear()
                    vertex_cache[key]=concrete_vertices(inputs)
                vertices.append(vertex_cache[key])
            for terms in product(*vertices):
                graph=nx.cycle_graph(size)
                for i,edge in enumerate(edges):
                    graph.edges[i,(i+1)%size]["particle"]=edge[0].label
                leaves=[];port_options=[]
                heavy={field_label_parts(f.label)[0] for f in fields}
                for i,term in enumerate(terms):
                    graph.nodes[i]["interaction"]=i
                    ports=defaultdict(list)
                    for ti,tensor in enumerate(term.tensors):
                        if not isinstance(tensor,IndexedField):
                            continue
                        base=field_label_parts(tensor.label)[0]
                        if base in heavy:
                            ports[base].append(ti)
                        else:
                            node=size+len(leaves)
                            leaves.append(Leaf(tensor,node))
                            graph.add_edge(i,node,particle=tensor.label)
                    incident=defaultdict(list)
                    for edge_index in ((i-1)%size,i):
                        base=field_label_parts(edges[edge_index][0].label)[0]
                        incident[base].append(tuple(sorted((edge_index,(edge_index+1)%size))))
                    options=[]
                    expected={tuple(sorted(((i-1)%size,i))):edges[i-1][1],
                              tuple(sorted((i,(i+1)%size))):edges[i][0]}
                    for assignments in product(*(permutations(incident[base]) for base in sorted(ports))):
                        mapped=[dict(tensor=ti,edge=edge) for base,assigned in zip(sorted(ports),assignments)
                                for ti,edge in zip(ports[base],assigned)]
                        if all(state_key(term.tensors[p["tensor"]]) == state_key(expected[p["edge"]])
                               for p in mapped):
                            options.append(mapped)
                    port_options.append(options)
                for ports in product(*port_options):
                    branch=graph.copy()
                    for i,items in enumerate(ports):
                        branch.nodes[i]["heavy_ports"]=items
                    used={f for edge in edges for f in edge[:2]}
                    yield Completion(EFF_OPERATORS["1"],tuple(leaves),branch,used,list(terms),
                        topology="T1" if size==4 else "T3"),tuple(edge[2] for edge in edges)


def candidate_loop_witness(candidate):
    """Return a fully contracted nonzero heavy-loop template, or exhaust them.

    The fermions form an open chain between the two external L fields. With
    zero external momentum, zero or two kinetic numerators give respectively
    a nonzero mass product or k² times a mass product; odd powers integrate to
    zero. The scalar denominators have generic nonzero masses. Gauge and
    identical-field cancellation are checked on the sewn tensor itself.
    """
    available=species_descriptors(candidate)
    for record,numerators in loop_records(candidate.exotics):
        amplitude=induced_amplitude_operator(record).safe_simplify()
        if amplitude == 0:
            continue
        required=species_descriptors(record)
        if any(available.get(label)!=descriptor for label,descriptor in required.items()):
            raise ValueError("ONE_LOOP_WITNESS_CHANGED_CANDIDATE_SPECIES")
        return dict(code="CANDIDATE_HEAVY_ONE_LOOP_WEINBERG",verified=True,
            scale_tev=ONE_LOOP_SCALE_TEV,required_species=tuple(sorted(required.values(),key=repr)),
            species_assignment=dict(sorted(required.items())),record=dumps_completion(record),
            numerator_kinds=numerators,amplitude_witness=str(amplitude),
            template_domain="heavy T1/T3; generic allowed couplings")
    return None
