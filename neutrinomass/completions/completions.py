#!/usr/bin/env python3

"""Functions to generate completions of operators with explicit SU(2) structure."""

from neutrinomass.tensormethod.core import (
    Index,
    Field,
    IndexedField,
    eps,
    delta,
    is_invariant_symbol,
    Operator,
    get_dynkin,
    D,
)
from neutrinomass.tensormethod.contract import (
    lorentz_singlets,
    colour_singlets,
    invariants,
    contract_su2,
)
from neutrinomass.tensormethod.lorentz import (
    LORENTZ_TYPES,
    LorentzBasis,
    LorentzContraction,
    lorentz_contraction,
    lorentz_field_key,
    lorentz_field_port_layout,
    primitive_coordinates,
)

from neutrinomass.utils import timeit
from neutrinomass.tensormethod.utils import safe_nocoeff
from neutrinomass.completions.equivalence import (
    equivalent_lagrangians,
    field_label_parts,
    remapped_field_label,
)
from neutrinomass.completions.utils import (
    flatten,
    chunks,
    factors,
    allowed_lor_dyn,
)
from neutrinomass.utils.functions import remove_equivalent
from neutrinomass.completions.core import (
    Completion,
    DerivativeRoute,
    LorentzProjection,
    MultiDerivativeProjection,
    PropagatorContribution,
    Model,
    FailedCompletion,
    EffectiveOperator,
    cons_completion_field,
    FieldType,
    VectorLikeDiracFermion,
    MajoranaFermion,
    ComplexScalar,
    RealScalar,
)
from neutrinomass.completions.topologies import get_topology_data, Leaf
from neutrinomass.utils import pmatch
from neutrinomass.utils.functions import stringify_qns, conjugate_term

from typing import Tuple, List, Dict, Union, Iterable
import networkx as nx
import networkx.algorithms.isomorphism as iso
from copy import copy, deepcopy
from alive_progress import alive_bar

from collections import Counter, defaultdict
from itertools import (
    permutations,
    combinations,
    combinations_with_replacement,
    product,
)
from sympy.tensor.tensor import Tensor
from sympy import Matrix, Rational, prime

from functools import lru_cache, reduce
from dataclasses import dataclass
import os


def replace(data, to_replace, replace_with, found=False) -> Tuple[tuple, bool]:
    """Replace first occurance of ``to_replace`` with ``replace_with`` in
    ``data``.

    Example:
        >>> replace((("F", 18), ("S", 45), ...), "F", L('u1 i1'))
        ((L(u1, i1), 18), ("S", 45), ...), True

    """
    if found:
        return data, found

    if isinstance(data, tuple):
        new_data = []
        for datai in data:
            new_datai, found = replace(datai, to_replace, replace_with, found)
            new_data.append(new_datai)

        f = lambda x: Leaf(*x) if isinstance(data, Leaf) else tuple(x)
        return f(new_data), found

    if data == to_replace:
        return replace_with, True

    return data, found


def replace_fields(fields: List[IndexedField], partition):
    """Takes the fields and puts them in place of the strings in the partition
    template.

        >>> replace_fields([H('i0_'), H('i1_'), L('u0_ i2_'), L('u1_ i3_')], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        ((L(u0_, i2_), 18), (H(i0_), 162), ((L(u1_, i3_), 6), (H(i1_), 54)))

    """
    for field in fields:
        char = "S" if field.is_boson else "F"
        partition, _ = replace(data=partition, to_replace=char, replace_with=field)

    return partition


def quick_remove_equivalent_partitions(partitions):
    """Just remove double ups. (For now.)

    This is also a good place to remove partitions that you know will be
    filtered out.

    """
    return list(dict.fromkeys(partitions))


def distribute_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template in every possible way.

        >>> distribute_fields([H('i0_'), H('i1_'), L('u0_ i2_'), L('u1_ i3_')], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        [((L(u0_, i2_), 18), (H(i0_), 162), ...), ((L(u1_, i3_), 18), (H(i0_), 162), ...), ...]

    Returns lots of double ups.

    """
    perms = permutations(fields)
    parts = (replace_fields(permutation, partition) for permutation in perms)
    return quick_remove_equivalent_partitions(parts)


def node_dictionary(
    partition: tuple, field_dict: Dict[IndexedField, int]
) -> Dict[int, str]:
    """Returns a dictionary mapping node to indexed field label.

    Example:
        >>> node_dictionary((((Q(u364_, c210_, i369_), 6), (L(u362_, i367_), 18)),
                            ((L(u361_, i366_), 54), (Q(u363_, c209_, i368_), 162)),
                            ((db(u368_, -c214_), 486), (db(u366_, -c212_), 1458))))
        {6: 'Q', 18: 'L', ...}

    """
    flat_data = list(flatten(partition))
    tuples = chunks(flat_data, 2)
    reversed_data = list(map(reversed, tuples))
    return {k: {"particle": v.label + str(field_dict[v])} for k, v in reversed_data}


def set_external_fields(
    partition: tuple, graph: nx.Graph, field_dict: Dict[IndexedField, int]
) -> nx.Graph:
    """Add indexed fields as edge attributes on graph through side effect."""
    g = deepcopy(graph)
    node_attrs = node_dictionary(partition, field_dict)

    edge_attrs = {}
    for edge in graph.edges:
        for n, field_dict in node_attrs.items():
            if n in edge:
                edge_attrs[edge] = field_dict

    nx.set_edge_attributes(g, edge_attrs)
    return g


def indexed_fields_with_counters(op: Operator) -> Dict[IndexedField, int]:
    """Return a dictionary mapping indexed fields to an integer labelling distinct
    fields to help with isomorphism filtering.

    TODO Need to rewrite this to include colour indices! Need to then move
    position of call to include operator with colour structure!

    """
    # idxs are the pairs of contracted isospin indices
    counts = defaultdict(list)
    idxs = []
    for f in op.tensors:
        if isinstance(f, IndexedField):
            counts[f.label].append(f)
        else:
            idxs.append(f.indices)

    labelled_counts = {k: [[f, i] for i, f in enumerate(v)] for k, v in counts.items()}
    for k, v in labelled_counts.items():
        for (f1, i1), (f2, i2) in combinations(v, 2):
            if not f1.indices_by_type["Isospin"]:
                # fields are interchangeable, replace
                f2_idx = labelled_counts[k].index([f2, i2])
                labelled_counts[k][f2_idx] = [f2, i1]
                continue

            iso1 = f1.indices_by_type["Isospin"][0]
            iso2 = f2.indices_by_type["Isospin"][0]
            if [-iso1, -iso2] in idxs or [-iso2, -iso1] in idxs:
                # combination of indices match an epsilon-index pair. In this
                # case, need to replace i2 with i1
                f2_idx = labelled_counts[k].index([f2, i2])
                labelled_counts[k][f2_idx] = [f2, i1]
            else:
                continue

    flat = reduce(lambda x, y: x + y, labelled_counts.values())
    return dict(flat)


def partitions(operator: EffectiveOperator, verbose=False) -> List[dict]:
    """Returns a list of operator partitions, epsilons and graphs of the form:

    {"fields": ((L(u0, I_0), 18), ...)
    "epsilons": (...),
    "graph": ...}

    from the partitions of the fields in the operator. This is all of the
    information required to find the completion.

    """
    topology_data_list = get_topology_data(**operator.topology_type)

    colour_ops = colour_singlets([operator.operator], overcomplete=True)
    colour_ops = [EffectiveOperator(operator.name, op) for op in colour_ops]

    if verbose:
        print(
            f"Finding partitions of {operator.name}. "
            + f"There are {len(colour_ops)} colour structures and "
            + f"{len(topology_data_list)} relevant topologies."
        )

    out = []
    counter = 1
    fields_and_counters = indexed_fields_with_counters(operator.operator)
    fields = list(fields_and_counters)
    for topology_data in topology_data_list:
        if verbose:
            print(f"Furnishing topology {counter}...")
            counter += 1

        perms = distribute_fields(fields, topology_data["partition"])
        for op in colour_ops:
            # col_out = []
            epsilons = op.operator.epsilons

            for perm in perms:
                g = topology_data["graph"]
                g = set_external_fields(perm, g, fields_and_counters)

                data = {
                    "operator": op,
                    "partition": perm,
                    "epsilons": epsilons,
                    "graph": g,
                    "topology": topology_data["topology"],
                    "canonical_topology": topology_data["canonical_topology"],
                }
                out.append(data)

            # if remove_isomorphic_diagrams:
            #     col_out = remove_isomorphic(col_out)

            # out += col_out

    return out


def are_equivalent_partitions(a, b):
    """Checks for partition equivalence by checking if the graphs are isomorphic."""
    ga = a["graph"]
    gb = b["graph"]

    if not iso.faster_could_be_isomorphic(ga, gb):
        return False

    em = iso.categorical_edge_match("particle", "exotic")
    return nx.is_isomorphic(ga, gb, edge_match=em)


def graph_fingerprint(part):
    g = part["graph"]
    degree = dict(g.degree())
    return sorted(degree.values())


def remove_isomorphic(partitions: List[dict]) -> List[dict]:
    """Same algorithm as removeIsomorphic in ``wolfram/`` directory. Remove
    isomorphic graphs (by side effect) to reduce double-ups of completions.

    """
    retained = []
    by_coarse_key = defaultdict(list)
    for partition in partitions:
        graph = partition["graph"]
        edge_particles = Counter(
            nx.get_edge_attributes(graph, "particle").values()
        )
        coarse_key = (
            partition.get("canonical_topology", partition.get("topology")),
            tuple(sorted(dict(graph.degree()).values())),
            tuple(sorted(edge_particles.items())),
        )
        if any(
            are_equivalent_partitions(partition, known)
            for known in by_coarse_key[coarse_key]
        ):
            continue
        retained.append(partition)
        by_coarse_key[coarse_key].append(partition)
    return retained


# The approach to finding the completions is the following: contract off fields
# and find corresponding exotic and term. Replace the fields by the exotic and
# keep track of the available epsilons and the terms by mutation. The pipeline is
#
# contract: returns exotic field, new gauge epsilons (fewer) and new lorentz
# epsilons (more)
#
# replace_and_mutate: returns a Leaf structure that enters the partition in
# place of the contracted fields, mutates terms, edge_dict of graph,
# gauge_epsilons and lorentz_epsilons
#
# reduce_partition: applies replace_and_mutate to a partition until last vertex.


def all_scalars(fields: List[Field]) -> bool:
    """Checks if all fields are scalars."""
    boolean = True
    for f in fields:
        boolean = boolean and f.is_boson

    return boolean


def all_fermions(fields: List[Field]) -> bool:
    """Checks if all fields are fermions."""
    boolean = True
    for f in fields:
        boolean = boolean and f.is_fermion

    return boolean


def drop_scalar(fields: List[Field]) -> List[Field]:
    """Given a list of fields with one scalar, return a list of only the
    fermions, i.e. remove the scalar.

    """
    scalars, fermions = [], []
    for f in fields:
        if f.is_boson:
            scalars.append(f)
        elif f.is_fermion:
            fermions.append(f)
    assert len(scalars) == 1
    return fermions


def get_lorentz_epsilons(fields: Tuple[IndexedField]) -> Tuple[bool, List[Tensor]]:
    """Takes a list of two or three fields (possibly with derivatives) and returns
    the lorentz epsilons that contract the fields to as low a Lorentz irrep as
    possible as well as a boolean indicating whether the contraction is allowed.

    """

    deriv_structure = [f.derivs for f in fields]
    n_derivs = sum(deriv_structure)

    if n_derivs > 2:
        raise Exception(
            f"Not currently supporting {n_derivs} derivatives in an operator."
        )

    if not n_derivs and len(fields) == 4:
        return True, []

    if not n_derivs and len(fields) == 3:
        if all_scalars(fields):
            return True, []

        elif all_fermions(fields):
            return False, []

        return get_lorentz_epsilons(drop_scalar(fields))

    if n_derivs == 2 and len(fields) == 3:
        fields = sorted(fields, key=lambda f: -f.derivs)

    prod = reduce(lambda x, y: x * y, fields)
    undotted, dotted, _, _, _, = prod.indices_by_type.values()

    # Reject vector contraction
    if len(undotted) == 1 and len(dotted) == 1:
        return False, []

    epsilons = []
    for indices in [undotted, dotted]:
        # skip single indices (fermion, scalar) contraction
        if len(indices) == 1:
            continue

        # pair up all even indices; if odd, leave last index
        if len(indices) % 2 != 0:
            indices.pop(-1)

        for i, j in chunks(indices, 2):
            epsilons.append(eps(f"-{i} -{j}"))

    return True, epsilons


def is_vector_fermion_contraction(fields: Tuple[IndexedField]) -> bool:
    """Return whether ``fields`` form a dotted-undotted fermion current.

    Such a current is not a renormalisable scalar interaction by itself.  In a
    derivative completion it can nevertheless be saturated by the momentum
    numerator of an arrow-preserving internal fermion propagator.
    """

    if not any(isinstance(f, FieldType) and f.is_fermion for f in fields):
        return False

    prod = reduce(lambda x, y: x * y, fields)
    undotted, dotted, _, _, _ = prod.indices_by_type.values()
    return len(undotted) == 1 and len(dotted) == 1


def differentiate_indexed_fermion(field: IndexedField) -> IndexedField:
    """Act a slashed derivative on ``field`` while preserving its gauge indices."""

    derivative = D(field.field, allowed_lor_dyn(field))
    fresh = derivative.fresh_indices()
    undotted, dotted, _, _, _ = fresh.indices_by_type.values()
    indices = (*undotted, *dotted, *field.gauge_indices)
    return derivative(" ".join(str(i) for i in indices))


MAX_DERIVATIVE_ROUTE_CHOICES = 2
PROJECTED_LORENTZ_OPERATORS = frozenset(
    {
        "D6a",
        "D6b",
        "D8a",
        "D8b",
        "D8c",
        "D8d",
        "D8e",
        "D8f",
        "D8g",
        "D8h",
        "D8i",
        "D9a",
        "D9b",
        "D12a",
        "D12b",
        "D14a",
        "D14b",
        "D14c",
        "D16a",
        "D16b",
        "D16c",
        "D17",
    }
)

HISTORICAL_IBP_RELATION = (
    "IBP defines the derivative-placement orbit; routed momentum uses the "
    "lexicographically first side of the cut without quotienting placements"
)
HISTORICAL_EOM_RELATION = "No equation-of-motion reduction is applied"
UNIQUE_MULTI_DERIVATIVE_IBP_RELATION = (
    "The named derivative placement is retained; each propagator momentum is "
    "the lexicographically first side of its cut"
)
UNREDUCED_SECOND_DERIVATIVE_IBP_RELATION = (
    "All explicit second-derivative placements generated from the "
    "lexicographically first side of each propagator cut are retained"
)
_MOMENTUM_FIELD_LABEL = "Pmomentum"


def _stripped_occurrence_key(field):
    if not field.derivs:
        return lorentz_field_key(field)
    stripped_field = field.strip_derivs()
    undotted, dotted, _, _, _ = (
        stripped_field.fresh_indices().indices_by_type.values()
    )
    indices = (*undotted, *dotted, *field.gauge_indices)
    stripped = stripped_field(" ".join(map(str, indices)))
    return lorentz_field_key(stripped)


def _normalised_pair(left, right):
    if left <= right:
        return (left, right), Rational(1)
    return (right, left), Rational(-1)


def _box_derivative_field(field):
    """Return the unreduced scalar :math:`D^2` component acting on ``field``."""

    base = field.strip_derivs() if field.derivs else field
    stripped = base.stripped
    if stripped is None:
        stripped = {
            "label": base.label,
            "dynkin": base.dynkin,
            "symmetry": base.symmetry,
            "charges": base.charges,
            "latex": base.latex,
        }
    boxed = Field(
        "D2" + base.label,
        dynkin=base.dynkin,
        charges=base.charges,
        comm=base.comm,
        is_conj=base.is_conj,
        nf=base.nf,
        derivs=2,
        stripped=stripped,
    )
    boxed.latex = f"(D^2 {base.get_latex()})"
    return boxed


def _embed_multi_derivative_contraction(
    contraction,
    derivative_operator,
    ambient_operator,
    momentum_assignments,
):
    """Embed explicit derivative irreps into labelled momentum space.

    ``momentum_assignments`` maps each stripped external-field occurrence to
    the momentum-spurion labels acting on it.  One label represents an
    ordinary derivative.  Two labels represent the unreduced scalar
    :math:`D^2` component.
    """

    derivative_layout, _ = lorentz_field_port_layout(derivative_operator)
    ambient_layout, ambient_counts = lorentz_field_port_layout(ambient_operator)
    momentum_labels = {
        label
        for labels in momentum_assignments.values()
        for label in labels
    }
    ambient_fields = {
        lorentz_field_key(field): ports
        for field, ports in ambient_layout
        if field.label not in momentum_labels
    }
    momentum_ports = {
        field.label: ports
        for field, ports in ambient_layout
        if field.label in momentum_labels
    }

    port_mapping = {short_type: {} for short_type in LORENTZ_TYPES}
    added_pairs = {short_type: [] for short_type in LORENTZ_TYPES}

    for field, derivative_ports in derivative_layout:
        occurrence_key = _stripped_occurrence_key(field)
        ambient_ports = ambient_fields[occurrence_key]
        assigned = tuple(momentum_assignments.get(occurrence_key, ()))
        if field.derivs != len(assigned):
            raise ValueError(
                "Derivative count does not match the momentum assignment"
            )
        for short_type in LORENTZ_TYPES:
            new_ports = derivative_ports[short_type]
            old_ports = ambient_ports[short_type]
            common = min(len(new_ports), len(old_ports))
            port_mapping[short_type].update(
                zip(new_ports[:common], old_ports[:common])
            )
            if not field.derivs:
                if len(new_ports) != len(old_ports):
                    raise ValueError("Non-derivative Lorentz ports changed")
                continue
            if field.derivs == 2:
                if len(new_ports) != len(old_ports):
                    raise ValueError(
                        "The retained D^2 component must preserve Lorentz ports"
                    )
                left, right = assigned
                added_pairs[short_type].append(
                    (
                        momentum_ports[left][short_type][0],
                        momentum_ports[right][short_type][0],
                    )
                )
                continue
            if field.derivs != 1:
                raise ValueError("Only one or two explicit derivatives are supported")
            momentum = momentum_ports[assigned[0]][short_type][0]
            if len(new_ports) == len(old_ports) + 1:
                port_mapping[short_type][new_ports[-1]] = momentum
            elif len(old_ports) == len(new_ports) + 1:
                added_pairs[short_type].append((momentum, old_ports[-1]))
            elif len(new_ports) != len(old_ports):
                raise ValueError("Derivative changes more than one Lorentz port")

    coefficient = Rational(contraction.coefficient)
    pairings = {short_type: [] for short_type in LORENTZ_TYPES}
    for short_type, source_pairs in contraction.pairings:
        for left, right in source_pairs:
            pair, sign = _normalised_pair(
                port_mapping[short_type][left],
                port_mapping[short_type][right],
            )
            coefficient *= sign
            pairings[short_type].append(pair)
    for short_type, source_pairs in added_pairs.items():
        for left, right in source_pairs:
            pair, sign = _normalised_pair(left, right)
            coefficient *= sign
            pairings[short_type].append(pair)

    for short_type, count in ambient_counts:
        used = sorted(port for pair in pairings[short_type] for port in pair)
        if used != list(range(count)):
            raise ValueError("Embedded derivative contraction leaves open ports")

    return LorentzContraction(
        coefficient=coefficient,
        pairings=tuple(
            (short_type, tuple(sorted(pairings[short_type])))
            for short_type in LORENTZ_TYPES
        ),
        port_counts=ambient_counts,
    )


def _embed_derivative_contraction(
    contraction, derivative_operator, ambient_operator, momentum_label
):
    """Embed one historical derivative contraction into momentum space."""

    derivative = next(
        field for field in derivative_operator.indexed_fields if field.derivs
    )
    return _embed_multi_derivative_contraction(
        contraction,
        derivative_operator,
        ambient_operator,
        {_stripped_occurrence_key(derivative): (momentum_label,)},
    )


@dataclass(frozen=True)
class HistoricalDerivativePlacement:
    index: int
    label: str
    occurrence_key: tuple
    operator: EffectiveOperator
    basis: LorentzBasis
    coordinate_start: int
    ambient_matrix: Matrix

    @property
    def dimension(self):
        return self.basis.dimension

    def project_ambient(self, contraction):
        target = Matrix(contraction.evaluation_vector())
        coordinates, parameters = self.ambient_matrix.gauss_jordan_solve(target)
        if parameters.rows:
            raise ValueError("Historical ambient basis is not independent")
        return tuple(Rational(value) for value in coordinates[: self.dimension])


@dataclass(frozen=True)
class HistoricalDerivativeBasis:
    requested_operator: EffectiveOperator
    stripped_operator: EffectiveOperator
    momentum_field: Field
    placements: tuple
    labels: tuple

    @classmethod
    def from_operator(cls, requested_operator):
        fields, epsilons, n_derivs = operator_strip_derivs(
            requested_operator.operator
        ).values()
        if n_derivs != 1:
            raise ValueError("Historical routing requires exactly one derivative")
        stripped_operator = EffectiveOperator(
            requested_operator.name, construct_operator(fields, epsilons)
        )
        momentum_field = Field(
            _MOMENTUM_FIELD_LABEL,
            dynkin="11000",
            charges={"y": 0, "3b": 0},
            latex="p",
        )
        ambient_operator = Operator(
            *stripped_operator.operator.tensors,
            momentum_field.fresh_indices(),
        )
        _, ambient_port_counts = lorentz_field_port_layout(ambient_operator)
        ambient_basis = LorentzBasis.from_port_counts(ambient_port_counts)

        placements = []
        labels = []
        coordinate_start = 0
        for spec in derivative_placement_combinations(requested_operator):
            basis = LorentzBasis.from_operator(spec.operator.operator)
            singlets = {
                lorentz_contraction(singlet, normalise=True).label: singlet
                for singlet in lorentz_singlets(spec.operator.operator)
                if singlet.safe_simplify() != 0
            }
            embedded_vectors = [
                _embed_derivative_contraction(
                    lorentz_contraction(singlets[label], normalise=True),
                    spec.operator.operator,
                    ambient_operator,
                    momentum_field.label,
                ).evaluation_vector()
                for label in basis.labels
            ]
            extension = list(embedded_vectors)
            rank = Matrix.hstack(*(Matrix(vector) for vector in extension)).rank()
            for vector in ambient_basis.vectors:
                trial = Matrix.hstack(
                    *(Matrix(item) for item in (*extension, vector))
                )
                trial_rank = trial.rank()
                if trial_rank == rank:
                    continue
                extension.append(vector)
                rank = trial_rank
            if rank != ambient_basis.dimension:
                raise ValueError("Historical placement does not extend to a basis")

            placement_label = f"p{spec.index}:{spec.field_label}"
            labels.extend(
                f"{placement_label}|{label}" for label in basis.labels
            )
            placements.append(
                HistoricalDerivativePlacement(
                    index=spec.index,
                    label=placement_label,
                    occurrence_key=spec.occurrence_key,
                    operator=spec.operator,
                    basis=basis,
                    coordinate_start=coordinate_start,
                    ambient_matrix=Matrix.hstack(
                        *(Matrix(vector) for vector in extension)
                    ),
                )
            )
            coordinate_start += basis.dimension

        return cls(
            requested_operator=requested_operator,
            stripped_operator=stripped_operator,
            momentum_field=momentum_field,
            placements=tuple(placements),
            labels=tuple(labels),
        )

    def _metadata(self, coordinates):
        derivative = next(
            field
            for field in self.requested_operator.operator.indexed_fields
            if field.derivs
        )
        return LorentzProjection(
            basis_labels=self.labels,
            coordinates=tuple(str(value) for value in coordinates),
            derivative_field=derivative.label,
            ibp_relation=HISTORICAL_IBP_RELATION,
            eom_relation=HISTORICAL_EOM_RELATION,
        )

    def project_local(self, placement, operator):
        local = placement.basis.project(operator, primitive=False)
        if local is None:
            raise ValueError("Local contraction is outside its Lorentz basis")
        coordinates = [Rational(0)] * len(self.labels)
        start = placement.coordinate_start
        coordinates[start : start + placement.dimension] = local
        return self._metadata(primitive_coordinates(coordinates))

    def placement_for_operator(self, operator):
        derivative_fields = [field for field in operator.indexed_fields if field.derivs]
        if len(derivative_fields) != 1:
            raise ValueError("Historical local record must contain one derivative")
        occurrence_key = _stripped_occurrence_key(derivative_fields[0])
        return next(
            placement
            for placement in self.placements
            if placement.occurrence_key == occurrence_key
        )

    def project_existing_local(self, operator):
        return self.project_local(self.placement_for_operator(operator), operator)

    def _closed_routed_contraction(self, explicit):
        free_lorentz = {}
        for short_type in LORENTZ_TYPES:
            index_type = Index.get_index_types()[short_type]
            indices = [
                index
                for index in explicit.free_indices
                if index.index_type == index_type
            ]
            if len(indices) != 1:
                raise ValueError("Routed contraction must expose one momentum port")
            free_lorentz[short_type] = indices[0]

        momentum_indices = []
        closing_epsilons = []
        for short_type in LORENTZ_TYPES:
            free_index = free_lorentz[short_type]
            if not free_index.is_up:
                momentum_indices.append(-free_index)
                continue
            fresh = Index.fresh(short_type)
            momentum_indices.append(fresh)
            closing_epsilons.append(eps(f"{-free_index} {-fresh}"))
        momentum = self.momentum_field(" ".join(map(str, momentum_indices)))
        closed = Operator(*explicit.tensors, momentum, *closing_epsilons)
        return lorentz_contraction(closed)

    def project_routed(self, explicit, partition, graph, route):
        contraction = self._closed_routed_contraction(explicit)
        cut_graph = graph.copy()
        cut_graph.remove_edge(*route.edge)
        external = {leaf.node: leaf.field for leaf in partition_leaves(partition)}
        sides = []
        for component in nx.connected_components(cut_graph):
            nodes = tuple(sorted(node for node in component if node in external))
            if nodes:
                sides.append(nodes)
        if len(sides) != 2:
            raise ValueError("A routed propagator must split external fields in two")
        side = min(sides)

        by_occurrence = {
            placement.occurrence_key: placement for placement in self.placements
        }
        coordinates = [Rational(0)] * len(self.labels)
        for node in side:
            placement = by_occurrence.get(_stripped_occurrence_key(external[node]))
            if placement is None:
                continue
            local = placement.project_ambient(contraction)
            start = placement.coordinate_start
            for offset, value in enumerate(local):
                coordinates[start + offset] += value
        coordinates = primitive_coordinates(coordinates)
        if not any(coordinates):
            return None
        return self._metadata(coordinates)


def _derivative_irrep_dynkins(field):
    """Return every Lorentz irrep in ``(1, 1) x field``.

    The completion database retains the irrep chosen by ``allowed_lor_dyn``.
    The remaining irreps provide the representation-theoretic complement
    needed to extract that component without imposing equations of motion.
    """

    undotted, dotted = map(int, field.dynkin[:2])
    candidates = []
    for delta_undotted in (-1, 1):
        for delta_dotted in (-1, 1):
            new_undotted = undotted + delta_undotted
            new_dotted = dotted + delta_dotted
            if new_undotted < 0 or new_dotted < 0:
                continue
            candidates.append(f"{new_undotted}{new_dotted}")
    allowed = allowed_lor_dyn(field)
    return tuple(sorted(candidates, key=lambda item: (item != allowed, item)))


def _indexed_stripped_occurrence(field, gauge_indices):
    undotted, dotted, _, _, _ = field.fresh_indices().indices_by_type.values()
    lorentz = " ".join(map(str, (*undotted, *dotted)))
    indexed = field(
        " ".join(part for part in (lorentz, gauge_indices) if part)
    )
    return lorentz_field_key(indexed)


def _derivative_irrep_operator(fields, epsilons, assignments):
    structure = []
    for index, (field, indices) in enumerate(fields):
        assignment = assignments.get(index)
        if assignment == "box":
            placed = _box_derivative_field(field)
        elif assignment is not None:
            placed = D(field, assignment)
        else:
            placed = field
        structure.append((placed, indices))
    return construct_operator(structure, epsilons)


def _basis_singlets(operator, basis):
    singlets = {
        lorentz_contraction(singlet, normalise=True).label: singlet
        for singlet in lorentz_singlets(operator)
        if singlet.safe_simplify() != 0
    }
    return tuple(singlets[label] for label in basis.labels)


@dataclass(frozen=True)
class SecondDerivativePlacement:
    indices: tuple
    occurrence_keys: tuple
    label: str
    operator: EffectiveOperator
    basis: LorentzBasis
    coordinate_start: int

    @property
    def dimension(self):
        return self.basis.dimension


@dataclass
class UnreducedSecondDerivativeBasis:
    """Canonical placement basis for two derivatives with no EOM reduction."""

    requested_operator: EffectiveOperator
    stripped_operator: EffectiveOperator
    fields: tuple
    epsilons: tuple
    occurrence_keys: tuple
    placements: tuple
    placement_by_indices: dict
    labels: tuple
    momentum_fields: tuple
    ambient_operator: Operator
    cross_matrices: dict

    @classmethod
    def from_operator(cls, requested_operator):
        fields, epsilons, n_derivs = operator_strip_derivs(
            requested_operator.operator
        ).values()
        if n_derivs != 2:
            raise ValueError("The unreduced placement basis requires two derivatives")
        fields = tuple(fields)
        epsilons = tuple(epsilons)
        stripped_operator = EffectiveOperator(
            requested_operator.name, construct_operator(fields, epsilons)
        )
        occurrence_keys = tuple(
            _indexed_stripped_occurrence(field, indices)
            for field, indices in fields
        )

        placements = []
        labels = []
        coordinate_start = 0
        for left, right in combinations_with_replacement(range(len(fields)), 2):
            assignments = (
                {left: "box"}
                if left == right
                else {
                    left: allowed_lor_dyn(fields[left][0]),
                    right: allowed_lor_dyn(fields[right][0]),
                }
            )
            derivative_operator = _derivative_irrep_operator(
                fields, epsilons, assignments
            )
            if derivative_operator.safe_simplify() == 0:
                continue
            try:
                basis = LorentzBasis.from_operator(derivative_operator)
            except ValueError:
                continue
            field_labels = tuple(
                fields[index][0].label_with_dagger for index in (left, right)
            )
            placement_label = (
                f"p{left},{right}:" + ",".join(field_labels)
            )
            labels.extend(
                f"{placement_label}|{label}" for label in basis.labels
            )
            placements.append(
                SecondDerivativePlacement(
                    indices=(left, right),
                    occurrence_keys=(
                        occurrence_keys[left],
                        occurrence_keys[right],
                    ),
                    label=placement_label,
                    operator=EffectiveOperator(
                        requested_operator.name, derivative_operator
                    ),
                    basis=basis,
                    coordinate_start=coordinate_start,
                )
            )
            coordinate_start += basis.dimension

        momentum_fields = tuple(
            Field(
                f"{_MOMENTUM_FIELD_LABEL}{index}",
                dynkin="11000",
                charges={"y": 0, "3b": 0},
                latex=f"p_{index}",
            )
            for index in range(2)
        )
        ambient_operator = Operator(
            *stripped_operator.operator.tensors,
            *(momentum.fresh_indices() for momentum in momentum_fields),
        )
        _, ambient_port_counts = lorentz_field_port_layout(ambient_operator)
        ambient_dimension = LorentzBasis.from_port_counts(
            ambient_port_counts
        ).dimension

        placement_by_indices = {
            placement.indices: placement for placement in placements
        }
        cross_matrices = {}
        momentum_labels = tuple(field.label for field in momentum_fields)
        for placement in placements:
            left, right = placement.indices
            if left == right:
                continue
            allowed_dynkins = (
                allowed_lor_dyn(fields[left][0]),
                allowed_lor_dyn(fields[right][0]),
            )
            irrep_choices = [
                allowed_dynkins,
                *(
                    choice
                    for choice in product(
                        _derivative_irrep_dynkins(fields[left][0]),
                        _derivative_irrep_dynkins(fields[right][0]),
                    )
                    if choice != allowed_dynkins
                ),
            ]
            for ordered_indices in ((left, right), (right, left)):
                assignments = {
                    occurrence_keys[ordered_indices[0]]: (momentum_labels[0],),
                    occurrence_keys[ordered_indices[1]]: (momentum_labels[1],),
                }
                extension = []
                rank = 0
                target_dimension = 0
                for choice_index, dynkins in enumerate(irrep_choices):
                    derivative_operator = _derivative_irrep_operator(
                        fields,
                        epsilons,
                        {left: dynkins[0], right: dynkins[1]},
                    )
                    if derivative_operator.safe_simplify() == 0:
                        if choice_index == 0:
                            raise ValueError(
                                "The retained derivative placement vanished"
                            )
                        continue
                    try:
                        basis = LorentzBasis.from_operator(derivative_operator)
                    except ValueError:
                        continue
                    vectors = [
                        _embed_multi_derivative_contraction(
                            lorentz_contraction(singlet, normalise=True),
                            derivative_operator,
                            ambient_operator,
                            assignments,
                        ).evaluation_vector()
                        for singlet in _basis_singlets(
                            derivative_operator, basis
                        )
                    ]
                    for vector in vectors:
                        trial = Matrix.hstack(
                            *(Matrix(item) for item in (*extension, vector))
                        )
                        trial_rank = trial.rank()
                        if trial_rank == rank:
                            if choice_index == 0:
                                raise ValueError(
                                    "Retained derivative basis is not independent"
                                )
                            continue
                        extension.append(vector)
                        rank = trial_rank
                    if choice_index == 0:
                        target_dimension = len(extension)
                        if target_dimension != placement.dimension:
                            raise ValueError(
                                "Retained derivative basis changed on embedding"
                            )
                    if rank == ambient_dimension:
                        break
                if rank != ambient_dimension:
                    raise ValueError(
                        "Derivative irreps do not span the momentum space"
                    )
                cross_matrices[ordered_indices] = Matrix.hstack(
                    *(Matrix(vector) for vector in extension)
                )

        return cls(
            requested_operator=requested_operator,
            stripped_operator=stripped_operator,
            fields=fields,
            epsilons=epsilons,
            occurrence_keys=occurrence_keys,
            placements=tuple(placements),
            placement_by_indices=placement_by_indices,
            labels=tuple(labels),
            momentum_fields=momentum_fields,
            ambient_operator=ambient_operator,
            cross_matrices=cross_matrices,
        )

    def _metadata(self, coordinates):
        derivative_fields = tuple(
            field.label
            for field in self.requested_operator.operator.indexed_fields
            if field.derivs
        )
        return MultiDerivativeProjection(
            basis_labels=self.labels,
            coordinates=tuple(str(value) for value in coordinates),
            derivative_fields=derivative_fields,
            ibp_relation=UNREDUCED_SECOND_DERIVATIVE_IBP_RELATION,
            eom_relation=HISTORICAL_EOM_RELATION,
        )

    def _placement_coordinates(self, placement, operator):
        local = placement.basis.project(operator, primitive=False)
        if local is None:
            raise ValueError("Contraction is outside its derivative-placement basis")
        coordinates = [Rational(0)] * len(self.labels)
        start = placement.coordinate_start
        coordinates[start : start + placement.dimension] = local
        return coordinates

    def placement_for_operator(self, operator):
        derivative_fields = [
            field for field in operator.indexed_fields if field.derivs
        ]
        if sum(field.derivs for field in derivative_fields) != 2:
            raise ValueError("Local record must contain two derivatives")
        indices = tuple(
            sorted(
                self.occurrence_keys.index(_stripped_occurrence_key(field))
                for field in derivative_fields
                for _ in range(field.derivs)
            )
        )
        return self.placement_by_indices[indices]

    def project_local(self, operator):
        placement = self.placement_for_operator(operator)
        coordinates = self._placement_coordinates(placement, operator)
        return self._metadata(primitive_coordinates(coordinates))

    def project_existing_local(self, operator):
        return self.project_local(operator)

    def _closed_momentum_square(self, explicit):
        momenta = tuple(
            momentum.fresh_indices() for momentum in self.momentum_fields
        )
        closing_epsilons = []
        for short_type in LORENTZ_TYPES:
            ports = [
                momentum.indices_by_type[Index.get_index_types()[short_type]][0]
                for momentum in momenta
            ]
            closing_epsilons.append(eps(f"{-ports[0]} {-ports[1]}"))
        closed = Operator(*explicit.tensors, *momenta, *closing_epsilons)
        return lorentz_contraction(closed)

    def _project_box(self, explicit, placement):
        occurrence_key = placement.occurrence_keys[0]
        target = next(
            field
            for field in explicit.indexed_fields
            if _stripped_occurrence_key(field) == occurrence_key
        )
        boxed = _box_derivative_field(target.field)(
            " ".join(map(str, target.indices))
        )
        boxed_operator = Operator(
            *(boxed if tensor is target else tensor for tensor in explicit.tensors)
        )
        local = placement.basis.project(boxed_operator, primitive=False)
        if local is None:
            raise ValueError("D^2 contraction is outside its placement basis")
        return local

    def _project_cross(self, contraction, ordered_indices):
        placement = self.placement_by_indices[tuple(sorted(ordered_indices))]
        matrix = self.cross_matrices[ordered_indices]
        target = Matrix(contraction.evaluation_vector())
        coordinates, parameters = matrix.gauss_jordan_solve(target)
        if parameters.rows:
            raise ValueError("Derivative-irrep decomposition is not independent")
        return tuple(
            Rational(value) for value in coordinates[: placement.dimension]
        )

    def project_denominator(self, explicit, partition, contribution):
        if contribution.denominator_order != 1:
            raise ValueError("Second-derivative projection requires one p^2 term")
        external = {
            leaf.node: self.occurrence_keys.index(
                _stripped_occurrence_key(leaf.field)
            )
            for leaf in partition_leaves(partition)
        }
        coordinates = [Rational(0)] * len(self.labels)
        momentum_square = self._closed_momentum_square(explicit)
        for left_node in contribution.cut_side:
            for right_node in contribution.cut_side:
                ordered_indices = (external[left_node], external[right_node])
                placement = self.placement_by_indices.get(
                    tuple(sorted(ordered_indices))
                )
                if placement is None:
                    continue
                if ordered_indices[0] == ordered_indices[1]:
                    local = self._project_box(explicit, placement)
                else:
                    local = self._project_cross(
                        momentum_square, ordered_indices
                    )
                start = placement.coordinate_start
                for offset, value in enumerate(local):
                    coordinates[start + offset] += value
        coordinates = primitive_coordinates(coordinates)
        if not any(coordinates):
            return None
        return self._metadata(coordinates)


def derivative_route_candidates(fields):
    """Return a deterministic ordering of differentiable fermion fields."""

    candidates = [field for field in fields if field.is_fermion and not field.derivs]
    return sorted(
        candidates,
        key=lambda field: (
            isinstance(field, FieldType),
            field.label,
            field.dynkin,
            tuple(str(index) for index in field.indices),
        ),
    )


def derivative_route_alternatives(fields):
    """Return every admissible heavy-fermion numerator at a vector current."""

    fields = tuple(fields)
    if not is_vector_fermion_contraction(fields):
        return ()

    internal_fermions = tuple(
        field
        for field in fields
        if isinstance(field, FieldType) and field.is_fermion
    )
    alternatives = []
    selected_numerators = set()
    for candidate in derivative_route_candidates(fields):
        if isinstance(candidate, FieldType):
            numerator = candidate
        elif len(internal_fermions) == 1:
            numerator = internal_fermions[0]
        else:
            continue

        if numerator in selected_numerators:
            continue

        differentiated = differentiate_indexed_fermion(candidate)
        routed_fields = tuple(
            differentiated if field is candidate else field for field in fields
        )
        allowed, lorentz_epsilons = get_lorentz_epsilons(routed_fields)
        if not allowed:
            continue

        selected_numerators.add(numerator)
        alternatives.append(
            (routed_fields, lorentz_epsilons, (numerator, candidate))
        )

    assert len(alternatives) <= MAX_DERIVATIVE_ROUTE_CHOICES
    return tuple(alternatives)


def route_derivative_to_internal_fermion(fields, derivative_state):
    """Select the next explicitly requested heavy-fermion momentum branch."""

    if derivative_state is None or derivative_state["remaining"] <= 0:
        return None

    alternatives = derivative_route_alternatives(fields)
    route_index = derivative_state.get("consumed", 0)
    route_choices = derivative_state.get(
        "route_choices", (derivative_state.get("route_choice", 0),)
    )
    if route_index >= len(route_choices):
        return None
    route_choice = route_choices[route_index]
    if route_choice >= len(alternatives):
        return None
    return alternatives[route_choice]


def prepare_lorentz_contraction(fields, derivative_state=None):
    """Return fields and epsilon tensors for a local Lorentz contraction."""

    fields = tuple(fields)
    allowed, lorentz_epsilons = get_lorentz_epsilons(fields)
    if allowed:
        return fields, lorentz_epsilons, None

    routed = route_derivative_to_internal_fermion(fields, derivative_state)
    if routed is not None:
        return routed

    return None


def consume_routed_derivative(derivative_state, derivative_route):
    """Consume a routed derivative after its interaction has been validated."""

    if derivative_route is None:
        return

    if derivative_state["pending_route"] is not None:
        raise ValueError("A routed derivative was not assigned to an edge")
    derivative_state["remaining"] -= 1
    derivative_state["consumed"] += 1
    derivative_state["pending_route"] = derivative_route


def record_routed_derivative_edges(derivative_state, edge_dict):
    """Resolve one selected heavy fermion to its graph edge and metadata."""

    pending_route = derivative_state["pending_route"]
    if pending_route is None:
        return

    numerator, differentiated_field = pending_route
    edge = edge_dict.get(numerator)
    assert edge is not None
    derivative_state["routes"].append(
        DerivativeRoute(
            edge=edge,
            numerator_field=numerator.label,
            numerator_lorentz=numerator.dynkin[:2],
            differentiated_field=differentiated_field.label,
            differentiated_lorentz=differentiated_field.dynkin[:2],
        )
    )
    derivative_state["pending_route"] = None


def is_contracted_epsilon(eps: Tensor, indices: List[Index]) -> bool:
    """Return True if two indices on epsilon are contracted, False otherwise."""
    i, j, *k = eps.indices

    # deal with su2 epsilon first
    if not k:
        if -i in indices and -j in indices:
            return True
        return False

    # su3 epsilon has three indices
    to_remove, free = [], []
    for idx in eps.indices:
        if -idx in indices:
            to_remove.append(-idx)
        else:
            free.append(idx)

    # is contracted epsilon
    if len(to_remove) == 2:
        assert len(free) == 1
        indices.append(free[0])
        return True

    if len(to_remove) == 3:
        return True

    return False


def separate_gauge_epsilons(
    fields: List[IndexedField], epsilons: List[Tensor]
) -> Tuple[List[Tensor], List[Tensor]]:
    """Return a 2-tuple with the spectator epsilon tensors carrying gauge indices
    and those that are contracted with the fields passed in.

    """

    prod_fields = reduce(lambda x, y: x * y, fields)
    free_indices = prod_fields.get_free_indices()
    contracted_epsilons, spectator_epsilons = [], []

    # contracted epsilons are those all of whose indices are contracted on the
    # set of fields passed in. Spectator epsilons are all of the others
    for epsilon in epsilons:
        if is_contracted_epsilon(epsilon, free_indices):
            contracted_epsilons.append(epsilon)
        else:
            spectator_epsilons.append(epsilon)

    return spectator_epsilons, contracted_epsilons


def check_charges(operator: Operator, ignore=[]) -> None:
    """Make sure sum of charges vanishes, i.e. term is a U(1) singlet."""
    fields = [f for f in operator.tensors if isinstance(f, IndexedField)]
    for q in fields[0].charges:
        if q in ignore:
            continue
        assert not sum(f.charges[q] for f in fields)


def check_singlet(operator: Operator, ignore=["3b"]) -> None:
    """Make sure operator is a SM and Lorentz singlet."""
    check_charges(operator, ignore=ignore)
    for free in operator.free_indices:
        assert free.index_type == "Generation"


def is_singlet(operator: Operator, ignore=("3b",)) -> bool:
    """Return whether an operator is a Lorentz and SM singlet."""

    fields = [field for field in operator.tensors if isinstance(field, IndexedField)]
    for charge in fields[0].charges:
        if charge not in ignore and sum(field.charges[charge] for field in fields):
            return False
    return all(index.index_type == "Generation" for index in operator.free_indices)


def is_vanishing_interaction(term: Operator) -> bool:
    """Return whether tensor symmetries make an interaction vanish."""

    return term.safe_simplify() == 0


def exotic_field_and_term(
    op: Operator, symbols: Dict[str, List[str]], field_dict: Dict[tuple, str]
) -> Tuple[IndexedField, IndexedField, Union[Operator, str]]:
    """Returns exotic field, partner (that couples in Lagrangian) and Lagrangian
    term. Mutates field_dict with exotic field's symbol.

    The last item returned may be a string explaining why the contraction
    failed.

    """

    fields = [f for f in op.tensors if isinstance(f, IndexedField)]
    exotic_charges = {}
    pairs = list(map(lambda x: x.charges.items(), fields))

    # ensure no fields have missing or extra charges
    fst, *rst = pairs
    for pair in rst:
        assert len(pair) == len(fst)

    for n_pairs in zip(*pairs):
        # fish out key
        (k, _), *rst = n_pairs
        # sum over values
        exotic_charges[k] = sum(map(lambda x: x[1], n_pairs))

    indices_by_type = op.indices_by_type.values()
    exotic_undotted, exotic_dotted, exotic_colour, exotic_isospin, _ = map(
        sorted, indices_by_type
    )

    exotic_indices = " ".join(
        str(i)
        for i in [*exotic_undotted, *exotic_dotted, *exotic_colour, *exotic_isospin]
    )

    # establish fermion or boson for symbols
    lorentz_dynkin = get_dynkin(exotic_indices)[:2]
    if lorentz_dynkin in ("10", "01"):
        symbols_to_use = symbols["fermion"]
    else:
        symbols_to_use = symbols["boson"]

    fs = tuple(sorted([f.field for f in fields], key=lambda x: x.label_with_dagger))
    field_key = (
        fs,
        get_dynkin(exotic_indices),
        tuple(sorted(exotic_charges.items())),
    )
    if field_key in field_dict:
        symbol = field_dict[field_key]
    else:
        symbol = symbols_to_use.pop(0)
        # field_dict mutated here!
        field_dict[field_key] = symbol

    # for Dirac and Majorana fermions, always keep plain symbol left handed
    to_conj = False
    if exotic_indices and exotic_indices[0] == "d":
        exotic_indices = Index.conj_index_string(exotic_indices)
        to_conj = True

    exotic_indexed_field = IndexedField(
        label=symbol, indices=exotic_indices, charges=exotic_charges
    )

    # construct MajoranaFermion, VectorLikeDiracFermion, ...
    # `partner` is in the Lagrangian, `exotic_field` transforms like contracted pair
    exotic_field = cons_completion_field(exotic_indexed_field)
    exotic_field = exotic_field.conj_indices if to_conj else exotic_field

    partner = exotic_field
    if isinstance(exotic_field, ComplexScalar):
        partner = exotic_field.conj
    elif isinstance(exotic_field, RealScalar):
        partner = exotic_field.swap_colour_indices()
    elif isinstance(exotic_field, VectorLikeDiracFermion):
        partner = exotic_field.dirac_partner()
    elif isinstance(exotic_field, MajoranaFermion):
        partner = exotic_field.majorana_partner()

    # Need additional su2 epsilons to fix su2 indices (since not working with
    # lowered indices at all). Won't need to do this if removing a derivative in
    # the process
    partner, fix_su2_epsilons = partner.lower_su2()
    term = reduce(lambda x, y: x * y, fix_su2_epsilons, op * partner)

    # construct term and check to see if vanishes. This is a very costly step,
    # check first whether there are any doubled up fields in the term and only
    # run on those
    if is_vanishing_interaction(term):
        return exotic_field, partner, f"Vanishing coupling at {term}"

    check_singlet(term)

    return exotic_field, partner, term


def process_derivative_term(op: Operator) -> Union[Operator, str]:
    """Process term containing derivatives, return corresponding term that would
    appear in the Lagrangian.

    """
    deriv_structure = [f.derivs for f in op.fields]
    n_derivs = sum(deriv_structure)

    if n_derivs == 0:
        return op

    # Remove derivatives and lorentz epsilons and call contract_su2 on Lorentz
    # structure.
    #
    # There are a number of cases here
    # 1. (DH)(DH)SS
    # 2. (DH)(DH)S
    # 3. (DH)ψF -> in all but this case, affected fields are clear
    # 4. S(Dψ)F
    # 5. (Dψ)(Dψ)S
    # 6. S(Dψ)ψ
    scalars, fermions, exotic_fermions, epsilons = [], [], [], []
    for t in op.tensors:
        if isinstance(t, IndexedField) and t.derivs > 0 and t.is_boson:
            no_deriv_field = t.strip_derivs_with_indices()
            scalars.append(no_deriv_field)
        elif isinstance(t, IndexedField) and t.derivs > 0 and t.is_fermion:
            no_deriv_field = t.strip_derivs_with_indices()
            fermions.append(no_deriv_field)
        elif isinstance(t, VectorLikeDiracFermion):
            fixed_field = t.dirac_partner().conj
            exotic_fermions.append(fixed_field)
        elif isinstance(t, MajoranaFermion):
            fixed_field = t.conj
            exotic_fermions.append(fixed_field)
        elif isinstance(t, IndexedField) and t.is_fermion:
            fermions.append(t)
        elif isinstance(t, IndexedField) and t.is_scalar:
            scalars.append(t)
        elif not isinstance(t, IndexedField):
            # is epsilon, keep gauge ones, not lorentz
            if t.indices[0].index_type in ("Undotted", "Dotted"):
                continue
            else:
                epsilons.append(t)

    # cases 1 and 2
    if len(scalars) > 2:
        term = reduce(lambda x, y: x * y, scalars + epsilons)
        if term.safe_simplify() == 0:
            return "Vanishing structure"
        return term
    # case 6
    if len(fermions) == 2 and n_derivs == 1:
        return "Not allowed contraction"
    # case 5
    if len(fermions) == 2 and n_derivs == 2:
        left, right = fermions
    # cases 3 and 4
    if len(exotic_fermions) == 1:
        assert len(fermions) == 1
        left, right = exotic_fermions[0], fermions[0]
    if len(exotic_fermions) == 2:
        left, right = exotic_fermions

    # include scalars and epsilons in su2 contraction
    right = reduce(lambda x, y: x * y, scalars + epsilons, right)
    lu, ld, _, _, _ = left.indices_by_type.values()
    ru, rd, _, _, _ = right.indices_by_type.values()

    # if the indices are equal after taking the conj, then there will be an
    # error. In this case, you can just lower one of them
    if lu == ru and ld == rd:
        partner, fix_su2_epsilons = left.lower_su2(skip=["Isospin"])
        assert len(fix_su2_epsilons) == 1
        return right * partner * fix_su2_epsilons[0]
    if lu:
        if not (len(lu) == 1 and len(ru) == 1):
            return "Not allowed contraction"
        index_str = " ".join(str(-i) for i in lu + ru)
    else:
        if not (len(ld) == 1 and len(rd) == 1):
            return "Not allowed contraction"
        index_str = " ".join(str(-i) for i in ld + rd)

    return right * left * eps(index_str)


def contract(
    fields: Tuple[IndexedField],
    symbols: Dict[str, List[str]],
    gauge_epsilons: list,
    field_dict: Dict[tuple, str],
    derivative_state=None,
) -> Union[Tuple[FieldType, Operator, List[Tensor], List[Tensor]], str]:
    """Takes two or three indexed fields and the epsilons [epsilons and deltas of
    SU(2) and SU(3) from the operator] and returns a new indexed field
    transforming in the same way as $x \\otimes y$.

    Gauge epsilons (and deltas) are going to be potentially used up in this
    process, while epsilons carrying Lorentz indices will be introduced
    enforcing the contractions between dotted and undotted indices in the
    generated operator.

    Returns a tuple with the field transforming like the product of `fields`,
    the term, and the new gauge and lorentz epsilons.

    If the contraction fails, returns a string with the reason it failed.

    Example:
        >>> field, term, gauge_epsilons, lorentz_epsilons = contract((H('i0'), H('i1')), [], {"fermion": [], "boson": ["S"]}, {})
        >>> field
        S(i0, i1)
        >>> field.y
        1

    """
    if len(fields) != 2 and len(fields) != 3:
        raise Exception("Too many fields passed to contract.")

    prepared = prepare_lorentz_contraction(fields, derivative_state)
    if prepared is None:
        # Bad lorentz contraction
        return "Bad Lorentz contraction."

    fields, lorentz_epsilons, derivative_route = prepared

    # some gauge epsilons will be removed in the contraction, the others will
    # just watch
    spectator_gauge_eps, eps_to_remove = separate_gauge_epsilons(fields, gauge_epsilons)
    # contracted_fields is the fields (with the derivatives still present) with
    # lorentz indices contracted
    contracted_fields = reduce(
        lambda x, y: x * y, (*fields, *lorentz_epsilons, *eps_to_remove)
    )

    exotic, partner, maybe_term = exotic_field_and_term(
        contracted_fields, symbols, field_dict
    )

    if isinstance(maybe_term, str):
        # return the reason
        return maybe_term

    check_singlet(maybe_term)

    # Check to see if there are any derivatives present, if there are process the term
    deriv_structure = [f.derivs for f in maybe_term.fields]
    n_derivs = sum(deriv_structure)

    if n_derivs == 0:
        no_deriv_maybe_term = maybe_term

        if isinstance(no_deriv_maybe_term, str):
            return no_deriv_maybe_term

    else:
        no_deriv_maybe_term = process_derivative_term(maybe_term)

        if isinstance(no_deriv_maybe_term, str):
            return no_deriv_maybe_term

        if no_deriv_maybe_term.safe_simplify() == 0:
            return f"Vanishing coupling at {maybe_term} after derivative processing."

    if not is_singlet(no_deriv_maybe_term):
        return f"Non-singlet coupling at {no_deriv_maybe_term}"
    check_singlet(no_deriv_maybe_term)

    consume_routed_derivative(derivative_state, derivative_route)

    return exotic, no_deriv_maybe_term, spectator_gauge_eps, lorentz_epsilons


def get_connecting_edge(graph: nx.Graph, nodes: List[int]) -> Tuple[int, int]:
    """Returns an edge that connects to nodes in ``nodes``.

    Taking ``graph`` to be:

            4
            |
            |
            3
           / \
          /   \
         1     2

    Example:
        >>> get_connecting_edge(graph, (1, 2))
        (3, 4)

    """
    neighbours = {}
    for node in nodes:
        neighbours[node] = set(graph.neighbors(node))

    fst, *rst = list(neighbours.values())
    intersection = fst.intersection(*rst)
    assert len(intersection) == 1
    connecting_node = list(intersection)[0]

    other_nodes = list(graph.neighbors(connecting_node))
    for node in nodes:
        other_nodes.remove(node)

    assert len(other_nodes) == 1
    return (connecting_node, other_nodes[0])


def replace_and_mutate(
    leaves: Tuple[Tuple[IndexedField, int]],
    symbols: List[str],
    gauge_epsilons: list,
    lorentz_epsilons: list,
    terms: list,
    edge_dict: Dict[FieldType, Tuple[int, int]],
    field_dict: Dict[tuple, str],
    graph: nx.Graph,
    derivative_state=None,
) -> Leaf:
    """Returns a Leaf structure that enters the partition in place of the contracted
    fields. Mutates major state of completion: terms, edge_dict of graph,
    gauge_epsilons and lorentz_epsilons. Mutation of the field_dict happens in
    `exotic_field_and_term` through `contract`.

    For a failed completion, keep reason in first element of leaf-tuple.

    """
    fields, nodes = [], []
    for leaf in leaves:
        if leaf[1] == None:
            return leaf

        field, node = leaf
        fields.append(field)
        nodes.append(node)

    # if only one field, SM field at last vertex
    if len(fields) == 1:
        return Leaf(field, node)

    # field_dict is updated in this call
    maybe_contract = contract(
        fields,
        symbols,
        gauge_epsilons,
        field_dict,
        derivative_state=derivative_state,
    )

    # For a failed completion, keep reason in first element of leaf-tuple.
    if isinstance(maybe_contract, str):
        return Leaf(maybe_contract, None)

    # mutate gauge_epsilons immediately
    exotic_field, term, gauge_epsilons, new_lorentz_epsilons = maybe_contract

    # mutate lorentz_epsilons, terms
    lorentz_epsilons += new_lorentz_epsilons

    check_singlet(term)
    terms.append(term)

    # update edge_dict
    exotic_edge = get_connecting_edge(graph, nodes)
    edge_dict[exotic_field] = exotic_edge

    if derivative_state is not None and derivative_state["pending_route"] is not None:
        record_routed_derivative_edges(derivative_state, edge_dict)

    return Leaf(exotic_field, exotic_edge[0])


def contains_only_leaves(xs: tuple) -> bool:
    if not isinstance(xs, tuple):
        return False

    for x in xs:
        if not isinstance(x, Leaf):
            return False

    return True


def reduced_row(row, func):
    """Helper function to apply recursive call until you reach leaves."""
    if isinstance(row, Leaf):
        return row

    # row is a tuple
    if contains_only_leaves(row):
        return func(row)

    return func(tuple(map(lambda a: reduced_row(a, func), row)))


def partition_leaves(partition):
    """Return the external leaves in a recursive partition."""

    if isinstance(partition, Leaf):
        return (partition,)
    return tuple(
        leaf for branch in partition for leaf in partition_leaves(branch)
    )


def canonical_rooted_partitions(partition, graph):
    """Reroot a recursive tree at every deterministic internal vertex."""

    leaves_by_node = {leaf.node: leaf for leaf in partition_leaves(partition)}

    def rooted_branch(node, parent):
        if node in leaves_by_node:
            return leaves_by_node[node]
        children = sorted(
            neighbour for neighbour in graph.neighbors(node) if neighbour != parent
        )
        return tuple(rooted_branch(child, node) for child in children)

    roots = sorted(node for node in graph if node not in leaves_by_node)
    if not roots:
        raise ValueError("A completion tree must have an internal vertex")
    return tuple(rooted_branch(root, None) for root in roots)


def canonical_rooted_partition(partition, graph):
    """Return the first canonical rooting for single-result compatibility."""

    return canonical_rooted_partitions(partition, graph)[0]


def derivative_route_choice_count(rooted_partition):
    """Bound the possible heavy-fermion numerator choices in a rooted tree.

    An internal tree edge carries a fermion only when the corresponding subtree
    contains an odd number of external fermions.  At a given vertex, the number
    of such child edges therefore bounds the number of distinct internal
    fermions that can supply a routed momentum numerator.  This structural
    preflight deliberately ignores Lorentz and gauge constraints: it only
    removes branches which cannot possibly consume the derivative.
    """

    def subtree_data(branch):
        if isinstance(branch, Leaf):
            return int(branch.field.is_fermion), 0

        children = [subtree_data(child) for child in branch]
        fermion_count = sum(count for count, _ in children)
        internal_fermions = sum(
            count % 2
            for child, (count, _) in zip(branch, children)
            if not isinstance(child, Leaf)
        )
        child_maximum = max((maximum for _, maximum in children), default=0)
        return fermion_count, max(internal_fermions, child_maximum)

    _, maximum = subtree_data(rooted_partition)
    return min(maximum, MAX_DERIVATIVE_ROUTE_CHOICES)


def construct_completion(
    partition,
    gauge_epsilons,
    graph,
    derivative_count=0,
    derivative_route_choice=0,
    derivative_route_choices=None,
) -> Union[str, tuple]:
    """Returns arguments needed to pass into Completion object contructor, or a
    string with the reason the completion failed.

    """
    lorentz_epsilons, terms, edge_dict, field_dict = [], [], {}, {}
    if derivative_route_choices is None:
        derivative_route_choices = (derivative_route_choice,) * derivative_count
    derivative_route_choices = tuple(derivative_route_choices)
    if len(derivative_route_choices) != derivative_count:
        raise ValueError(
            "The route-choice sequence must match the routed derivative count"
        )
    derivative_state = {
        "remaining": derivative_count,
        "consumed": 0,
        "pending_route": None,
        "routes": [],
        "route_choices": derivative_route_choices,
    }
    more_fermion_symbols = ["f" + str(i) for i in range(10)]
    more_scalar_symbols = ["S" + str(i) for i in range(10)]
    symbols = {
        "fermion": ["ψ", "χ", "f", "ζ", "θ"] + more_fermion_symbols,
        "boson": ["φ", "η", "s", "ω", "σ"] + more_scalar_symbols,
    }

    func = lambda leaves: replace_and_mutate(
        leaves=leaves,
        symbols=symbols,
        gauge_epsilons=gauge_epsilons,
        lorentz_epsilons=lorentz_epsilons,
        terms=terms,
        edge_dict=edge_dict,
        field_dict=field_dict,
        graph=graph,
        derivative_state=derivative_state,
    )

    reduced_partition = [reduced_row(row, func) for row in partition]

    # construct final interaction term and add to terms
    prod = None
    for i in reduced_partition:
        f = i.field
        if isinstance(f, str):
            return f

        if prod is None:
            prod = f
        else:
            prod *= f

    fields = [f for f in prod.tensors if isinstance(f, IndexedField)]

    prepared = prepare_lorentz_contraction(fields, derivative_state)
    if prepared is None:
        return "Bad Lorentz contraction."

    fields, new_lorentz_epsilons, derivative_route = prepared
    prod = reduce(lambda x, y: x * y, fields)

    _, eps_to_remove = separate_gauge_epsilons(fields, gauge_epsilons)

    for e in [*new_lorentz_epsilons, *eps_to_remove]:
        prod *= e

    # mutate Lorentz epsilons with last contraction
    lorentz_epsilons += new_lorentz_epsilons

    # Check to see if there are any derivatives present, if there are process the term
    deriv_structure = [f.derivs for f in prod.fields]
    n_derivs = sum(deriv_structure)

    if n_derivs == 0:
        proc_term = prod

        if isinstance(proc_term, str):
            return proc_term

    else:
        proc_term = process_derivative_term(prod)

        if isinstance(proc_term, str):
            return proc_term

        if proc_term.safe_simplify() == 0:
            return f"Vanishing coupling at {prod} after derivative processing."

    if not is_singlet(proc_term):
        return f"Non-singlet coupling at {proc_term}"

    # make sure the term is a singlet
    check_singlet(proc_term)

    # append the processed term to terms
    terms.append(proc_term)

    consume_routed_derivative(derivative_state, derivative_route)
    record_routed_derivative_edges(derivative_state, edge_dict)
    if derivative_state["remaining"]:
        return "Unresolved derivative insertion."
    if len(derivative_state["routes"]) != derivative_count:
        return "Incorrect number of routed derivative insertions."

    if n_derivs == 0 and is_vanishing_interaction(proc_term):
        return f"Vanishing coupling at {proc_term}"

    return (
        terms,
        edge_dict,
        field_dict,
        lorentz_epsilons,
        tuple(derivative_state["routes"]),
    )


def restore_routed_derivative_operator(
    requested_operator: Operator,
    stripped_operator: Operator,
    lorentz_epsilons: list,
) -> Operator:
    """Restore the requested derivative on the UV-induced stripped contraction.

    The heavy-fermion numerator leaves one undotted and one dotted external
    Lorentz index.  Attaching those ports to the originally differentiated field
    chooses the representative used by the derivative-operator archive.
    """

    explicit = reduce(
        lambda left, right: left * right, lorentz_epsilons, stripped_operator
    )
    derivative = next(
        field for field in requested_operator.indexed_fields if field.derivs
    )
    stripped_field = derivative.strip_derivs()
    target = next(
        field
        for field in explicit.indexed_fields
        if field.field == stripped_field
        and tuple(map(str, field.gauge_indices))
        == tuple(map(str, derivative.gauge_indices))
    )

    free_lorentz = {
        short_type: [
            index
            for index in explicit.free_indices
            if index.index_type == Index.get_index_types()[short_type]
        ]
        for short_type in ("u", "d")
    }
    if any(len(indices) != 1 for indices in free_lorentz.values()):
        raise ValueError(
            "A routed derivative must expose one undotted and one dotted port"
        )

    derivative_indices = []
    restoring_epsilons = []
    for short_type in ("u", "d"):
        free_index = free_lorentz[short_type][0]
        if not free_index.is_up:
            derivative_indices.append(-free_index)
            continue
        fresh = Index.fresh(short_type)
        derivative_indices.append(fresh)
        restoring_epsilons.append(eps(f"{-free_index} {-fresh}"))

    restored = derivative.field(
        " ".join(map(str, (*derivative_indices, *target.gauge_indices)))
    )
    return Operator(
        *(restored if tensor is target else tensor for tensor in explicit.tensors),
        *restoring_epsilons,
    )


def partition_completion(partition) -> Union[Completion, FailedCompletion]:
    """Return the completion object associated with a partition."""
    part = partition["partition"]
    gauge_epsilons = partition["epsilons"]
    graph = partition["graph"]
    derivative_count = partition.get("derivative_count", 0)
    if derivative_count and not partition.get("derivative_partition_is_canonical"):
        part = canonical_rooted_partition(part, graph)
    op = partition["operator"]
    topo = partition["topology"]
    canonical_topo = partition.get("canonical_topology", topo)

    # if args is a string, then it's the reason the completion failed
    args = construct_completion(
        part,
        gauge_epsilons,
        graph,
        derivative_count=derivative_count,
        derivative_route_choice=partition.get("derivative_route_choice", 0),
        derivative_route_choices=partition.get("derivative_route_choices"),
    )
    if not isinstance(args, str):
        terms, edge_dict, field_dict, lorentz_epsilons, derivative_routes = args
    else:
        return FailedCompletion(args)

    completion_operator = partition.get("completion_operator")
    historical_basis = partition.get("historical_derivative_basis")
    lorentz_projection = None
    if historical_basis is not None:
        requested_operator = partition["requested_operator"]
        explicit = reduce(
            lambda left, right: left * right, lorentz_epsilons, op.operator
        )
        explicit_op = restore_routed_derivative_operator(
            requested_operator.operator, op.operator, lorentz_epsilons
        )
        route = derivative_routes[0]
        lorentz_projection = historical_basis.project_routed(
            explicit, part, graph, route
        )
        if lorentz_projection is None:
            return FailedCompletion(
                "Zero projection in the historical derivative-placement basis."
            )
    elif completion_operator is None:
        explicit_op = reduce(lambda a, b: a * b, lorentz_epsilons, op.operator)
    else:
        explicit_op = completion_operator.operator
    exotics = set(f for f in edge_dict.keys())
    eff_operator = EffectiveOperator(op.name, explicit_op)

    new_edge_attrs = {v: {"particle": k.label} for k, v in edge_dict.items()}
    nx.set_edge_attributes(graph, new_edge_attrs)

    return Completion(
        operator=eff_operator,
        partition=part,
        graph=graph,
        exotics=exotics,
        terms=terms,
        topology=topo,
        canonical_topology=canonical_topo,
        derivative_routes=derivative_routes,
        lorentz_projection=lorentz_projection,
    )


def operator_completions(
    operator: EffectiveOperator, verbose=False, canonical_partitions=False
) -> List[Completion]:
    """Return a list of the completions of an effective operator."""

    parts = partitions(operator, verbose=verbose)
    if canonical_partitions:
        parts = remove_isomorphic(parts)
    if verbose:
        mode = "canonical" if canonical_partitions else "raw"
        print(f"Starting with {len(parts)} {mode} partitions...")

    if verbose:
        print(f"Finding completions of {len(parts)} partitions...")
        with alive_bar(len(parts)) as bar:
            for p in parts:
                # completions.append(partition_completion(p))
                comp = partition_completion(p)
                if not isinstance(comp, FailedCompletion):
                    yield comp
                bar()
    else:
        for p in parts:
            comp = partition_completion(p)
            if not isinstance(comp, FailedCompletion):
                yield comp

        # completions = [partition_completion(p) for p in parts]

    # good_completions = [c for c in completions if not isinstance(c, FailedCompletion)]
    # return good_completions


def check_remapping_on_terms(terms1, terms2, remapping):
    """Return the remapping on the field labels in the terms that would get you from
    one to the other, i.e. return the isomorphism if one exists, otherwise
    return the empty dictionary.

    """
    if equivalent_lagrangians(terms1, terms2, remapping):
        return remapping

    # otherwise, no equivalence, return empty dict
    return {}


def base_exotic_label(label: str) -> str:
    """Return the particle-species label without conjugate/Dirac suffixes."""

    return field_label_parts(label)[0]


def exotic_species_kind(field: FieldType) -> str:
    """Return the mass/particle nature omitted from the gauge quantum numbers."""

    if isinstance(field, VectorLikeDiracFermion):
        return "dirac_fermion"
    if isinstance(field, MajoranaFermion):
        return "majorana_fermion"
    if isinstance(field, RealScalar):
        return "real_scalar"
    if isinstance(field, ComplexScalar):
        return "complex_scalar"
    raise TypeError(f"Unrecognised exotic species {type(field).__name__}")


def exotic_species(completion: Completion) -> Dict[str, tuple]:
    """Map each distinct particle species to its physical descriptor.

    A base label denotes one physical species. Repeated occurrences of that label
    are interaction-edge occurrences, whereas different labels remain distinct even
    when their representations coincide. Democratic filtering deliberately applies
    the separate, coarser convention of collapsing equal representations.
    """

    species = {}
    for field, quantum_numbers in completion.exotic_info().items():
        label = base_exotic_label(field.label)
        descriptor = (exotic_species_kind(field), quantum_numbers)
        known_descriptor = species.get(label)
        if known_descriptor is not None and known_descriptor != descriptor:
            raise ValueError(f"Inconsistent quantum numbers for exotic species {label}")
        species[label] = descriptor
    return species


def exotic_field_occurrences(completion: Completion):
    """Return field-factor signatures grouped by particle-species label."""

    occurrences = defaultdict(Counter)
    species_labels = exotic_species(completion)
    for term in completion.terms:
        for field in term.indexed_fields:
            label, is_conjugate, is_dirac_partner = field_label_parts(field.label)
            if label not in species_labels:
                continue
            signature = (
                is_conjugate,
                is_dirac_partner,
                field.dynkin,
                tuple(
                    sorted(
                        (name, str(value))
                        for name, value in field.charges.items()
                    )
                ),
                field.comm,
                field.derivs,
            )
            occurrences[label][signature] += 1
    return occurrences


def oriented_exotic_mappings(comp1, comp2, remapping):
    """Yield suffix orientations compatible with the physical field factors."""

    occurrences1 = exotic_field_occurrences(comp1)
    occurrences2 = exotic_field_occurrences(comp2)
    species1 = exotic_species(comp1)
    orientation_groups = []
    source_labels = sorted(remapping)
    for source_label in source_labels:
        target_label = remapping[source_label]
        is_dirac = species1[source_label][0] == "dirac_fermion"
        options = [(target_label, False, False)]
        for conjugate_flip in (False, True):
            for dirac_flip in ((False, True) if is_dirac else (False,)):
                remapped_occurrences = Counter()
                for signature, multiplicity in occurrences1[source_label].items():
                    is_conjugate, is_dirac_partner, *field_data = signature
                    remapped_signature = (
                        is_conjugate ^ conjugate_flip,
                        is_dirac_partner ^ dirac_flip,
                        *field_data,
                    )
                    remapped_occurrences[remapped_signature] += multiplicity
                if remapped_occurrences == occurrences2[target_label]:
                    option = (target_label, conjugate_flip, dirac_flip)
                    if option not in options:
                        options.append(option)
        orientation_groups.append(options)

    for orientations in product(*orientation_groups):
        yield dict(zip(source_labels, orientations))


def exotic_label_bijections(comp1: Completion, comp2: Completion):
    """Yield all representation-preserving species relabellings."""

    species1 = exotic_species(comp1)
    species2 = exotic_species(comp2)
    quantum_numbers1 = Counter(species1.values())
    quantum_numbers2 = Counter(species2.values())
    if quantum_numbers1 != quantum_numbers2:
        return

    labels1_by_quantum_numbers = defaultdict(list)
    labels2_by_quantum_numbers = defaultdict(list)
    for label, quantum_numbers in species1.items():
        labels1_by_quantum_numbers[quantum_numbers].append(label)
    for label, quantum_numbers in species2.items():
        labels2_by_quantum_numbers[quantum_numbers].append(label)

    quantum_number_classes = sorted(labels1_by_quantum_numbers, key=repr)
    permutation_groups = []
    source_groups = []
    for quantum_numbers in quantum_number_classes:
        source_labels = sorted(labels1_by_quantum_numbers[quantum_numbers])
        target_labels = sorted(labels2_by_quantum_numbers[quantum_numbers])
        source_groups.append(source_labels)
        permutation_groups.append(tuple(permutations(target_labels)))

    for target_groups in product(*permutation_groups):
        mapping = {}
        for source_labels, target_labels in zip(source_groups, target_groups):
            mapping.update(zip(source_labels, target_labels))
        yield mapping


def momentum_contributions_equivalent(comp1, comp2, oriented_remapping):
    """Compare propagator-expansion terms under an exotic-field relabelling."""

    contributions1 = getattr(comp1, "momentum_contributions", ())
    contributions2 = getattr(comp2, "momentum_contributions", ())
    if not contributions1 or not contributions2:
        return True
    if len(contributions1) != len(contributions2):
        return False

    legacy_numerators = all(
        contribution.numerator_kind == "momentum"
        and contribution.denominator_order == 0
        for contribution in (*contributions1, *contributions2)
    )

    def contribution_signature(contribution, label_mapping):
        mapped_label = remapped_field_label(
            contribution.particle, label_mapping
        )
        if legacy_numerators:
            return field_label_parts(mapped_label)[0]
        return (
            field_label_parts(mapped_label)[0],
            contribution.numerator_kind,
            contribution.denominator_order,
            tuple(sorted(contribution.cut_side)),
        )

    signatures1 = Counter(
        contribution_signature(contribution, oriented_remapping)
        for contribution in contributions1
    )
    signatures2 = Counter(
        contribution_signature(contribution, {})
        for contribution in contributions2
    )
    return signatures1 == signatures2


def derivative_routes_equivalent(comp1, comp2, oriented_remapping):
    """Backwards-compatible alias for propagator-contribution equivalence."""

    return momentum_contributions_equivalent(
        comp1, comp2, oriented_remapping
    )


def lorentz_projections_equivalent(comp1, comp2):
    """Keep distinct effective Lorentz components as distinct exact classes."""

    projection1 = getattr(comp1, "lorentz_projection", None)
    projection2 = getattr(comp2, "lorentz_projection", None)
    if projection1 is None or projection2 is None:
        return projection1 is projection2
    return (
        projection1.basis_labels == projection2.basis_labels
        and projection1.coordinates == projection2.coordinates
        and projection1.derivative_field == projection2.derivative_field
    )


def compare_terms(comp1: Completion, comp2: Completion) -> Dict[str, str]:
    """Returns a dictionary representing the field relabellings that would need to
    be applied to the terms of comp1 to make it equivalent to comp2. This
    includes the identity remapping. That is, if the terms of comp1 are the same
    as the terms in comp2, the function returns a dictionary like

       {"φ": "φ", "η": "η", ...}

    """
    # cannot be equivalent
    if len(comp1.terms) != len(comp2.terms):
        return {}

    for remapping in exotic_label_bijections(comp1, comp2):
        for oriented_remapping in oriented_exotic_mappings(comp1, comp2, remapping):
            if equivalent_lagrangians(
                comp1.terms, comp2.terms, oriented_remapping
            ) and momentum_contributions_equivalent(
                comp1, comp2, oriented_remapping
            ) and lorentz_projections_equivalent(comp1, comp2):
                return remapping

    return {}


def are_equivalent_completions(comp1: Completion, comp2: Completion) -> bool:
    """Checks to see if the Lagrangian terms describing two completions are
    equivalent.

    Two completions are equivalent if their exact contraction graphs are the same
    up to representation-preserving field relabellings.

    """
    return bool(compare_terms(comp1, comp2))


def slow_remove_equivalent_completions(
    comps: List[Completion], verbose: bool = False
) -> List[Completion]:
    """Compares completions by comparing Lagrangian terms. Removes duplicates and
    returns copied list.

    """
    remove_equivalent(comps, are_equivalent_completions)


def collect_completions(
    completions: Iterable[Completion], key=None
) -> Dict[tuple, List[Completion]]:
    """Return dictionary mapping field content to list of completions.

    `key` is a function that takes a completion and returns a dictionary mapping
    FieldType to a tuple of numbers representing that field. This defaults to
    the `exotic_info` method.

    Not for general user interface.

    """
    out = defaultdict(list)

    if key is None:
        key = lambda x: x.exotic_info()

    for completion in completions:
        model_key = tuple(sorted(set(key(completion).values())))
        out[model_key].append(completion)

    return dict(out)


def prime_registry(sieve: Dict[tuple, List[Completion]]) -> Dict[tuple, int]:
    """Ascribe a unique prime number to each exotic appearing in `sieve`.

    `sieve` is a dictionary mapping a tuple of field information to a list of
    completions.

    """
    reg = {}
    counter = 1
    for k, v in sieve.items():
        for field in k:
            if field not in reg:
                reg[field] = prime(counter)
                counter += 1
    return reg


def model_registry(completions, registry) -> Dict[tuple, int]:
    """Assigns an unique integer to every model by multiplying primes of fields."""
    reg = {}
    for k in completions:
        prod = 1
        for field in k:
            prod *= registry[field]

        reg[k] = prod

    return reg


def filter_completions(
    completions: Dict[tuple, List[Completion]], sieve: Dict[tuple, List[Completion]]
) -> Dict[tuple, List[Completion]]:
    # establish prime registry
    registry = prime_registry({**sieve, **completions})

    # construct dictionaries mapping tuples of field info to integers (products
    # of primes)
    completions_model_registry = model_registry(completions, registry)
    sieve_model_registry = model_registry(sieve, registry)

    unique = {}
    for k, v in completions_model_registry.items():
        factors_ = factors(v)
        for ref_val in sieve_model_registry.values():
            if ref_val in factors_:
                break
        else:  # no break => unique model
            unique[k] = completions[k]

    return unique


def operator_strip_derivs(op: Operator) -> List[Operator]:
    """Removes the derivatives from the operator and returns a dictionary of the
    fields, epsilons and number of derivatives. The fields output is an
    association list between Field and list of guage indices (including
    generational indices) for the field in the operator.

    """
    tensors = op.tensors
    new_fields = []
    epsilons = []
    n_derivs = 0
    for field in tensors:
        if isinstance(field, Field):
            if field.derivs:
                n_derivs += field.derivs
                new_field = field.strip_derivs()
                indices = field.gauge_indices
                new_fields.append((new_field, " ".join(str(i) for i in indices)))
                # new_fields.append((new_field, indices))
            else:
                indices = field.gauge_indices
                new_fields.append((field.field, " ".join(str(i) for i in indices)))
                # new_fields.append((field.field, indices))
        else:
            epsilons.append(field)

    return {"fields": new_fields, "epsilons": epsilons, "n_derivs": n_derivs}


def construct_operator(
    fields: List[Tuple[Field, str]], epsilons: List[Tensor]
) -> Operator:
    """Helper function to construct operator."""
    tensors = []
    for field, index_string in fields:
        u, d, _, _, _ = field.fresh_indices().indices_by_type.values()
        lor_idx_str = " ".join(str(i) for i in u + d)
        tensors.append(field(lor_idx_str + " " + index_string))

    return reduce(lambda x, y: x * y, tensors + list(epsilons))


def unique_lorentz_completion_operator(operator: EffectiveOperator):
    """Return the unique nonzero Lorentz singlet associated with ``operator``.

    Multi-dimensional spaces are handled separately by explicit placement
    projectors.
    """

    singlets = {}
    for singlet in lorentz_singlets(operator.operator):
        simple = singlet.safe_simplify()
        if simple != 0:
            singlets[str(safe_nocoeff(simple))] = singlet

    if len(singlets) != 1:
        return None

    return EffectiveOperator(operator.name, next(iter(singlets.values())))


def unique_multi_derivative_projection(operator):
    """Return the normalised coordinate of a unique multi-derivative singlet."""

    _, _, n_derivs = operator_strip_derivs(operator.operator).values()
    if n_derivs < 2:
        return None
    basis = LorentzBasis.from_operator(operator.operator)
    if basis.dimension != 1:
        return None
    derivative_fields = tuple(
        field.label
        for field in operator.operator.indexed_fields
        if field.derivs
    )
    return MultiDerivativeProjection(
        basis_labels=basis.labels,
        coordinates=("1",),
        derivative_fields=derivative_fields,
        ibp_relation=UNIQUE_MULTI_DERIVATIVE_IBP_RELATION,
        eom_relation=HISTORICAL_EOM_RELATION,
    )


def canonical_propagator_cut(partition, graph, edge):
    """Return a deterministic external-node side of a cut internal edge."""

    cut_graph = graph.copy()
    cut_graph.remove_edge(*edge)
    external_nodes = {leaf.node for leaf in partition_leaves(partition)}
    sides = []
    for component in nx.connected_components(cut_graph):
        side = tuple(sorted(component & external_nodes))
        if side:
            sides.append(side)
    if len(sides) != 2:
        raise ValueError("A propagator edge must split the external fields in two")
    return min(sides)


def weak_compositions(total, length):
    """Yield deterministic weak compositions of ``total`` into ``length`` parts."""

    if length == 0:
        if total == 0:
            yield ()
        return
    for first in range(total + 1):
        for rest in weak_compositions(total - first, length - 1):
            yield (first, *rest)


def _completion_with_momentum_contributions(completion, contributions):
    """Copy a completion while replacing its propagator-expansion metadata."""

    return Completion(
        operator=completion.operator,
        partition=completion.partition,
        graph=deepcopy(completion.graph),
        exotics=completion.exotics,
        terms=completion.terms,
        topology=completion.topology,
        canonical_topology=completion.canonical_topology,
        momentum_contributions=tuple(contributions),
        lorentz_projection=completion.lorentz_projection,
    )


def _edge_particle_field(completion, edge):
    particle = completion.graph.edges[edge]["particle"]
    exact_matches = [
        field for field in completion.exotics if field.label == particle
    ]
    matches = exact_matches or [
        field
        for field in completion.exotics
        if base_exotic_label(field.label) == base_exotic_label(particle)
    ]
    if not matches or len({field.is_fermion for field in matches}) != 1:
        labels = tuple(sorted(field.label for field in completion.exotics))
        raise ValueError(
            f"Cannot identify particle {particle!r} on propagator edge "
            f"{edge}; completion exotics are {labels}"
        )
    return particle, matches[0]


def expand_propagator_denominators(completion, denominator_order):
    """Attach every propagator-denominator term of a fixed total order.

    ``denominator_order`` counts powers of :math:`p^2/M^2`, so it supplies
    twice as many EFT derivatives.  Existing fermion momentum numerators are
    retained and may carry additional denominator powers on the same edge.
    """

    if denominator_order < 0:
        raise ValueError("Propagator denominator order must be nonnegative")
    exotic_labels = {
        base_exotic_label(field.label) for field in completion.exotics
    }
    internal_edges = tuple(
        sorted(
            tuple(sorted(edge))
            for edge, particle in nx.get_edge_attributes(
                completion.graph, "particle"
            ).items()
            if base_exotic_label(particle) in exotic_labels
        )
    )
    if not internal_edges:
        return [] if denominator_order else [completion]

    canonical_contributions = [
        contribution._replace(
            edge=tuple(sorted(contribution.edge)),
            cut_side=canonical_propagator_cut(
                completion.partition,
                completion.graph,
                contribution.edge,
            ),
        )
        for contribution in completion.momentum_contributions
    ]
    expanded = []
    for allocation in weak_compositions(denominator_order, len(internal_edges)):
        contributions = list(canonical_contributions)
        by_edge = {
            contribution.edge: index
            for index, contribution in enumerate(contributions)
        }
        for edge, order in zip(internal_edges, allocation):
            if not order:
                continue
            if edge in by_edge:
                index = by_edge[edge]
                contribution = contributions[index]
                contributions[index] = contribution._replace(
                    denominator_order=contribution.denominator_order + order
                )
                continue

            particle, field = _edge_particle_field(completion, edge)
            contributions.append(
                PropagatorContribution(
                    edge=edge,
                    particle=particle,
                    numerator_kind="mass" if field.is_fermion else "scalar",
                    denominator_order=order,
                    numerator_lorentz="",
                    differentiated_field="",
                    differentiated_lorentz="",
                    cut_side=canonical_propagator_cut(
                        completion.partition, completion.graph, edge
                    ),
                )
            )
        contributions.sort(
            key=lambda item: (
                item.edge,
                item.numerator_kind,
                item.denominator_order,
                item.particle,
            )
        )
        expanded.append(
            _completion_with_momentum_contributions(completion, contributions)
        )
    return expanded


def momentum_routed_completion_stream(
    operator: EffectiveOperator,
    verbose=False,
    canonical_partitions=False,
    historical_basis=None,
    second_derivative_basis=None,
) -> Iterable[Completion]:
    """Yield completions using the low-momentum propagator expansion.

    External derivative labels are removed while furnishing the graph.
    Fermion momentum numerators supply odd derivative degree and denominator
    corrections supply even degree.  Multi-derivative routing is currently
    enabled either when the requested Lorentz structure is unique or when an
    explicit unreduced placement projector is available.  No EOM reduction is
    applied.
    """

    fields, epsilons, n_derivs = operator_strip_derivs(operator.operator).values()
    if not n_derivs:
        return

    completion_operator = unique_lorentz_completion_operator(operator)
    unique_projection = unique_multi_derivative_projection(operator)
    if completion_operator is None:
        if n_derivs == 2:
            if second_derivative_basis is None:
                second_derivative_basis = (
                    UnreducedSecondDerivativeBasis.from_operator(operator)
                )
        elif n_derivs != 1 or operator.name not in PROJECTED_LORENTZ_OPERATORS:
            return
        elif historical_basis is None:
            historical_basis = HistoricalDerivativeBasis.from_operator(operator)

    stripped_operator = EffectiveOperator(
        operator.name, construct_operator(fields, epsilons)
    )
    routed_partitions = partitions(stripped_operator, verbose=verbose)
    if canonical_partitions:
        routed_partitions = remove_isomorphic(routed_partitions)
    for partition in routed_partitions:
        partition_candidates = []
        for numerator_count in range(n_derivs % 2, n_derivs + 1, 2):
            denominator_order = (n_derivs - numerator_count) // 2
            if numerator_count:
                rooted_partitions = canonical_rooted_partitions(
                    partition["partition"], partition["graph"]
                )
            else:
                rooted_partitions = (partition["partition"],)

            for rooted_partition in rooted_partitions:
                if numerator_count:
                    route_choice_count = derivative_route_choice_count(
                        rooted_partition
                    )
                    route_choice_sequences = product(
                        range(route_choice_count), repeat=numerator_count
                    )
                else:
                    route_choice_sequences = ((),)

                for route_choices in route_choice_sequences:
                    branch = dict(partition)
                    branch["partition"] = rooted_partition
                    branch["graph"] = deepcopy(partition["graph"])
                    branch["completion_operator"] = completion_operator
                    if historical_basis is not None:
                        branch["historical_derivative_basis"] = historical_basis
                        branch["requested_operator"] = operator
                    branch["derivative_count"] = numerator_count
                    branch["derivative_route_choices"] = route_choices
                    branch["derivative_partition_is_canonical"] = True
                    completion = partition_completion(branch)
                    if isinstance(completion, FailedCompletion):
                        continue
                    if unique_projection is not None:
                        completion.lorentz_projection = unique_projection
                    expanded = expand_propagator_denominators(
                        completion, denominator_order
                    )
                    if second_derivative_basis is not None:
                        if numerator_count:
                            raise ValueError(
                                "An unreduced two-numerator branch requires "
                                "an explicit spinor-momentum projection"
                            )
                        projected = []
                        for candidate in expanded:
                            contribution = next(
                                contribution
                                for contribution in candidate.momentum_contributions
                                if contribution.denominator_order
                            )
                            projection = second_derivative_basis.project_denominator(
                                candidate.operator.operator,
                                candidate.partition,
                                contribution,
                            )
                            if projection is None:
                                continue
                            candidate.operator = operator
                            candidate.lorentz_projection = projection
                            projected.append(candidate)
                        expanded = projected
                    partition_candidates.extend(expanded)
        partition_branches = []
        append_unique_completions(partition_branches, partition_candidates)
        yield from partition_branches


def momentum_routed_completions(
    operator: EffectiveOperator,
    verbose=False,
    canonical_partitions=False,
    historical_basis=None,
    second_derivative_basis=None,
) -> List[Completion]:
    """Return the routed completion stream as a compatibility list."""

    return list(
        momentum_routed_completion_stream(
            operator,
            verbose=verbose,
            canonical_partitions=canonical_partitions,
            historical_basis=historical_basis,
            second_derivative_basis=second_derivative_basis,
        )
    )


def exact_completion_bucket_key(completion):
    """Return a necessary physical-equivalence key for exact deduplication.

    This deliberately retains particle kind and species multiplicity.  The
    democratic filtering key is coarser and is not safe for exact classes.
    """

    species = Counter(exotic_species(completion).values())
    projection = getattr(completion, "lorentz_projection", None)
    projection_key = None
    if projection is not None:
        projection_key = (
            projection.basis_labels,
            projection.coordinates,
            projection.derivative_field,
        )
    return (
        completion.operator.name,
        len(completion.terms),
        tuple(sorted(species.items(), key=repr)),
        projection_key,
    )


def append_unique_completions(completions, candidates):
    """Append candidates whose Lagrangians are not already represented."""

    by_model = defaultdict(list)

    for completion in completions:
        by_model[exact_completion_bucket_key(completion)].append(completion)

    for candidate in candidates:
        key = exact_completion_bucket_key(candidate)
        if any(
            are_equivalent_completions(candidate, known)
            for known in by_model[key]
        ):
            continue

        completions.append(candidate)
        by_model[key].append(candidate)


@dataclass(frozen=True)
class DerivativePlacementSpec:
    index: int
    field_label: str
    occurrence_key: tuple
    operator: EffectiveOperator


def derivative_placement_combinations(op: EffectiveOperator):
    """Return nonzero historical single-derivative placements with labels."""

    fields, epsilons, n_derivs = operator_strip_derivs(op.operator).values()
    if n_derivs != 1:
        raise ValueError("Historical placement basis requires one derivative")

    placements = []
    for placement_index, (target_field, target_indices) in enumerate(fields):
        structure = []
        for field_index, (field, indices) in enumerate(fields):
            placed = (
                D(field, allowed_lor_dyn(field))
                if field_index == placement_index
                else field
            )
            structure.append((placed, indices))
        new_operator = construct_operator(structure, epsilons)
        if not new_operator.safe_simplify():
            continue

        undotted, dotted, _, _, _ = (
            target_field.fresh_indices().indices_by_type.values()
        )
        lorentz_indices = " ".join(map(str, (*undotted, *dotted)))
        original = target_field(
            " ".join(part for part in (lorentz_indices, target_indices) if part)
        )
        placements.append(
            DerivativePlacementSpec(
                index=placement_index,
                field_label=target_field.label
                + ("†" if target_field.is_conj else ""),
                occurrence_key=lorentz_field_key(original),
                operator=EffectiveOperator(op.name, new_operator),
            )
        )
    return placements


def derivative_combinations(
    op: Union[Operator, EffectiveOperator]
) -> Union[List[Operator], List[EffectiveOperator]]:
    """Takes an operator with derivatives and returns a list of operators with
    equivalent SU2 structure with the derivative acted in all possible ways.

    Function expects a specific kind of input: no double derivatives on a single
    field.

    """
    eff_op = None
    if isinstance(op, EffectiveOperator):
        eff_op = op
        op = op.operator

    fields, epsilons, n_derivs = operator_strip_derivs(op).values()
    if n_derivs == 1:
        temporary = eff_op or EffectiveOperator("__derivative_placements__", op)
        placements = derivative_placement_combinations(temporary)
        if eff_op is not None:
            return [placement.operator for placement in placements]
        return [placement.operator.operator for placement in placements]

    deriv_id_func = lambda x: x
    act_deriv = lambda field: D(field, allowed_lor_dyn(field))
    deriv_tuple = [act_deriv for _ in range(n_derivs)] + [
        deriv_id_func for _ in range(len(fields) - n_derivs)
    ]

    structs = []
    for perm in permutations(deriv_tuple):
        structs.append(
            [
                (derivative_action(field), indices)
                for (field, indices), derivative_action in zip(fields, perm)
            ]
        )
    remove_equivalent(structs, eq_func=lambda left, right: left == right)

    out = []
    for struct in structs:
        new_operator = construct_operator(struct, epsilons)
        if new_operator.safe_simplify():
            out.append(new_operator)
    if eff_op is not None:
        return [EffectiveOperator(eff_op.name, item) for item in out]
    return out


def deriv_operator_completion_stream(
    operator: EffectiveOperator, verbose=False, canonical_partitions=False
) -> Iterable[Completion]:
    """Yield completions of every supported derivative placement.

    Unlike :func:`deriv_operator_completions`, this function does not retain
    generated local completions or compare routed candidates against them.
    Downstream artifact generation can therefore write each candidate as soon
    as it is furnished and leave exact uniqueness to the disk-backed
    deduplicator.  Two-derivative multidimensional operators retain both
    distributed derivatives and explicit :math:`D^2` placements without EOM
    reduction.

    """
    historical_basis = None
    unique_projection = unique_multi_derivative_projection(operator)
    second_derivative_basis = None
    _, _, derivative_count = operator_strip_derivs(operator.operator).values()
    if (
        operator.name in PROJECTED_LORENTZ_OPERATORS
        and unique_lorentz_completion_operator(operator) is None
    ):
        historical_basis = HistoricalDerivativeBasis.from_operator(operator)
        placements = historical_basis.placements
        deriv_combos = [placement.operator for placement in placements]
    elif unique_projection is None and derivative_count == 2:
        second_derivative_basis = UnreducedSecondDerivativeBasis.from_operator(
            operator
        )
        placements = second_derivative_basis.placements
        deriv_combos = [placement.operator for placement in placements]
    else:
        placements = None
        deriv_combos = derivative_combinations(operator)

    if verbose:
        print(f"Finding completions of {len(deriv_combos)} IBP-related operators...")

    for combo_number, combo in enumerate(deriv_combos):
        if combo.operator.simplify() == 0:
            continue
        for completion in operator_completions(
            combo,
            verbose=verbose,
            canonical_partitions=canonical_partitions,
        ):
            if historical_basis is not None:
                placement = placements[combo_number]
                completion.lorentz_projection = historical_basis.project_local(
                    placement, completion.operator.operator
                )
            elif unique_projection is not None:
                completion.lorentz_projection = unique_projection
            elif second_derivative_basis is not None:
                completion.lorentz_projection = (
                    second_derivative_basis.project_local(
                        completion.operator.operator
                    )
                )
            yield completion

    yield from momentum_routed_completion_stream(
        operator,
        verbose=verbose,
        canonical_partitions=canonical_partitions,
        historical_basis=historical_basis,
        second_derivative_basis=second_derivative_basis,
    )


def deriv_operator_completions(
    operator: EffectiveOperator, verbose=False, canonical_partitions=False
) -> List[Completion]:
    """Return the derivative completion stream as a compatibility list.

    Routed candidates retain the historical API invariant of being exactly
    distinct from the generated local candidates and from earlier routes.
    """

    local = []
    routed = []
    for completion in deriv_operator_completion_stream(
        operator,
        verbose=verbose,
        canonical_partitions=canonical_partitions,
    ):
        if completion.momentum_contributions:
            routed.append(completion)
        else:
            local.append(completion)
    append_unique_completions(local, routed)

    return local


def exact_completions(operator: EffectiveOperator, verbose=False) -> List[Completion]:
    """Return exact UV-Lagrangian classes after symbolic construction.

    The canonical-partition preflight is safe for ordinary and one-derivative
    operators.  For several derivatives, propagator numerator orientations can
    distinguish partitions only after furnishing, so the complete raw set is
    constructed before applying exact interaction-graph equivalence.
    """

    if any(field.derivs for field in operator.fields):
        derivative_count = operator_strip_derivs(operator.operator)["n_derivs"]
        candidates = deriv_operator_completions(
            operator,
            verbose=verbose,
            canonical_partitions=derivative_count <= 1,
        )
    else:
        candidates = operator_completions(
            operator,
            verbose=verbose,
            canonical_partitions=True,
        )
    exact = []
    append_unique_completions(exact, candidates)
    return exact


def completions(operator: EffectiveOperator, *args, **kwargs):
    """General dispatch function for completions"""
    if any(field.derivs for field in operator.fields):
        return deriv_operator_completions(operator, *args, **kwargs)
    return operator_completions(operator, *args, **kwargs)


def collect_models(comps):
    """Group models by particle content.

    A bit cumbersome to use. Should be refactored out of tests at some point.
    """
    collected = collect_completions(comps)
    return [Model(cs) for _, cs in list(collected.items())]


def cons_term_prime_dict(completions: List[Completion]) -> Dict[tuple, int]:
    # begin by generating term prime dictionary
    term_dict = {}
    counter = 1
    for comp in completions:
        n_terms = len(comp.terms)
        for term in comp.terms:
            # sort all of the terms by side effect
            new_term = tuple(sorted(stringify_qns(f) for f in term.fields))
            if new_term not in term_dict:
                # add conj term first so that when filtering you keep
                # the unconjugated term
                term_dict[conjugate_term(new_term)] = prime(counter)
                term_dict[new_term] = prime(counter)
                counter += 1

    return term_dict


def completion_characteristic_number(
    comp: Completion, prime_dict: Dict[tuple, int]
) -> int:
    prod = 1
    for term in comp.terms:
        new_term = tuple(sorted(stringify_qns(f) for f in term.fields))
        prod *= prime_dict[new_term]
    return prod


def clean_completions(completions: List[Completion]) -> List[Completion]:
    """A fast way of removing equivalent completions using prime label method on terms.

    """
    completions = list(completions)
    prime_dict = cons_term_prime_dict(completions)

    comp_dict = {}
    for comp in completions:
        num = completion_characteristic_number(comp, prime_dict)
        comp_dict[(num, comp.topology)] = comp

    return sorted((v for k, v in comp_dict.items()), key=lambda x: x.topology)
