#!/usr/bin/env python3

from neutrinomass.completions.core import Completion, cons_completion_field
from neutrinomass.completions.topologies import Leaf
from neutrinomass.tensormethod.core import Field, IndexedField, eps, delta, Operator
from neutrinomass.completions.core import FieldType
from neutrinomass.utils.functions import stringify_qns

import networkx as nx


def export_tensor(tensor):
    indices = " ".join(str(i) for i in tensor.indices)
    if isinstance(tensor, FieldType):
        # exotic field
        charges = {k: str(v) for k, v in tensor.charges.items()}
        kwargs = (
            f"label='{tensor.label}'",
            f"indices='{indices}'",
            f"symmetry={tensor.symmetry}",
            f"charges={charges}",
            f"nf={tensor.nf}",
            f"dynkin='{tensor.dynkin}'",
            f"comm='{tensor.comm}'",
            f"latex='{tensor.latex}'",
            f"is_conj={tensor.is_conj}",
        )

        # For VectorLikeDiracFermion, keep track of (un)barred boolean
        is_unbarred = None
        if hasattr(tensor, "is_unbarred"):
            is_unbarred = tensor.is_unbarred
        kwargs += (f"is_unbarred={is_unbarred}",)

        kwargs_str = ", ".join(kwargs)
        indexed_field = f"ExoticField(IndexedField({kwargs_str}))"
        return indexed_field

    if isinstance(tensor, Field):
        # SM field

        if tensor.derivs:
            assert tensor.derivs == 1
            a, b = tensor.lorentz_irrep
            dynkin_label = str(a) + str(b)

            field = tensor.strip_derivs()
            label = field.label + (".conj" if tensor.is_conj else "")
            return f"D({label}, '{dynkin_label}')('{indices}')"

        label = tensor.field.label + (".conj" if tensor.is_conj else "")
        return f"{label}('{indices}')"

    if str(tensor).startswith("metric") or str(tensor).startswith("Eps"):
        return f"eps('{indices}')"

    if str(tensor).startswith("KD"):
        return f"delta('{indices}')"

    else:
        raise ValueError(f"Unrecognised tensor: {tensor}")


def export_operator(op: Operator):
    tensors = op.tensors
    tensor_strings = [export_tensor(t) for t in tensors]
    return "*".join(s for s in tensor_strings)


def proc_dicts(expr):
    if isinstance(expr, (int, str)):
        return expr

    # Used to have field as value, now using field.label, keep this here just in
    # case it gets changed back
    if isinstance(expr, Field):
        return export_tensor(expr)

    return {k: proc_dicts(v) for k, v in expr.items()}


def export_graph(g: nx.Graph):
    dict_of_dicts = nx.to_dict_of_dicts(g)
    return "Graph(" + str(proc_dicts(dict_of_dicts)).replace('"', "") + ")"


def proc_partition(expr):
    if isinstance(expr, Leaf):
        return Leaf(export_tensor(expr.field), expr.node)

    return tuple([proc_partition(i) for i in expr])


def export_partition(partition):
    return "Partition" + str(proc_partition(partition)).replace('"', "")


def export_terms(terms):
    guts = ", ".join(export_operator(op) for op in terms)
    return f"[{guts}]"


def export_exotics(exotics: set):
    return str(set([export_tensor(f) for f in exotics])).replace('"', "")


def export_completion(c: Completion, lazy=True):

    name = c.operator.name
    op = export_operator(c.operator.operator)
    eff_op = f"EffectiveOperator(name='{name}', operator={op})"

    graph = export_graph(c.graph)
    part = export_partition(c.partition)
    terms = export_terms(c.terms)
    exotics = export_exotics(c.exotics)
    topo = c.topology

    quantum_numbers = []
    # lorentz irrep (string), colour dynkins, isospin dynkin, 3 * B, hypercharge
    for l, cu, cd, i, (bl, b), (yl, y) in c.exotic_info().values():
        # quantum_numbers.append((l, (cu, cd), i, (bl, b), (yl, str(y))))
        quantum_numbers.append(f"{l},{cu}{cd},{i},{str(y)},{b}")
    quantum_numbers.sort()

    plain_terms = []
    for t in c.terms:
        plain_term = []
        for f in t.fields:
            plain_term.append(stringify_qns(f))
        plain_terms.append(tuple(plain_term))

    head = (
        "{'operator_name': '%s', 'quantum_numbers': %s, 'terms': %s, 'topology': '%s'}"
        % (name, str(quantum_numbers), str(plain_terms), topo)
    )
    completion_string = f"Completion(operator={eff_op}, partition={part}, graph={graph}, exotics={exotics}, terms={terms}, topology='{topo}')"
    export_string = f"""LazyCompletion(head={head}, tail="{completion_string}")"""

    return export_string if lazy else completion_string
