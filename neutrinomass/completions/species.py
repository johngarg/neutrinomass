"""Finite candidate-specific species assignments before vertex rejection."""

from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy

_assignment = ContextVar("completion_species_assignment",default=None)


class SpeciesChoiceRequired(Exception):
    def __init__(self, choices):
        self.choices = choices


class SpeciesUnavailable(Exception):
    pass


def species_orbit(field):
    from neutrinomass.completions.core import VectorLikeDiracFermion
    states = [field,field.conj]
    if isinstance(field,VectorLikeDiracFermion):
        states += [field.dirac_partner(),field.dirac_partner().conj]
    unique = {}
    for state in states:
        key = (state.label,state.dynkin,tuple(sorted(state.charges.items())),
               state.is_conj,getattr(state,"is_unbarred",None))
        unique[key] = state
    return tuple(unique[key] for key in sorted(unique,key=repr))


def assign_candidate_species(field):
    selector = _assignment.get()
    return field if selector is None else selector(field)


class CandidateAssignment:
    def __init__(self, fields, choices):
        self.states = tuple(state for field in fields for state in species_orbit(field))
        self.choices = choices
        self.cursor = 0

    def __call__(self, generated):
        from neutrinomass.database.serialization import tensor_from_data, tensor_to_data
        from neutrinomass.completions.completions import exotic_species_kind
        allowed = [state for state in self.states
                   if exotic_species_kind(state) == exotic_species_kind(generated)
                   and state.dynkin == generated.dynkin
                   and state.charges == generated.charges]
        if not allowed:
            raise SpeciesUnavailable()
        if self.cursor == len(self.choices):
            raise SpeciesChoiceRequired(len(allowed))
        selected = self.choices[self.cursor]
        self.cursor += 1
        if selected >= len(allowed):
            raise ValueError("Non-reproducible candidate species choice")
        data = tensor_to_data(allowed[selected])
        data["indices"] = [str(index) for index in generated.indices]
        return tensor_from_data(data)


@contextmanager
def candidate_assignment(fields, choices):
    token = _assignment.set(CandidateAssignment(fields,choices))
    try:
        yield
    finally:
        _assignment.reset(token)


def candidate_partition_completions(partition, fields):
    """Replay a small deterministic choice tree; never add an admitted copy."""
    from neutrinomass.completions.completions import partition_completion
    from neutrinomass.completions.core import FailedCompletion
    pending = [()]
    while pending:
        choices = pending.pop()
        try:
            with candidate_assignment(fields,choices):
                # SymPy's generic deepcopy reconstructs IndexedField without
                # its charge/statistics kwargs. External tensors are immutable;
                # only the furnishing graph needs a fresh mutable copy.
                branch = dict(partition)
                branch["graph"] = deepcopy(partition["graph"])
                branch["epsilons"] = list(partition["epsilons"])
                completion = partition_completion(branch)
        except SpeciesChoiceRequired as required:
            pending.extend(choices+(i,) for i in reversed(range(required.choices)))
        except SpeciesUnavailable:
            continue
        else:
            if not isinstance(completion,FailedCompletion):
                yield completion


def candidate_completion_stream(operator, fields, *, topology=None):
    """Use the ordinary complete tree enumeration with fixed admitted species.

    The selector acts within exotic-field construction, before the first
    symmetry/vanishing check. Standard completion entry points are unchanged
    when no candidate context is active.
    """
    from neutrinomass.completions.completions import partition_stream
    for partition in partition_stream(operator):
        if topology is None or partition["canonical_topology"] == topology:
            yield from candidate_partition_completions(partition,fields)
