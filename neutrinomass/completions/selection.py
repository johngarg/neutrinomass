"""Stable generation selection, applied after bare-graph canonicalisation."""

from contextlib import contextmanager
from contextvars import ContextVar

_selection = ContextVar("generation_selection", default=None)


@contextmanager
def generation_selection(selection, boundary=None):
    token = _selection.set((dict(selection), boundary))
    try:
        yield
    finally:
        _selection.reset(token)


def selected_partition_family(topology, colour):
    active = _selection.get()
    if active is None:
        return True
    selection, _ = active
    return ((selection.get("topology") is None or selection["topology"] == topology)
            and (selection.get("colour") is None or selection["colour"] == colour))


def selected_assignment(partition):
    active = _selection.get()
    if active is None:
        return True
    selection, _ = active
    ordinal = partition["assignment"]
    return (ordinal >= selection.get("assignment_start", 0)
            and (selection.get("assignment_stop") is None
                 or ordinal < selection["assignment_stop"]))


def selected_producer(producer, placement=None):
    active = _selection.get()
    if active is None:
        return True
    selection, _ = active
    return ((selection.get("producer") is None or selection["producer"] == producer)
            and (producer != "local" or selection.get("placement") is None
                 or selection["placement"] == placement))


def assignment_boundary(partition):
    active = _selection.get()
    if active is not None and active[1] is not None:
        active[1](partition)
