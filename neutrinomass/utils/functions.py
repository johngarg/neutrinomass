#!/usr/bin/env python3

import time
from typing import List, Callable, TypeVar, Set

T = TypeVar("T")


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def timeit(method):
    """Decorator for timing function execution.

    """

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result

    return timed


def remove_equivalent(
    l: List[T], eq_func: Callable[[T, T], bool], verbose: bool = False
) -> None:
    """Iterates through the list `l` and removes duplicate items according to
    `eq_func`.

    """

    i = 0
    while i < len(l) - 1:

        j = i + 1
        while j <= len(l) - 1:
            if eq_func(l[i], l[j]):
                l.pop(j)
            else:
                j += 1
        i += 1


def remove_equivalent_nopop(
    l: List[T], eq_func: Callable[[T, T], bool], verbose: bool = False
) -> List[T]:
    """Iterates through the list `l` and removes duplicate items according to
    `eq_func`.

    """

    to_remove: Set[int] = set()
    i = 0
    while i < len(l) - 1:

        if i in to_remove:
            i += 1
            continue

        j = i + 1
        while j <= len(l) - 1:
            if j in to_remove:
                j += 1
                continue

            if eq_func(l[i], l[j]):
                to_remove.add(j)
            else:
                j += 1

        i += 1

    return dropped_vals(l, to_remove)


def dropped_vals(l: list, vals: set):
    out = [None] * (len(l) - len(vals))
    i = 0
    for idx, item in enumerate(l):
        if idx in vals:
            continue

        out[i] = item
        i += 1

    return out
