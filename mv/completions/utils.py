#!/usr/bin/env python3

import re
from typing import List, Callable, TypeVar
from copy import deepcopy

from functools import reduce

T = TypeVar("T")


def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def foldr(func, xs, acc=None):
    if acc is None:
        return reduce(lambda x, y: func(y, x), xs[::-1])

    return reduce(lambda x, y: func(y, x), xs[::-1], acc)


def fixedpoint(f, x):
    """Applies f o f o f ... f(x) until output is unchanged."""
    app = f(x)
    while application != x:
        return fixedpoint(f, app)
    return x


def factors(n):
    return set(
        reduce(
            list.__add__,
            ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0),
        )
    )


def remove_equivalent(l: List[T], eq_func: Callable[[T, T], bool]) -> List[T]:
    """Iterates through the list `l` and removes duplicate items according to
    `eq_func`. Returns a copy of `l` with duplicates removed.

    """
    list_copy = deepcopy(l)

    i = 0
    while i < len(list_copy) - 1:
        j = i + 1
        while j <= len(list_copy) - 1:
            if eq_func(list_copy[i], list_copy[j]):
                list_copy.pop(j)
            else:
                j += 1
        i += 1

    return list_copy


def multiple_replace(mapping: dict, text: str) -> str:
    regex = re.compile("|".join(map(re.escape, mapping.keys())))
    return regex.sub(lambda mo: mapping[mo.group(0)], text)
