#!/usr/bin/env python3

from functools import reduce


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
