#!/usr/bin/env python3


def _extend_if_consistent(pattern, data, frame):
    if pattern in frame:
        return pattern_match(frame[pattern], data, frame)
    else:
        return {**frame, **{pattern: data}}


def _match_on_collection(pattern, data, frame):
    try:
        return pattern_match(
            pattern[1:], data[1:], pattern_match(pattern[0], data[0], frame)
        )
    except:
        return False


def is_pattern_var(symbol):
    return isinstance(symbol, str) and symbol[0] == "?"


# TODO write matching tests
# sympy expression, python expression -> dictionary or bool
def pattern_match(pattern, data, frame):
    """A simple pattern matcher for use in constraining completions. Returns a
    dictionary mapping sympy wilds to the values that make `pattern` match
    `data`. If no such values exist the function returns False, while a trivial
    success results in the empty dictionary. Pattern symbols are bound to their
    matches in the `frame` dictionary.

    """
    if frame == False:
        return False
    elif is_pattern_var(pattern):
        return _extend_if_consistent(pattern, data, frame)
    # elif isinstance(pattern, sym.Add) or isinstance(pattern, sym.Mul):
    #     symbol = list(pattern.free_symbols)[0]
    #     value = sym.solve(pattern - data, symbol)[0]
    #     # return _extend_if_consistent(symbol, value, frame)
    #     return pattern_match(symbol, value, frame)
    # elif isinstance(pattern, Field) and isinstance(data, Field):
    #     pattern, data = _field_match_form(pattern), _field_match_form(data)
    #     return pattern_match(pattern, data, frame)
    elif pattern == data:
        return frame
    elif isinstance(pattern, (list, tuple)) and isinstance(data, (list, tuple)):
        return _match_on_collection(pattern, data, frame)
    else:
        return False


# sympy expression, python expression -> bool
def is_match(pattern, data):
    """Returns True if pattern matches data, False otherwise.

    """
    return pattern_match(pattern, data, {}) != False


def pmatch(pattern, data):
    return pattern_match(pattern, data, {})
