# coding: utf-8

r"""Memoization."""


def memoize(func):
    """Memoization decorator for functions taking args and kwargs."""
    # cache = func.cache = {}  # causes an issue with Cython
    cache = {}

    # @wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]

    return memoized_func
