def cache(func):
    internal_cache = None

    def memoized(self):
        nonlocal internal_cache
        if internal_cache is None:
            internal_cache = func(self)

        return internal_cache

    return memoized
