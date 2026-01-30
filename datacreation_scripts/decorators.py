#!/usr/bin/env python


import time
from datetime import timedelta
from functools import wraps

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.time()
        total_time = timedelta(minutes=time.perf_counter()-start_time)
        print("\n"+"-"*10)
        print(f"Function {func.__name__} took {total_time}")
        print("-"*10+"\n")
        return result
    return timeit_wrapper

