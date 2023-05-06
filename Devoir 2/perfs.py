import logging
from time import perf_counter
from typing import Any, Callable


def perfs(func: Callable[..., Any]) -> Callable[..., Any]:
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = perf_counter()
        value = func(*args, **kwargs)
        end_time = perf_counter()
        run_time = end_time - start_time
        logging.info(print(f"Execution of {func.__name__}took {run_time: }seconds."))
        #logging.info(f"Inputs are {args,kwargs}")
        return value
    return wrapper

logging.basicConfig(level=logging.INFO)