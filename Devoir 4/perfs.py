def perfs(func: Callable[..., Any]) -> Callable[..., Any]:
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = perf_counter()
        value = func(*args, **kwargs)
        end_time = perf_counter()
        run_time = end_time - start_time
        logging.info(prt(f"Execution of {func.__name__}took {run_time: }seconds.",color="BLUE"))
        #logging.info(f"Inputs are {args,kwargs}")
        return value
    return wrapper
logging.basicConfig(level=logging.INFO)
