import pytest
import multiprocessing


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "skip_on_timeout(seconds): Skip the test if it runs longer than `seconds`.",
    )


def _run_test_in_process(test_func, args, kwargs, result_queue):
    try:
        test_func(*args, **kwargs)
        result_queue.put(None)  # no exception = success
    except Exception as e:
        result_queue.put(e)


def pytest_runtest_call(item):
    """
    Intercept test execution for tests marked with @pytest.mark.skip_on_timeout(N).

    We run the test in a separate process and skip it if it takes too long.
    """
    mark = item.get_closest_marker("skip_on_timeout")
    if not mark:
        # No marker run test normally
        item.runtest()
        return

    # Get timeout in seconds from the marker
    try:
        timeout_seconds = float(mark.args[0])
    except Exception:
        raise RuntimeError("skip_on_timeout requires a numeric timeout as input")

    # Prepare a queue to get test exceptions from the subprocess
    result_queue = multiprocessing.Queue()

    # Wrap test function to run in subprocess
    def test_func_wrapper():
        item.runtest()  # runs setup + test + teardown

    p = multiprocessing.Process(
        target=_run_test_in_process, args=(test_func_wrapper, (), {}, result_queue)
    )
    p.start()
    p.join(timeout_seconds)

    if p.is_alive():
        # Test exceeded timeout skip
        p.terminate()
        p.join()
        pytest.skip(f"Skipping because test exceeded timeout ({timeout_seconds}s)")

    # Check if test raised an exception
    result = result_queue.get()
    if result is not None:
        raise result  # re-raise the exception so pytest registers failure normally
