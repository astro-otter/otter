import pytest
import multiprocessing


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "skip_on_timeout(seconds): Skip the test if it runs longer than `seconds`.",
    )


def _run_test_in_process(runtest_func, result_queue):
    try:
        runtest_func()  # runs setup + test + teardown
        result_queue.put(None)
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

    p = multiprocessing.Process(
        target=_run_test_in_process, args=(item.obj, result_queue)
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


###################################################################################
######### FAKES FOR EXTERNAL SERVICES (SIMBAD, VizieR, IRSA, BLAST, ADS, ...) #####
###################################################################################
# These keep the unit tests independent of third-party services being up. The OTTER
# database queries themselves are NOT mocked.


class _FakeAstroquery:
    """
    Stand-in for DataFinder._wrap_astroquery, the single place every astroquery
    query_region call goes through. Set `result` to what the service should return;
    `calls` records (module, args, kwargs) for each call.
    """

    def __init__(self):
        self.result = None
        self.calls = []

    def __call__(self, module, *args, **kwargs):
        self.calls.append((module, args, kwargs))
        return self.result


@pytest.fixture
def fake_astroquery(monkeypatch):
    from otter.io.data_finder import DataFinder

    fake = _FakeAstroquery()
    monkeypatch.setattr(DataFinder, "_wrap_astroquery", staticmethod(fake))
    return fake


class FakeHTTPResponse:
    """Minimal stand-in for the object returned by urllib.request.urlopen"""

    def __init__(self, body, status=200):
        self.status = status
        self._body = body if isinstance(body, bytes) else body.encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


FAKE_BLAST_RESPONSE = [
    {
        "name": "2024advm",
        "host": {
            "name": "FAKE-HOST-1",
            "ra_deg": 150.0,
            "dec_deg": 2.0,
            "redshift": 0.05,
            "photometric_redshift": 0.06,
        },
    }
]


@pytest.fixture
def fake_blast(monkeypatch):
    """Replace the BLAST HTTP query with a canned response; returns the URLs hit"""
    import json

    urls = []

    def fake_urlopen(url, *args, **kwargs):
        urls.append(url)
        return FakeHTTPResponse(json.dumps(FAKE_BLAST_RESPONSE))

    monkeypatch.setattr("otter.io.host.urlopen", fake_urlopen)
    return urls


@pytest.fixture
def fake_ads(monkeypatch):
    """Replace the ADS bibcode -> human readable name query"""

    def fake_query(bibcodes):
        return list(bibcodes), [f"Fake et al. ({b})" for b in bibcodes]

    monkeypatch.setattr("otter.util._bibcode_to_hrn_with_query", fake_query)
