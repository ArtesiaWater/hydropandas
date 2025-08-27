# %%
import os

import requests

import hydropandas as hpd


def test_groundwaterstations_vitens():
    lizard_GWS_endpoint = "https://vitens.lizard.net/api/v4/groundwaterstations/"
    code = "27BP0003"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}?code={code}"
    r = requests.get(url_groundwaterstation_code)

    r.raise_for_status()
    assert r.json()["count"] > 0, "request returned empty json"


def test_timeseries_vitens():
    lizard_GWS_endpoint = "https://vitens.lizard.net/api/v4/timeseries/"
    uuid = "44d5b8f9-5680-4fb5-8010-30c85632d106"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}{uuid}"
    r = requests.get(url_groundwaterstation_code)

    r.raise_for_status()


def test_groundwaterstations_rotterdam():
    lizard_GWS_endpoint = "https://rotterdam.lizard.net/api/v4/groundwaterstations/"
    auth = ("__key__", os.environ["LIZARD_ROTTERDAM_API_KEY"])
    code = "GMW000000036819"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}?code={code}"
    r = requests.get(url_groundwaterstation_code, auth=auth)

    r.raise_for_status()
    assert r.json()["count"] > 0, "request returned empty json"


def test_timeseries_rotterdam():
    lizard_GWS_endpoint = "https://rotterdam.lizard.net/api/v4/timeseries/"
    auth = ("__key__", os.environ["LIZARD_ROTTERDAM_API_KEY"])
    uuid = "286b0c4a-6a6a-4c5a-9c81-9e088ba6a5e7"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}{uuid}"
    r = requests.get(url_groundwaterstation_code, auth=auth)

    r.raise_for_status()


def test_single_observation_rotterdam():
    code = "GMW000000036819"
    auth = ("__key__", os.environ["LIZARD_ROTTERDAM_API_KEY"])
    o = hpd.GroundwaterObs.from_lizard(code, organisation="rotterdam", auth=auth)
    assert o.tube_nr == 1


def test_extent_rotterdam():
    extent = [68_500, 69_500, 443_500, 444_500]
    auth = ("__key__", os.environ["LIZARD_ROTTERDAM_API_KEY"])
    oc = hpd.read_lizard(extent, organisation="rotterdam", auth=auth)
    assert not oc.empty


@pytest.mark.slow
def test_codes_rotterdam():
    auth = ("__key__", os.environ["LIZARD_ROTTERDAM_API_KEY"])
    oc = hpd.read_lizard(
        codes=["GMW000000036819", "GMW000000037933"],
        type_timeseries="merge",
        organisation="rotterdam",
        auth=auth,
    )
    assert not oc.empty
