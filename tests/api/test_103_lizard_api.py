# %%

import requests


def test_groundwaterstations():
    lizard_GWS_endpoint = "https://vitens.lizard.net/api/v4/groundwaterstations/"
    code = "27BP0003"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}?code={code}"
    r = requests.get(url_groundwaterstation_code)

    r.raise_for_status()


def test_timeseries():
    lizard_GWS_endpoint = "https://vitens.lizard.net/api/v4/timeseries/"
    uuid = "44d5b8f9-5680-4fb5-8010-30c85632d106"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}{uuid}"
    r = requests.get(url_groundwaterstation_code)

    r.raise_for_status()
