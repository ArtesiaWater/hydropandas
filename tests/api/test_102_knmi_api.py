# %%
import logging

from hydropandas.io import knmi

logging.basicConfig(level=logging.DEBUG)


def test_daily_rainfall_url():
    stn = 550
    stn_name = "DE-BILT"
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie/"
        f"gegevens/monv_reeksen/neerslaggeg_{stn_name}_{stn}.zip"
    )

    f = knmi.request_url(url)
    assert "20200101" in f.read(), "date not found in response"


def test_daily_meteo_url():
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie"
        "/gegevens/daggegevens/etmgeg_260.zip"
    )

    f = knmi.request_url(url)
    assert "20200101" in f.read(), "date not found in response"


def test_daily_rainfall_api():
    params = {"vars": "RD", "stns": "550", "start": "20200101", "end": "20200131"}
    f = knmi.request_api(knmi.URL_DAILY_PREC, params)
    assert params["start"] in f.read(), "Start date not found in response"


def test_daily_meteo_api():
    params = {"vars": "RH", "stns": "260", "start": "20200101", "end": "20200131"}
    f = knmi.request_api(knmi.URL_DAILY_METEO, params)
    assert params["start"] in f.read(), "Start date not found in response"


def test_hourly_meteo_api():
    params = {"vars": "RH", "stns": "260", "start": "2020010224", "end": "2020010301"}
    f = knmi.request_api(knmi.URL_HOURLY_METEO, params)
    assert params["end"][:-2] in f.read(), "End date not found in response"
