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

    knmi.request_url(url)  # , fname='test.txt')


def test_daily_meteo_url():
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie"
        "/gegevens/daggegevens/etmgeg_260.zip"
    )

    knmi.request_url(url)


def test_daily_rainfall_api():
    params = {"vars": "RD", "stns": "550", "start": "20200101", "end": "20200131"}
    knmi.request_api(knmi.URL_DAILY_PREC, params)


def test_daily_meteo_api():
    params = {"vars": "RH", "stns": "260", "start": "20200101", "end": "20200131"}
    knmi.request_api(knmi.URL_DAILY_METEO, params)


def test_hourly_meteo_api():
    params = {"vars": "RH", "stns": "260", "start": "2020010224", "end": "2020010301"}

    knmi.request_api(knmi.URL_HOURLY_METEO, params)
