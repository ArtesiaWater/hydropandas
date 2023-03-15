# -*- coding: utf-8 -*-
"""Created on Sat Dec 21 09:24:29 2019

@author: oebbe
"""

import numpy as np
import pandas as pd
import pytest
import hydropandas as hpd
from hydropandas.io import knmi

import logging

logging.basicConfig(level=logging.DEBUG)


def test_read_knmi_file1():
    # De Bilt neerslagstation
    fname = "./tests/data/2023-KNMI-test/neerslaggeg_ESBEEK_831.txt"
    knmi.read_knmi_daily_rainfall_file(fname)

    return


def test_read_knmi_file2():
    # De Bilt neerslagstation
    fname = "./tests/data/2023-KNMI-test/neerslaggeg_VILSTEREN_342.txt"
    knmi.read_knmi_daily_rainfall_file(fname)

    return


def test_read_knmi_file3():
    # De Bilt neerslagstation
    fname = "./tests/data/2023-KNMI-test/precipitation_st_anthonis.txt"
    knmi.read_knmi_daily_rainfall_file(fname)

    return


def test_get_knmi_precip_neerslagstation():
    # De Bilt neerslagstation
    knmi.get_knmi_timeseries_stn(
        "550", "RD", start="2010-1-1", end="2010-1-10", settings=None
    )
    return


def test_read_knmi_meteostation_file():

    fname = "./tests/data/2023-KNMI-test/etmgeg_260.txt"
    knmi.read_knmi_daily_meteo_file(fname, "EV24")

    return


def test_get_knmi_precip_meteostation_fill_missing():

    # De Bilt meteostation
    knmi.get_knmi_timeseries_stn(
        260, "RH", start="2010-1-1", end="2010-1-10", settings=None
    )
    return


def test_get_knmi_precip_meteostation_hourly():

    # De Bilt meteostation uurlijks
    knmi.get_knmi_timeseries_stn(
        260,
        "RH",
        start="2010-1-1",
        end="2010-1-10",
        settings={"interval": "hourly", "fill_missing_obs": False},
    )
    return


def test_get_pressure_hourly():

    # De Bilt meteostation uurlijks
    knmi.get_knmi_timeseries_stn(
        310,
        "P",
        start="2010-1-1",
        end="2010-1-10",
        settings={"interval": "hourly", "fill_missing_obs": False},
    )
    return


def test_get_knmi_precip_neerslagstation_no_api():
    # De Bilt neerslagstation
    knmi.get_knmi_timeseries_stn(
        "550", "RD", start="2010-1-1", end="2010-1-10", settings={"use_api": False}
    )
    return


def test_get_knmi_precip_meteostation_no_api():
    # De Bilt neerslagstation
    knmi.get_knmi_timeseries_stn(
        260, "RH", start="2010-1-1", end="2010-1-10", settings={"use_api": False}
    )
    return


def test_get_knmi_precip_neerslagstation_fill_missing():
    settings = knmi._get_default_settings(None)
    knmi.download_knmi_data(
        "550", meteo_var="RD", start="1952", end=None, settings=settings
    )
    return


def test_download_rd_550_no_api():
    settings = knmi._get_default_settings(None)
    settings["use_api"] = False
    knmi.download_knmi_data(
        "550",
        stn_name="DE-BILT",
        meteo_var="RD",
        start="1952",
        end=None,
        settings=settings,
    )
    return


def test_download_rd_12():
    settings = knmi._get_default_settings(None)
    knmi.download_knmi_data(
        "12", meteo_var="RD", start="2010", end=None, settings=settings
    )
    return


def test_download_without_data():
    settings = knmi._get_default_settings(None)
    try:
        knmi.download_knmi_data(
            "324", meteo_var="RD", start="2018", end="2020", settings=settings
        )
    except ValueError:
        pass

    return


def test_download_without_data_no_error():
    settings = knmi._get_default_settings(settings={"raise_exceptions": False})

    knmi_df, variables, stations = knmi.download_knmi_data(
        "324", meteo_var="RD", start="2018", end="2020", settings=settings
    )

    assert (knmi_df.empty, variables == {}, stations.empty) == (True, True, True)

    return


def test_download_ev24_210():
    settings = knmi._get_default_settings(None)

    knmi.download_knmi_data(
        210, meteo_var="EV24", start="1952", end=None, settings=settings
    )
    return


def test_download_ev24_210_no_api():
    settings = knmi._get_default_settings({"use_api": False})
    knmi.download_knmi_data(
        210, meteo_var="EV24", start="1952", end=None, settings=settings
    )
    return


def test_get_knmi_daily_meteo_ev24_empty():
    start, end = knmi._start_end_to_datetime("1959", "1963")
    knmi.get_knmi_daily_meteo_api(knmi.URL_DAILY_METEO, 265, "EV24", start, end, False)

    return


def test_fill_missing_measurements_ev24_278():
    knmi.fill_missing_measurements(278, meteo_var="EV24", start="1959", end="1963")
    return


def test_fill_missing_measurements_rh_273():
    knmi.fill_missing_measurements(273, meteo_var="RH", start="1995", end="2000")
    return


def test_fill_missing_measurements_rh_278():
    settings = knmi._get_default_settings({"exclude_neerslag_stn": True})
    knmi.fill_missing_measurements(
        278, meteo_var="RH", start="1990", end=None, settings=settings
    )
    return


def test_obslist_from_grid():
    xy = [[x, y] for x in [104150.0, 104550.0] for y in [510150.0, 510550.0]]

    knmi.get_knmi_obslist(
        xy=xy,
        meteo_vars=["RH", "EV24"],
        starts=["2010", "2010"],
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )

    return


def test_obslist_from_locations():
    locations = pd.DataFrame(data={"x": [100000], "y": [350000]})
    knmi.get_knmi_obslist(
        locations, meteo_vars=["EV24"], ObsClasses=[hpd.EvaporationObs]
    )

    return


def test_obslist_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    knmi.get_knmi_obslist(
        stns=stns,
        meteo_vars=["RH", "EV24"],
        starts=["2010", "2010"],
        ends=["2015", "2015"],
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )

    return


def test_obslist_from_stns_single_startdate():
    stns = [344, 260]  # Rotterdam en de Bilt
    knmi.get_knmi_obslist(
        stns=stns,
        meteo_vars=["RH", "EV24"],
        starts="2010",
        ends="2015",
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )

    return
