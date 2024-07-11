# %%
import logging

import pandas as pd

import hydropandas as hpd
from hydropandas.io import knmi

logging.basicConfig(level=logging.DEBUG)


# %% test observations


def test_read_knmi_files():
    # neerslagstation
    knmi.get_knmi_obs(
        fname="./tests/data/2023-KNMI-test/neerslaggeg_ESBEEK_831.txt",
        start="2010-1-1",
        end="2010-1-10",
    )

    # neerslagstation
    knmi.get_knmi_obs(fname="./tests/data/2023-KNMI-test/neerslaggeg_VILSTEREN_342.txt")

    # neerslagstation
    knmi.get_knmi_obs(fname="./tests/data/2023-KNMI-test/precipitation_st_anthonis.txt")

    # neerslagstation
    knmi.get_knmi_obs(
        fname="./tests/data/2023-KNMI-test/etmgeg_260.txt", meteo_var="EV24"
    )

    # neerslagstation
    knmi.get_knmi_obs(
        fname="./tests/data/2023-KNMI-test/uurgeg_260_2001-2010.txt",
        meteo_var="RH",
        interval="hourly",
    )


def test_download_knmi_de_bilt():
    # De Bilt precipitation station daily api
    knmi.get_knmi_obs(
        stn=550,
        meteo_var="RD",
        start="2010-1-1",
        end="2010-1-10",
    )

    # De Bilt precipitation from meteo station daily api
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="RH",
        start="2010-1-1",
        end="2010-1-10",
    )

    # De Bilt evaporation from meteo station daily api
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="EV24",
        start="2010-1-1",
        end="2010-1-10",
    )

    # De Bilt precipitation meteostation hourly api
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="RH",
        interval="hourly",
        start=pd.Timestamp("2010-1-5 01:00"),
        end=pd.Timestamp("2010-1-6 23:00"),
    )

    # De Bilt precipitation station daily no api
    knmi.get_knmi_obs(
        stn=550, meteo_var="RD", use_api=False, start="2010-1-1", end="2010-1-10"
    )

    # De Bilt precipitation meteo station daily no api
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="RH",
        use_api=False,
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    # De Bilt evaporation station daily no api
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="EV24",
        use_api=False,
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )


def test_xy():
    # empty dataframe because nearest station has no data in time frame
    _ = knmi.get_knmi_obs(
        xy=(100000, 330000),
        meteo_var="EV24",
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    # fill missing obs does work
    _ = knmi.get_knmi_obs(
        xy=(100000, 330000),
        meteo_var="EV24",
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
        fill_missing_obs=True,
    )

    _ = knmi.get_knmi_obs(
        xy=(150000, 330000),
        meteo_var="RD",
        start=pd.Timestamp("1953-1-1"),
        end=pd.Timestamp("1953-1-10"),
    )


def test_calculate_evaporation():
    knmi.get_knmi_obs(
        stn=260,
        meteo_var="makkink",
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    knmi.get_knmi_obs(
        stn=260,
        meteo_var="penman",
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    knmi.get_knmi_obs(
        stn=260,
        meteo_var="hargreaves",
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )


def test_download_knmi_xy():
    df1, _ = knmi.get_knmi_obs(meteo_var="RH", stn=344)
    df2, _ = knmi.get_knmi_obs(meteo_var="RH", xy=(90600, 442800))

    assert df1.equals(df2), "Dataframes should be identical"


def test_schiermonnikoog_precipitation_station():
    knmi.get_knmi_obs(12, meteo_var="RD", start=pd.Timestamp("2023-1-1"), end=None)


def test_download_without_data():
    dfrd, _ = knmi.get_knmi_obs(
        324,
        meteo_var="RD",
        start=pd.Timestamp("2018"),
        end=pd.Timestamp("2020"),
        raise_exceptions=True,
    )
    assert dfrd.empty, "expected empty DataFrame"

    dfev, _ = knmi.get_knmi_obs(
        265,
        meteo_var="EV24",
        start=pd.Timestamp("1959"),
        end=pd.Timestamp("1963"),
        raise_exceptions=True,
    )
    assert dfev.empty, "expected empty DataFrame"


# %%
def test_fill_missing_measurements():
    settings = knmi._get_default_settings({"fill_missing_obs": True})

    # nothing is missing
    knmi.get_knmi_timeseries_stn(
        260,
        "RH",
        settings=settings,
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    # missing all data
    _, meta = knmi.get_knmi_timeseries_stn(
        265,
        meteo_var="EV24",
        settings=settings,
        start=pd.Timestamp("1959-1-1"),
        end=pd.Timestamp("1959-1-10"),
    )

    assert meta["station"] == 260, "expected metadata from different station"

    # missing some data
    _, meta = knmi.get_knmi_timeseries_stn(
        273,
        meteo_var="RH",
        settings=settings,
        start=pd.Timestamp("1998-9-1"),
        end=pd.Timestamp("1998-9-10"),
    )

    # no data at all (test is disabled because of too many requests)
    # df, meta = knmi.get_knmi_timeseries_stn(
    #     260,
    #     meteo_var="EV24",
    #     settings=settings,
    #     start=pd.Timestamp("1951-1-1"),
    #     end=pd.Timestamp("1951-1-10"),
    # )

    # assert df.empty, "expected empty dataframe"


# %% obs collections


def test_obslist_from_grid():
    xy = [[x, y] for x in [104150.0, 104550.0] for y in [510150.0, 510550.0]]

    knmi.get_knmi_obslist(
        xy=xy,
        meteo_vars=["RH", "EV24"],
        starts=["2010", "2010"],
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )


def test_obslist_from_locations():
    locations = pd.DataFrame(data={"x": [100000], "y": [350000]})
    knmi.get_knmi_obslist(
        locations, meteo_vars=["EV24"], ObsClasses=[hpd.EvaporationObs]
    )


def test_obslist_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    knmi.get_knmi_obslist(
        stns=stns,
        meteo_vars=["RH", "EV24"],
        starts=["2010", "2010"],
        ends=["2015", "2015"],
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )


def test_obslist_from_stns_single_startdate():
    stns = [344, 260]  # Rotterdam en de Bilt
    knmi.get_knmi_obslist(
        stns=stns,
        meteo_vars=["RH", "EV24"],
        starts="2010",
        ends="2015",
        ObsClasses=[hpd.PrecipitationObs, hpd.EvaporationObs],
    )
