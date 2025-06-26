# %%
import logging

import pandas as pd

import hydropandas as hpd
from hydropandas.io import knmi

logging.basicConfig(level=logging.DEBUG)


# compare api calls with pre-downloaded files
def test_knmi_meteo_station_hourly_api_values():
    stn = 260
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2001-01-01")
    df, meta = knmi.get_hourly_meteo_api(stn=stn, meteo_var="RH", start=start, end=end)
    df2, _ = knmi.interpret_knmi_file(
        df,
        meta=meta,
        meteo_var="RH",
        start=start,
        end=end,
        add_day=False,
        add_hour=True,
    )
    truth, _ = knmi.parse_data("./tests/data/2023-KNMI-test/uurgeg_260_1991-2000.txt")

    # check raw data
    pd.testing.assert_series_equal(
        df["RH"].loc["2000-01-01 01:00:00":"2001-01-01 00:00:00"], truth["RH"]
    )
    # check after interpretation, since interpretation converts to UTC+1,
    # values have shifted 1h
    assert (
        df2.loc["2000-01-01 06:00:00", "RH"] * 1e4
        == truth.loc["2000-01-01 05:00:00", "RH"]
    )


def test_knmi_meteo_station_daily_api_values():
    stn = 260
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2001-01-01")
    df, meta = knmi.get_daily_meteo_api(stn=stn, meteo_var="RH", start=start, end=end)
    df3, _ = knmi.interpret_knmi_file(
        df,
        meta,
        meteo_var="RH",
        start=start,
        end=end,
        add_day=True,
        add_hour=True,
    )
    truth, _ = knmi.parse_data("./tests/data/2023-KNMI-test/etmgeg_260.txt")
    # check raw data
    pd.testing.assert_series_equal(
        df["RH"].loc["2000"], truth["RH"].loc["2000"], check_dtype=False
    )
    # check after interpretation, registration moved to end of day
    assert df3.loc["2000-01-02 01:00:00", "RH"] * 1e4 == truth.loc["2000-01-01", "RH"]


def test_knmi_meteo_station_daily_url_values():
    stn = 260
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2001-01-01")
    # daily data from meteorological stations
    df, meta = knmi.get_daily_meteo_url(stn=stn)
    df2, _ = knmi.interpret_knmi_file(
        df,
        meta,
        meteo_var="RH",
        start=start,
        end=end,
        add_day=True,
        add_hour=True,
    )
    truth, _ = knmi.parse_data("./tests/data/2023-KNMI-test/etmgeg_260.txt")
    # check raw data
    pd.testing.assert_series_equal(
        df["RH"].loc["2000"], truth["RH"].loc["2000"], check_dtype=False
    )
    # check after interpretation, registration moved to end of day
    assert df2.loc["2000-01-02 01:00:00", "RH"] * 1e4 == truth.loc["2000-01-01", "RH"]

    # also check if equal to api result
    df_api, meta_api = knmi.get_daily_meteo_api(
        stn=stn, meteo_var="RH", start=start, end=end
    )
    df_api, _ = knmi.interpret_knmi_file(
        df_api,
        meta_api,
        meteo_var="RH",
        start=start,
        end=end,
        add_day=True,
        add_hour=True,
    )
    pd.testing.assert_frame_equal(df2, df_api)


def test_knmi_daily_rainfall_api_values():
    stn = 550
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2000-12-31")
    df, meta = knmi.get_daily_rainfall_api(stn=stn, start=start, end=end)
    df2, _ = knmi.interpret_knmi_file(
        df,
        meta,
        meteo_var="RD",
        start=start,
        end=end,
        add_day=False,
        add_hour=True,
    )
    truth, _ = knmi.parse_data(
        "./tests/data/2023-KNMI-test/neerslaggeg_DE-BILT_550.txt"
    )
    # check raw data
    pd.testing.assert_series_equal(df.loc[start:end, "RD"], truth["RD"])
    # check after interpretation
    pd.testing.assert_series_equal(
        df2["RD"] * 1e4,
        truth["RD"],
        check_index=False,
        check_dtype=False,
        atol=1e-8,
        rtol=1e-8,
    )


def test_knmi_daily_rainfall_url_values():
    stn = 550
    stn_name = knmi.get_station_name(stn=stn)
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2000-12-31")
    # daily data from rainfall-stations
    df, meta = knmi.get_daily_rainfall_url(stn, stn_name)

    df2, _ = knmi.interpret_knmi_file(
        df,
        meta,
        meteo_var="RD",
        start=start,
        end=end,
        add_day=False,
        add_hour=True,
    )

    truth, _ = knmi.parse_data(
        "./tests/data/2023-KNMI-test/neerslaggeg_DE-BILT_550.txt"
    )
    # check raw data
    pd.testing.assert_series_equal(
        df.loc[start:end, "RD"], truth["RD"], check_dtype=False
    )
    # check after interpretation
    pd.testing.assert_series_equal(
        df2["RD"] * 1e4,
        truth["RD"],
        check_index=False,
        check_dtype=False,
        atol=1e-8,
        rtol=1e-8,
    )


# test reading other files
def test_read_daily_rainfall():
    # neerslagstation
    knmi.get_knmi_obs(
        fname="./tests/data/2023-KNMI-test/neerslaggeg_ESBEEK_831.txt",
        start="2010-1-1",
        end="2010-1-10",
    )


def test_read_daily_rainfall2():
    # neerslagstation
    knmi.get_knmi_obs(fname="./tests/data/2023-KNMI-test/neerslaggeg_VILSTEREN_342.txt")


def test_read_daily_rainfall3():
    # neerslagstation
    knmi.get_knmi_obs(fname="./tests/data/2023-KNMI-test/precipitation_st_anthonis.txt")


def test_read_hourly_meteo():
    # neerslagstation
    knmi.get_knmi_obs(
        fname="./tests/data/2023-KNMI-test/uurgeg_260_2001-2010.txt",
        meteo_var="RH",
        interval="hourly",
    )


def test_xy():
    # get only results from 1951-7-1
    _, meta1 = knmi.get_knmi_obs(
        xy=(150000, 330000),
        meteo_var="RD",
        start=pd.Timestamp("1951-1-1"),
        end=pd.Timestamp("1952-1-1"),
    )

    assert meta1["station"] == 976

    # fill results before 1951-7-1 with nearest station
    _ = knmi.get_knmi_obs(
        xy=(150000, 330000),
        meteo_var="RD",
        start=pd.Timestamp("1951-1-1"),
        end=pd.Timestamp("1952-1-1"),
        fill_missing_obs=True,
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


def test_download_without_data():
    dfrd, _ = knmi.get_knmi_obs(
        324, meteo_var="RD", start=pd.Timestamp("2018"), end=pd.Timestamp("2020")
    )
    assert dfrd.empty, "expected empty DataFrame"

    dfev, _ = knmi.get_knmi_obs(
        265, meteo_var="EV24", start=pd.Timestamp("1959"), end=pd.Timestamp("1963")
    )
    assert dfev.empty, "expected empty DataFrame"


def test_fill_missing_measurements():
    settings = knmi._get_default_settings({"fill_missing_obs": True})

    # nothing is missing
    knmi.get_timeseries_stn(
        260,
        "RH",
        settings=settings,
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    # missing all data
    df, meta = knmi.get_timeseries_stn(
        265,
        meteo_var="EV24",
        settings=settings,
        start=pd.Timestamp("1959-1-1"),
        end=pd.Timestamp("1959-1-10"),
    )

    assert meta["station"] == 265, "expected metadata from different station"
    assert (df["station"] == "260").all(), "expected metadata from different station"

    # missing some data
    _, meta = knmi.get_timeseries_stn(
        273,
        meteo_var="RH",
        settings=settings,
        start=pd.Timestamp("1998-9-8"),
        end=pd.Timestamp("1998-9-10"),
    )

    # no data at all (test is disabled because of too many requests)
    # df, meta = knmi.get_timeseries_stn(
    #     260,
    #     meteo_var="EV24",
    #     settings=settings,
    #     start=pd.Timestamp("1951-1-1"),
    #     end=pd.Timestamp("1951-1-10"),
    # )
    # assert df.empty, "expected empty dataframe"

    # # maximum fill meteostation den Bosch
    # # loops through almost all meteo stations because there is no data in part of 1945
    # df, meta = knmi.get_timeseries_stn(
    #     265,
    #     meteo_var="RH",
    #     settings=settings,
    #     start=pd.Timestamp("1895-1-1"),
    #     end=pd.Timestamp.now(),
    # )
    # assert not df.empty, "expected filled df"


def test_fill_missing_measurements_neerslag():
    settings = knmi._get_default_settings({"fill_missing_obs": True})

    # nothing is missing RD neerslagstation De Bilt
    knmi.get_timeseries_stn(
        550,
        "RD",
        settings=settings,
        start=pd.Timestamp("2010-1-1"),
        end=pd.Timestamp("2010-1-10"),
    )

    # missing all data RD neerslagstation De Bilt
    knmi.get_timeseries_stn(
        550,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp("1853-1-1"),
        end=pd.Timestamp("1853-1-10"),
    )

    # missing some data at the start RD meteostation 72
    df, meta = knmi.get_timeseries_stn(
        72,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp("1986-03-25"),
        end=pd.Timestamp("1986-03-29"),
    )
    assert df["station"].astype(int).isin([72, 78]).all(), (
        "expected data from these two stations"
    )

    # missing some data at the end RD meteostation 57
    df, meta = knmi.get_timeseries_stn(
        57,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp("1942-11-29"),
        end=pd.Timestamp("1942-12-3"),
    )
    assert df["station"].astype(int).isin([57, 79]).all(), (
        "expected data from these two stations"
    )

    # no historical data at any neerslagstation for this time period
    df, meta = knmi.get_timeseries_stn(
        550,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp("1800-1-1"),
        end=pd.Timestamp("1800-1-10"),
    )
    assert df.empty, "expected empty dataframe"

    # no recent data at any neerslagstation for this time period
    df, meta = knmi.get_timeseries_stn(
        550,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp.now(),
        end=pd.Timestamp.now() + pd.Timedelta(days=2),
    )
    assert df.empty, "expected empty dataframe"

    # no data for neerslagstation de Bilt for this time period
    df, meta = knmi.get_timeseries_stn(
        550,
        meteo_var="RD",
        settings=settings,
        start=pd.Timestamp("1896-1-1"),
        end=pd.Timestamp("1896-1-10"),
    )
    assert not df.empty, "expected filled df"

    # # maximum fill neerslagstation den Bosch
    # # loops through all neerslagstations because there is not measurement at 29-10-1885
    # df, meta = knmi.get_timeseries_stn(
    #     872,
    #     meteo_var="RD",
    #     settings=settings,
    #     start=pd.Timestamp("1800-1-1"),
    #     end=pd.Timestamp.now(),
    # )
    # assert not df.empty, "expected filled df"


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


def test_knmi_daily_rainfall():
    stn = 550
    meteo_var = "RD"
    stn_name = knmi.get_station_name(stn=stn)
    start = pd.Timestamp("2000-01-01")
    end = pd.Timestamp("2000-12-31")

    # daily data from rainfall-stations
    df, meta = knmi.get_daily_rainfall_url(stn, stn_name)

    df2, _ = knmi.interpret_knmi_file(
        df,
        meta,
        meteo_var=meteo_var,
        start=start,
        end=end,
        add_day=False,
        add_hour=True,
    )

    truth, _ = knmi.parse_data(
        "./tests/data/2023-KNMI-test/neerslaggeg_DE-BILT_550.txt"
    )

    # check raw data
    pd.testing.assert_series_equal(
        df.loc[start:end, "RD"], truth["RD"], check_dtype=False
    )
    # check after interpretation
    pd.testing.assert_series_equal(
        df2["RD"] * 1e4,
        truth["RD"],
        check_index=False,
        check_dtype=False,
        atol=1e-8,
        rtol=1e-8,
    )


def test_meteo_station_methods():
    """Test if 3 different ways (api, non-api and .txt file) of obtaining data from a meteo
    station yield the same results"""

    start = pd.Timestamp("2000-1-1")
    end = pd.Timestamp("2000-1-10")

    # daily meteo station (default)
    precip1 = hpd.PrecipitationObs.from_knmi(stn=260, start=start, end=end)

    # daily meteo station without api
    precip2 = hpd.PrecipitationObs.from_knmi(
        stn=260,
        start=start,
        end=end,
        use_api=False,
    )

    # daily meteo station from file
    precip3 = hpd.PrecipitationObs.from_knmi(
        fname=r"./tests/data/2023-KNMI-test/etmgeg_260.txt",
        start=start,
        end=end,
    )

    assert precip1.equals(precip2)
    assert precip1.equals(precip3)


def test_rainfall_station_methods():
    """Test if 3 different ways (api, non-api and .txt file) of obtaining data from a rainfall
    station yield the same results"""

    start = pd.Timestamp("2000-1-1")
    end = pd.Timestamp("2000-1-10")

    # daily rainfall station
    precip1 = hpd.PrecipitationObs.from_knmi(
        stn=550, meteo_var="RD", start=start, end=end
    )

    # daily rainfall station without api
    precip2 = hpd.PrecipitationObs.from_knmi(
        meteo_var="RD",
        stn=550,
        start=start,
        end=end,
        use_api=False,
    )

    # daily rainfall station from file
    precip3 = hpd.PrecipitationObs.from_knmi(
        fname=r"./tests/data/2023-KNMI-test/neerslaggeg_DE-BILT_550.txt",
        start=start,
        end=end,
    )

    assert precip1.equals(precip2)
    assert precip1.equals(precip3)
