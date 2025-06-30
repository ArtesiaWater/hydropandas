import datetime as dt
import os

import pandas as pd
import pastastore as pst
import pytest
from requests.exceptions import ConnectionError, HTTPError

import hydropandas as hpd

# %% BRO


def test_bro_gld():
    # single observation
    bro_id = "GLD000000012893"
    hpd.GroundwaterObs.from_bro(bro_id)


def test_bro_gmn():
    # single observation
    bro_id = "GMN000000000001"  # 6 objects as per 2025-4-17
    hpd.read_bro(bro_id=bro_id, only_metadata=True)


def test_bro_extent():
    # extent = (210260, 213550, 459890, 473920)  # extent more than 1000 observations
    # extent = (213260, 213550, 473890, 473920)  # extent skip duplicates
    extent = (102395, 103121, 434331, 434750)  # 4 observations within extent

    hpd.read_bro(extent=extent, only_metadata=True)


def test_bro_extent_too_big():
    extent = (102395, 213550, 334331, 473920)  # too many observations in extent

    with pytest.raises(HTTPError):
        hpd.read_bro(extent=extent, only_metadata=True)


# %% DINO
dinozip = "./tests/data/2019-Dino-test/dino.zip"


def test_observation_gwq():
    # single observation
    path = "./tests/data/2019-Dino-test/Grondwatersamenstellingen_Put/B52C0057.txt"
    hpd.WaterQualityObs.from_dino(path)


def test_observation_wl():
    path = "./tests/data/2019-Dino-test/Peilschaal/P58A0001.csv"
    hpd.WaterlvlObs.from_dino(path)


def observation_gw_dino_old():
    path = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv"
    o = hpd.GroundwaterObs.from_dino(path=path)
    return o


def observation_gw_dino_new():
    path = "./tests/data/2024-Dino-test/DINO_Grondwaterstanden/B02H0090001.csv"
    o = hpd.GroundwaterObs.from_dino(path=path)
    return o


def test_observation_gw():
    observation_gw_dino_old()
    observation_gw_dino_new()


def test_obscollection_from_directory_old_school():
    dino_gw = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    obs_list = [o for o in dino_gw.obs.values]
    hpd.ObsCollection(obs_list)


def test_obscollection_from_directory_new_school():
    dinozip = "./tests/data/2024-Dino-test/dino.zip"
    hpd.read_dino(dirname=dinozip, ObsClass=hpd.GroundwaterObs, keep_all_obs=True)


def test_obscollection_from_df():
    df = pd.DataFrame(index=["pb1", "pb2"], data={"tube_nr": [1, 1]})

    hpd.ObsCollection(df)


def test_obscollection_empty():
    hpd.ObsCollection()


# read dino directories


def obscollection_dinozip_gw():
    # groundwater quantity
    oc = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=False,
    )
    return oc


def obscollection_dinozip_gw_keep_all_obs():
    # do not delete empty dataframes
    oc = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    return oc


def obscollection_dinozip_wl():
    # surface water
    oc = hpd.read_dino(
        dirname=dinozip, ObsClass=hpd.WaterlvlObs, subdir="Peilschaal", suffix=".csv"
    )

    return oc


def test_obscollection_dinozip_gw():
    # groundwater quantity
    obscollection_dinozip_gw()


def test_obscollection_dinozip_gw_keep_all_obs():
    obscollection_dinozip_gw_keep_all_obs()


def test_obscollection_dinozip_wl():
    obscollection_dinozip_wl()


def test_obscollection_dinozip_gwq():
    # groundwater quality
    hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.WaterQualityObs,
        subdir="Grondwatersamenstellingen_Put",
        suffix=".txt",
    )


# %% FEWS


def obscollection_fews_lowmemory():
    oc = hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=None,
        low_memory=True,
    )
    return oc


def test_obscollection_fews_highmemory():
    hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        to_mnap=False,
        remove_nan=False,
        low_memory=False,
    )


def test_obscollection_fews_lowmemory():
    obscollection_fews_lowmemory()


def test_obscollection_fews_selection():
    hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=("MPN-N-2",),
    )


# %% WISKI
@pytest.mark.slow
def test_observation_wiskicsv_gw():
    hpd.GroundwaterObs.from_wiski(
        "./tests/data/2019-WISKI-test/1016_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        dayfirst=True,
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )


@pytest.mark.slow
def test_obscollection_wiskizip_gw():
    hpd.read_wiski(
        r"./tests/data/2019-WISKI-test/1016_PBF.zip",
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
        sep=r"\s+",
        header_sep=":",
        dayfirst=True,
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
    )


# %% PASTASTORE
def test_to_pastastore():
    dino_gw = obscollection_dinozip_gw()
    # drop duplicate
    dino_gw.drop("B22D0155-001", inplace=True)
    pstore = dino_gw.to_pastastore()
    # export to zip for read test
    pstore.to_zip("test_pastastore.zip")


def test_from_pastastore():
    pstore = pst.PastaStore.from_zip(
        "test_pastastore.zip", conn=pst.DictConnector("pastas_db")
    )
    _ = hpd.read_pastastore(pstore, "oseries")
    os.remove("test_pastastore.zip")


# %% excel


def test_to_excel():
    oc = hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=None,
        low_memory=True,
    )

    oc.to_excel("tests/data/excel/test.xlsx")


def test_from_excel():
    hpd.read_excel("tests/data/excel/test.xlsx")


# %% Meteo


def test_pressure_obs_from_stn():
    hpd.MeteoObs.from_knmi(
        stn=310, meteo_var="P", interval="hourly", start="2022-1-1", end="2023-1-1"
    )


def test_pressure_read_knmi():
    hpd.read_knmi(
        stns=(310,),
        meteo_vars=("P",),
        interval="hourly",
        starts="2022-1-1",
        ends="2023-1-1",
    )


# %% Evaporation


def test_evap_obs_from_file():
    path = "./tests/data/2023-KNMI-test/etmgeg_260.txt"
    hpd.EvaporationObs.from_knmi(fname=path)


def test_evap_obs_from_stn():
    hpd.EvaporationObs.from_knmi(stn=260, meteo_var="EV24")


def test_evap_obs_from_stn_makkink():
    hpd.EvaporationObs.from_knmi(stn=260, meteo_var="makkink")


def test_evap_obs_from_stn_penman():
    hpd.EvaporationObs.from_knmi(stn=260, meteo_var="penman")


def test_evap_obs_from_stn_hargreaves():
    hpd.EvaporationObs.from_knmi(stn=260, meteo_var="hargreaves")


# %% Precipitation


def test_precip_obs_from_file():
    path = "./tests/data/2023-KNMI-test/neerslaggeg_ESBEEK_831.txt"
    hpd.PrecipitationObs.from_knmi(fname=path)


def test_precip_obs_from_stn():
    hpd.PrecipitationObs.from_knmi(stn=233, meteo_var="RD")


def test_knmi_obs_from_stn_no_api():
    hpd.PrecipitationObs.from_knmi(stn=233, meteo_var="RD", use_api=False)


def test_knmi_obs_from_stn_without_any_data():
    hpd.EvaporationObs.from_knmi(
        stn=210, start="19500101", end="19600101", fill_missing_obs=False
    )


def test_knmi_obs_from_stn_with_missing_data_in_time_period():
    hpd.PrecipitationObs.from_knmi(stn=441, meteo_var="RD", start="2010-1-2")


def test_knmi_obs_from_xy():
    hpd.PrecipitationObs.from_knmi(xy=(100000, 350000))


# @pytest.xfail(
#     "Station HEIBLOEM 967 not available. See issue"
#     " https://github.com/ArtesiaWater/hydropandas/issues/103"
# )
def test_knmi_collection_from_locations():
    obsc = obscollection_dinozip_gw()
    try:
        hpd.read_knmi(
            locations=obsc,
            meteo_vars=["EV24", "RD"],
            starts="2010",
            ends="2015",
            raise_exceptions=True,
        )
    except ConnectionError:
        pass


def test_knmi_collection_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    hpd.read_knmi(
        stns=stns,
        meteo_vars=["EV24", "RH"],
        starts=["2010", "2010"],
        ends=["2015", "2015"],
    )


def test_knmi_collection_from_grid():
    # somewhere in Noord-Holland (near Castricum)
    xy = [[104150.0, 510150.0], [104550.0, 510550.0]]
    hpd.read_knmi(
        xy=xy,
        meteo_vars=["RH"],
        starts=["2010"],
        ends=["2015"],
    )


# %% WATERINFO


def test_waterinfo_from_dir():
    path = "./tests/data/2023-waterinfo-test"
    hpd.read_waterinfo(file_or_dir=path)


def test_waterinfo_ddlpy():
    grootheid_code = "WATHTE"
    locatie = "SCHOONHVN"
    tmin = dt.datetime(2020, 1, 1)
    tmax = dt.datetime(2020, 1, 5)
    hpd.WaterlvlObs.from_waterinfo(
        grootheid_code=grootheid_code, locatie=locatie, tmin=tmin, tmax=tmax
    )


def test_waterinfo_ddlpy_extent():
    grootheid_code = "WATHTE"
    tmin = dt.datetime(2020, 1, 1)
    tmax = dt.datetime(2020, 1, 2)
    extent = (110000, 125000, 429550, 449900)
    oc = hpd.read_waterinfo(
        extent=extent, grootheid_code=grootheid_code, tmin=tmin, tmax=tmax
    )
    assert not oc.empty


# %% MENYANTHES


def test_obscollection_menyanthes():
    fname = "./tests/data/2023-MEN-test/test.men"
    hpd.read_menyanthes(fname, ObsClass=hpd.GroundwaterObs)
