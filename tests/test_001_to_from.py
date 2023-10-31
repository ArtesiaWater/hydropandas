import pandas as pd
import pytest
from requests.exceptions import ConnectionError

import hydropandas as hpd

# %% BRO


def test_bro_gld():
    # single observation
    bro_id = "GLD000000012893"
    hpd.GroundwaterObs.from_bro(bro_id)
    return


def test_bro_gmn():
    # single observation
    bro_id = "GMN000000000163"
    hpd.read_bro(bro_id=bro_id, only_metadata=True)
    return


def test_bro_extent():
    # extent = (210260, 213550, 459890, 473920)  # extent more than 1000 observations
    # extent = (213260, 213550, 473890, 473920)  # extent skip duplicates
    extent = (102395, 103121, 434331, 434750)  # 4 observations within extent

    hpd.read_bro(extent=extent, only_metadata=True)
    return


def test_bro_extent_too_big():
    extent = (102395, 213550, 334331, 473920)  # to many observations in extent

    with pytest.raises(RuntimeError):
        hpd.read_bro(extent=extent, only_metadata=True)


# %% DINO

dinozip = "./tests/data/2019-Dino-test/dino.zip"


def test_observation_gwq():
    # single observation
    path = "./tests/data/2019-Dino-test/Grondwatersamenstellingen_Put/B52C0057.txt"
    hpd.WaterQualityObs.from_dino(path)
    return


def test_observation_wl():
    path = "./tests/data/2019-Dino-test/Peilschaal/P58A0001.csv"
    hpd.WaterlvlObs.from_dino(path)
    return


def test_observation_gw():
    path = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv"
    o = hpd.GroundwaterObs.from_dino(path=path)
    return o


def test_obscollection_from_list():
    dino_gw = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    obs_list = [o for o in dino_gw.obs.values]
    hpd.ObsCollection(obs_list)
    return


def test_obscollection_from_df():
    df = pd.DataFrame(index=["pb1", "pb2"], data={"tube_nr": [1, 1]})

    hpd.ObsCollection(df)


def test_obscollection_empty():
    hpd.ObsCollection()


# read dino directories


def test_obscollection_dinozip_gw():
    # groundwater quantity
    oc = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=False,
    )
    return oc


def test_obscollection_dinozip_gw_keep_all_obs():
    # do not delete empty dataframes
    oc = hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    return oc


def test_obscollection_dinozip_wl():
    # surface water
    oc = hpd.read_dino(
        dirname=dinozip, ObsClass=hpd.WaterlvlObs, subdir="Peilschaal", suffix=".csv"
    )

    return oc


def test_obscollection_dinozip_gwq():
    # groundwater quality
    hpd.read_dino(
        dirname=dinozip,
        ObsClass=hpd.WaterQualityObs,
        subdir="Grondwatersamenstellingen_Put",
        suffix=".txt",
    )
    return


# %% FEWS
def test_obscollection_fews_highmemory():
    hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        to_mnap=False,
        remove_nan=False,
        low_memory=False,
    )
    return


def test_obscollection_fews_lowmemory():
    oc = hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=None,
        low_memory=True,
    )
    return oc


def test_obscollection_fews_selection():
    hpd.read_fews(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=("MPN-N-2",),
    )
    return


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
    dino_gw = test_obscollection_dinozip_gw()
    # drop duplicate
    dino_gw.drop("B22D0155-001", inplace=True)
    dino_gw.to_pastastore()


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
    obsc = test_obscollection_dinozip_gw()
    try:
        hpd.read_knmi(
            locations=obsc, meteo_vars=["EV24", "RD"], starts="2010", ends="2015"
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
    hpd.read_waterinfo(path)


# %% MENYANTHES


def test_obscollection_menyanthes():
    fname = "./tests/data/2023-MEN-test/test.men"
    hpd.read_menyanthes(fname, ObsClass=hpd.GroundwaterObs)
