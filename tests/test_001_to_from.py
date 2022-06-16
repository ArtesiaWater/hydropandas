# import os
import numpy as np
import pandas as pd
import pytest
from hydropandas import obs_collection as oc
from hydropandas import observation as obs

# %% DINO

dinozip = "./tests/data/2019-Dino-test/dino.zip"


def test_observation_gwq():
    # single observation
    fname = "./tests/data/2019-Dino-test/Grondwatersamenstellingen_Put/B52C0057.txt"
    ogq = obs.GroundwaterQualityObs.from_dino(fname)
    return ogq


def test_observation_wl():
    fname = "./tests/data/2019-Dino-test/Peilschaal/P58A0001.csv"
    wl = obs.WaterlvlObs.from_dino(fname)
    return wl


def test_observation_gw():
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv"
    gw = obs.GroundwaterObs.from_dino(fname=fname)
    return gw


def test_obscollection_from_list():
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    obs_list = [o for o in dino_gw.obs.values]
    oc_list = oc.ObsCollection.from_list(obs_list)
    return oc_list


def test_obscollection_from_df():
    df = pd.DataFrame(index=["pb1", "pb2"], data={"filternr": [1, 1]})

    df_oc = oc.ObsCollection.from_dataframe(df)

    return df_oc


# read dino directories
def test_obscollection_dinozip_gw():
    # groundwater quantity
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=False,
    )
    return dino_gw


def test_obscollection_dinozip_gw_keep_all_obs():
    # do not delete empty dataframes
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        keep_all_obs=True,
    )
    return dino_gw


def test_obscollection_dinozip_wl():
    # surface water
    dino_ps = oc.ObsCollection.from_dino(
        dirname=dinozip, ObsClass=obs.WaterlvlObs, subdir="Peilschaal", suffix=".csv"
    )
    return dino_ps


def test_obscollection_dinozip_gwq():
    # groundwater quality
    dino_gwq = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterQualityObs,
        subdir="Grondwatersamenstellingen_Put",
        suffix=".txt",
    )
    return dino_gwq


def test_obscollection_dino_download_bbox_empty():
    # download DINO from bbox
    bbox = [88596.63500000164, 407224.8449999988, 89623.4149999991, 407804.27800000086]
    dino_gw_bbox = oc.ObsCollection.from_dino(bbox=bbox, ObsClass=obs.GroundwaterObs)
    return dino_gw_bbox


# %% FEWS
def test_obscollection_fews_highmemory():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        translate_dic={"locationId": "locatie"},
        to_mnap=False,
        remove_nan=False,
        low_memory=False,
    )
    return fews_gw_prod


def test_obscollection_fews_lowmemory():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=None,
        low_memory=True,
    )
    return fews_gw_prod


def test_obscollection_fews_selection():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=("MPN-N-2",),
    )
    return fews_gw_prod


# %% WISKI
@pytest.mark.slow
def test_observation_wiskicsv_gw():
    wiski_gw = obs.GroundwaterObs.from_wiski(
        "./tests/data/2019-WISKI-test/1016_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )

    return wiski_gw


@pytest.mark.slow
def test_obscollection_wiskizip_gw():
    wiski_col = oc.ObsCollection.from_wiski(
        r"./tests/data/2019-WISKI-test/1016_PBF.zip",
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
        sep=r"\s+",
        header_sep=":",
        dayfirst=True,
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
    )

    return wiski_col


# %% PASTASTORE
def test_to_pastastore():

    dino_gw = test_obscollection_dinozip_gw()
    pstore = dino_gw.to_pastastore()

    return pstore


#%% Evaporation

def test_evap_obs_from_stn():
    return obs.EvaporationObs.from_knmi(260, et_type='EV24')

def test_evap_obs_from_stn_makkink():
    return obs.EvaporationObs.from_knmi(260, et_type='makkink')

def test_evap_obs_from_stn_penman():
    return obs.EvaporationObs.from_knmi(260, et_type='penman')

def test_evap_obs_from_stn_hargreaves():
    return obs.EvaporationObs.from_knmi(260, et_type='hargreaves')

#%% Precipitation

def test_precip_obs_from_stn():
    return obs.PrecipitationObs.from_knmi(233, "precipitation")


def test_knmi_obs_from_stn_no_api():
    return obs.PrecipitationObs.from_knmi(233, "precipitation", use_api=False)


def test_knmi_obs_from_stn_without_any_data():

    obs.EvaporationObs.from_knmi(
        210, startdate="19500101", enddate="19600101", fill_missing_obs=False
    )

    return 1


def test_knmi_obs_from_stn_with_missing_data_in_time_period():
    return obs.PrecipitationObs.from_knmi("441", "precipitation", startdate="2010-1-2")


def test_knmi_obs_from_xy():
    return obs.PrecipitationObs.from_nearest_xy(100000, 350000)


def test_knmi_obs_from_obs():
    pb = test_observation_gw()
    return obs.PrecipitationObs.from_obs(pb, fill_missing_obs=False)


def test_knmi_collection_from_locations():
    obsc = test_obscollection_dinozip_gw()
    oc_knmi = oc.ObsCollection.from_knmi(
        locations=obsc, meteo_vars=["EV24", "RD"], start="2010", end="2015", cache=False
    )
    return oc_knmi


def test_knmi_collection_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    oc_knmi = oc.ObsCollection.from_knmi(
        stns=stns,
        meteo_vars=["EV24", "RH"],
        start=["2010", "2010"],
        end=["2015", "2015"],
        ObsClass=[obs.EvaporationObs, obs.PrecipitationObs],
    )
    return oc_knmi


def test_knmi_collection_from_grid():
    # somewhere in Noord-Holland (near Castricum)
    x = np.array([104150.0, 104550.0])
    y = np.array([510150.0, 510550.0])
    oc_knmi = oc.ObsCollection.from_knmi(
        x=x,
        y=y,
        meteo_vars=["RH"],
        start=["2010", "2010"],
        end=["2015", "2015"],
    )
    return oc_knmi


# %% WATERINFO


def test_waterinfo_from_dir():
    path = "./tests/data/waterinfo-test"
    wi = oc.ObsCollection.from_waterinfo(path)
    return wi


# %% MENYANTHES (still need a small menyanthes file to do the test)

# def test_obscollection_menyanthes():
#
#    fname = r'export_from_ADI.men'
#    obsc = oc.ObsCollection.from_menyanthes(fname)
#
#    return obsc
