# import os
import numpy as np
import pandas as pd
import pytest
from hydropandas import obs_collection as oc
from hydropandas import observation as obs

# %% DINO

dinozip = './tests/data/2019-Dino-test/dino.zip'


def test_observation_gwq():
    # single observation
    fname = './tests/data/2019-Dino-test/Grondwatersamenstellingen_Put/B52C0057.txt'
    ogq = obs.GroundwaterQualityObs.from_dino(fname, verbose=True)
    return ogq


def test_observation_wl():
    fname = './tests/data/2019-Dino-test/Peilschaal/P58A0001.csv'
    wl = obs.WaterlvlObs.from_dino(fname, verbose=True)
    return wl


def test_observation_gw():
    fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
    gw = obs.GroundwaterObs.from_dino(fname=fname, verbose=True)
    return gw


def test_obscollection_fieldlogger():
    # collection of observations
    fl = oc.ObsCollection.from_fieldlogger(
        './tests/data/2019-Dino-test/fieldlogger/locations.csv')
    return fl


def test_obscollection_from_list():
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        keep_all_obs=True,
        verbose=False)
    obs_list = [o for o in dino_gw.obs.values]
    oc_list = oc.ObsCollection.from_list(obs_list)
    return oc_list


def test_obscollection_from_df():
    df = pd.DataFrame(index=['pb1', 'pb2'],
                      data={'filternr': [1, 1]})

    df_oc = oc.ObsCollection.from_dataframe(df)

    return df_oc


# read dino directories
def test_obscollection_dinozip_gw():
    # groundwater quantity
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        keep_all_obs=False,
        verbose=False)
    return dino_gw


def test_obscollection_dinozip_gw_keep_all_obs():
    # do not delete empty dataframes
    dino_gw = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        keep_all_obs=True,
        verbose=False)
    return dino_gw


def test_obscollection_dinozip_wl():
    # surface water
    dino_ps = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.WaterlvlObs,
        subdir='Peilschaal',
        suffix='.csv',
        verbose=True)
    return dino_ps


def test_obscollection_dinozip_gwq():
    # groundwater quality
    dino_gwq = oc.ObsCollection.from_dino(
        dirname=dinozip,
        ObsClass=obs.GroundwaterQualityObs,
        subdir='Grondwatersamenstellingen_Put',
        suffix='.txt',
        verbose=True)
    return dino_gwq


def test_obscollection_dino_download_bbox_empty():
    # download DINO from bbox
    bbox = [88596.63500000164, 407224.8449999988,
            89623.4149999991, 407804.27800000086]
    dino_gw_bbox = oc.ObsCollection.from_dino(
        bbox=bbox, ObsClass=obs.GroundwaterObs, verbose=True)
    return dino_gw_bbox


# collection methods
def test_obscollection_to_fieldlogger():
    dino_gw = test_obscollection_dinozip_gw()
    fdf = dino_gw.to_fieldlogger(
        './tests/data/2019-Dino-test/fieldlogger/locations.csv', verbose=True)
    return fdf


# %% FEWS
def test_obscollection_fews_highmemory():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        './tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip',
        translate_dic={'locationId': 'locatie'},
        verbose=True,
        to_mnap=False,
        remove_nan=False,
        low_memory=False)
    return fews_gw_prod


def test_obscollection_fews_lowmemory():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        './tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip',
        verbose=True,
        locations=None,
        low_memory=True)
    return fews_gw_prod


def test_obscollection_fews_selection():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        './tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip',
        verbose=True,
        locations=("MPN-N-2",)
    )
    return fews_gw_prod


# %% WISKI
@pytest.mark.slow
def test_observation_wiskicsv_gw():
    wiski_gw = obs.GroundwaterObs.from_wiski(
        "./tests/data/2019-WISKI-test/1016_PBF.csv",
        sep=r'\s+',
        header_sep=':',
        header_identifier=':',
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        translate_dic={
            'name': 'Station Number',
            'x': 'GlobalX',
            'y': 'GlobalY'},
        verbose=True)

    return wiski_gw


@pytest.mark.slow
def test_obscollection_wiskizip_gw():
    wiski_col = oc.ObsCollection.from_wiski(
        r"./tests/data/2019-WISKI-test/1016_PBF.zip",
        translate_dic={
            'name': 'Station Number',
            'x': 'GlobalX',
            'y': 'GlobalY'},
        sep=r'\s+',
        header_sep=':',
        dayfirst=True,
        header_identifier=':',
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        verbose=True)

    return wiski_col


# %% PASTAS PROJECTS AND PASTASTORE
@pytest.mark.skip(reason="needs installation pastastore")
def test_to_pastas_project():

    dino_gw = test_obscollection_dinozip_gw()
    pr = dino_gw.to_pastas_project(verbose=True)

    return pr


@pytest.mark.skip(reason="needs installation pastastore")
def test_to_pastastore():

    dino_gw = test_obscollection_dinozip_gw()
    pstore = dino_gw.to_pastastore(verbose=True)

    return pstore


@pytest.mark.skip(reason="needs installation pastastore")
def test_from_pastas_project():

    pr = test_to_pastas_project()
    pr_oc = oc.ObsCollection.from_pastas_project(pr)

    return pr_oc


# %% PYSTORE

def test_obscollection_to_pystore():
    obsc = test_obscollection_fews_lowmemory()
    obsc.to_pystore("test_pystore", "./tests/data/2019-Pystore-test",
                    groupby="locatie", overwrite=True)


def test_obscollection_from_pystore():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore", "./tests/data/2019-Pystore-test")
    return obsc


def test_obscollection_pystore_only_metadata():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./tests/data/2019-Pystore-test",
                                         read_series=False)
    return obsc


def test_obscollection_pystore_extent():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./tests/data/2019-Pystore-test",
                                         extent=[115534, 115539, 0, 10000000]
                                         )
    return obsc


def test_obscollection_pystore_item_names():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./tests/data/2019-Pystore-test",
                                         item_names=['MPN-N-2']
                                         )
    return obsc


def test_obs_from_pystore_item():
    import pystore
    pystore.set_path("./tests/data/2019-Pystore-test")
    store = pystore.store("test_pystore")
    coll = store.collection(store.collections[0])
    item = coll.item(list(coll.list_items())[0])
    o = obs.GroundwaterObs.from_pystore_item(item)
    return o


# %% KNMI
def test_knmi_obs_from_stn():
    return obs.KnmiObs.from_knmi(233, "RD", verbose=True)

def test_knmi_obs_from_stn_no_api():
    return obs.KnmiObs.from_knmi(233, "RD", verbose=True,
                                 use_api=False)


def test_knmi_obs_from_stn_without_any_data():
    try:
        obs.KnmiObs.from_knmi(210, "RD",
                              verbose=True)
    except (ValueError, KeyError):
        pass

    return 1


def test_knmi_obs_from_stn_with_missing_data_in_time_period():
    return obs.KnmiObs.from_knmi(441, "RD", startdate='2010-1-2',
                                 verbose=True)


def test_knmi_obs_from_xy():
    return obs.KnmiObs.from_nearest_xy(100000, 350000, "RD")


def test_knmi_obs_from_obs():
    pb = test_observation_gw()
    return obs.KnmiObs.from_obs(pb, "EV24", fill_missing_obs=False)


def test_knmi_collection_from_locations():
    obsc = test_obscollection_dinozip_gw()
    oc_knmi = oc.ObsCollection.from_knmi(locations=obsc,
                                         meteo_vars=["EV24", "RD"],
                                         start=['2010', '2010'],
                                         end=['2015', '2015'],
                                         verbose=True, cache=False)
    return oc_knmi


def test_knmi_collection_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    oc_knmi = oc.ObsCollection.from_knmi(stns=stns,
                                         meteo_vars=["EV24", "RH"],
                                         start=['2010', '2010'],
                                         end=['2015', '2015'],
                                         verbose=True)
    return oc_knmi


def test_knmi_collection_from_grid():
    # somewhere in Noord-Holland (near Castricum)
    xmid = np.array([104150., 104550.])
    ymid = np.array([510150., 510550.])
    oc_knmi = oc.ObsCollection.from_knmi(xmid=xmid, ymid=ymid,
                                         meteo_vars=["RD"],
                                         start=['2010', '2010'],
                                         end=['2015', '2015'],
                                         verbose=True)
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
#    obsc = oc.ObsCollection.from_menyanthes(fname, verbose=True)
#
#    return obsc
