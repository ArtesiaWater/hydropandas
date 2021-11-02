# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 09:24:29 2019

@author: oebbe
"""

import numpy as np
import pandas as pd
import pytest
from hydropandas import observation as obs
from hydropandas.io import io_knmi

import logging
logging.basicConfig(level=logging.DEBUG)


def test_get_knmi_precip_neerslagstation():
    # De Bilt neerslagstation
    ts, meta = io_knmi.get_knmi_timeseries_stn(
        '550_neerslag_station',
        "RH",
        start='2010-1-1',
        end='2010-1-10',
        settings=None)
    return ts, meta

def test_get_knmi_precip_meteostation_fill_missing():
    
    # De Bilt meteostation
    ts2, meta2 = io_knmi.get_knmi_timeseries_stn(
        260,
        "RH",
        start='2010-1-1',
        end='2010-1-10',
        settings=None)
    return ts2, meta2

@pytest.mark.skip(reason="KNMI API is down")
def test_get_knmi_precip_meteostation_hourly():
    
    # De Bilt meteostation uurlijks
    ts3, meta3 = io_knmi.get_knmi_timeseries_stn(
        260,
        "RH",
        start='2010-1-1',
        end='2010-1-10',
        settings={'interval':'hourly',
                  'fill_missing_obs':False})
    return ts3, meta3

def test_get_knmi_precip_neerslagstation_no_api():
    # De Bilt neerslagstation
    ts4, meta4 = io_knmi.get_knmi_timeseries_stn(
        '550_neerslag_station',
        "RH",
        start='2010-1-1',
        end='2010-1-10',
        settings={'use_api': False})
    return ts4, meta4


def test_get_knmi_precip_meteostation_no_api():
    # De Bilt neerslagstation
    ts5, meta5 = io_knmi.get_knmi_timeseries_stn(
        260,
        "RH",
        start='2010-1-1',
        end='2010-1-10',
        settings={'use_api': False})
    return ts5, meta5


def test_get_knmi_precip_neerslagstation_fill_missing():
    settings = io_knmi._get_default_settings(None)
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        '550_neerslag_station',
        meteo_var='RH',
        start='1952',
        end=None,
        settings=settings)
    return knmi_df, variables, stations


def test_download_rd_550_no_api():
    settings = io_knmi._get_default_settings(None)
    settings['use_api'] = False
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        '550_neerslag_station',
        meteo_var='RH',
        start='1952',
        end=None,
        settings=settings)
    return knmi_df, variables, stations


def test_download_rd_12():
    settings = io_knmi._get_default_settings(None)
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        '12_neerslag_station',
        meteo_var='RH',
        start='2010',
        end=None,
        settings=settings)
    return knmi_df, variables, stations


def test_download_without_data():
    settings = io_knmi._get_default_settings(None)
    try:
        knmi_df, variables, stations = io_knmi.download_knmi_data(
            '324_neerslag_station',
            meteo_var='RH',
            start='2018',
            end='2020',
            settings=settings)
    except ValueError:
        pass

    return 1


def test_download_without_data_no_error():
    settings = io_knmi._get_default_settings(settings={'raise_exceptions':False})
    
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        '324_neerslag_station',
        meteo_var='RH',
        start='2018',
        end='2020',
        settings=settings)

    assert (knmi_df.empty, variables == {},
            stations.empty) == (True, True, True)

    return knmi_df, variables, stations


def test_download_ev24_210():
    settings = io_knmi._get_default_settings(None)
    
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        210,
        meteo_var='EV24',
        start='1952',
        end=None,
        settings=settings)
    return knmi_df, variables, stations


def test_download_ev24_210_no_api():
    settings = io_knmi._get_default_settings({'use_api':False})
    knmi_df, variables, stations = io_knmi.download_knmi_data(
        210,
        meteo_var='EV24',
        start='1952',
        end=None,
        settings=settings)
    return knmi_df, variables, stations

@pytest.mark.skip(reason="KNMI API is down")
def test_get_knmi_daily_meteo_ev24_empty():
    start, end = io_knmi._start_end_to_datetime('1959', '1963')
    knmi_df, variables, stations = io_knmi.get_knmi_daily_meteo_api(
                    io_knmi.URL_DAILY_METEO, 
                    265, 'EV24', start, end, False)
            
    return knmi_df, variables, stations



def test_fill_missing_measurements_ev24_278():
    knmi_df, variables, stations = io_knmi.fill_missing_measurements(
        278,
        meteo_var='EV24',
        start='1959',
        end='1963')
    return knmi_df, variables, stations


def test_fill_missing_measurements_rh_273():
    knmi_df, variables, stations = io_knmi.fill_missing_measurements(
        273,
        meteo_var='RH',
        start='1995',
        end='2000')
    return knmi_df, variables, stations

def test_fill_missing_measurements_rh_278():
    settings = io_knmi._get_default_settings({'exclude_neerslag_stn':True})
    knmi_df, variables, stations = io_knmi.fill_missing_measurements(
        278,
        meteo_var='RH',
        start='1990',
        end=None,
        settings=settings)
    return knmi_df, variables, stations


def test_obslist_from_grid():
    xmid = np.array([104150., 104550.])
    ymid = np.array([510150., 510550.])
    obs_list = io_knmi.get_knmi_obslist(xmid=xmid, ymid=ymid,
                                        meteo_vars=['RH', 'EV24'],
                                        start=['2010', '2010'],
                                        end=[None, None],
                                        ObsClass=[obs.PrecipitationObs, 
                                                  obs.EvaporationObs])

    return obs_list


def test_obslist_from_locations():
    locations = pd.DataFrame(data={'x': [100000], 'y': [350000]})
    obs_list = io_knmi.get_knmi_obslist(locations, meteo_vars=['EV24'],
                                        start=[None, None],
                                        ObsClass=[obs.EvaporationObs],
                                        end=[None, None])

    return obs_list


def test_obslist_from_stns():
    stns = [344, 260]  # Rotterdam en de Bilt
    obs_list = io_knmi.get_knmi_obslist(stns=stns, meteo_vars=['RH', 'EV24'],
                                        start=['2010', '2010'],
                                        end=['2015', '2015'],
                                        ObsClass=[obs.PrecipitationObs, 
                                                  obs.EvaporationObs])

    return obs_list

def test_obslist_from_stns_single_startdate():
    stns = [344, 260]  # Rotterdam en de Bilt
    obs_list = io_knmi.get_knmi_obslist(stns=stns, meteo_vars=['RH', 'EV24'],
                                        start='2010',
                                        end='2015',
                                        ObsClass=[obs.PrecipitationObs, 
                                                  obs.EvaporationObs])

    return obs_list
