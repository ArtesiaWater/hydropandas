import os

import numpy as np
import pandas as pd

import pastas as ps

from .util import _import_art_tools


def get_stations(variable='RD'):
    """[summary]

    Parameters
    ----------
    variable : str, optional
        [description], by default 'RD'

    Returns
    -------
    [type]
        [description]
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))

    if variable == "RD":
        fname = "../data/knmi_neerslagstation.json"
    else:
        fname = "../data/knmi_meteostation.json"

    stations = pd.read_json(os.path.join(dir_path, fname))

    if variable == 'PG':
        # in Ell wordt geen luchtdruk gemeten
        stations.drop(377, inplace=True)
    if variable == 'EV24':
        # in Woensdrecht wordt geen verdamping gemeten
        stations.drop(340, inplace=True)

    return stations


def get_nearest_stations_xy(x, y, variable, n=1, stations=None, ignore=None):

    if stations is None:
        stations = get_stations(variable=variable)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
        if stations.empty:
            return None

    d = np.sqrt((stations.x - x)**2 + (stations.y - y)**2)

    return d.nsmallest(n).index.to_list()


def get_nearest_station_df(locations, stations=None, variable="RD", ignore=None):

    if stations is None:
        stations = get_stations(variable=variable)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
        if stations.empty:
            return None

    xo = pd.to_numeric(locations.x)
    xt = pd.to_numeric(stations.x)
    yo = pd.to_numeric(locations.y)
    yt = pd.to_numeric(stations.y)

    xh, xi = np.meshgrid(xt, xo)
    yh, yi = np.meshgrid(yt, yo)

    distances = pd.DataFrame(np.sqrt((xh - xi) ** 2 + (yh - yi) ** 2),
                             index=locations.index,
                             columns=stations.index)

    stns = distances.idxmin(axis=1).unique()
    if np.any(np.isnan(stns)):
        stns = stns[~np.isnan(stns)].astype(int)

    return stns


def _download_knmi_data(stn, var, start, end, verbose):

    # load art_tools module
    art = _import_art_tools()

    if start is None:
        if var == "EV24":
            start = "1957-7"
        else:
            start = "1900"

    if end is None:
        end = pd.datetime.today()

    # download data
    knmi = art.fill_missing_measurements(start=start, stns=stn, vars=var,
                                         end_date=end, verbose=verbose)
    return knmi


def get_knmi_timeseries_xy(x, y, var, start, end, normalize_index=True,
                           verbose=False):

    # get data
    stations = get_stations(variable=var)
    stn = get_nearest_stations_xy(x, y, var, stations=stations)

    # download data
    knmi = _download_knmi_data(stn, var, start, end, verbose)

    # set metadata
    name = var + ' ' + stations.loc[stn[0], 'naam']
    x = stations.loc[stn[0], 'x']
    y = stations.loc[stn[0], 'y']
    meta = {'x': x, 'y': y, 'station': str(stn[0]), 'name': name}

    v = knmi.data.loc[:, var]

    # normalize index
    if normalize_index:
        v.index = v.index.normalize()

    return v, meta


def get_knmi_timeseries_stn(stn, var, start, end, normalize_index=True,
                            verbose=False):
    # get data
    stations = get_stations(variable=var)

    # download data
    knmi = _download_knmi_data(stn, var, start, end, verbose)

    # set metadata
    name = var + ' ' + stations.loc[stn, 'naam']
    x = stations.loc[stn, 'x']
    y = stations.loc[stn, 'y']
    meta = {'x': x, 'y': y, 'station': str(stn), 'name': name}

    v = knmi.data.loc[:, var]

    # normalize index
    if normalize_index:
        v.index = v.index.normalize()

    return v, meta
