import json
import logging
from io import StringIO
from pathlib import Path
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from pyproj import Proj, Transformer
from shapely.geometry import Point, box
from tqdm import tqdm

from ..util import EPSG_28992

URL = "https://noos.matroos.rws.nl/direct/get_series.php?"

logger = logging.getLogger(__name__)


def get_obs_list_from_extent(
    ObsClass,
    extent=None,
    locations=None,
    sources=None,
    units=None,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=False,
    force_download=False,
    **kwargs,
):
    """Get observations within a specific extent.

    Parameters
    ----------
    ObsClass : type
        class of the observations, e.g. GroundwaterObs
    extent : list, tuple, numpy-array or None, optional
        get measurements within this extent
        [xmin, xmax, ymin, ymax]
    locations :
    sources :
    units :

    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : bool, optional
        if False, only observations with measurements are kept. The default
        is True.
    location_gdf : GeoDataFrame, optional
        geodataframe with the locations of the water drill holes you want to include.
    update : bool, optional
        if True new locations are downloaded and stored locally (slow) otherwise a
        cached version of the locations is used. By default False
    force_download : bool, optional
        if True an attempt is made to download data even if the params dic indicates
        there is no data for that combination of location, source and unit.
    **kwargs
        additional keyword arguments are passed to the ObsClass.from_waterconnect()
        method

    Returns
    -------
    obs_list : list
        list with Obs objects

    """

    if only_metadata and not keep_all_obs:
        logger.warning(
            "Option with only_metadata is True and keep_all_obs is False not supported "
            "setting keep_all_obs=True"
        )
        keep_all_obs = True

    params_dic = get_loc_source_units()

    # parse tmin and tmax
    if isinstance(tmin, str):
        tmin = pd.to_datetime(tmin)
    elif tmin is None:
        tmin = (pd.Timestamp.now() - pd.Timedelta(10, unit="D")).strftime("%Y%m%d%H%M")

    if isinstance(tmax, str):
        tmax = pd.to_datetime(tmax)
    elif tmax is None:
        tmax = pd.Timestamp.now()

    if extent is not None:
        if locations is None:
            locations = []
        locations += get_locations_within_extent(params_dic, extent)

    if len(locations) < 1:
        raise ValueError(f"no locations found within {extent=}")

    download_pars = select_by_params(
        locations=locations,
        sources=sources,
        units=units,
        force=force_download,
        params_dic=params_dic,
        keep_coords=False,
    )

    nobs = sum([len(s) for units in download_pars.values() for s in units.values()])
    logger.info(f"downloading {nobs} observations from {len(download_pars)} locations")
    logger.debug(f"download parameters {download_pars}")

    obs_list = []
    for location, units in tqdm(
        download_pars.items(), total=len(download_pars), desc="location"
    ):
        for unit, sources in units.items():
            for source in sources:
                logger.debug(
                    f"downloading matroos measurements for {location=}, {source=}, {unit=}, between {tmin=}, {tmax=}"
                )
                o = ObsClass.from_matroos(
                    location,
                    source,
                    unit,
                    tmin=tmin,
                    tmax=tmax,
                    only_metadata=only_metadata,
                    **kwargs,
                )
                if o.empty and not keep_all_obs:
                    logger.info(
                        f"no measurements found for {location=}, {source=}, {unit=}, between {tmin=}, {tmax=} not adding to obs_list"
                    )
                    continue

                obs_list.append(o)

    return obs_list


def _params_dic_to_df(params_dic):
    """convert params_dic to a pandas dataframe

    Parameters
    ----------
    params_dic : dict
        dictionary with locations, sources and units

    Returns
    -------
    pd.DataFrame
        dataframe with columns location, unit, source, x, y
    """
    df = pd.DataFrame(params_dic).T
    return df


def _params_dic_to_gdf(params_dic):
    """convert params_dic to a pandas dataframe

    Parameters
    ----------
    params_dic : dict
        dictionary with locations, sources and units

    Returns
    -------
    gpd.GeoDataFrame
        geodataframe with columns location, unit, source, x, y
    """
    df = _params_dic_to_df(params_dic)

    if "coords" not in df.columns:
        raise ValueError(
            "no coordinates found in params_dic, please use keep_coords=True in select_by_params function"
        )

    gdf = gpd.GeoDataFrame(
        df, geometry=df["coords"].apply(lambda x: Point(x[0], x[1])), crs="EPSG:28992"
    )
    return gdf


def select_by_params(
    locations=None,
    sources=None,
    units=None,
    force=False,
    params_dic=None,
    keep_coords=True,
    astype: Literal["dict", "dataframe", "geodataframe"] = "dict",
):
    # make sure locations, sources and units are lists, tuples or None
    if locations is not None and not isinstance(locations, (list, tuple)):
        locations = [locations]
    if sources is not None and not isinstance(sources, (list, tuple)):
        sources = [sources]
    if units is not None and not isinstance(units, (list, tuple)):
        units = [units]

    # get all locations, sources and units
    if params_dic is None:
        params_dic = get_loc_source_units()

    if locations is None:
        locations = params_dic.keys()

    if astype in ["geodataframe"]:
        keep_coords = True

    # create parameter dictionary
    selection_dic = {}
    for loc in locations:
        possible_units = params_dic[loc]["units"].keys()
        unit_dic = {}
        if units is None:
            dunits = possible_units
        elif force:
            dunits = units
        else:
            dunits = set(possible_units) & set(units)

        for unit in dunits:
            possible_sources = params_dic[loc]["units"][unit]
            if sources is None:
                dsources = possible_sources
            elif force:
                dsources = sources
            else:
                dsources = set(possible_sources) & set(sources)

            if dsources:
                unit_dic[unit] = dsources
        if unit_dic:
            selection_dic[loc] = unit_dic
            if keep_coords:
                selection_dic[loc]["coords"] = params_dic[loc]["coords"]

    if astype == "dict":
        return selection_dic

    df = pd.DataFrame(selection_dic).T
    if astype == "dataframe":
        return df
    elif astype == "geodataframe":
        gdf = gpd.GeoDataFrame(
            df,
            geometry=df["coords"].apply(lambda x: Point(x[0], x[1])),
            crs="EPSG:28992",
        )
        return gdf
    else:
        raise TypeError('astype should be one of "dict", "dataframe" or "geodataframe"')


def get_locations_within_extent(params_dic, extent):
    """get locations within a certain extent

    Parameters
    ----------
    params_dic : dict
        locations and coordinates
    extent :  list, tuple or numpy-array
        get locations within this extent
        [xmin, xmax, ymin, ymax]

    Returns
    -------
    list
        locations
    """

    pol_extent = box(*tuple(np.asarray(extent)[[0, 2, 1, 3]]))

    locations = [
        loc
        for loc, item in params_dic.items()
        if pol_extent.contains(Point(item["coords"]))
    ]

    return locations


def get_loc_source_units():
    """get the possible values for location, source and unit + the coordinates of each location

    Returns
    -------
    dictionary
    {location:
        {'units':
            {unit:
                [source_1, source_2, ..., source_n],
            },
        'coords':
            (x,y)
        }
    }
    """

    dir_path = Path(__file__).parent

    with open(dir_path / "../data/matroos_params.json", "r") as f:
        d = json.load(f)
    return d


def request_api(location, source, unit, tmin, tmax, timeout=600, fname=None):
    """make a request to the Matroos api for a specific location, source and unit

    Parameters
    ----------
    location : str
        location e.g. 'krimpen a/d lek'
    source : str
        source e.g. 'observed'
    unit : str
        unit e.g. 'waterlevel'
    tmin : str
        start of requested time series, format '%Y%m%d%H%M' e.g. '202505100110'
    tmax : str
        end of requested time series, format '%Y%m%d%H%M' e.g. '202505100110'
    timeout : int, optional
        timeout for request in seconds, by default 600
    fname : str, optional
        filename to save result of request in text format, by default None

    Returns
    -------
    StringIO
    """
    params = {
        "loc": location,
        "source": source,
        "unit": unit,
        "tstart": tmin,
        "tstop": tmax,
    }

    r = requests.get(URL, params=params, timeout=timeout)

    r.raise_for_status()

    result_str = r.text

    if fname is not None:
        with open(fname, "w") as f:
            f.write(result_str)

    return StringIO(result_str)


def get_matroos_obs(
    location, source, unit, tmin=None, tmax=None, only_metadata=False, **kwargs
):
    """get observations for a certain location, source and unit between
    tmin and tmax.

    Parameters
    ----------
    location : str
        location e.g. 'krimpen a/d lek'
    source : str
        source e.g. 'observed'
    unit : str
        unit e.g. 'waterlevel'
    tmin : pd.Timestamp, str or None, optional
        start of time series if None tmin is 10 days ago, by default None
    tmax : pd.Timestamp, str or None, optional
        start of time series if None tmin is today, by default None
    only_metadata : bool, optional
        if True download only metadata, slightly faster. The default
        is False.
    **kwargs are passed to request_api function

    Returns
    -------
    pd.DataFrame: measurements

    dict: metadata
    """
    # check if location, source and unit are valid
    params_dic = get_loc_source_units()
    if location not in params_dic.keys():
        msg = f"{location} not listed in possible locations, please select location from {params_dic.keys()}"
        logger.warning(msg)

    elif unit not in params_dic[location]["units"].keys():
        msg = f"{unit} not available for {location=}, please select unit from {params_dic[location]['units'].keys()}"
        logger.warning(msg)

    elif source not in params_dic[location]["units"][unit]:
        msg = f"{source} not available for {location=} and {unit=}, please select source from {params_dic[location]['units'][unit]}"
        logger.warning(msg)

    # parse tmin and tmax
    if isinstance(tmin, str):
        tmin = pd.to_datetime(tmin)
    elif tmin is None:
        tmin = (pd.Timestamp.now() - pd.Timedelta(10, unit="D")).strftime("%Y%m%d%H%M")
    if isinstance(tmin, pd.Timestamp):
        tmin = tmin.strftime("%Y%m%d%H%M")

    if isinstance(tmax, str):
        tmax = pd.to_datetime(tmax)
    elif tmax is None:
        tmax = pd.Timestamp.now()
    if isinstance(tmax, pd.Timestamp):
        tmax = tmax.strftime("%Y%m%d%H%M")

    # request api
    f = request_api(location, source, unit, tmin, tmax, **kwargs)

    # read first 4 lines:
    breakline = f.readline()
    line = f.readline()
    while line != breakline:
        line = f.readline()

    # read metadata
    line = f.readline()

    meta = {
        "location": location,
        "source": "Matroos",
        "name": f"{location.replace('/', '').replace(' ', '_')}_{source}_{unit}",
    }
    while line != breakline:
        key, item = line.strip("#").split(":")
        if "Position" in key:
            lon, lat = (float(a) for a in item.strip()[1:-1].split(","))
            proj_from = Proj("EPSG:4326")
            proj_to = Proj(EPSG_28992)
            transformer = Transformer.from_proj(proj_from, proj_to)
            xy = transformer.transform(lat, lon)
            meta["x"] = xy[0]
            meta["y"] = xy[1]
        elif "Analyse time" in key:
            if "*** no data found ***" in key:
                msg = f"no measurement data found for {location=}, {source=}, {unit=}, {tmax=}, {tmin=}. Only returning metadata"
                logger.warning(msg)
                only_metadata = True
        elif ("Location" in key) or ("Source" in key) or ("Unit" in key):
            pass
        else:
            meta[key.strip()] = item.strip()
        line = f.readline()

    if only_metadata:
        return pd.DataFrame(), meta

    # read measurements
    df = pd.read_csv(
        f,
        header=None,
        sep=r"\s+",
        index_col=0,
        names=["datetime", "waterlevel"],
        parse_dates=[0],
    )

    return df, meta
