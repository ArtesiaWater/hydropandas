"""
Module with functions to read or download time series with observations from knmi.

function levels:
1. get_knmi_obs_list: list of observations
    2. get_knmi_obs: single observation
        3. get_timeseries_stn, get_timeseries_from_file: get time series from station or file
            4. fill_missing_measurements, get_evaporation: aggregate time series data
                5. download_knmi_data: download a single timeseries
                    6a. get_hourly_meteo_api, get_daily_meteo_api, get_daily_rainfall_api,
                        get_daily_meteo_url, get_daily_rainfall_url
                        7a. request_api, request_url
                        7b. parse_data
                    6b. interpret_knmi_file
"""

import datetime as dt
import logging
import os
import warnings
from functools import lru_cache
from io import BytesIO, StringIO
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
from zipfile import ZipFile

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

URL_DAILY_PREC = "https://www.daggegevens.knmi.nl/klimatologie/monv/reeksen"
URL_DAILY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/daggegevens"
URL_HOURLY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/uurgegevens"


def get_knmi_obs(
    stn: Union[int, None] = None,
    fname: Union[str, None] = None,
    xy: Union[List[float], Tuple[float], None] = None,
    meteo_var: Union[str, None] = None,
    start: Union[pd.Timestamp, str, None] = None,
    end: Union[pd.Timestamp, str, None] = None,
    **kwargs,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """get knmi observation from stn, fname or nearest xy coordinates.

    Parameters
    ----------
    stn : int, str or None, optional
        measurement station e.g. 829. The default is None.
    fname : str, path object, file-like object or None, optional
        filename of a knmi file. The default is None.
    xy : list, tuple or None, optional
        RD coördinates of a location in the Netherlands. The station nearest
        to this location used. The Default is None.
    meteo_var : str or None, optional
        meteo variable e.g. "RH" or "EV24". See list with all options in the
        hydropandas documentation.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    **kwargs:
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False. Note: if the given stn has no data between start and
            end the data from nearby stations is used. In this case the metadata of the
            Observation is the metadata from the nearest station that has any
            measurement in the given period.
        interval : str, optional
            desired time interval for observations. Options are 'daily' and
            'hourly'. The default is 'daily'.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary folder and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

    Returns
    -------
    pd.DataFrame, measurements
    dict, metadata

    Raises
    ------
    ValueError
        if no meteo_var is given or stn, fname and xy are all None.
    """
    if meteo_var is None and fname is None:
        raise ValueError("To get knmi data a meteo_var should be specified")
    elif meteo_var is not None:
        if not isinstance(meteo_var, str):
            raise (TypeError(f"meteo var should be string not {type(meteo_var)}"))

    settings = _get_default_settings(kwargs)

    start = start if start is None else pd.to_datetime(start)
    end = end if end is None else pd.to_datetime(end)

    if (stn in (913, 967)) and (meteo_var == "RD") and settings["use_api"]:
        msg = (
            f"precipitation data not available for station {stn} via the API. "
            "setting use_api to False, more info here: "
            "https://github.com/ArtesiaWater/hydropandas/issues/245"
        )
        logger.warning(msg)
        settings["use_api"] = False

    if stn is not None:
        stn = int(stn)

        start_str = str(start).replace(" 00:00:00", "")
        end_str = str(end).replace(" 00:00:00", "")

        logger.info(
            f"get data from station {stn} and variable {meteo_var} "
            f"from {start_str} to {end_str}"
        )
        ts, meta = get_timeseries_stn(
            stn=stn, meteo_var=meteo_var, settings=settings, start=start, end=end
        )
    elif fname is not None:
        logger.info(f"get KNMI data from file {fname} and meteo variable {meteo_var}")
        ts, meta = get_timeseries_from_file(
            fname=str(fname),
            meteo_var=meteo_var,
            settings=settings,
            start=start,
            end=end,
        )
    elif xy is not None:
        logger.info(
            f"get KNMI data from station nearest to coordinates {xy} and meteo"
            f"variable {meteo_var}"
        )
        stns = get_n_nearest_stations_xy(
            xy=xy,
            meteo_var=meteo_var,
            start=start,
            end=end,
            n=1,
            stations=None,
            ignore=None,
        )
        ts, meta = get_timeseries_stn(
            stn=stns[0], meteo_var=meteo_var, settings=settings, start=start, end=end
        )
    else:
        raise ValueError(
            "specify KNMI station (stn), filename (fname) or coordinates (xy)"
        )

    return ts, meta


def get_knmi_timeseries_fname(
    fname: str,
    meteo_var: str,
    settings: Dict[str, Any],
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Get a knmi time series and metadata from a file.

    .. deprecated:: 0.13.3
        `get_knmi_timeseries_fname` will be removed in hydropandas 1.0.0, it is replaced
        by `get_timeseries_from_file`.

    Parameters
    ----------
    fname : str
        filename of the knmi file.
    meteo_var : str
        observation type e.g. "RH" or "EV24". See list with all options in the
        hpd.read_knmi function.
    settings : dict
        settings for obtaining the right time series, see _get_default_settings
        for more information
    start : pd.Timestamp
        start date of observations.
    end : pd.Timestamp
        end date of observations.

    Returns
    -------
    ts_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """
    warnings.warn(
        "the function 'get_knmi_timeseries_fname' is deprecated and will eventually be "
        "removed, please use 'get_timeseries_from_file'.",
        DeprecationWarning,
    )

    return get_timeseries_from_file(fname, meteo_var, settings, start, end)


def get_timeseries_from_file(
    fname: str,
    meteo_var: str,
    settings: Dict[str, Any],
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Get a knmi time series and metadata from a file.

    Parameters
    ----------
    fname : str
        filename of the knmi file.
    meteo_var : str
        observation type e.g. "RH" or "EV24". See list with all options in the
        hpd.read_knmi function.
    settings : dict
        settings for obtaining the right time series, see _get_default_settings
        for more information
    start : pd.Timestamp
        start date of observations.
    end : pd.Timestamp
        end date of observations.

    Returns
    -------
    ts_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """

    df, meta = parse_data(fname)
    # first try to figure out filetype by it's name
    if "neerslaggeg" in fname:
        # neerslagstation
        meteo_var = "RD"
        add_day = False
    elif "etmgeg" in fname:
        # meteo station
        add_day = True
    elif "uurgeg" in fname:
        add_day = False
        # hourly station
    # if that doesn't work try to figure out by the meteo_var and settings
    elif meteo_var is None or meteo_var == "RD":
        # neerslagstation
        meteo_var = "RD"
        add_day = False
    elif settings["interval"] == "daily":
        # meteo station
        add_day = True
    elif settings["interval"] == "hourly":
        # uurlijks station
        add_day = False
    else:
        raise ValueError(
            "please indicate how to read the file by specifying a meteo_var and"
            " an interval"
        )
    if df.empty:
        logger.warning(
            f"No data for {meteo_var=} in {fname=} between{start=} and {end=}."
        )
    else:
        ts, meta = interpret_knmi_file(
            df=df,
            meta=meta,
            meteo_var=meteo_var,
            start=start,
            end=end,
            add_day=add_day,
            add_hour=True,
        )

    stn = meta["station"]
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn=stn, stations=stations)

    # set metadata
    meta.update(
        {
            "x": stations.loc[stn, "x"],
            "y": stations.loc[stn, "y"],
            "name": f"{meteo_var}_{stn_name}_{stn}",
            "location": stn_name,
            "source": "KNMI",
            "filename": fname,
        }
    )

    return ts, meta


def _get_default_settings(settings=None) -> Dict[str, Any]:
    """adds the default settings to a dictinary with settings. If settings
    is None all the settings are default. If there are already settings given
    only the non-existing settings are added with their default value.

    The default settings are:
    fill_missing_obs = False
        nan values in time series are filled with nearby time series.
    interval = 'daily'
        desired time interval for observations. Can be 'daily' or 'hourly'.
        'hourly' is only for precipitation ('RH') data from meteo stations.
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
        https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        if False a text file is downloaded into a temporary folder and the data is read
        from that file
        Default is True since the api is back online (July 2021).
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is False.

    Parameters
    ----------
    settings : dict or None
        settings dictionary without default values.

    Returns
    -------
    settings : dict
        settings dictionary with default values of no values was given in
        the input dictionary.
    """

    default_settings = {
        "fill_missing_obs": False,
        "interval": "daily",
        "use_api": True,
        "raise_exceptions": True,
    }

    if settings is None:
        settings = {}

    if "fill_missing_obs" in settings.keys():
        if "raise_exceptions" in settings.keys():
            if settings["fill_missing_obs"] and settings["raise_exceptions"]:
                logger.debug(
                    "set raise_exceptions=False because fill_missing_obs is True"
                )
                settings["raise_exceptions"] = False
        else:
            settings["raise_exceptions"] = False

    for key, value in default_settings.items():
        if key not in settings.keys():
            settings[key] = value

    return settings


def get_knmi_timeseries_stn(
    stn: int,
    meteo_var: str,
    settings: Dict[str, Any],
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Get a knmi time series and metadata.

    .. deprecated:: 0.13.3
        `get_knmi_timeseries_stn` will be removed in hydropandas 1.0.0, it is replaced
        by `get_timeseries_stn`.

    Parameters
    ----------
    stn : int
        measurement station e.g. 829.
    meteo_var : str
        observation type e.g. "RH" or "EV24". See list with all options in the
        hpd.read_knmi function.
    settings : dict
        settings for obtaining the right time series, see _get_default_settings
        for more information
    start : pd.TimeStamp or None, optional
        start date of observations. The default is None.
    end : pd.TimeStamp or None, optional
        end date of observations. The default is None.

    Returns
    -------
    ts_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """
    warnings.warn(
        "the function 'get_knmi_timeseries_stn' is deprecated and will eventually be "
        "removed, please use 'get_timeseries_stn'.",
        DeprecationWarning,
    )

    return get_timeseries_stn(stn, meteo_var, settings, start, end)


def get_timeseries_stn(
    stn: int,
    meteo_var: str,
    settings: Dict[str, Any],
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Get a knmi time series and metadata.

    Parameters
    ----------
    stn : int
        measurement station e.g. 829.
    meteo_var : str
        observation type e.g. "RH" or "EV24". See list with all options in the
        hpd.read_knmi function.
    settings : dict
        settings for obtaining the right time series, see _get_default_settings
        for more information
    start : pd.TimeStamp or None, optional
        start date of observations. The default is None.
    end : pd.TimeStamp or None, optional
        end date of observations. The default is None.

    Returns
    -------
    ts_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn=stn, stations=stations)

    # raise error if hourly neerslag station data is requested
    if (meteo_var == "RD") and settings["interval"].startswith("hour"):
        message = (
            "No hourly precipitation data available at precipitation station, "
            "use meteo_var 'RH' to obtain hourly precipitation data from "
            "meteo stations."
        )
        raise ValueError(message)
    elif (
        (meteo_var == "EV24")
        and settings["interval"].startswith("hour")
        and settings["use_api"]
    ):
        message = (
            "No hourly evaporation data available through the api, set use_api=False."
        )
        raise ValueError(message)
    elif settings["fill_missing_obs"]:
        # download data
        ts_df, meta = fill_missing_measurements(
            stn=stn,
            meteo_var=meteo_var,
            start=start,
            end=end,
            settings=settings,
            stn_name=stn_name,
        )

        return ts_df, meta
    elif meteo_var in ["penman", "makkink", "hargreaves"]:
        # compute evaporation from data
        ts_df, meta = get_evaporation(
            meteo_var=meteo_var, stn=stn, start=start, end=end, settings=settings
        )
    else:
        ts_df, variables, station_meta = download_knmi_data(
            stn=stn,
            meteo_var=meteo_var,
            start=start,
            end=end,
            settings=settings,
            stn_name=stn_name,
        )
        if ts_df.empty:
            logger.warning(
                f"No data for {meteo_var=} at {stn=} between{start=} and {end=}."
            )
        if str(stn) in station_meta.index:
            meta = station_meta.loc[f"{stn}"].to_dict()
        else:
            meta = {}

        # set metadata
        x = stations.at[stn, "x"] if stn in stations.index else np.nan
        y = stations.at[stn, "y"] if stn in stations.index else np.nan
        meta.update(
            {
                "x": x,
                "y": y,
                "station": stn,
                "name": f"{meteo_var}_{stn_name}_{stn}",
                "location": stn_name,
                "source": "KNMI",
            }
        )
        meta.update(variables)

    return ts_df, meta


def get_stations(
    meteo_var: str,
    start: Union[pd.Timestamp, str, None] = None,
    end: Union[pd.Timestamp, str, None] = None,
) -> pd.DataFrame:
    """get knmi stations from json files according to variable.

    Parameters
    ----------
    meteo_var : str, optional
        type of meteodata, by default 'RH'
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.

    Returns
    -------
    pandas DataFrame with stations, names and coordinates (Lat/Lon & RD)
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))

    mstations = pd.read_json(os.path.join(dir_path, "../data/knmi_meteostation.json"))
    pstations = pd.read_json(
        os.path.join(dir_path, "../data/knmi_neerslagstation.json")
    )

    stations = pd.concat([mstations, pstations], axis=0)
    stations = stations.where(~stations.isna(), False)
    if meteo_var in ("makkink", "penman", "hargreaves"):
        meteo_var = "EV24"

    # select only stations with meteo_var
    if meteo_var == slice(None):
        meteo_mask = stations.loc[:, meteo_var].any(axis=1)
    else:
        meteo_mask = stations.loc[:, meteo_var]
    stations = stations.loc[
        meteo_mask, ["lon", "lat", "name", "x", "y", "altitude", "tmin", "tmax"]
    ]

    # select only stations with measurement
    if start is not None or end is not None:
        stations = _get_stations_tmin_tmax(stations, start, end)

    return stations


def _get_stations_tmin_tmax(stations_df, start, end):
    """select stations within period defined by start and end.

    Parameters
    ----------
    stations_df : pd.DataFrame
        stations
    start : datetime or None, optional
        start date of observations.
    end : datetime or None, optional
        end date of observations.

    Returns
    -------
    DataFrame with all station with measurements in selected period.

    Notes
    -----
    Does not work on a DataFrames with duplicate indices
    """
    if stations_df.index.duplicated().any():
        raise IndexError("function does not work for dataframe with duplicated index")

    if end is None:
        tmin_stns = set(stations_df.index)
    else:
        # keep stations where tmin is unknown (=False)
        stns_unknown_tmin = set(stations_df.loc[stations_df["tmin"] == False].index)
        tmin_available = stations_df.loc[stations_df["tmin"] != False, "tmin"]
        tmin_within_range = pd.to_datetime(tmin_available) < end
        tmin_stns = set(tmin_available.loc[tmin_within_range].index) | stns_unknown_tmin

    if start is None:
        tmax_stns = set(stations_df.index)
    else:
        stns_unknown_tmax = set(stations_df.loc[stations_df["tmax"] == False].index)
        tmax_available = stations_df.loc[stations_df["tmax"] != False, "tmax"]
        tmax_available.loc[tmax_available.isnull()] = dt.datetime.now().date()
        tmax_within_range = pd.to_datetime(tmax_available) > start
        tmax_stns = set(tmax_available.loc[tmax_within_range].index) | stns_unknown_tmax

    return stations_df.loc[list(tmin_stns & tmax_stns)]


def get_station_name(stn: int, stations: Union[pd.DataFrame, None] = None) -> str:
    """Returns the station name from a KNMI station.

    Modifies the station name in such a way that a valid url can be obtained.

    Parameters
    ----------
    stn : int
        station number.
    stations : pandas DataFrame or None
        DataFrame with the station metadata.

    Returns
    -------
    str or None
        Name of the station or None if station is not found.
    """
    if stations is None:
        stations = get_stations(meteo_var=slice(None))

    if stn not in stations.index:
        logger.warning(f"station {stn} not found")
        return None

    stn_name = stations.at[stn, "name"]
    if isinstance(stn_name, pd.Series):
        raise ValueError(
            f'station {stn} is a meteo- and a precipitation station, please indicate which one you want to use using a "meteo_var"'
        )

    stn_name = stn_name.upper().replace(" ", "-").replace("(", "").replace(")", "")
    return stn_name


def fill_missing_measurements(
    stn: int,
    meteo_var: str,
    start: pd.Timestamp,
    end: pd.Timestamp,
    settings: Dict[str, Any],
    stn_name: Union[str, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any], pd.DataFrame]:
    """fill missing measurements in knmi data.

    Parameters
    ----------
    stn : int
        measurement station.
    meteo_var : str
        observation type.
    start : pd.TimeStamp
        start date of observations.
    end : pd.TimeStamp
        end date of observations.
    settings : dict
        settings for obtaining data.
    stn_name : str, optional
        station name. If None the name is obtained form te stn number. The
        default is None.

    Returns
    -------
    knmi_df : pandas DataFrame
        data from one station from one type of observation, with additional
        column to see which station is used to fill the value
    meta : dictionary
        metadata from the originally requested station even if this station
        has no data
    """
    if settings["interval"] == "hourly":
        raise NotImplementedError("cannot yet fill missing values in hourly data")

    if (start is None) or (end is None):
        raise ValueError(
            "a start and end date have to be specified if fill_missing_obs is True"
        )

    # 1. get stations
    stations = get_stations(meteo_var=meteo_var)
    if stn not in stations.index:
        logger.error(f"station {stn} does not exists or does not measure {meteo_var}")

    stations_period = get_stations(meteo_var=meteo_var, start=start, end=end)
    if stn_name is None:
        stn_name = get_station_name(stn=stn, stations=stations)

    meta = stations.loc[stn].to_dict()
    meta.update(
        {
            "station": stn,
            "name": f"{meteo_var}_{stn_name}_{stn}",
            "location": stn_name,
            "source": "KNMI",
        }
    )

    # 2. Download the data from the given station
    if stn not in stations_period.index:
        # no measurements in given period, continue with empty dataframe
        ts_df = pd.DataFrame()

        # add location of station without data to dataframe
        stations_period = pd.concat([stations_period, stations.loc[[stn]]])
    else:
        # download data from station if it has data between start-end
        ts_df, variables, _station_meta = download_knmi_data(
            stn, meteo_var, start, end, settings, stn_name
        )

    # 3. Change start date
    # NOTE: Assuming there is no data before the 'tmin' dates in the json files.
    tmin_first_meas = pd.Timestamp(stations["tmin"].min())
    if end < tmin_first_meas:
        msg = (
            f"There are no measurements available between {start} and {end} "
            f"for variable {meteo_var}."
        )
        logger.warning(msg)

        return ts_df, meta

    elif start < tmin_first_meas:
        msg = (
            f"There are no measurements available for {meteo_var} before "
            f"{tmin_first_meas} and a start of {start} was requested. Changing "
            f"start to {tmin_first_meas}"
        )
        logger.info(msg)
        start = tmin_first_meas

    # 4. Change end date
    # NOTE: Assuming there is no data after the last measurement available in de Bilt.
    stn_de_bilt = 550 if meteo_var == "RD" else 260
    first_meas_de_bilt = pd.Timestamp(stations.loc[stn_de_bilt, "tmin"])

    # only change end if dataframe does not have measurements at the end date
    if end < first_meas_de_bilt:
        msg = (
            f"knmi station De Bilt has no measurements before {first_meas_de_bilt}"
            "no need to change the end date"
        )
        logger.debug(msg)
    elif ts_df.empty or (end > ts_df.index[-1]):
        # get measurements at de Bilt
        start_de_bilt = start if ts_df.empty else ts_df.index[-1]
        ts_de_bilt, _, _ = download_knmi_data(
            stn_de_bilt,
            meteo_var,
            start=start_de_bilt,
            end=end,
            settings=settings,
            stn_name="De Bilt",
        )

        if ts_de_bilt.empty:
            msg = (
                f"knmi station De Bilt has no measurements between {start_de_bilt} and"
                f" {end} for variable {meteo_var}. It is assumed no data is available"
                f" for this period"
            )
            logger.warning(msg)
            return ts_df, meta

        new_end = ts_de_bilt.index[-1]
        if new_end < end:
            msg = (
                f"knmi station De Bilt has no measurements for {meteo_var} after "
                f"{new_end} and an end date of {end} was requested. Changing end to "
                f"{new_end}"
            )
            logger.info(msg)
            end = new_end

    # 5. get starting dataframe with measurements
    # if there are no measurements that fit the requirements at the given stn
    # read the nearby station
    ignore = [stn]
    while ts_df.empty:
        stn_lst = get_nearest_station_df(
            stations_period.loc[[ignore[0]]],
            meteo_var=meteo_var,
            start=start,
            end=end,
            stations=stations_period,
            ignore=ignore,
        )
        if stn_lst is None:
            logger.error(
                f"there is no station with measurements of {meteo_var} between "
                f"{start} and {end}"
            )
            return pd.DataFrame(), meta

        msg = (
            f"station {stn} has no measurements between {start} and {end} "
            f"trying to get measurements from nearest station -> {stn_lst[0]}"
        )
        logger.info(msg)

        stn = stn_lst[0]
        stn_name = get_station_name(stn=stn, stations=stations_period)
        ts_df, variables, _station_meta = download_knmi_data(
            stn, meteo_var, start, end, settings, stn_name
        )
        # do not use station number for metadata
        if "station" in variables:
            variables.pop("station")

        # ignore this station
        ignore.append(stn)

    # 6. find and fill missing values
    ts_df = _add_missing_indices(ts_df, stn, start, end)
    missing = ts_df[meteo_var].isna()
    logger.debug(f"station {stn} has {missing.sum()} missing measurements")
    ts_df.loc[~missing, "station"] = str(stn)

    while np.any(missing) and not np.all(missing):
        stn_comp = get_nearest_station_df(
            stations_period.loc[[stn]],
            meteo_var=meteo_var,
            start=start,
            end=end,
            stations=stations_period,
            ignore=ignore,
        )

        if stn_comp is None:
            logger.warning(
                "Could not fill all missing measurements as there are "
                "no stations left to check!"
            )
            missing[:] = False
            break
        else:
            stn_comp = stn_comp[0]

        stn_name_comp = get_station_name(stn_comp, stations_period)
        logger.debug(
            f"Trying to fill {missing.sum()} missing measurements with "
            f"station {stn_comp} {stn_name_comp}"
        )

        ts_df_comp, _, __ = download_knmi_data(
            stn_comp, meteo_var, start, end, settings, stn_name_comp
        )

        if ts_df_comp.empty:
            logger.debug(f"No data available for station {stn_comp}")
        else:
            # dropnans from new data
            ts_df_comp = ts_df_comp.loc[~ts_df_comp[meteo_var].isna(), :]
            # get index of missing data in original timeseries
            missing_idx = missing.loc[missing].index
            # if any missing are in the new data, update
            if missing_idx.isin(ts_df_comp.index).any():
                # index for missing but in newly downloaded data
                ix_idx = missing_idx.intersection(ts_df_comp.index)
                # update missing data
                ts_df.loc[ix_idx, meteo_var] = ts_df_comp.loc[ix_idx, meteo_var]
                # add source station number
                ts_df.loc[ix_idx, "station"] = str(stn_comp)
                logger.info(
                    f"Filled {ix_idx.size} observations from station {stn_comp} "
                    f"{stn_name_comp} -> {stn} {stn_name}"
                )
            else:
                logger.debug(
                    f"No new data available from {stn_comp} {stn_name_comp} "
                    f"for filling missing measurements"
                )
        missing = ts_df[meteo_var].isna()
        ignore.append(stn_comp)

    meta.update(variables)

    return ts_df, meta


def download_knmi_data(
    stn: int,
    meteo_var: str,
    start: pd.Timestamp,
    end: pd.Timestamp,
    settings: Dict[str, Any],
    stn_name: Union[str, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any], pd.DataFrame]:
    """download knmi data of a measurements station for certain observation
    type.

    Parameters
    ----------
    stn : int
        measurement station.
    meteo_var : str
        observation type.
    start : pd.TimeStamp
        start date of observations.
    end : pd.TimeStamp
        end date of observations.
    settings : dict
        settings for obtaining data
    stn_name : str, optional
        station name. If None the name is obtained form te stn number. The
        default is None.

    Raises
    ------
    ValueError
        if the data from knmi cannot not be read a ValueError is raised.
        Unless raise_exceptions is False

    Returns
    -------
    ts_df : pandas DataFrame
        data from one station from one type of observation
    variables : dictionary
        information about the observed variables
    stations : pandas DataFrame
        information about the measurement station
    """
    msg = f"{stn}-{stn_name} between {start} and {end}"
    logger.debug(f"download KNMI {meteo_var} data from station " + msg)

    # define variables
    ts_df = pd.DataFrame()
    variables = {}

    # download and read data
    try:
        try:
            if settings["use_api"]:
                if settings["interval"].startswith("hour"):
                    # hourly data from meteorological stations
                    df, meta = get_hourly_meteo_api(
                        stn=stn, meteo_var=meteo_var, start=start, end=end
                    )
                    add_day = False
                elif meteo_var == "RD":
                    # daily data from rainfall-stations
                    df, meta = get_daily_rainfall_api(stn=stn, start=start, end=end)
                    add_day = False
                else:
                    # daily data from meteorological stations
                    df, meta = get_daily_meteo_api(
                        stn=stn, meteo_var=meteo_var, start=start, end=end
                    )
                    add_day = True
                if not df.empty:
                    ts_df, variables = interpret_knmi_file(
                        df=df,
                        meta=meta,
                        meteo_var=meteo_var,
                        start=start,
                        end=end,
                        add_day=add_day,
                        add_hour=True,
                    )

        except (RuntimeError, requests.ConnectionError) as e:
            logger.warning(
                "KNMI API failed, try setting the 'use_api' argument to 'False'"
            )
            if settings["raise_exceptions"]:
                raise e
            logger.info("Try non api method")
            settings = settings.copy()
            settings["use_api"] = False

        if not settings["use_api"]:
            if settings["interval"].startswith("hour"):
                # hourly data from meteorological stations
                raise NotImplementedError()
            elif meteo_var == "RD":
                # daily data from rainfall-stations
                df, meta = get_daily_rainfall_url(stn, stn_name)
                add_day = False
            else:
                # daily data from meteorological stations
                df, meta = get_daily_meteo_url(stn=stn)
                add_day = True
            if not df.empty:
                ts_df, variables = interpret_knmi_file(
                    df=df,
                    meta=meta,
                    meteo_var=meteo_var,
                    start=start,
                    end=end,
                    add_day=add_day,
                    add_hour=True,
                )
    except (ValueError, KeyError, pd.errors.EmptyDataError) as e:
        logger.error(f"{e} {msg}")
        if settings["raise_exceptions"]:
            raise e

    stations = get_stations(meteo_var=meteo_var).loc[[stn], :]

    return ts_df, variables, stations


def get_knmi_daily_rainfall_api(
    stn: int,
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily rainfall.

    .. deprecated:: 0.13.3
        `get_knmi_daily_rainfall_api` will be removed in hydropandas 1.0.0, it is replaced
        by `get_daily_rainfall_api`.

    Parameters
    ----------
    stn : int
        station number.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.

    Raises
    ------
    ValueError
        if there is no data for the provided stn an error is raised.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    warnings.warn(
        "the function 'get_knmi_daily_rainfall_api' is deprecated and will eventually be "
        "removed, please use 'get_daily_rainfall_api'.",
        DeprecationWarning,
    )

    return get_daily_rainfall_api(stn, start, end)


@lru_cache()
def get_daily_rainfall_api(
    stn: int,
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily rainfall.

    Parameters
    ----------
    stn : int
        station number.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.

    Raises
    ------
    ValueError
        if there is no data for the provided stn an error is raised.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    meteo_var = "RD"
    url = URL_DAILY_PREC
    params = {"vars": meteo_var, "stns": str(stn)}

    if start is not None:
        params["start"] = (start - pd.Timedelta(days=1)).strftime("%Y%m%d")
    if end is not None:
        params["end"] = end.strftime("%Y%m%d")
    strio = request_api(url, params)

    return parse_data(strio)


def request_url(url: str, fname=None) -> StringIO:
    """download and read knmi daily rainfall data.

    Parameters
    ----------
    stn : int
        station number.
    fname : str or None, optional
        filename to save the data to, only used if not None. The default is None.

    Returns
    -------
    StringIO
    """

    # request zipfile
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()

    # unpack data
    with ZipFile(BytesIO(r.content), mode="r") as zf:
        zipfile = zf.read(zf.namelist()[0])
        result_str = zipfile.decode().replace("\r\n", "\n")

    if fname is not None:
        with open(fname, "w") as f:
            f.write(result_str)

    return StringIO(result_str)


def get_knmi_daily_rainfall_url(
    stn: int,
    stn_name: str,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily rainfall.

    .. deprecated:: 0.13.3
        `get_knmi_daily_rainfall_url` will be removed in hydropandas 1.0.0, it is replaced
        by `get_daily_rainfall_url`.

    Parameters
    ----------
    stn : int
        station number.
    stn_name : str
        the name of the station in capital letters, can be tricky

    Raises
    ------
    ValueError
        if there is no data for the provided stn an error is raised.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    warnings.warn(
        "the function 'get_knmi_daily_rainfall_url' is deprecated and will eventually be "
        "removed, please use 'get_daily_rainfall_url'.",
        DeprecationWarning,
    )

    get_daily_rainfall_url(stn, stn_name)


@lru_cache()
def get_daily_rainfall_url(
    stn: int,
    stn_name: str,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily rainfall.

    Parameters
    ----------
    stn : int
        station number.
    stn_name : str
        the name of the station in capital letters, can be tricky

    Raises
    ------
    ValueError
        if there is no data for the provided stn an error is raised.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    stn = f"{stn:03d}"  # make sure there are leading zeros
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie/"
        f"gegevens/monv_reeksen/neerslaggeg_{stn_name}_{stn}.zip"
    )

    strio = request_url(url)
    return parse_data(strio)


def _transform_variables(
    df: pd.DataFrame, variables: Dict[str, Any]
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Transforms the timeseries to default units and settings.

    Does 3 things:
        1. all values equal to -1 are converted to zero
        2. the units are changed from 0.1 mm to 1 mm.
        3. the units are changed from mm to m.

    Parameters
    ----------
    df : pandas DataFrame
        time series.
    variables : dictionary
        description of variables in time series.

    Raises
    ------
    NameError
        if there are columns in the DataFrame and no matching key in the
        variables dictionary.

    Returns
    -------
    df : pandas DataFrame
        time series.
    variables : dictionary
        description of variables in time series.
    """
    add_m_unit = False
    for key, value in variables.items():
        # test if key existst in data
        if key not in df.columns:
            if key == "YYYYMMDD" or key == "HH":
                pass
            elif key == "T10N":
                variables.pop(key)
                key = "T10"
            else:
                raise NameError(key + " does not exist in data")

        if "(-1 voor <0.05 mm)" in value:
            # remove -1 for precipitation smaller than <0.05 mm
            df.loc[df.loc[:, key] == -1, key] = 0.0
            value = value.replace("(-1 voor <0.05 mm)", "(0 voor <0.05mm)").replace(
                "(-1 for <0.05 mm)", "(0 for <0.05mm)"
            )
            logger.debug(f"transform {key}, {value} lower than 0.05 mm to 0 mm")

        if "0.1 " in value:
            logger.debug(f"transform {key}, {value} from 0.1 to 1")

            df[key] = df[key] * 0.1
            value = value.replace("0.1 ", "")
        if " tiende " in value:
            logger.debug(f"transform {key}, {value} from 0.1 to 1")

            df[key] = df[key] * 0.1
            value = value.replace(" tiende ", " ")
        if " mm" in value:
            logger.debug(f"transform {key}, {value} from mm to m")

            df[key] = df[key] * 0.001
            value = value.replace(" mm", " m")
            add_m_unit = True
        if " millimeters" in value:
            logger.debug(f"transform {key}, {value} from mm to m")
            df[key] = df[key] * 0.001
            add_m_unit = True
        if "08.00 UTC" in value:
            logger.debug(f"transform {key}, {value} from UTC to UTC+1")

            # over the period 08.00 preceding day - 08.00 UTC present day
            df.index = df.index + pd.to_timedelta(8, unit="h")

            value = value.replace("08.00", "09.00").replace("UTC", "UTC+1")

        # Store new variable
        variables[key] = value

    if add_m_unit:
        variables["unit"] = "m"
    else:
        variables["unit"] = ""
    return df, variables


def request_api(url: str, params: Dict[str, str], fname=None) -> StringIO:
    """Download KNMI data from the API

    Parameters
    ----------
    url : str
        URL to parse the request to
    params : Dict[str, str]
        Dictionary with parameters that are parsed to the request get
    fname : str or None, optional
        filename to save the data to only used if not None, by default None

    Returns
    -------
    StringIO

    """
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()

    result_str = r.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down\n" + result_str)

    if fname is not None:
        with open(fname, "w") as f:
            f.write(result_str)

    return StringIO(result_str)


def get_knmi_daily_meteo_api(
    stn, start=None, end=None, meteo_var: Union[str, None] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily meteo data.

    Parameters
    ----------
    stn : int
        station number.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    stations : pandas DataFrame
        additional data about the measurement station
    """

    warnings.warn(
        "the function 'get_knmi_daily_meteo_api' is deprecated and will eventually be "
        "removed, please use 'get_daily_meteo_api'.",
        DeprecationWarning,
    )

    return get_daily_meteo_api(stn, start, end, meteo_var)


@lru_cache()
def get_daily_meteo_api(
    stn, start=None, end=None, meteo_var: Union[str, None] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily meteo data.

    Parameters
    ----------
    stn : int
        station number.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    stations : pandas DataFrame
        additional data about the measurement station
    """

    url = URL_DAILY_METEO

    params = {
        "stns": str(stn),
    }
    if meteo_var is not None:
        params["vars"] = meteo_var

    # modify start and end date for the api call because the dates are later
    # transformed, see example notebook 02_knmi_observations.
    if start is not None:
        start = pd.Timestamp(start) - dt.timedelta(days=1)
        params["start"] = start.strftime("%Y%m%d")
    if end is not None:
        end = pd.Timestamp(end) - dt.timedelta(days=1)
        params["end"] = end.strftime("%Y%m%d")

    strio = request_api(url, params)

    return parse_data(strio)


def get_knmi_daily_meteo_url(stn: int) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily meteo data.

    .. deprecated:: 0.13.3
        `get_knmi_daily_meteo_url` will be removed in hydropandas 1.0.0, it is replaced
        by `get_daily_meteo_url`.

    Parameters
    ----------
    stn : int
        station number.

    Returns
    -------
    pandas DataFrame
        measurements.
    meta : dictionary
        additional information about the variables
    """
    warnings.warn(
        "the function 'get_knmi_daily_meteo_url' is deprecated and will eventually be "
        "removed, please use 'get_daily_meteo_url'.",
        DeprecationWarning,
    )

    return get_daily_meteo_url(stn)


@lru_cache()
def get_daily_meteo_url(stn: int) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """download and read knmi daily meteo data.

    Parameters
    ----------
    stn : int
        station number.

    Returns
    -------
    pandas DataFrame
        measurements.
    meta : dictionary
        additional information about the variables
    """
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie"
        f"/gegevens/daggegevens/etmgeg_{stn}.zip"
    )

    strio = request_url(url)
    return parse_data(strio)


def parse_data(
    path: Union[str, Path, StringIO],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """read knmi daily meteo data from a file

    Parameters
    ----------
    path : str or StringIO
        file path of .txt file or str representation.

    Returns
    -------
    pandas DataFrame
        measurements.
    meta : dictionary
        additional information about the variables
    """
    meta_id1 = " = "
    meta_id2 = " : "
    data_id = "STN,"
    meta = {}
    df = None
    # with open(path, "r") as f:
    if isinstance(path, str):
        f = open(path, "r")
    elif isinstance(path, StringIO):
        f = path
    else:
        raise TypeError("Expected a path as str or StringIO object")

    line = f.readline()
    while line:
        if meta_id1 in line or meta_id2 in line:
            key, item = (
                line.split(meta_id1) if meta_id1 in line else line.split(meta_id2)
            )
            meta[key.strip("# ")] = item.strip()

        elif data_id in line:
            try:
                columns = line.strip("# ").strip("\n").split(",")
                columns = [x.strip(" ") for x in columns]
                df = pd.read_csv(
                    f,
                    header=None,
                    na_values="     ",
                    delimiter=",",
                    skip_blank_lines=True,
                )
                if columns[-2:] == ["HH", "YYYYMMDD"]:
                    columns = columns[:-2]
                df.columns = columns
                hour_col = "HH" if "HH" in columns else None
                if hour_col is not None:
                    df.loc[df[hour_col] == 24, hour_col] = 0
                    datetime = pd.to_datetime(
                        df.YYYYMMDD.astype(str) + df[hour_col].astype(str).str.zfill(2),
                        format="%Y%m%d%H",
                    )
                    datetime.loc[datetime.dt.hour == 0] = datetime + dt.timedelta(
                        days=1
                    )
                    df = df.drop(columns=["YYYYMMDD", hour_col]).loc[
                        df.index.notnull(), :
                    ]
                else:
                    datetime = pd.to_datetime(df.YYYYMMDD, format="%Y%m%d")
                    df = df.drop(columns=["YYYYMMDD"]).loc[df.index.notnull(), :]

                df = df.set_index(datetime)
            except pd.errors.EmptyDataError as e:
                logger.warning(f"{str(e)}. Returning empty DataFrame.")
                df = pd.DataFrame()
            f.close()

            return df, meta

        line = f.readline()

    f.close()

    if df is None:
        raise ValueError(f"Could not read: {path}")


def interpret_knmi_file(
    df: pd.DataFrame,
    meta: Dict[str, Any],
    meteo_var: str,
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
    add_day: bool = False,
    add_hour: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """interpret data from knmi by selecting meteo_var data and meta
    and transforming the variables

    Parameters
    ----------
    df: DataFrame
        dataframe with meteo_var as column
    meta: dictionary
        dictionary with meteo_var as key
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.
    add_day : boolean, optional
        add 1 day so that the timestamp is at the end of the period the data describes,
        default is False, and has to be set per type of file.
    add_hour : boolean, optional
        add 1 hour to convert from UT to UT+1 (standard-time in the Netherlands),
        default is True as this is usually the case.


    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    station : int
        meteostation
    """

    variables = {meteo_var: meta[meteo_var]}
    stn = None
    if not df.empty:
        unique_stn = df["STN"].unique()
        if len(unique_stn) > 1:
            raise ValueError(
                f"Cannot handle multiple stations {unique_stn} in single file"
            )
        stn = unique_stn[0]

        if add_day or add_hour:
            if add_day and add_hour:
                timedelta = pd.Timedelta(1, "d") + pd.Timedelta(1, "h")
            elif add_hour:
                timedelta = pd.Timedelta(1, "h")
            else:
                timedelta = pd.Timedelta(1, "d")

            df = df.copy()
            df.index = df.index + timedelta

        if df.index.has_duplicates:
            df = df.loc[~df.index.duplicated(keep="first")]
            logger.info("duplicate indices removed from RD measurements")

        istart = (
            df.index.get_indexer([start], method="backfill")[0]
            if start is not None
            else 0
        )
        iend = (
            df.index.get_indexer([end], method="backfill")[0] if end is not None else -1
        )
        iend = len(df) if iend == -1 else iend + 1
        icol = df.columns.get_indexer([meteo_var])
        meteo_df = df.iloc[istart:iend, icol].dropna()

        if not meteo_df.empty:
            mdf, var = _transform_variables(meteo_df, variables)
            variables["station"] = stn
            return mdf, var

    return pd.DataFrame(), variables


def get_knmi_hourly_meteo_api(
    stn: int, start: pd.Timestamp, end: pd.Timestamp, meteo_var: Union[str, None] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Retrieve hourly meteorological data from the KNMI API.

    .. deprecated:: 0.13.3
        `get_knmi_hourly_meteo_api` will be removed in hydropandas 1.0.0, it is replaced
        by `get_hourly_meteo_api`.

    Parameters
    ----------
    stn : int
        The station number.
    start : pd.Timestamp or None
        The start date and time of the data.
    end : pd.Timestamp or None
        The end date and time of the data.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the requested meteorological variable data.
    dict
        A dictionary containing information about the variables in the DataFrame.

    Raises
    ------
    requests.ConnectionError
        If there is a connection error while accessing the KNMI API.
    RuntimeError
        If the KNMI API is down or returns unexpected data.
    """
    warnings.warn(
        "the function 'get_knmi_hourly_meteo_api' is deprecated and will eventually be "
        "removed, please use 'get_hourly_meteo_api'.",
        DeprecationWarning,
    )

    return get_hourly_meteo_api(stn, start, end, meteo_var)


@lru_cache()
def get_hourly_meteo_api(
    stn: int, start: pd.Timestamp, end: pd.Timestamp, meteo_var: Union[str, None] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Retrieve hourly meteorological data from the KNMI API.

    Parameters
    ----------
    stn : int
        The station number.
    start : pd.Timestamp or None
        The start date and time of the data.
    end : pd.Timestamp or None
        The end date and time of the data.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the requested meteorological variable data.
    dict
        A dictionary containing information about the variables in the DataFrame.

    Raises
    ------
    requests.ConnectionError
        If there is a connection error while accessing the KNMI API.
    RuntimeError
        If the KNMI API is down or returns unexpected data.
    """
    url = URL_HOURLY_METEO

    params = {
        "stns": str(stn),
    }
    if meteo_var is not None:
        params["vars"] = meteo_var

    # It looks like the API has an option to select the start and end hour that
    # you want but in reality it does not have this option.
    if start is None:
        raise ValueError("A start date is required when using hourly interval")
    if end is None:
        raise ValueError("An end date is required when using hourly interval")

    if (end - start).days > 365 * 10:
        raise ValueError("time span for hourly data cannot be greater than 10 years")
    if (end - start).days < 1:
        raise ValueError("time span should be more than 1 day")

    params["end"] = end.strftime("%Y%m%d") + "24"

    s = start - pd.Timedelta(1, unit="h")
    params["start"] = s.strftime("%Y%m%d") + "01"

    strio = request_api(url=url, params=params)

    return parse_data(strio)


def get_nearest_station_df(
    locations: pd.DataFrame,
    xcol: str = "x",
    ycol: str = "y",
    stations: Union[pd.DataFrame, None] = None,
    meteo_var: str = "RH",
    start: Union[pd.Timestamp, str, None] = None,
    end: Union[pd.Timestamp, str, None] = None,
    ignore: Union[List[str], None] = None,
) -> list[int]:
    """Find the KNMI stations that measure 'meteo_var' closest to the
    coordinates in 'locations'.

    Parameters
    ----------
    locations : pandas.DataFrame
        DataFrame containing x and y coordinates
    xcol : str
        name of the column in the locations dataframe with the x values
    ycol : str
        name of the column in the locations dataframe with the y values
    stations : pandas DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    meteo_var : str
        measurement variable e.g. 'RH' or 'EV24'
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    ignore : list, optional
        list of stations to ignore. The default is None.

    Returns
    -------
    stns : list
        station numbers.
    """
    if stations is None:
        stations = get_stations(meteo_var=meteo_var, start=start, end=end)
    if ignore is not None:
        stations = stations.drop(ignore)
        if stations.empty:
            return None

    xo = pd.to_numeric(locations[xcol])
    xt = pd.to_numeric(stations.x)
    yo = pd.to_numeric(locations[ycol])
    yt = pd.to_numeric(stations.y)

    xh, xi = np.meshgrid(xt, xo)
    yh, yi = np.meshgrid(yt, yo)

    distances = pd.DataFrame(
        np.sqrt((xh - xi) ** 2 + (yh - yi) ** 2),
        index=locations.index,
        columns=stations.index,
    )

    stns = distances.idxmin(axis=1).to_list()

    return stns


def get_nearest_station_xy(
    xy: List[List[float]],
    stations: Union[pd.DataFrame, None] = None,
    meteo_var: str = "RH",
    ignore: Union[List[str], None] = None,
) -> List[int]:
    """find the KNMI stations that measure 'meteo_var' closest to the given
    x and y coordinates.

    Parameters
    ----------
    xy : list or numpy array, optional
        xy coordinates of the locations. e.g. [[10,25], [5,25]]
    stations : pandas DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    meteo_var : str
        measurement variable e.g. 'RH' or 'EV24'
    ignore : list, optional
        list of stations to ignore. The default is None.

    Returns
    -------
    stns : list
        station numbers.

    Notes
    -----
    assumes you have a structured rectangular grid.
    """

    locations = pd.DataFrame(data=xy, columns=["x", "y"])

    stns = get_nearest_station_df(
        locations,
        xcol="x",
        ycol="y",
        stations=stations,
        meteo_var=meteo_var,
        ignore=ignore,
    )

    return stns


def get_n_nearest_stations_xy(
    xy: List[List[float]],
    meteo_var: str,
    start: Union[pd.Timestamp, str, None] = None,
    end: Union[pd.Timestamp, str, None] = None,
    n: int = 1,
    stations: Union[pd.DataFrame, None] = None,
    ignore: Union[List[str], None] = None,
) -> List[int]:
    """Find the N nearest KNMI stations that measure variable 'meteo_var' to
    the x, y coordinates.

    Parameters
    ----------
    xy : list, tuple or numpy.array of int or float
        sinlge pair of xy coordinates. e.g. (150_000., 400_000.)
    meteo_var : str
        measurement variable e.g. 'RH' or 'EV24'
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    n : int, optional
        number of stations you want to return. The default is 1.
    stations : pandas DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    ignore : list, optional
        list of stations to ignore. The default is None.

    Returns
    -------
    list
        station numbers.
    """

    if stations is None:
        stations = get_stations(meteo_var=meteo_var, start=start, end=end)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
        if stations.empty:
            return None

    distance = np.sqrt((stations.x - xy[0]) ** 2 + (stations.y - xy[1]) ** 2)
    stns = distance.nsmallest(n).index.to_list()

    return stns


def _add_missing_indices(
    ts_df: pd.DataFrame, stn: int, start: pd.Timestamp, end: pd.Timestamp
) -> pd.DataFrame:
    """When downloading KNMI data you don't always get a DataFrame with the
    periods that you provided in your request. Thus the index does not cover
    the complete period that you are interested in. This function adds the
    missing period to the index of the DataFrame.

    Parameters
    ----------
    ts_df : pandas DataFrame
        data from one station from one type of observation, with additional
        column to see which station is used to fill the value
    stn : int or str
        measurement station.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.

    Returns
    -------
    ts_df : pandas DataFrame
        data from one station from one type of observation
    """
    # check if given dates are more or less similar than measurement dates

    if (ts_df.index[0] - start).days < 1:
        new_start = ts_df.index[0]
    else:
        new_start = pd.Timestamp(
            year=start.year,
            month=start.month,
            day=start.day,
            hour=ts_df.index[0].hour,
            minute=ts_df.index[0].minute,
            second=ts_df.index[0].second,
        )
        logger.info(f"station {stn} has no measurements before {ts_df.index[0]}")

    if (end - ts_df.index[-1]).days < 0:
        new_end = ts_df.index[-1]
    else:
        new_end = pd.Timestamp(
            year=end.year,
            month=end.month,
            day=end.day,
            hour=ts_df.index[-1].hour,
            minute=ts_df.index[-1].minute,
            second=ts_df.index[-1].second,
        )
        logger.info(f"station {stn} has no measurements after {ts_df.index[-1]}")

    # add missing indices
    new_index = pd.date_range(new_start, new_end, freq="D")
    ts_df = ts_df.reindex(new_index)

    return ts_df


def get_knmi_obslist(
    locations: Union[pd.DataFrame, None] = None,
    stns: Union[List[int], None] = None,
    xy: Union[List[List[float]], None] = None,
    meteo_vars: Tuple[str] = ("RH",),
    starts: Union[pd.Timestamp, List[pd.Timestamp], None] = None,
    ends: Union[pd.Timestamp, List[pd.Timestamp], None] = None,
    ObsClasses: List[Any] = None,
    **kwargs,
) -> List[Any]:
    """Get a list of observations of knmi stations. Either specify a list of
    knmi stations (stns) or a dataframe with x, y coordinates (locations).

    Parameters
    ----------
    locations : pandas DataFrame or None
        dataframe with x and y coordinates. The default is None
    stns : list of int or None
        list of knmi stations. The default is None
    xy : list or numpy array, optional
        xy coordinates of the locations. e.g. [[10,25], [5,25]]
    meteo_vars : list or tuple of str
        meteo variables e.g. ["RH", "EV24"]. The default is ("RH")
    starts : None, str, datetime or list, optional
        start date of observations per meteo variable. The start date is
        included in the time series.
        If start is None the start date will be January 1st of the
        previous year.
        if start is str it will be converted to datetime
        if start is a list it should be the same length as meteo_vars and
        the start time for each variable. The default is None
    ends : list of str, datetime or None
        end date of observations per meteo variable. The end date is
        included in the time series.
        If end is None the start date will be January 1st of the
        previous year.
        if end is a str it will be converted to datetime
        if end is a list it should be the same length as meteo_vars and
        the end time for each meteo variable. The default is None
    ObsClasses : list of type or None
        class of the observations, can be PrecipitationObs or
        EvaporationObs. The default is None.
    **kwargs:
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False. Note: if the given stn has no data between start and
            end the data from nearby stations is used. In this case the metadata of the
            Observation is the metadata from the nearest station that has any
            measurement in the given period.
        interval : str, optional
            desired time interval for observations. Options are 'daily' and
            'hourly'. The default is 'daily'.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary folder and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

    Returns
    -------
    obs_list : list of observation objects
        collection of multiple point observations
    """

    settings = _get_default_settings(kwargs)

    if isinstance(meteo_vars, str):
        meteo_vars = [meteo_vars]

    if starts is None:
        starts = [None] * len(meteo_vars)
    elif isinstance(starts, (str, dt.datetime)):
        starts = [starts] * len(meteo_vars)
    elif isinstance(starts, list):
        pass
    else:
        raise TypeError("must be None, str, dt.datetime or list")

    if ends is None:
        ends = [None] * len(meteo_vars)
    elif isinstance(ends, (str, dt.datetime)):
        ends = [ends] * len(meteo_vars)
    elif isinstance(ends, list):
        pass
    else:
        raise TypeError("must be None, str, dt.datetime or list")

    if not isinstance(xy, (tuple, list, np.ndarray, type(None))):
        raise TypeError("must be tuple, list or numpy.ndarray")

    obs_list = []
    for meteo_var, start, end, ObsClass in zip(meteo_vars, starts, ends, ObsClasses):
        start = start if start is None else pd.to_datetime(start)
        end = end if end is None else pd.to_datetime(end)

        # get stations
        if stns is None:
            stations = get_stations(meteo_var=meteo_var, start=start, end=end)
            if (locations is None) and (xy is not None):
                _stns = get_nearest_station_xy(
                    xy, stations=stations, meteo_var=meteo_var
                )
            elif locations is not None:
                _stns = get_nearest_station_df(
                    locations,
                    stations=stations,
                    meteo_var=meteo_var,
                    start=start,
                    end=end,
                )
            else:
                raise ValueError(
                    "stns, location and xy are all None. "
                    "Please specify one of these arguments."
                )
            _stns = np.unique(_stns)

            if locations is not None:
                xy = locations[["x", "y"]].values
        else:
            _stns = stns

        for stn in _stns:
            o = ObsClass.from_knmi(
                meteo_var=meteo_var,
                stn=stn,
                start=start,
                end=end,
                **settings,
            )

            obs_list.append(o)

    return obs_list


def get_evaporation(
    meteo_var: str,
    stn: int = 260,
    start: Union[pd.Timestamp, None] = None,
    end: Union[pd.Timestamp, None] = None,
    settings: Union[Dict[str, Any], None] = None,
) -> pd.DataFrame:
    """Collect different types of (reference) evaporation
    from KNMI weather stations

    Parameters
    ----------
    meteo_var : str
        Choice between 'penman', 'makkink' or 'hargraves'.
    stn : str
        station number, defaults to 260 De Bilt
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    settings : dict or None, optional
        settings for the time series

    Returns
    ------
    pandas.DataFrame

    """
    if meteo_var == "hargreaves":
        if not settings["use_api"]:
            raise NotImplementedError(
                "cannot use hargreaves without the api because hargreaves needs "
                "the lattitude"
            )
        d = {}
        mvs = ["TG", "TN", "TX"]
        for mv in mvs:
            ts, meta = get_timeseries_stn(
                stn, meteo_var=mv, start=start, end=end, settings=settings
            )
            d[mv] = ts.squeeze()
        et = hargreaves(
            d["TG"],
            d["TN"],
            d["TX"],
            d["TG"].index,
            meta["lat"] if "lat" in meta else 52.1,
        ).to_frame(name=meteo_var)
    elif meteo_var == "makkink":
        d = {}
        mvs = ["TG", "Q"]
        for mv in mvs:
            ts, meta = get_timeseries_stn(
                stn, meteo_var=mv, start=start, end=end, settings=settings
            )
            d[mv] = ts.squeeze()
        et = makkink(d["TG"], d["Q"]).to_frame(name=meteo_var)
    elif meteo_var == "penman":
        if not settings["use_api"]:
            raise NotImplementedError(
                "cannot use penman without the api because penman needs the "
                "lattitude and height of the meteo station"
            )
        d = {}
        mvs = ["TG", "TN", "TX", "Q", "FG", "UG"]
        for mv in mvs:
            ts, meta = get_timeseries_stn(
                stn, meteo_var=mv, start=start, end=end, settings=settings
            )
            d[mv] = ts.squeeze()
        et = penman(
            d["TG"],
            d["TN"],
            d["TX"],
            d["Q"],
            d["FG"],
            d["UG"],
            d["TG"].index,
            meta["lat"] if "lat" in meta else 52.1,
            meta["altitude"] if "altitude" in meta else 0.0,
        ).to_frame(name=meteo_var)
    else:
        raise ValueError(
            "Provide valid argument for meteo_var -> 'hargreaves', 'makkink' or "
            "'penman'"
        )

    stn_name = get_station_name(meta["station"], stations=get_stations(meteo_var))
    meta["name"] = f"{meteo_var}_{stn_name}_{stn}"
    meta["location"] = stn_name
    meta["unit"] = "m"

    return et, meta


def makkink(tmean: pd.Series, K: pd.Series) -> pd.Series:
    """Estimate of Makkink reference evaporation
    according to KNMI.

    Parameters
    ----------
    tmean : tmean : pandas.Series
        Daily mean temperature in Celsius
    K : pandas.Series
        Global radiation estimate in J/cm2

    Returns
    -------
    pandas.Series

    """

    a = 0.646 + 0.0006 * tmean
    b = 1 + tmean / 237.3
    c = 7.5 * np.log(10) * 6.107 * 10 ** (7.5 * (1 - 1 / b))
    et = 0.0065 * (1 - a / (c / (237.3 * b * b) + a)) / (2501 - 2.38 * tmean) * K
    return et


def penman(
    tmean: pd.Series,
    tmin: pd.Series,
    tmax: pd.Series,
    K: pd.Series,
    wind: pd.Series,
    rh: pd.Series,
    dates: pd.Series,
    z: float = 1.0,
    lat: float = 52.1,
    G: float = 0.0,
    wh: float = 10.0,
    tdew: Union[pd.Series, None] = None,
) -> pd.Series:
    """Estimate of Penman reference evaporation
    according to Allen et al 1990.

    Parameters
    ----------
    tmean : pandas.Series
        Daily mean temperature in Celsius
    tmin : pandas.Series
        Daily minimum temperature in Celsius
    tmax : pandas.Series
        Daily maximum temperature in Celsius
    K : pandas.Series
        Global radiation estimate in J/cm2
    wind : pandas.Series
        Daily mean wind speed in m/s
    rh : pandas.Series
        Relative humidity in %
    dates : pandas.Series
        Dates
    z : float, optional
        Elevation of station in m, by default 1.0
    lat : float, optional
        Latitude of station, by default 52.1
    G : float, optional
        Ground flux in MJ/m2, by default 0.0
    wh : float, optional
        Height of wind measurement in m, by default 10.0
    tdew : pandas.Series
        Dew point temperature in C, by default None

    Returns
    -------
    pandas.Series

    """
    K = K * 1e4 * 0.000012  # J/cm2 to W/m2/d
    P = 101.3 * ((293 - 0.0065 * z) / 293) ** 5.26  # kPa
    gamma = 0.665e-3 * P  # kPa/C
    tg = (tmax - tmin) / 2  # C
    s = 4098 * (0.6108 * np.exp(17.27 * tg / (tg + 237.3)) / (tg + 237.3) ** 2)  # kPa/C
    es0 = 0.6108 * np.exp(17.27 * tmean / (tmean + 237.3))  # kPa
    esmax = 0.6108 * np.exp(17.27 * tmax / (tmax + 237.3))  # kPa
    esmin = 0.6108 * np.exp(17.27 * tmin / (tmin + 237.3))  # kPa
    es = (esmax + esmin) / 2  # kpa
    if tdew:
        ea = 0.6108 * np.exp(17.27 * tdew / (tdew + 237.3))  # kPa
    else:
        ea = es0 * rh / 100  # kPa
    u2 = wind * 4.87 / np.log(67.8 * wh - 5.42)  # m/s
    J = pd.Series([int(date.strftime("%j")) for date in dates], index=dates)
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)  # -
    delt = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39)  # rad
    phi = lat * np.pi / 180  # rad
    ws = np.arccos(-np.tan(phi) * np.tan(delt))  # rad
    Kext = (
        1366
        * dr
        / np.pi
        * (ws * np.sin(phi) * np.sin(delt) + np.cos(phi) * np.cos(delt) * np.sin(ws))
    )  # W/m2/d
    K0 = Kext * (0.75 + 2e-5 * z)  # W/m2/d
    Lout_in = (
        5.67e-8
        * (((tmax + 273) ** 4 + (tmin + 273) ** 4) / 2)
        * (0.34 - 0.14 * np.sqrt(ea))
        * (1.35 * K / K0 - 0.35)
    )  # W/m2/d
    Q = (1 - 0.23) * K + Lout_in  # W/m2/d
    Q = Q * 0.0864  # W/m2/d to MJ/m2
    straling = 0.408 * s * (Q - G)
    windterm = gamma * 900 / (tmean + 273) * u2 * (es - ea)
    denom = s + gamma * (1 + 0.34 * u2)
    et = (straling + windterm) / denom * 1e-3
    et[et < 0] = 0
    return et


def hargreaves(
    tmean: pd.Series,
    tmin: pd.Series,
    tmax: pd.Series,
    dates: pd.Series,
    lat: float = 52.1,
    x: Union[List[float], None] = None,
) -> pd.Series:
    """Estimate of Hargraves potential evaporation
    according to Allen et al. 1990.

    Parameters
    ----------
    tmean : pandas.Series
        Daily mean temperature
    tmin : pandas.Series
        Daily minimum temperature
    tmax : pandas.Series
        Daily maximum temperature
    dates : pandas.Series.index
        Dates
    lat : float, optional
        Latitude of station, by default 52.1
    x : _type_, optional
        Optional parameter to scale evaporation estimate, by default None

    Returns
    -------
    pandas.Series

    """
    J = pd.Series([int(date.strftime("%j")) for date in dates], index=dates)
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
    delt = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39)
    phi = lat * np.pi / 180
    ws = np.arccos(-np.tan(phi) * np.tan(delt))
    Kext = (
        12
        * 60
        * 0.0820
        * dr
        / np.pi
        * (ws * np.sin(phi) * np.sin(delt) + np.cos(phi) * np.cos(delt) * np.sin(ws))
    )  # MJ/m2/d
    Kext = 0.408 * Kext  # MJ/m2/d to mm/d
    et = 0.0023 * (tmean + 273 + 17.8) / np.sqrt(tmax - tmin) * Kext * 1e-3
    if x:
        et = x[0] + x[1] * et
    return et
