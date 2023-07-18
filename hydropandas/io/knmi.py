"""Module with functions to read or download time series with observations from knmi.

The main function to download time series are:

- get_knmi_timeseries_xy: obtain a knmi station based on the xy location
- get_knmi_timeseries_stn: obtain a knmi station based on the station number

The main function to read time series txt files is:

 - read_knmi_timeseries_file: read a knmi txt file

"""
import datetime as dt
import logging
import os
import re
import tempfile
from functools import lru_cache
from io import StringIO

import numpy as np
import pandas as pd
import requests

from .. import util

logger = logging.getLogger(__name__)

URL_DAILY_PREC = "https://www.daggegevens.knmi.nl/klimatologie/monv/reeksen"
URL_DAILY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/daggegevens"
URL_HOURLY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/uurgegevens"
LOOK_BACK_DAYS = 365


def get_knmi_obs(
    stn=None,
    fname=None,
    xy=None,
    meteo_var=None,
    start=None,
    end=None,
    **kwargs,
):
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
        meteo variable e.g. "RH" or "EV24". See list with al options in the
        hydropandas documentation.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    **kwargs:
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False.
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

    start, end = _start_end_to_datetime(start, end)

    if stn is not None:
        stn = int(stn)

        logger.info(
            f"get KNMI data from station {stn} and meteo variable {meteo_var} from {start} to {end}"
        )
        ts, meta = get_knmi_timeseries_stn(stn, meteo_var, settings, start, end)
    elif fname is not None:
        logger.info(f"get KNMI data from file {fname} and meteo variable {meteo_var}")
        ts, meta = get_knmi_timeseries_fname(fname, meteo_var, settings, start, end)

    elif xy is not None:
        logger.info(
            f"get KNMI data from station nearest to coordinates {xy} and meteo variable {meteo_var}"
        )
        stns = get_n_nearest_stations_xy(xy, meteo_var, n=1)
        ts, meta = get_knmi_timeseries_stn(stns[0], meteo_var, settings, start, end)
    else:
        raise ValueError(
            "specify KNMI station (stn), filename (fname) or coordinates (xy)"
        )

    return ts, meta


def get_knmi_timeseries_fname(fname, meteo_var, settings, start, end):
    # first try to figure out filetype by it's name
    if "neerslaggeg" in fname:
        # neerslagstation
        ts, meta = read_knmi_daily_rainfall_file(fname, start, end)
        meteo_var = "RD"
    elif "etmgeg" in fname:
        # meteo station
        ts, meta = read_knmi_daily_meteo_file(fname, meteo_var, start, end)
    elif "uurgeg" in fname:
        # hourly station
        ts, meta = read_knmi_hourly(fname, meteo_var, start, end)
    # if that doesn't work try to figure out by the meteo_var and settings
    elif meteo_var is None or meteo_var == "RD":
        # neerslagstation
        ts, meta = read_knmi_daily_rainfall_file(fname, start, end)
        meteo_var = "RD"
    elif settings["interval"] == "daily":
        # meteo station
        ts, meta = read_knmi_daily_meteo_file(fname, meteo_var, start, end)
    elif settings["interval"] == "hourly":
        # uurlijks station
        ts, meta = read_knmi_hourly(fname, meteo_var, start, end)
    else:
        raise ValueError(
            "please indicate how to read the file by specifying a meteo_var and an interval"
        )
    stn = meta["station"]
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn, stations)

    # set metadata
    meta.update(
        {
            "x": stations.loc[stn, "x"],
            "y": stations.loc[stn, "y"],
            "name": f"{meteo_var}_{stn_name}",
            "source": "KNMI",
            "filename": fname,
        }
    )

    return ts, meta


def _get_default_settings(settings=None):
    """adds the default settings to a dictinary with settings. If settings
    is None all the settings are default. If there are already settings given
    only the non-existing settings are added with their default value.

    The default settings are:
    fill_missing_obs = True
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
            logger.debug("set raise_exceptions=False because fill_missing_obs is True")
            if settings["fill_missing_obs"] and settings["raise_exceptions"]:
                settings["raise_exceptions"] = False
        else:
            settings["raise_exceptions"] = False

    for key, value in default_settings.items():
        if key not in settings.keys():
            settings[key] = value

    return settings


def _start_end_to_datetime(start, end):
    """convert start and endtime to datetime.

    Parameters
    ----------
    start : str, datetime, None
        start time
    end : str, datetime, None
        start time

    Returns
    -------
    start : pd.TimeStamp or None
        start time
    end : pd.TimeStamp or None
        end time
    """

    if start is not None:
        start = pd.to_datetime(start)

    if end is not None:
        end = pd.to_datetime(end)

    return start, end


def get_knmi_timeseries_stn(stn, meteo_var, settings, start=None, end=None):
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
    knmi_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn, stations)

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
            "No hourly evaporation data available through the api, "
            "set use_api=False."
        )
        raise ValueError(message)

    # download data
    elif settings["fill_missing_obs"]:
        knmi_df, meta = fill_missing_measurements(
            stn, meteo_var, start, end, settings, stn_name
        )
    else:
        knmi_df, variables, station_meta = download_knmi_data(
            stn, meteo_var, start, end, settings, stn_name
        )
        meta = station_meta.to_dict()

        # set metadata
        meta.update(
            {
                "x": stations.loc[stn, "x"],
                "y": stations.loc[stn, "y"],
                "station": stn,
                "name": f"{meteo_var}_{stn_name}",
                "source": "KNMI",
            }
        )

        meta.update(variables)

    return knmi_df, meta


def get_stations(meteo_var):
    """get knmi stations from json files according to variable.

    Parameters
    ----------
    meteo_var : str, optional
        type of meteodata, by default 'RH'

    Returns
    -------
    pandas DataFrame with stations, names and coordinates (Lat/Lon & RD)
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))

    if meteo_var == "RD":
        # read precipitation station data only
        fname = "../data/knmi_neerslagstation.json"
        stations = pd.read_json(os.path.join(dir_path, fname))
    else:
        fname = "../data/knmi_meteostation.json"
        stations = pd.read_json(os.path.join(dir_path, fname))

    return stations


def get_station_name(stn, stations=None):
    """Returns the station name from a KNMI station.

    Modifies the station name in such a way that a valid url can be obtained.

    Parameters
    ----------
    stn : int
        station number.
    stations : pandas DataFrame
        DataFrame with the station metadata.

    Returns
    -------
    str or None
        Name of the station or None if station is not found.
    """
    if stations is None:
        stations = pd.concat([get_stations("RD"), get_stations("EV24")], axis=0)

    if stn not in stations.index:
        logger.warning("station {stn} not found")
        return None

    stn_name = stations.loc[stn, "naam"]
    stn_name = stn_name.upper().replace(" ", "-")
    stn_name = stn_name.replace("(", "").replace(")", "")

    return stn_name


def fill_missing_measurements(stn, meteo_var, start, end, settings, stn_name=None):
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
    variables : dictionary
        information about the observerd variables
    stations : pandas DataFrame
        information about the measurement station.
    """
    if settings["interval"] == "hourly":
        raise NotImplementedError("cannot yet fill missing values in hourly data")

    if (start is None) or (end is None):
        raise ValueError(
            "a start and end date have to be specified if fill_missing_obs is True"
        )

    # get the location of the stations
    stations = get_stations(meteo_var=meteo_var)
    if stn_name is None:
        stn_name = get_station_name(stn, stations)

    # check latest date at which measurements are available at De Bilt
    if (meteo_var in ["RD", "RH"]) and (
        end > (dt.datetime.now() - pd.Timedelta(LOOK_BACK_DAYS, unit="D"))
    ):
        end = min(
            end,
            _check_latest_measurement_date_RD_debilt(
                meteo_var, use_api=settings["use_api"]
            ),
        )
        logger.info(f'changing end_date to {end.strftime("%Y-%m-%d")}')

    # download data from station
    knmi_df, variables, station_meta = download_knmi_data(
        stn, meteo_var, start, end, settings, stn_name
    )

    # if the first station cannot be read, read another station as the first
    ignore = [stn]
    while knmi_df.empty:
        logger.info(f"station {stn} has no measurements between {start} and {end}")
        logger.info("trying to get measurements from nearest station")

        stn_lst = get_nearest_station_df(
            stations.loc[[stn]], meteo_var=meteo_var, ignore=ignore
        )
        if stn_lst is None:
            logger.warning(
                "there is no station with measurements of {meteo_var} between {start} and {end}"
            )
            return pd.DataFrame(), dict()

        stn = stn_lst[0]
        stn_name = get_station_name(stn, stations)
        knmi_df, variables, station_meta = download_knmi_data(
            stn, meteo_var, start, end, settings, stn_name
        )
        ignore.append(stn)

    # find missing values
    knmi_df = _add_missing_indices(knmi_df, stn, start, end)

    missing = knmi_df[meteo_var].isna()
    logger.info(f"station {stn} has {missing.sum()} missing measurements")

    knmi_df.loc[~missing, "station_opvulwaarde"] = str(stn)

    # fill missing values
    while np.any(missing) and not np.all(missing):
        stn_comp = get_nearest_station_df(
            stations.loc[[stn]], meteo_var=meteo_var, ignore=ignore
        )

        logger.info(
            f"trying to fill {missing.sum()} " f"measurements with station {stn_comp}"
        )

        if stn_comp is None:
            logger.info(
                "could not fill all missing measurements there are "
                "no stations left to check"
            )

            missing[:] = False
            break
        else:
            stn_comp = stn_comp[0]
        stn_name_comp = get_station_name(stn_comp, stations)
        knmi_df_comp, _, __ = download_knmi_data(
            stn_comp, meteo_var, start, end, settings, stn_name_comp
        )

        if knmi_df_comp.empty:
            logger.warning(f"station {stn_comp} cannot be downloaded")

        else:
            # dropnans from new data
            knmi_df_comp = knmi_df_comp.loc[~knmi_df_comp[meteo_var].isna(), :]
            # get index of missing data in original timeseries
            missing_idx = missing.loc[missing].index
            # if any missing are in the new data, update
            if missing_idx.isin(knmi_df_comp.index).any():
                # index for missing but in newly downloaded data
                ix_idx = missing_idx.intersection(knmi_df_comp.index)
                # update missing data
                knmi_df.loc[ix_idx, meteo_var] = knmi_df_comp.loc[ix_idx, meteo_var]
                # add source station number
                knmi_df.loc[ix_idx, "station_opvulwaarde"] = str(stn_comp)

        missing = knmi_df[meteo_var].isna()
        ignore.append(stn_comp)

    meta = station_meta.to_dict()

    # set metadata
    meta.update(
        {
            "x": stations.loc[stn, "x"],
            "y": stations.loc[stn, "y"],
            "station": stn,
            "name": f"{meteo_var}_{stn_name}",
            "source": "KNMI",
        }
    )

    meta.update(variables)

    return knmi_df, meta


def download_knmi_data(stn, meteo_var, start, end, settings, stn_name=None):
    """download knmi data of a measurements station for certain observation
    type.

    Parameters
    ----------
    stn : int
        measurement station.
    meteo_var : str
        observation type.
    start : pd.TimeStamp or None
        start date of observations.
    end : pd.TimeStamp or None
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
    knmi_df : pandas DataFrame
        data from one station from one type of observation
    variables : dictionary
        information about the observerd variables
    stations : pandas DataFrame
        information about the measurement station.
    """

    logger.info(
        f"download knmi {meteo_var} data from station "
        f"{stn}-{stn_name} between {start} and {end}"
    )

    # define variables
    knmi_df = pd.DataFrame()
    variables = {}
    stations = pd.DataFrame()

    # download and read data
    try:
        try:
            if settings["use_api"]:
                if settings["interval"].startswith("hour"):
                    # hourly data from meteorological stations
                    knmi_df, variables = get_knmi_hourly_api(stn, meteo_var, start, end)
                elif meteo_var == "RD":
                    # daily data from rainfall-stations
                    knmi_df, variables = get_knmi_daily_rainfall_api(stn, start, end)
                else:
                    # daily data from meteorological stations
                    knmi_df, variables, stations = get_knmi_daily_meteo_api(
                        stn, meteo_var, start, end
                    )

        except (RuntimeError, requests.ConnectionError) as e:
            logger.warning("KNMI API failed")
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
                knmi_df, variables = get_knmi_daily_rainfall_url(
                    stn, stn_name, start, end
                )
            else:
                # daily data from meteorological stations
                knmi_df, variables = get_knmi_daily_meteo_url(
                    stn, meteo_var, start, end
                )
    except (ValueError, KeyError) as e:
        logger.error(e)
        if settings["raise_exceptions"]:
            raise ValueError(e)

    return knmi_df, variables, stations


@lru_cache()
def get_knmi_daily_rainfall_api(stn, start=None, end=None):
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
        params["start"] = start.strftime("%Y%m%d")
    if end is not None:
        params["end"] = end.strftime("%Y%m%d")

    result = requests.get(url, params=params)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)
    knmi_df, variables = read_knmi_daily_rainfall(f, meteo_var)
    if stn not in knmi_df.STN.unique():
        return pd.DataFrame(), variables

    return knmi_df[[meteo_var]], variables


@lru_cache()
def get_knmi_daily_rainfall_url(stn, stn_name, start=None, end=None):
    """download and read knmi daily rainfall.

    Parameters
    ----------
    stn : int
        station number.
    stn_name : str
        the name of the station in capital letters, can be tricky
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
    if end is not None:
        end = end + dt.timedelta(days=1)

    stn = str(stn)
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie/"
        f"gegevens/monv_reeksen/neerslaggeg_{stn_name}_{stn}.zip"
    )

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), "knmi")
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(), "knmi", f"neerslaggeg_{stn}.zip")
    fname_dir = os.path.join(tempfile.gettempdir(), "knmi", f"neerslaggeg_{stn}")
    fname_txt = os.path.join(fname_dir, f"neerslaggeg_{stn_name}_{stn}.txt")

    # download zip file
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        raise ValueError(f"invalid url {url} please check station name {stn_name}")
    with open(fname_zip, "wb") as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    # unzip file
    util.unzip_file(fname_zip, fname_dir, force=True, preserve_datetime=True)

    return read_knmi_daily_rainfall_file(fname_txt, start=start, end=end)


def read_knmi_daily_rainfall_file(fname_txt, start=None, end=None):
    """read a knmi file with daily rainfall data.

    Parameters
    ----------
    fname_txt : str
        file path of a knmi .txt file.
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

    """
    meteo_var = "RD"

    with open(fname_txt, "r") as f:
        line = f.readline()
        # get meteo var
        for _ in range(50):
            if "RD   " in line:  # RD komt te vaak voor vandaar de spaties
                key, item = line.split("=")
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        # get dataframe
        for _ in range(50):
            if "STN" in line:
                columns = line.strip("# ").strip("\n").split(",")
                columns = [x.strip(" ") for x in columns]
                values = f.readline()
                if values == "\n":
                    values = f.readline()
            line = f.readline()

        df = pd.read_csv(f, header=None, names=columns, na_values="     ")

        df.set_index(pd.to_datetime(df.YYYYMMDD, format="%Y%m%d"), inplace=True)
        df = df.drop("YYYYMMDD", axis=1)

        if df.index.duplicated().sum() > 0:
            df = df.loc[~df.index.duplicated(keep="first")]
            logger.info("duplicate indices removed from RD measurements")

        # sometimes the last row is empty, check for that and remove it
        if not df.empty:
            if df.iloc[-1].isna().any():
                logger.debug("duplicate indices removed from RD measurements")
                df = df.drop(index=df.index[-1])
                df.loc[:, meteo_var] = df[meteo_var].astype(float)

        df, variables = _transform_variables(df, variables)
        variables["unit"] = "m"

        # add station to variables
        if len(df["STN"].unique()) != 1:
            raise ValueError("multiple stations in single file")
        else:
            variables["station"] = df["STN"].iloc[0]

    return df.loc[start:end, [meteo_var]], variables


def _read_knmi_header(f, meteo_var):
    variables = {}
    line = f.readline()
    for iline in range(500):
        if " = " in line or " : " in line:
            line = line.lstrip(" #").strip("\n")
            if " = " in line:
                varDes = line.split(" = ")
            else:
                varDes = line.split(" : ")
            if varDes[0].strip() == meteo_var:
                variables[varDes[0].strip()] = varDes[1].strip()

        if "STN,YY" in line:
            header = line.replace("#", "").split(",")
            header = [item.lstrip().rstrip() for item in header]
            break

        line = f.readline()

    if iline > 498:
        raise ValueError("cannot read measurements from file")

    return f, variables, header


def _transform_variables(df, variables):
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
        if " millimeters" in value:
            logger.debug(f"transform {key}, {value} from mm to m")

            df[key] = df[key] * 0.001
            value = value.replace(" millimeters", " m")
        if "08.00 UTC" in value:
            logger.debug(f"transform {key}, {value} from UTC to UTC+1")

            # over the period 08.00 preceding day - 08.00 UTC present day
            df.index = df.index + pd.to_timedelta(8, unit="h")

            # from UT to UT+1 (standard-time in the Netherlands)
            df.index = df.index + pd.to_timedelta(1, unit="h")

            value = value.replace("08.00", "09.00").replace("UTC", "UTC+1")

        # Store new variable
        variables[key] = value

    return df, variables


def read_knmi_daily_rainfall(f, meteo_var):
    """
    Read daily rainfall data from a KNMI file.

    Parameters
    ----------
    f : file-like object
        The file object containing the KNMI data.
    meteo_var : str
        The meteorological variable to extract.

    Returns
    -------
    pandas.DataFrame
        The DataFrame containing the extracted daily rainfall data.
    dict
        A dictionary containing information about the variables in the DataFrame.

    Notes
    -----
    This function assumes that the file object `f` is already open and positioned at the start of the data.
    The file is expected to have a header with variable names and a corresponding data table.

    The DataFrame returned by this function has the following modifications:
    - The index is set to the datetime values derived from the 'YYYYMMDD' column.
    - The 'YYYYMMDD' column is dropped.
    - Duplicate indices are removed, keeping the first occurrence.
    - If the last row has missing values, it is removed.
    - The 'meteo_var' column is cast to float data type.
    - Variables are transformed using an internal function `_transform_variables`.
    - The unit of measurement for the 'meteo_var' variable is set to 'm'.

    """
    f, variables, header = _read_knmi_header(f, meteo_var)

    df = pd.read_csv(f, header=None, names=header, na_values="     ")
    f.close()

    df.set_index(pd.to_datetime(df.YYYYMMDD, format="%Y%m%d"), inplace=True)
    df = df.drop("YYYYMMDD", axis=1)

    if df.index.duplicated().sum() > 0:
        df = df.loc[~df.index.duplicated(keep="first")]
        logger.debug("duplicate indices removed from RD measurements")

    # sometimes the last row is messed up, check for that and remove it
    if not df.empty:
        if df.iloc[-1].isna().any():
            logger.debug("last row contains no data, remove last row")

            df = df.drop(index=df.index[-1])
            df.loc[:, meteo_var] = df[meteo_var].astype(float)

    df, variables = _transform_variables(df, variables)
    variables["unit"] = "m"

    return df, variables


def _read_station_location(f):
    stations = None

    line = f.readline()
    for _ in range(30):
        if "STN" in line:
            titels = line.strip("# ").split()
            titels = [x.replace("(", "_") for x in titels]
            titels = [x.replace(r")", "") for x in titels]

            values = f.readline().strip("# ").strip().replace(":", "")
            values = re.split(r"\s{2,}", values)

            # Create pd.DataFrame for station data
            stations = pd.DataFrame(columns=titels, data=[values])
            stations.set_index(["STN"], inplace=True)
            for col in stations.columns:
                try:
                    stations[col] = stations[col].astype(float)
                except ValueError:
                    pass

            break

        line = f.readline()

    if stations is None:
        logger.warning("could not find stations")

    return f, stations


# @lru_cache()
def get_knmi_daily_meteo_api(stn, meteo_var, start=None, end=None):
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
        "vars": meteo_var,
        "stns": str(stn),
    }

    # modify start and end date for the api call because the dates are later
    # transformed, see example notebook 02_knmi_observations.
    if start is not None:
        start = start - dt.timedelta(days=1)
        params["start"] = start.strftime("%Y%m%d")
    if end is not None:
        end = end - dt.timedelta(days=1)
        params["end"] = end.strftime("%Y%m%d")

    result = requests.get(url, params=params)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)
    knmi_df, variables, stations = read_knmi_daily_meteo(f, meteo_var)

    knmi_df.dropna(subset=[meteo_var], inplace=True)

    return knmi_df[[meteo_var]], variables, stations


@lru_cache()
def get_knmi_daily_meteo_url(stn, meteo_var, start=None, end=None):
    """download and read knmi daily meteo data.

    Parameters
    ----------
    stn : str
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
    url = (
        "https://cdn.knmi.nl/knmi/map/page/klimatologie"
        f"/gegevens/daggegevens/etmgeg_{stn}.zip"
    )

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), "knmi")
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(), "knmi", f"etmgeg_{stn}.zip")
    fname_dir = os.path.join(tempfile.gettempdir(), "knmi", f"etmgeg_{stn}")
    fname_txt = os.path.join(fname_dir, f"etmgeg_{stn}.txt")

    # download zip file
    r = requests.get(url, stream=True)
    with open(fname_zip, "wb") as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    # unzip file
    util.unzip_file(fname_zip, fname_dir, force=True, preserve_datetime=True)

    return read_knmi_daily_meteo_file(fname_txt, meteo_var, start, end)


def read_knmi_daily_meteo_file(fname, meteo_var, start=None, end=None):
    """read knmi daily meteo data from a file

    Parameters
    ----------
    fname : str
        file path.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp or None
        start time of observations.
    end : pd.TimeStamp or None
        end time of observations.

    Raises
    ------
    ValueError
        If the meteo var is not in the file.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    variables = None
    with open(fname, "r") as f:
        line = f.readline()
        # get meteo var
        for _ in range(50):
            if meteo_var in line:
                key, item = line.split("=")
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        if variables is None:
            raise ValueError(f"could not find {meteo_var} in file {fname}")

        # get dataframe
        for _ in range(50):
            if "STN" in line:
                columns = line.strip("# ").strip("\n").split(",")
                columns = [x.strip(" ") for x in columns]
                values = f.readline()
                if values == "\n":
                    values = f.readline()
                df = pd.read_csv(f, header=None, names=columns, na_values="     ")
                df.set_index(pd.to_datetime(df.YYYYMMDD, format="%Y%m%d"), inplace=True)
                df = df.drop("YYYYMMDD", axis=1)

                df = df.loc[df.index.notnull(), :]

                # add a full day for meteorological data, so that the
                # timestamp is at the end of the period in the data
                df.index = df.index + pd.to_timedelta(1, unit="d")

                # from UT to UT+1 (standard-time in the Netherlands)
                df.index = df.index + pd.to_timedelta(1, unit="h")

                # add station to variables
                if len(df["STN"].unique()) != 1:
                    raise ValueError("multiple stations in single file")
                else:
                    station = df["STN"].iloc[0]

                df = df.loc[start:end, [meteo_var]]
                df = df.dropna()
                df, variables = _transform_variables(df, variables)
                variables["unit"] = ""
                variables["station"] = station
                break

            line = f.readline()

    return df, variables


def read_knmi_daily_meteo(f, meteo_var):
    f, stations = _read_station_location(f)
    f, variables, header = _read_knmi_header(f, meteo_var)
    header[0] = header[0].lstrip("# ")
    df = pd.read_csv(f, header=None, names=header, na_values="     ")
    f.close()

    df.set_index(pd.to_datetime(df.YYYYMMDD, format="%Y%m%d"), inplace=True)
    df = df.drop("YYYYMMDD", axis=1)

    df = df.loc[df.index.notnull(), :]

    # add a full day for meteorological data, so that the
    # timestamp is at the end of the period in the data
    df.index = df.index + pd.to_timedelta(1, unit="d")

    # from UT to UT+1 (standard-time in the Netherlands)
    df.index = df.index + pd.to_timedelta(1, unit="h")

    df, variables = _transform_variables(df, variables)
    variables["unit"] = ""

    return df, variables, stations


# @lru_cache()
def get_knmi_hourly_api(stn, meteo_var, start=None, end=None):
    """Retrieve hourly meteorological data from the KNMI API.

    Parameters
    ----------
    stn : int
        The station number.
    meteo_var : str
        The meteorological variable to retrieve.
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
        "vars": meteo_var,
        "stns": str(stn),
    }

    # It looks like the API has an option to select the start and end hour that
    # you want but in reality it does not have this option.
    if start is not None:
        s = start - pd.Timedelta(1, unit="H")
        params["start"] = s.strftime("%Y%m%d") + "01"
    if end is not None:
        params["end"] = end.strftime("%Y%m%d") + "24"

    result = requests.get(url, params=params)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)

    return read_knmi_hourly(f, meteo_var, start, end)


def read_knmi_hourly(f, meteo_var, start=None, end=None):
    """Read hourly KNMI file.

    Parameters
    ----------
    f : str or filelike object
        path to file or filelike object

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing data
    variables : dict
        dictionary containing metadata about the variables
    """

    ispath = False
    if isinstance(f, str):
        f = open(f, "r")
        ispath = True

    f, variables, header = _read_knmi_header(f, meteo_var)
    df = pd.read_csv(f, header=None, names=header, na_values="     ")

    if ispath:
        f.close()

    # convert 24 to 0
    if "H" in df.columns:
        hour_col = "H"
    elif "HH" in df.columns:
        hour_col = "HH"
    else:
        raise ValueError("could not find column with hour")

    df.loc[df[hour_col] == 24, hour_col] = 0
    datetime = pd.to_datetime(
        df.YYYYMMDD.astype(str) + df[hour_col].astype(str).str.zfill(2),
        format="%Y%m%d%H",
    )
    datetime.loc[datetime.dt.hour == 0] = datetime + dt.timedelta(days=1)

    # add station to variables
    if len(df["STN"].unique()) != 1:
        raise ValueError("multiple stations in single file")
    else:
        station = df["STN"].iloc[0]

    # set index
    df.set_index(datetime, inplace=True)

    df = df[[meteo_var]]
    variables = variables

    df, variables = _transform_variables(df, variables)
    variables["unit"] = "m"
    variables["station"] = station

    return df.loc[start:end, [meteo_var]], variables


def _check_latest_measurement_date_RD_debilt(meteo_var, use_api=True):
    """According to the website of the knmi it can take up to 3 weeks before
    precipitation data is updated. If you use the fill_missing_measurements
    method to fill a time series untill today, it will keep looking at all
    stations to find the data of these last days that probably does not exist.

    To prevent this, this function will find the last day there are measure-
    ments at knmi station de Bilt. It is assumed that if de Bilt doesn't have
    recent measurements, no station will have measurements for these dates.

    website knmi: https://www.knmi.nl/nederland-nu/klimatologie/monv/reeksen

    Parameters
    ----------
    meteo_var : str
        meteo variable.
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        Default is True.

    Returns
    -------
    last_measurement_date_debilt : pd.TimeStamp
        last date with measurements at station de Bilt
    """

    start = dt.datetime.now() - pd.Timedelta(LOOK_BACK_DAYS, unit="D")
    end = dt.datetime.now() + pd.Timedelta(10, unit="D")
    if meteo_var == "RD":
        if use_api:
            try:
                knmi_df, _ = get_knmi_daily_rainfall_api(550, start, end)
            except (RuntimeError, requests.ConnectionError):
                logger.warning("KNMI API failed, switching to non-API method")
                knmi_df, _ = get_knmi_daily_rainfall_url(550, "DE-BILT", start, end)
        else:
            knmi_df, _ = get_knmi_daily_rainfall_url(550, "DE-BILT", start, end)
    else:
        if use_api:
            try:
                knmi_df, _, _ = get_knmi_daily_meteo_api(
                    260, meteo_var, start=start, end=end
                )
            except (RuntimeError, requests.ConnectionError):
                logger.warning("KNMI API failed, switching to non-API method")
                knmi_df, _, _ = get_knmi_daily_meteo_url(260, meteo_var, start, end)
        else:
            knmi_df, _, _ = get_knmi_daily_meteo_url(260, meteo_var, start, end)

    knmi_df = knmi_df.dropna()
    if knmi_df.empty:
        raise ValueError(
            "knmi station de Bilt has no measurements "
            f"in the past {LOOK_BACK_DAYS} days for variable {meteo_var}."
        )

    last_measurement_date_debilt = knmi_df.index[-1]

    logger.debug(
        f"last {meteo_var} measurement available at the Bilt is from"
        f' {last_measurement_date_debilt.strftime("%Y-%m-%d")}'
    )
    logger.debug(
        f"assuming no {meteo_var} measurements are available at "
        "other stations after this date"
    )

    return last_measurement_date_debilt


def get_nearest_station_df(
    locations, xcol="x", ycol="y", stations=None, meteo_var="RH", ignore=None
):
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
    ignore : list, optional
        list of stations to ignore. The default is None.

    Returns
    -------
    stns : list
        station numbers.
    """
    if stations is None:
        stations = get_stations(meteo_var=meteo_var)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
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


def get_nearest_station_xy(xy, stations=None, meteo_var="RH", ignore=None):
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


def get_n_nearest_stations_xy(xy, meteo_var, n=1, stations=None, ignore=None):
    """Find the N nearest KNMI stations that measure variable 'meteo_var' to
    the x, y coordinates.

    Parameters
    ----------
    xy : list, tuple or numpy.array of int or float
        sinlge pair of xy coordinates. e.g. (150_000., 400_000.)
    meteo_var : str
        measurement variable e.g. 'RH' or 'EV24'
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
        stations = get_stations(meteo_var=meteo_var)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
        if stations.empty:
            return None

    distance = np.sqrt((stations.x - xy[0]) ** 2 + (stations.y - xy[1]) ** 2)
    stns = distance.nsmallest(n).index.to_list()

    return stns


def _add_missing_indices(knmi_df, stn, start, end):
    """When downloading KNMI data you don't always get a DataFrame with the
    periods that you provided in your request. Thus the index does not cover
    the complete period that you are interested in. This function adds the
    missing period to the index of the DataFrame.

    Parameters
    ----------
    knmi_df : pandas DataFrame
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
    knmi_df : pandas DataFrame
        data from one station from one type of observation
    """
    # check if given dates are more or less similar than measurement dates
    if (knmi_df.index[0] - start).days < 2:
        new_start = knmi_df.index[0]
    else:
        new_start = pd.Timestamp(
            year=start.year,
            month=start.month,
            day=start.day,
            hour=knmi_df.index[0].hour,
            minute=knmi_df.index[0].minute,
            second=knmi_df.index[0].second,
        )
        logger.info(f"station {stn} has no measurements before {knmi_df.index[0]}")

    if (end - knmi_df.index[-1]).days < 2:
        new_end = knmi_df.index[-1]
    else:
        new_end = pd.Timestamp(
            year=end.year,
            month=end.month,
            day=end.day,
            hour=knmi_df.index[-1].hour,
            minute=knmi_df.index[-1].minute,
            second=knmi_df.index[-1].second,
        )
        logger.info(f"station {stn} has no measurements after {knmi_df.index[-1]}")

    # add missing indices
    new_index = pd.date_range(new_start, new_end, freq="D")
    knmi_df = knmi_df.reindex(new_index)

    return knmi_df


def get_knmi_obslist(
    locations=None,
    stns=None,
    xy=None,
    meteo_vars=("RH",),
    starts=None,
    ends=None,
    ObsClasses=None,
    **kwargs,
):
    """Get a list of observations of knmi stations. Either specify a list of
    knmi stations (stns) or a dataframe with x, y coordinates (locations).

    Parameters
    ----------
    locations : pandas DataFrame or None
        dataframe with x and y coordinates. The default is None
    stns : list of str or None
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
            The default is False.
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
    obs_list : list of obsevation objects
        collection of multiple point observations
    """

    settings = _get_default_settings(kwargs)

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
        start, end = _start_end_to_datetime(start, end)

        # get stations
        if stns is None:
            stations = get_stations(meteo_var=meteo_var)
            if (locations is None) and (xy is not None):
                _stns = get_nearest_station_xy(
                    xy, stations=stations, meteo_var=meteo_var
                )
            elif locations is not None:
                _stns = get_nearest_station_df(
                    locations, stations=stations, meteo_var=meteo_var
                )
            else:
                raise ValueError(
                    "stns, location and x are all None" "please specify one of these"
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
