"""Module with functions to read time series with observations from knmi.

The main function to obtain a time series are:

- get_knmi_timeseries_xy: obtain a knmi station based on the xy location
- get_knmi_timeseries_stn: obtain a knmi station based on the station number

There is some ugly code involved if you want to obtain precipitation data.
This code combines the data from the meteo_stations and the neerslagstations.
"""
from functools import lru_cache
import datetime as dt
import logging
import os
import re
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
from scipy.interpolate import RBFInterpolator
import requests

from .. import util

logger = logging.getLogger(__name__)

URL_DAILY_NEERSLAG = "https://www.daggegevens.knmi.nl/klimatologie/monv/reeksen"
URL_DAILY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/daggegevens"
URL_HOURLY_METEO = "https://www.daggegevens.knmi.nl/klimatologie/uurgegevens"
LOOK_BACK_DAYS = 365


def get_stations(meteo_var="RH"):
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

    fname = "../data/knmi_meteostation.json"

    stations = pd.read_json(os.path.join(dir_path, fname))

    if meteo_var == "RD":
        # read preciptiation station data only
        fname = "../data/knmi_neerslagstation.json"
        stations = pd.read_json(os.path.join(dir_path, fname))

    if meteo_var == "PG":
        # in Ell wordt geen luchtdruk gemeten
        stations.drop(377, inplace=True)

    elif meteo_var == "EV24":
        # in Woensdrecht wordt geen globale straling gemeten
        stations.drop(340, inplace=True)
        # op Vlieland wordt geen globale straling gemeten
        stations.drop(265, inplace=True)

    return stations


def get_station_name(stn, stations):
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
    str
        Name of the station.
    """

    stn_name = stations.loc[stn, "naam"]
    stn_name = stn_name.upper().replace(" ", "-")
    stn_name = stn_name.replace("(", "").replace(")", "")

    return stn_name


def get_nearest_stations_xy(xy, meteo_var, n=1, stations=None, ignore=None):
    """find the KNMI stations that measure 'variable' closest to the x, y
    coordinates.

    Parameters
    ----------
    xy : list. tuple or numpy array of int or floats
        xy coördinates. e.g.(0.1, 25.05)
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


def get_nearest_station_df(
    locations, xcol="x", ycol="y", stations=None, meteo_var="RH", ignore=None
):
    """find the KNMI stations that measure 'meteo_var' closest to the
    coordinates in 'locations'.

    Parameters
    ----------
    locations : pandas DataFrame
        x and y coordinates
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
    x and y coördinates.

    Parameters
    ----------
    xy : list or numpy array, optional
        xy coördinates of the locations. e.g. [[10,25], [5,25]]
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
    start : pd.TimeStamp
        start time
    end : pd.TimeStamp
        end time
    """

    if start is None:
        start = pd.Timestamp(pd.Timestamp.today().year - 1, 1, 1)
    else:
        start = pd.to_datetime(start)

    if end is None:
        end = pd.Timestamp.today() - pd.Timedelta(1, unit="D")
    else:
        end = pd.to_datetime(end)

    return start, end


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
        Default is False as the API is down (apr 2021).

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
                knmi_df, _ = get_knmi_daily_rainfall_api(
                    URL_DAILY_NEERSLAG, 550, "RD", start=start, end=end, inseason=False
                )
            except (RuntimeError, requests.ConnectionError):
                logger.info("KNMI API failed, switching to non-API method")
                knmi_df, _ = get_knmi_daily_rainfall_url(
                    550, "DE-BILT", "RD", start, end, inseason=False
                )
        else:
            knmi_df, _ = get_knmi_daily_rainfall_url(
                550, "DE-BILT", "RD", start, end, inseason=False
            )
    else:
        if use_api:
            try:
                knmi_df, _, _ = get_knmi_daily_meteo_api(
                    URL_DAILY_METEO,
                    260,
                    meteo_var,
                    start=start,
                    end=end,
                    inseason=False,
                )
            except (RuntimeError, requests.ConnectionError):
                logger.info("KNMI API failed, switching to non-API method")
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

    logger.info(
        f"last {meteo_var} measurement available at the Bilt is from"
        f' {last_measurement_date_debilt.strftime("%Y-%m-%d")}'
    )
    logger.info(
        f"assuming no {meteo_var} measurements are available at "
        "other stations after this date"
    )

    return last_measurement_date_debilt


def download_knmi_data(
    stn, stn_name=None, meteo_var="RH", start=None, end=None, settings=None
):
    """download knmi data of a measurements station for certain observation
    type.

    Parameters
    ----------
    stn : int or str
        number of measurement station
    stn_name : str
        Only necesary if meteo_var is RD and use_api is True.
    meteo_var : str, optional
        measurement type 'RH' or 'EV24'. The default is 'RH'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    settings

    Raises
    ------
    NotImplementedError
        different time intervals and inseason data is not yet working.
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

    # checks
    if settings["inseason"]:
        raise NotImplementedError("inseason stuff not implemented")

    start, end = _start_end_to_datetime(start, end)

    # raise error if hourly neerslag station data is requested
    if (meteo_var == "RD") and settings["interval"].startswith("hour"):
        message = 'No hourly precipitation data available at precipitation station, use meteo_var "RH" to obtain hourly precipitation data from meteo stations'
        raise (ValueError(message))

    # convert stn to string
    stn = str(stn)

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
                    end = end - dt.timedelta(days=1)
                    knmi_df, variables = get_knmi_hourly_api(
                        URL_HOURLY_METEO, stn, meteo_var, start, end
                    )
                elif meteo_var == "RD":
                    # daily data from rainfall-stations
                    knmi_df, variables = get_knmi_daily_rainfall_api(
                        URL_DAILY_NEERSLAG,
                        stn,
                        meteo_var,
                        start,
                        end,
                        settings["inseason"],
                    )
                else:
                    start = start - dt.timedelta(days=1)
                    end = end - dt.timedelta(days=1)
                    # daily data from meteorological stations
                    knmi_df, variables, stations = get_knmi_daily_meteo_api(
                        URL_DAILY_METEO,
                        stn,
                        meteo_var,
                        start,
                        end,
                        settings["inseason"],
                    )
        except (RuntimeError, requests.ConnectionError):
            settings["use_api"] = False
            logger.info("KNMI API failed, switching to non-API method")

        if not settings["use_api"]:
            end = end + dt.timedelta(days=1)
            if settings["interval"].startswith("hour"):
                # hourly data from meteorological stations
                raise NotImplementedError()
            elif meteo_var == "RD":
                # daily data from rainfall-stations
                knmi_df, variables = get_knmi_daily_rainfall_url(
                    stn, stn_name, meteo_var, start, end, settings["inseason"]
                )
            else:
                # daily data from meteorological stations
                knmi_df, variables, stations = get_knmi_daily_meteo_url(
                    stn, meteo_var, start, end, settings["inseason"]
                )
    except (ValueError, KeyError) as e:
        logger.error(e)
        if settings["raise_exceptions"]:
            raise ValueError(e)

    return knmi_df, variables, stations


@lru_cache()
def get_knmi_daily_rainfall_api(
    url, stn, meteo_var, start=None, end=None, inseason=False
):
    """download and read knmi daily rainfall.

    Parameters
    ----------
    url : str
        download url.
    stn : str
        station number.
    meteo_var : str
        must be 'RD'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    inseason : boolean
        flag to obtain inseason data.

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

    data = {"inseason": str(int(inseason)), "vars": meteo_var, "stns": stn}

    if start is not None:
        data["start"] = start.strftime("%Y%m%d")

    if end is not None:
        data["end"] = end.strftime("%Y%m%d")

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)
    knmi_df, variables = read_knmi_daily_rainfall(f, meteo_var)
    if int(stn) not in knmi_df.STN.unique():
        raise ValueError(
            f"KNMI station {stn} does not exists or has no "
            "measurements in the given period"
        )

    return knmi_df[[meteo_var]], variables


@lru_cache()
def get_knmi_daily_rainfall_url(
    stn, stn_name, meteo_var, start=None, end=None, inseason=False, use_cache=True
):
    """download and read knmi daily rainfall.

    Parameters
    ----------
    stn : str
        station number.
    stn_name : str
        the name of the station in capital letters, can be tricky
    meteo_var : str
        must be 'RD'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    inseason : boolean
        flag to obtain inseason data.

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
    stn = f"{int(stn):03}"
    url = f"https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/monv_reeksen/neerslaggeg_{stn_name}_{stn}.zip"

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), "knmi")
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(), "knmi", f"neerslaggeg_{stn}.zip")
    fname_dir = os.path.join(tempfile.gettempdir(), "knmi", f"neerslaggeg_{stn}")
    fname_txt = os.path.join(fname_dir, f"neerslaggeg_{stn_name}_{stn}.txt")

    # check if file should be downloaded and unzipped
    donwload_and_unzip = True
    if os.path.exists(fname_txt) and use_cache:
        fname_time_mod = dt.datetime.fromtimestamp(os.path.getmtime(fname_txt))
        if fname_time_mod > end:
            donwload_and_unzip = False

    if donwload_and_unzip:
        # download zip file
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            raise ValueError(f"invalid url {url} please check station name {stn_name}")
        with open(fname_zip, "wb") as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        # unzip file
        util.unzip_file(fname_zip, fname_dir, force=True, preserve_datetime=True)

    with open(fname_txt, "r") as f:
        line = f.readline()
        # get meteo var
        for iline in range(50):
            if "RD   " in line:  # RD komt te vaak voor vandaar de spaties
                key, item = line.split("=")
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        # get dataframe
        for iline in range(50):
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

    return df.loc[start:end, [meteo_var]], variables


def _read_knmi_header(f):

    variables = {}
    line = f.readline()
    if "DOCTYPE HTML PUBLIC" in line:
        logger.error(f.read())
        raise ValueError("Internal Server Error")
    for iline in range(500):
        if " = " in line or " : " in line:
            line = line.lstrip(" #").strip("\n")
            if " = " in line:
                varDes = line.split(" = ")
            else:
                varDes = line.split(" : ")
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
            logger.info(f"transform {key}, {value} lower than 0.05 mm to 0 mm")

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

    f, variables, header = _read_knmi_header(f)

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

    return df, variables


def _read_station_location(f):

    stations = None

    line = f.readline()
    for iline in range(30):
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
                    stations.loc[:, col] = stations[col].astype(float)
                except ValueError:
                    pass

            break

        line = f.readline()

    if stations is None:
        logger.warning("could not find stations")

    return f, stations


@lru_cache()
def get_knmi_daily_meteo_api(url, stn, meteo_var, start, end, inseason):
    """download and read knmi daily meteo data.

    Parameters
    ----------
    url : str
        download url.
    stn : str
        station number.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    inseason : boolean
        flag to obtain inseason data.

    Returns
    -------
    pandas DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    stations : pandas DataFrame
        additional data about the measurement station
    """
    data = {
        "start": start.strftime("%Y%m%d"),
        "end": end.strftime("%Y%m%d"),
        "inseason": str(int(inseason)),
        "vars": meteo_var,
        "stns": stn,
    }

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)
    knmi_df, variables, stations = read_knmi_daily_meteo(f)

    knmi_df.dropna(subset=[meteo_var], inplace=True)

    return knmi_df[[meteo_var]], variables, stations


@lru_cache()
def get_knmi_daily_meteo_url(stn, meteo_var, start, end, use_cache=True):
    """download and read knmi daily meteo data.

    Parameters
    ----------
    url : str
        download url.
    stn : str
        station number.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
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
    url = f"https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/daggegevens/etmgeg_{stn}.zip"

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), "knmi")
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(), "knmi", f"etmgeg_{stn}.zip")
    fname_dir = os.path.join(tempfile.gettempdir(), "knmi", f"etmgeg_{stn}")
    fname_txt = os.path.join(fname_dir, f"etmgeg_{stn}.txt")

    # check if file should be downloaded and unzipped
    donwload_and_unzip = True
    if os.path.exists(fname_txt) and use_cache:
        fname_time_mod = dt.datetime.fromtimestamp(os.path.getmtime(fname_txt))
        if fname_time_mod > end:
            donwload_and_unzip = False

    if donwload_and_unzip:
        # download zip file
        r = requests.get(url, stream=True)
        with open(fname_zip, "wb") as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        # unzip file
        util.unzip_file(fname_zip, fname_dir, force=True, preserve_datetime=True)

    variables = None
    with open(fname_txt, "r") as f:
        line = f.readline()
        # get meteo var
        for iline in range(50):
            if meteo_var in line:
                key, item = line.split("=")
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        if variables is None:
            raise ValueError(f"could not find {meteo_var} for station {stn}")

        # get dataframe
        for iline in range(50):
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
                df = df.loc[start:end, [meteo_var]]
                df = df.dropna()
                df, variables = _transform_variables(df, variables)
                break

            line = f.readline()

    return df, variables, None


def read_knmi_daily_meteo(f):

    f, stations = _read_station_location(f)
    f, variables, header = _read_knmi_header(f)
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

    return df, variables, stations


@lru_cache()
def get_knmi_hourly_api(url, stn, meteo_var, start, end):

    data = {
        "start": start.strftime("%Y%m%d") + "01",
        "end": end.strftime("%Y%m%d") + "24",
        "vars": meteo_var,
        "stns": stn,
    }

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    result_str = result.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("KNMI API down")

    f = StringIO(result_str)
    df, variables = read_knmi_hourly(f)

    return df[[meteo_var]], variables


def read_knmi_hourly(f):

    f, variables, header = _read_knmi_header(f)
    df = pd.read_csv(f, header=None, names=header, na_values="     ")

    # convert 24 to 0
    df.loc[df["H"] == 24, "H"] = 0
    datetime = pd.to_datetime(
        df.YYYYMMDD.astype(str) + df.H.astype(str).str.zfill(2), format="%Y%m%d%H"
    )
    datetime.loc[datetime.dt.hour == 0] = datetime + dt.timedelta(days=1)

    # set index
    df.set_index(datetime, inplace=True)

    df = df.drop(["YYYYMMDD", "H"], axis=1)

    df, variables = _transform_variables(df, variables)

    return df, variables


def get_knmi_timeseries_xy(xy, meteo_var, start, end, stn_name=None, settings=None):
    """Get timeseries with measurements from station closest to the given (x,y)
    coördinates.

    Parameters
    ----------
    xy : list. tuple or numpy array of int or floats
        xy coördinates. e.g. [10.1,25.05]
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    stn_name : str
        name of the KNMI station
    settings : dict or None
        settings to obtain data

    Returns
    -------
    knmi_df : pandas DataFrame
        time series.
    meta : dic
        dictionary with metadata.
    """
    settings = _get_default_settings(settings)

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn = get_nearest_stations_xy(xy, meteo_var, stations=stations)[0]

    knmi_df, meta = get_knmi_timeseries_stn(
        stn, meteo_var, start, end, settings=settings
    )

    return knmi_df, meta


def _get_default_settings(settings):
    """adds the default settings to a dictinary with settings. If settings
    is None all the settings are default. If there are already settings given
    only the non-existing settings are added with their default value.

    The default settings are:
    fill_missing_obs = True
        nan values in time series are filled with nearby time series.
    interval = 'daily'
        desired time interval for observations. Can be 'daily' or 'hourly'.
        'hourly' is only for precipitation ('RH') data from meteo stations.
    inseason = False
        flag to obtain inseason data. Not implemented very well yet.
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        if False a text file is downloaded into a temporary folder and the data is read from there.
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
        "fill_missing_obs": True,
        "interval": "daily",
        "inseason": False,
        "use_api": True,
        "raise_exceptions": True,
        "normalize_index": True,
    }

    if settings is None:
        settings = {}

    for key, value in default_settings.items():
        if key not in settings.keys():
            settings[key] = value

    return settings


def get_knmi_timeseries_stn(stn, meteo_var, start, end, settings=None):
    """Get a knmi time series and metadata.

    Parameters
    ----------
    stn : int or str
        measurement station e.g. 829.
    meteo_var : str, optional
        observation type e.g. "RH" or "EV24". See list with all options below.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    settings : dict or None, optional
        settings for obtaining the right time series, options are:
            fill_missing_obs : bool, optional
                if True nan values in time series are filled with nearby time series.
                The default is True.
            interval : str, optional
                desired time interval for observations. The default is 'daily'.
            inseason : boolean, optional
                flag to obtain inseason data. The default is False
            use_api : bool, optional
                if True the api is used to obtain the data, API documentation is here:
                    https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
                Default is True as the API since (July 2021).
            raise_exceptions : bool, optional
                if True you get errors when no data is returned. The default is True.

    List of possible variables:
        neerslagstations:
        RD    = de 24-uurs neerslagsom, gemeten van 0800 utc op de
        voorafgaande dag tot 0800 utc op de vermelde datum

        meteostations:
        DDVEC = Vectorgemiddelde windrichting in graden (360=noord,
        90=oost, 180=zuid, 270=west, 0=windstil/variabel). Zie
        http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
        / Vector mean wind direction in degrees (360=north, 90=east,
        180=south, 270=west, 0=calm/variable)
        FHVEC = Vectorgemiddelde windsnelheid (in 0.1 m/s). Zie
        http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
        / Vector mean windspeed (in 0.1 m/s)
        FG    = Etmaalgemiddelde windsnelheid (in 0.1 m/s) / Daily mean
        windspeed (in 0.1 m/s)
        FHX   = Hoogste uurgemiddelde windsnelheid (in 0.1 m/s) / Maximum
        hourly mean windspeed (in 0.1 m/s)
        FHXH  = Uurvak waarin FHX is gemeten / Hourly division in which
        FHX was measured
        FHN   = Laagste uurgemiddelde windsnelheid (in 0.1 m/s) / Minimum
        hourly mean windspeed (in 0.1 m/s)
        FHNH  = Uurvak waarin FHN is gemeten / Hourly division in which
        FHN was measured
        FXX   = Hoogste windstoot (in 0.1 m/s) / Maximum wind gust (in
        0.1 m/s)
        FXXH  = Uurvak waarin FXX is gemeten / Hourly division in which
        FXX was measured
        TG    = Etmaalgemiddelde temperatuur (in 0.1 graden Celsius) /
        Daily mean temperature in (0.1 degrees Celsius)
        TN    = Minimum temperatuur (in 0.1 graden Celsius) / Minimum
        temperature (in 0.1 degrees Celsius)
        TNH   = Uurvak waarin TN is gemeten / Hourly division in which TN
        was measured
        TX    = Maximum temperatuur (in 0.1 graden Celsius) / Maximum
        temperature (in 0.1 degrees Celsius)
        TXH   = Uurvak waarin TX is gemeten / Hourly division in which TX
        was measured
        T10N  = Minimum temperatuur op 10 cm hoogte (in 0.1 graden
        Celsius) / Minimum temperature at 10 cm above surface (in 0.1
        degrees Celsius)
        T10NH = 6-uurs tijdvak waarin T10N is gemeten / 6-hourly division
        in which T10N was measured; 6=0-6 UT, 12=6-12 UT, 18=12-18 UT,
        24=18-24 UT
        SQ    = Zonneschijnduur (in 0.1 uur) berekend uit de globale
        straling (-1 voor <0.05 uur) / Sunshine duration (in 0.1 hour)
        calculated from global radiation (-1 for <0.05 hour)
        SP    = Percentage van de langst mogelijke zonneschijnduur /
        Percentage of maximum potential sunshine duration
        Q     = Globale straling (in J/cm2) / Global radiation (in J/cm2)
        DR    = Duur van de neerslag (in 0.1 uur) / Precipitation duration
        (in 0.1 hour)
        RH    = Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) /
        Daily precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
        RHX   = Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05
        mm) / Maximum hourly precipitation amount (in 0.1 mm) (-1 for
        <0.05 mm)
        RHXH  = Uurvak waarin RHX is gemeten / Hourly division in which
        RHX was measured
        PG    = Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1
        hPa) berekend uit 24 uurwaarden / Daily mean sea level pressure
        (in 0.1 hPa) calculated from 24 hourly values
        PX    = Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau
        (in 0.1 hPa) / Maximum hourly sea level pressure (in 0.1 hPa)
        PXH   = Uurvak waarin PX is gemeten / Hourly division in which PX
        was measured
        PN    = Laagste uurwaarde van de luchtdruk herleid tot zeeniveau
        (in 0.1 hPa) / Minimum hourly sea level pressure (in 0.1 hPa)
        PNH   = Uurvak waarin PN is gemeten / Hourly division in which PN
        was measured
        VVN   = Minimum opgetreden zicht / Minimum visibility; 0: <100 m,
        1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
        56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
        ..., 89: >70 km)
        VVNH  = Uurvak waarin VVN is gemeten / Hourly division in which
        VVN was measured
        VVX   = Maximum opgetreden zicht / Maximum visibility; 0: <100 m,
        1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
        56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
        ..., 89: >70 km)
        VVXH  = Uurvak waarin VVX is gemeten / Hourly division in which
        VVX was measured
        NG    = Etmaalgemiddelde bewolking (bedekkingsgraad van de
        bovenlucht in achtsten, 9=bovenlucht onzichtbaar) / Mean daily
        cloud cover (in octants, 9=sky invisible)
        UG    = Etmaalgemiddelde relatieve vochtigheid (in procenten) /
        Daily mean relative atmospheric humidity (in percents)
        UX    = Maximale relatieve vochtigheid (in procenten) / Maximum
        relative atmospheric humidity (in percents)
        UXH   = Uurvak waarin UX is gemeten / Hourly division in which UX
        was measured
        UN    = Minimale relatieve vochtigheid (in procenten) / Minimum
        relative atmospheric humidity (in percents)
        UNH   = Uurvak waarin UN is gemeten / Hourly division in which UN
        was measured
        EV24  = Referentiegewasverdamping (Makkink) (in 0.1 mm) /
        Potential evapotranspiration (Makkink) (in 0.1 mm)

    Returns
    -------
    knmi_df : pandas DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """
    if isinstance(stn, str):
        stn = int(stn)

    settings = _get_default_settings(settings)

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn, stations)

    # download data
    if settings["fill_missing_obs"] and (settings["interval"] == "hourly"):
        raise NotImplementedError("cannot yet fill missing values in hourly data")

    elif settings["fill_missing_obs"]:
        knmi_df, variables, station_meta = fill_missing_measurements(
            stn, stn_name, meteo_var, start, end, settings
        )
    else:
        knmi_df, variables, station_meta = download_knmi_data(
            stn, stn_name, meteo_var, start, end, settings
        )

    if not station_meta is None:
        meta = station_meta.to_dict()
    else:
        meta = {}
    meta.update(variables)

    # set metadata
    x = stations.loc[stn, "x"]
    y = stations.loc[stn, "y"]
    meta.update({"x": x, "y": y, "station": stn, "name": f"{meteo_var}_{stn_name}"})

    return knmi_df, meta


def get_knmi_obslist(
    locations=None,
    stns=None,
    xy=None,
    meteo_vars=("RH",),
    starts=None,
    ends=None,
    ObsClasses=None,
    settings=None,
    raise_exceptions=False,
    method="nearest",
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
        xy coördinates of the locations. e.g. [[10,25], [5,25]]
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
    settings : dict or None, optional
        settings for obtaining the right time series, options are:
            fill_missing_obs : bool, optional
                if True nan values in time series are filled with nearby time series.
                The default is True.
            interval : str, optional
                desired time interval for observations. The default is 'daily'.
            inseason : boolean, optional
                flag to obtain inseason data. The default is False
            use_api : bool, optional
                if True the api is used to obtain the data, API documentation is here:
                    https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
                Default is True as the API since (July 2021).
            raise_exceptions : bool, optional
                if True you get errors when no data is returned. The default is True.
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.

    Returns
    -------
    obs_list : list of obsevation objects
        collection of multiple point observations
    """

    if settings is None:
        settings = _get_default_settings(settings)
    settings["raise_exceptions"] = raise_exceptions

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
        obs_kwargs = {}
        if meteo_var in ("RD", "RH", "EV24"):
            if meteo_var == "RD":
                obs_kwargs["stn_type"] = "precipitation"
            elif meteo_var == "RH":
                obs_kwargs["stn_type"] = "meteo"
            elif meteo_var == "EV24":
                obs_kwargs["et_type"] = "EV24"

        start, end = _start_end_to_datetime(start, end)

        # get stations
        if (stns is None) and (method == "nearest"):
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
                    "stns, location and x are both None" "please specify one of these"
                )
            _stns = np.unique(_stns)
        elif (stns is None) and (method == "interpolation"):
            _stns = get_stations(meteo_var=meteo_var)

            if locations is not None:
                xy = locations[["x", "y"]].values
        else:
            _stns = stns

        if method == "nearest":
            for stn in _stns:
                o = ObsClass.from_knmi(
                    stn,
                    startdate=start,
                    enddate=end,
                    fill_missing_obs=settings["fill_missing_obs"],
                    interval=settings["interval"],
                    inseason=settings["inseason"],
                    raise_exceptions=raise_exceptions,
                    **obs_kwargs,
                )

                if settings["normalize_index"]:
                    o.index = o.index.normalize()

                obs_list.append(o)

        elif method == "interpolation":
            settings["fill_missing_obs"] = False  # dont fill missing obs

            # get dataframe with data from all knmi stations
            df = pd.DataFrame(
                index=pd.date_range(start="1900", end=end, freq="H"),
                columns=_stns.index,
            )
            for stn in _stns.index:  # fill dataframe with measurements
                et, meta = get_evaporation(
                    stn, meteo_var, start=start, end=end, settings=settings
                )
                df[stn] = et

            ts = interpolate(
                np.asarray(xy), _stns, df.dropna(how="all"), meteo_var=meteo_var
            )
            ts[ts < 0] = 0

            for i, col in enumerate(ts.columns):
                meta = {
                    "station": "interpolation thin plate sline",
                    "x": xy[i][0],
                    "y": xy[i][1],
                    "name": col,
                    "meteo_var": meteo_var,
                }
                o = ObsClass(ts.loc[:, [col]].copy(), **meta)

                if settings["normalize_index"]:
                    o.index = o.index.normalize()
                obs_list.append(o)

    return obs_list


def add_missing_indices(knmi_df, stn, start, end):
    """when downloading KNMI data you don't always get a DataFrame with the
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


def fill_missing_measurements(
    stn, stn_name=None, meteo_var="RH", start=None, end=None, settings=None
):
    """fill missing measurements in knmi data.

    Parameters
    ----------
    stn : int or str
        measurement station.
    meteo_var : str, optional
        observation type. The default is 'RH'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    settings : dict, optional
        settings for obtaining data

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
    settings = _get_default_settings(settings)

    if not isinstance(meteo_var, str):
        raise (TypeError(f"meteo var should be string not {type(meteo_var)}"))

    # get the location of the stations
    stations = get_stations(meteo_var=meteo_var)
    if stn_name is None:
        stn_name = get_station_name(stn, stations)

    # get start and end date
    start, end = _start_end_to_datetime(start, end)

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
        stn, stn_name, meteo_var, start=start, end=end, settings=settings
    )

    # if the first station cannot be read, read another station as the first
    ignore = [stn]
    while knmi_df.empty:
        logger.info(f"station {stn} has no measurements between {start} and {end}")
        logger.info("trying to get measurements from nearest station")
        stn = get_nearest_station_df(
            stations.loc[[stn]], meteo_var=meteo_var, ignore=ignore
        )[0]
        stn_name = get_station_name(stn, stations)
        knmi_df, variables, station_meta = download_knmi_data(
            stn, stn_name, meteo_var, start=start, end=end, settings=settings
        )
        ignore.append(stn)

    # find missing values
    knmi_df = add_missing_indices(knmi_df, stn, start, end)

    missing = knmi_df[meteo_var].isna()
    logger.info(f"station {stn} has {missing.sum()} missing measurements")

    # fill missing values
    settings["raise_exceptions"] = False
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
            stn_comp, stn_name_comp, meteo_var, start=start, end=end, settings=settings
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

    return knmi_df, variables, station_meta


def makkink(tmean, K):
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
    tmean, tmin, tmax, K, wind, rh, dates, z=1.0, lat=52.1, G=0.0, wh=10.0, tdew=None
):
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


def hargreaves(tmean, tmin, tmax, dates, lat=52.1, x=None):
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


def get_evaporation(stn=260, et_type="EV24", start=None, end=None, settings=None):
    """Collect different types of (reference) evaporation
    from KNMI weather stations

    Parameters
    ----------
    stn : str
        station number, defaults to 260 De Bilt
    et_type : str
        Choice between 'EV24', 'penman', 'makkink' or 'hargraves'.
        Defaults to 'EV24' which collects the KNMI Makkink EV24 evaporation.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    settings : dict or None, optional
        settings for obtaining the right time series, options are:
    fill_missing_obs : bool, optional
        if True nan values in time series are filled with nearby time series.
        The default is True.
    interval : str, optional
        desired time interval for observations. The default is 'daily'.
    inseason : boolean, optional
        flag to obtain inseason data. The default is False
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        Default is True as the API since (July 2021).
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.


    Returns
    ------
    pandas DataFrame

    """
    if et_type == "EV24":
        et, meta = get_knmi_timeseries_stn(
            stn, meteo_var="EV24", start=start, end=end, settings=settings
        )
    elif et_type == "hargreaves":
        d = {}
        mvs = ["TG", "TN", "TX"]
        for mv in mvs:
            ts, meta = get_knmi_timeseries_stn(
                stn, meteo_var=mv, start=start, end=end, settings=settings
            )
            d[mv] = ts.squeeze()
        et = hargreaves(
            d["TG"],
            d["TN"],
            d["TX"],
            d["TG"].index,
            meta["LAT_north"][str(meta["station"])],
        ).to_frame(name=et_type)
    elif et_type == "makkink":
        d = {}
        mvs = ["TG", "Q"]
        for mv in mvs:
            ts, meta = get_knmi_timeseries_stn(
                stn, meteo_var=mv, start=start, end=end, settings=settings
            )
            d[mv] = ts.squeeze()
        et = makkink(d["TG"], d["Q"]).to_frame(name=et_type)
    elif et_type == "penman":
        d = {}
        mvs = ["TG", "TN", "TX", "Q", "FG", "UG"]
        for mv in mvs:
            ts, meta = get_knmi_timeseries_stn(
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
            meta["LAT_north"][str(meta["station"])],
            meta["ALT_m"][str(meta["station"])],
        ).to_frame(name=et_type)
    else:
        raise ValueError(
            "Provide valid argument evaporation type \
                         'hargreaves', 'makkink' or 'penman'"
        )

    return et, meta


def interpolate(
    xy,
    stations,
    obs,
    meteo_var="EV24",
    kernel="thin_plate_spline",
    kernel2="linear",
    epsilon=None,
):
    """Interpolation method using the Scipy radial basis function (RBF)

    Parameters
    ----------
    xy : list or numpy array, optional
        xy coördinates in m RD of the locations. e.g. [[10,25], [5,25]]
    stations : pandas DataFrame
        Dataframe containing the observation locations at the index
        and the x and y coordinates as a column 'x' and 'y' respectively.
    obs : pandas DataFrame
        Dataframe containing the observation locations as columns and
        the observations at a measurement time in each row.
    meteo_var : str, optional
        Kind of observation to put in the column name DataFrame
    kernel : str, optional
        Type of radial basis funtion, by default thin_plate_spline.
        Other options are linear, gaussian, inverse_quadratic,
        multiquadric, inverse_multiquadric, cubic or quintic
    kernel2 : str, optional
        Kernel in case there are not enough observations (3 or 6) for
        time step, by default linear. Other options are gaussian,
        inverse_quadratic, multiquadric, or inverse_multiquadric.
    epsilon : float, optional
        Shape parameter that scales the input to the RBF. If kernel is
        linear, thin_plate_spline, cubic, or quintic, this defaults to 1.
        Otherwise this must be specified.

    Returns
    -------
    pandas DataFrame
        DataFrame with the interpolated observations.
    """

    if (kernel == "thin_plate_spline") or (kernel == "cubic"):
        min_val = 3
    elif kernel == "quintic":
        min_val = 6
    else:
        min_val = len(obs.columns)

    xy_obs = stations.loc[obs.columns, ["x", "y"]]

    df = pd.DataFrame(
        index=obs.index,
        columns=[f"{meteo_var}_{int(_xy[0])}_{int(_xy[1])}" for _xy in xy],
    )

    for idx in obs.index:
        # get all stations with values for this date
        val = obs.loc[idx].dropna()
        # get stations for this date
        coor = xy_obs.loc[val.index]

        if len(val) >= min_val:
            kernel = kernel
        else:
            kernel = kernel2

        # create an scipy interpolator
        rbf = RBFInterpolator(
            coor.to_numpy(), val.to_numpy(), epsilon=epsilon, kernel=kernel
        )

        val_rbf = rbf(xy)
        df.loc[idx] = val_rbf

    return df.astype(float)
