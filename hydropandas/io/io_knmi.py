import datetime as dt
import os
import re
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
import requests

from .. import util

import logging
logger = logging.getLogger(__name__)

URL_DAILY_RD = 'https://www.daggegevens.knmi.nl/klimatologie/monv/reeksen'
URL_DAILY_METEO = 'https://www.daggegevens.knmi.nl/klimatologie/daggegevens'
URL_HOURLY_METEO = 'https://www.daggegevens.knmi.nl/klimatologie/uurgegevens'


def get_stations(meteo_var='RD'):
    """get knmi stations from json files according to variable.

    Parameters
    ----------
    meteo_var : str, optional
        type of meteodata, by default 'RD'

    Returns
    -------
    pandas DataFrame with stations, names and coordinates (Lat/Lon & RD)
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))

    if meteo_var == "RD":
        fname = "../data/knmi_neerslagstation.json"
    else:
        fname = "../data/knmi_meteostation.json"

    stations = pd.read_json(os.path.join(dir_path, fname))

    if meteo_var == 'PG':
        # in Ell wordt geen luchtdruk gemeten
        stations.drop(377, inplace=True)
    if meteo_var == 'EV24':
        # in Woensdrecht wordt geen verdamping gemeten
        stations.drop(340, inplace=True)

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
    stn_name = stations.loc[stn, 'naam']
    stn_name = stn_name.upper().replace(' ', '-')
    stn_name = stn_name.replace('(', '').replace(')', '')

    return stn_name


def get_nearest_stations_xy(x, y, meteo_var, n=1, stations=None, ignore=None):
    """find the KNMI stations that measure 'variable' closest to the x, y
    coordinates.

    Parameters
    ----------
    x : int or float
        x coordinate in RD
    y : int or float
        x coordinate in RD
    meteo_var : str
        measurement variable e.g. 'RD' or 'EV24'
    n : int, optional
        number of stations you want to return. The default is 1.
    stations : pd.DataFrame, optional
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

    d = np.sqrt((stations.x - x)**2 + (stations.y - y)**2)

    return d.nsmallest(n).index.to_list()


def get_nearest_station_df(locations, xcol='x', ycol='y', stations=None,
                           meteo_var="RD", ignore=None):
    """find the KNMI stations that measure 'meteo_var' closest to the
    coordinates in 'locations'.

    Parameters
    ----------
    locations : pd.DataFrame
        x and y coordinates
    xcol : str
        name of the column in the locations dataframe with the x values
    ycol : str
        name of teh column in the locations dataframe with the y values
    stations : pd.DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    meteo_var : str
        measurement variable e.g. 'RD' or 'EV24'
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

    distances = pd.DataFrame(np.sqrt((xh - xi) ** 2 + (yh - yi) ** 2),
                             index=locations.index,
                             columns=stations.index)

    stns = distances.idxmin(axis=1).unique()
    if np.any(np.isnan(stns)):
        stns = stns[~np.isnan(stns)].astype(int)

    return stns


def get_nearest_station_grid(xmid, ymid, stations=None, meteo_var="RD",
                             ignore=None):
    """find the KNMI stations that measure 'meteo_var' closest to all cells in
    a grid.

    Parameters
    ----------
    xmid : np.array
        x coördinates of the cell centers of your grid shape(ncol)
    ymid : np.array
        y coördinates of the cell centers of your grid shape(nrow)
    stations : pd.DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    meteo_var : str
        measurement variable e.g. 'RD' or 'EV24'
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
    mg = np.meshgrid(xmid, ymid)

    locations = pd.DataFrame(data={'x': mg[0].ravel(),
                                   'y': mg[1].ravel()})

    stns = get_nearest_station_df(locations, xcol='x', ycol='y',
                                  stations=stations, meteo_var=meteo_var,
                                  ignore=ignore)

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
    # start date one day before because later the datetime index is modified
    start = start - pd.Timedelta(1, 'D')
    if end is None:
        end = pd.Timestamp.today() - pd.Timedelta(1, unit='D')
    else:
        end = pd.to_datetime(end)

    return start, end


def _check_latest_measurement_date_RD_debilt(use_api=True):
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
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        Default is False as the API is down (apr 2021).

    Returns
    -------
    last_measurement_date_debilt : pd.TimeStamp
        last date with measurements at station de Bilt
    """

    look_back_days = 90
    start = dt.datetime.now() - pd.Timedelta(look_back_days, unit='D')
    end = dt.datetime.now() + pd.Timedelta(10, unit='D')
    if use_api:
        knmi_df, _ = get_knmi_daily_rainfall_api(URL_DAILY_RD, 550, "RD",
                                                 start=start,
                                                 end=end, inseason=False)
    else:
        knmi_df, _ = get_knmi_daily_rainfall_url(550, 'DE-BILT', 'RD', start,
                                                 end, inseason=False)

    knmi_df = knmi_df.dropna()
    if knmi_df.empty:
        raise ValueError(
            f'knmi station de Bilt has no RD measurements in the past {look_back_days} days.')

    last_measurement_date_debilt = knmi_df.index[-1]

    logger.info(f'last RD measurement available at the Bilt is from'
                f' {last_measurement_date_debilt.strftime("%Y-%m-%d")}')
    logger.info('assuming no measurements are available at other stations before this date')

    return last_measurement_date_debilt


def download_knmi_data(stn, stn_name=None,
                       meteo_var='RD', start=None, end=None,
                       interval='daily', inseason=False,
                       use_api=True, raise_exceptions=True):
    """download knmi data of a measurements station for certain observation
    type.

    Parameters
    ----------
    stn : int or str
        number of measurement station
    stn_name : str
        Only necesary if meteo_var is RD and use_api is True.
    meteo_var : str, optional
        measurement type 'RD' or 'EV24'. The default is 'RD'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    interval : str, optional
        time interval of observations. The default is 'daily'.
    inseason : bool, optional
        passed to the knmi api. The default is False.
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        Default is True as the API since (July 2021).
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.

    Raises
    ------
    NotImplementedError
        different time intervals and inseason data is not yet working.
    ValueError
        if the data from knmi cannot not be read a ValueError is raised.
        Unless raise_exceptions is False

    Returns
    -------
    knmi_df : pd.DataFrame
        data from one station from one type of observation
    variables : dictionary
        information about the observerd variables
    stations : pd.DataFrame
        information about the measurement station.
    """
    # checks
    if interval.startswith('hour') and meteo_var == 'RD':
        message = 'Interval cannot be hourly for rainfall-stations'
        raise (ValueError(message))

    if interval != 'daily':
        raise NotImplementedError('only daily intervals are working now')

    if inseason:
        raise NotImplementedError('season stuff not implemented')

    start, end = _start_end_to_datetime(start, end)

    # convert possible integer to string
    stn = str(stn)

    logger.info(f'download knmi {meteo_var} data from station {stn}-{stn_name} between {start} and {end}')

    # define variables
    knmi_df = pd.DataFrame()
    variables = {}
    stations = pd.DataFrame()

    # download and read data
    try:
        if use_api:
            if interval.startswith('hour'):
                # hourly data from meteorological stations
                knmi_df = get_knmi_hourly_api(URL_HOURLY_METEO, stn, meteo_var, start, end)

            elif meteo_var == 'RD':
                # daily data from rainfall-stations
                knmi_df, variables = get_knmi_daily_rainfall_api(
                    URL_DAILY_RD, stn, meteo_var, start, end, inseason)
            else:
                # daily data from meteorological stations
                knmi_df, variables, stations = get_knmi_daily_meteo_api(
                    URL_DAILY_METEO, stn, meteo_var, start, end, inseason)
        else:
            if interval.startswith('hour'):
                # hourly data from meteorological stations
                raise NotImplementedError()
            elif meteo_var == 'RD':
                # daily data from rainfall-stations
                if stn_name is None:
                    stations_df = get_stations(meteo_var=meteo_var)
                    stn_name = get_station_name(int(stn), stations_df)
                knmi_df, variables = get_knmi_daily_rainfall_url(
                    stn, stn_name, meteo_var, start, end, inseason)
            else:
                # daily data from meteorological stations
                knmi_df, variables, stations = get_knmi_daily_meteo_url(
                    stn, meteo_var, start, end, inseason)
    except (ValueError, KeyError) as e:
        logger.error(e)
        if raise_exceptions:
            raise ValueError(e)

    return knmi_df, variables, stations


def get_knmi_daily_rainfall_api(url, stn, meteo_var,
                                start=None, end=None, inseason=False):
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
    pd.DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """

    data = {
        'inseason': str(int(inseason)),
        'vars': meteo_var,
        'stns': stn}

    if start is not None:
        data['start'] = start.strftime('%Y%m%d')

    if end is not None:
        data['end'] = end.strftime('%Y%m%d')

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    f = StringIO(result.text)
    knmi_df, variables = read_knmi_daily_rainfall(
        f, meteo_var)
    if int(stn) not in knmi_df.STN.unique():
        raise ValueError(f'KNMI station {stn} does not exists or has no '
                         'measurements in the given period')

    return knmi_df[[meteo_var]], variables


def get_knmi_daily_rainfall_url(stn, stn_name, meteo_var,
                                start=None, end=None, inseason=False,
                                use_cache=True):
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
    pd.DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    """
    stn = f'{int(stn):03}'
    url = f"https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/monv_reeksen/neerslaggeg_{stn_name}_{stn}.zip"

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), 'knmi')
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(),
                             'knmi', f'neerslaggeg_{stn}.zip')
    fname_dir = os.path.join(tempfile.gettempdir(),
                             'knmi', f'neerslaggeg_{stn}')
    fname_txt = os.path.join(fname_dir, f'neerslaggeg_{stn_name}_{stn}.txt')

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
            raise ValueError(
                f'invalid url {url} please check station name {stn_name}')
        with open(fname_zip, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        # unzip file
        util.unzip_file(fname_zip, fname_dir, force=True,
                        preserve_datetime=True)

    with open(fname_txt, 'r') as f:
        line = f.readline()
        # get meteo var
        for iline in range(50):
            if 'RD   ' in line:  # RD komt te vaak voor vandaar de spaties
                key, item = line.split('=')
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        # get dataframe
        for iline in range(50):
            if 'STN' in line:
                columns = line.strip('# ').strip('\n').split(',')
                columns = [x.strip(' ') for x in columns]
                values = f.readline()
                if values == '\n':
                    values = f.readline()
            line = f.readline()

        df = pd.read_csv(f, header=None, names=columns, na_values='     ')

        df.set_index(pd.to_datetime(df.YYYYMMDD, format='%Y%m%d'),
                     inplace=True)
        df = df.drop('YYYYMMDD', axis=1)

        if df.index.duplicated().sum() > 0:
            df = df.loc[~df.index.duplicated(keep='first')]
            logger.info('duplicate indices removed from RD measurements')

        # sometimes the last row is empty, check for that and remove it
        if not df.empty:
            if df.iloc[-1].isna().any():
                logger.debug('duplicate indices removed from RD measurements')
                df = df.drop(index=df.index[-1])
                df.loc[:, meteo_var] = df[meteo_var].astype(float)

        # daily precipitation amount in 0.1 mm over the period 08.00
        # preceding day - 08.00 UTC present day
        df.index = df.index + pd.to_timedelta(8, unit='h')

        # from UT to UT+1 (standard-time in the Netherlands)
        df.index = df.index + pd.to_timedelta(1, unit='h')

        df, variables = _transform_variables(df, variables)

    return df.loc[start:end, [meteo_var]], variables


def _read_knmi_header(f):

    variables = {}
    line = f.readline()
    if 'DOCTYPE HTML PUBLIC' in line:
        logger.error(f.read())
        raise ValueError('Internal Server Error')
    for iline in range(500):
        if ' = ' in line or ' : ' in line:
            line = line.lstrip(' #').strip('\n')
            if ' = ' in line:
                varDes = line.split(' = ')
            else:
                varDes = line.split(' : ')
            variables[varDes[0].strip()] = varDes[1].strip()

        if 'STN,YY' in line:
            header = line.replace('#', '').split(',')
            header = [item.lstrip().rstrip() for item in header]
            break

        line = f.readline()

    if iline > 498:
        raise ValueError('cannot read measurements from file')

    return f, variables, header


def _transform_variables(df, variables):

    for key, value in variables.items():
        # test if key existst in data
        if key not in df.keys():
            if key == 'YYYYMMDD' or key == 'HH':
                pass
            elif key == 'T10N':
                variables.pop(key)
                key = 'T10'
            else:
                raise NameError(key + ' does not exist in data')

        if '0.1 ' in value:
            logger.debug(f'transform {key}, {value} from 0.1 to 1')

            df[key] = df[key] * 0.1
            value = value.replace('0.1 ', '')
        if ' tiende ' in value:
            logger.debug(f'transform {key}, {value} from 0.1 to 1')

            df[key] = df[key] * 0.1
            value = value.replace(' tiende ', ' ')
        if ' mm' in value:
            logger.debug(f'transform {key}, {value} from mm to m')

            df[key] = df[key] * 0.001
            value = value.replace(' mm', ' m')
        if ' millimeters' in value:
            logger.debug(f'transform {key}, {value} from mm to m')

            df[key] = df[key] * 0.001
            value = value.replace(' millimeters', ' m')
        # Store new variable
        variables[key] = value

    return df, variables


def read_knmi_daily_rainfall(f, meteo_var):

    f, variables, header = _read_knmi_header(f)

    df = pd.read_csv(f, header=None, names=header, na_values='     ')
    f.close()

    df.set_index(pd.to_datetime(df.YYYYMMDD, format='%Y%m%d'),
                 inplace=True)
    df = df.drop('YYYYMMDD', axis=1)

    if df.index.duplicated().sum() > 0:
        df = df.loc[~df.index.duplicated(keep='first')]
        logger.debug('duplicate indices removed from RD measurements')

    # sometimes the last row is messed up, check for that and remove it
    if not df.empty:
        if df.iloc[-1].isna().any():
            logger.debug('last row contains no data, remove last row')

            df = df.drop(index=df.index[-1])
            df.loc[:, meteo_var] = df[meteo_var].astype(float)

    # daily precipitation amount in 0.1 mm over the period 08.00
    # preceding day - 08.00 UTC present day
    df.index = df.index + pd.to_timedelta(8, unit='h')

    # from UT to UT+1 (standard-time in the Netherlands)
    df.index = df.index + pd.to_timedelta(1, unit='h')

    df, variables = _transform_variables(df, variables)

    return df, variables


def _read_station_location(f):

    stations = None

    line = f.readline()
    for iline in range(30):
        if 'STN' in line:
            titels = line.strip('# ').split()
            titels = [x.replace('(', '_') for x in titels]
            titels = [x.replace(r')', '') for x in titels]

            values = f.readline().strip('# ').strip().replace(':', '')
            values = re.split(r'\s{2,}', values)

            # Create pd.DataFrame for station data
            stations = pd.DataFrame(columns=titels, data=[values])
            stations.set_index(['STN'], inplace=True)
            for col in stations.columns:
                try:
                    stations.loc[:, col] = stations[col].astype(float)
                except ValueError:
                    pass

            break

        line = f.readline()

    if stations is None:
        logger.warning('could not find stations')

    return f, stations


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
    pd.DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    stations : pd.DataFrame
        additional data about the measurement station
    """
    data = {
        'start': start.strftime('%Y%m%d'),
        'end': end.strftime('%Y%m%d'),
        'inseason': str(int(inseason)),
        'vars': meteo_var,
        'stns': stn,
    }

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    f = StringIO(result.text)
    knmi_df, variables, stations = read_knmi_daily_meteo(f)

    knmi_df.dropna(subset=[meteo_var], inplace=True)

    return knmi_df[[meteo_var]], variables, stations


def get_knmi_daily_meteo_url(stn, meteo_var, start, end,
                             use_cache=True):
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
    pd.DataFrame
        measurements.
    variables : dictionary
        additional information about the variables
    stations : pd.DataFrame
        additional data about the measurement station
    """
    url = f"https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/daggegevens/etmgeg_{stn}.zip"

    # get file name
    basedir = os.path.join(tempfile.gettempdir(), 'knmi')
    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    fname_zip = os.path.join(tempfile.gettempdir(),
                             'knmi', f'etmgeg_{stn}.zip')
    fname_dir = os.path.join(tempfile.gettempdir(), 'knmi', f'etmgeg_{stn}')
    fname_txt = os.path.join(fname_dir, f'etmgeg_{stn}.txt')

    # check if file should be downloaded and unzipped
    donwload_and_unzip = True
    if os.path.exists(fname_txt) and use_cache:
        fname_time_mod = dt.datetime.fromtimestamp(os.path.getmtime(fname_txt))
        if fname_time_mod > end:
            donwload_and_unzip = False

    if donwload_and_unzip:
        # download zip file
        r = requests.get(url, stream=True)
        with open(fname_zip, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        # unzip file
        util.unzip_file(fname_zip, fname_dir, force=True,
                        preserve_datetime=True)

    variables = None
    with open(fname_txt, 'r') as f:
        line = f.readline()
        # get meteo var
        for iline in range(50):
            if meteo_var in line:
                key, item = line.split('=')
                variables = {key.strip(): item.strip()}
                break
            line = f.readline()

        if variables is None:
            raise ValueError(f'could not find {meteo_var} for station {stn}')

        # get dataframe
        for iline in range(50):
            if 'STN' in line:
                columns = line.strip('# ').strip('\n').split(',')
                columns = [x.strip(' ') for x in columns]
                values = f.readline()
                if values == '\n':
                    values = f.readline()
                df = pd.read_csv(f, header=None, names=columns,
                                 na_values='     ')
                df.set_index(pd.to_datetime(df.YYYYMMDD, format='%Y%m%d'),
                             inplace=True)
                df = df.drop('YYYYMMDD', axis=1)

                df = df.loc[df.index.notnull(), :]

                # add a full day for meteorological data, so that the
                # timestamp is at the end of the period in the data
                df.index = df.index + pd.to_timedelta(1, unit='d')

                # from UT to UT+1 (standard-time in the Netherlands)
                df.index = df.index + pd.to_timedelta(1, unit='h')
                df = df.loc[start:end, [meteo_var]]
                df = df.dropna()
                df, variables = _transform_variables(df, variables)
                break

            line = f.readline()

    return df, variables, None


def read_knmi_daily_meteo(f):

    f, stations = _read_station_location(f)
    f, variables, header = _read_knmi_header(f)
    header[0] = header[0].lstrip('# ')
    df = pd.read_csv(f, header=None, names=header, na_values='     ')
    f.close()

    df.set_index(pd.to_datetime(df.YYYYMMDD, format='%Y%m%d'),
                 inplace=True)
    df = df.drop('YYYYMMDD', axis=1)

    df = df.loc[df.index.notnull(), :]

    # add a full day for meteorological data, so that the
    # timestamp is at the end of the period in the data
    df.index = df.index + pd.to_timedelta(1, unit='d')

    # from UT to UT+1 (standard-time in the Netherlands)
    df.index = df.index + pd.to_timedelta(1, unit='h')

    df, variables = _transform_variables(df, variables)

    return df, variables, stations

   # return knmi_series


def get_knmi_hourly_api(url, stn, meteo_var, start, end):

    data = {
        'start': start.strftime('%Y%m%d') + '01',
        'end': end.strftime('%Y%m%d') + '24',
        'vars': meteo_var,
        'stns': stn,
    }

    result = requests.get(url, params=data)

    if result.status_code != 200:
        raise requests.ConnectionError(f"Cannot connect to {url}")

    f = StringIO(result.text)
    knmi_series = read_knmi_hourly(f)

    return knmi_series


def read_knmi_hourly(f):
    raise NotImplementedError('work in progress')
    knmi_df = pd.DataFrame()
    return knmi_df


def get_knmi_timeseries_xy(x, y, meteo_var, start, end,
                           stn_name=None,
                           fill_missing_obs=True,
                           interval='daily',
                           inseason=False,
                           raise_exceptions=False):
    """Get timeseries with measurements from station closest to the given (x,y)
    coördinates.

    Parameters
    ----------
    x : int or float
        x coördinate in m RD.
    y : int or float
        y coördinate in m RD.
    meteo_var : str
        e.g. 'EV24'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    fill_missing_obs : bool
        if True missing observations are filled with values of next closest
        KNMI station
    interval : str, optional
        desired time interval for observations. The default is 'daily'.
    inseason : boolean, optional
        flag to obtain inseason data. The default is False
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is False.

    Returns
    -------
    knmi_df : pd.DataFrame
        time series.
    meta : dic
        dictionary with metadata.
    """

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn = get_nearest_stations_xy(x, y, meteo_var, stations=stations)[0]

    # download data
    if fill_missing_obs:
        knmi_df, variables, station_meta = \
            fill_missing_measurements(stn, stn_name, meteo_var, start, end,
                                      interval, raise_exceptions)
    else:
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, stn_name, meteo_var, start, end,
                               interval, inseason, raise_exceptions)

    if not station_meta is None:
        meta = station_meta.to_dict()
    else:
        meta = {}
    meta.update(variables)

    # set metadata
    name = meteo_var + ' ' + stations.loc[stn, 'naam']
    x = stations.loc[stn, 'x']
    y = stations.loc[stn, 'y']
    meta.update({'x': x, 'y': y, 'station': stn, 'name': name})

    return knmi_df, meta


def get_knmi_timeseries_stn(stn, meteo_var, start, end,
                            fill_missing_obs=True, interval='daily',
                            inseason=False, use_api=True, raise_exceptions=False):
    """Get a knmi time series and metadata.

    Parameters
    ----------
    stn : int or str
        measurement station e.g. 829.
    meteo_var : str, optional
        observation type e.g. "RD" or "EV24".
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
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
    -------
    knmi_df : pandas.DataFrame
        time series with measurements.
    meta : dictionary
        metadata from the measurement station.
    """

    # get station
    stations = get_stations(meteo_var=meteo_var)
    stn_name = get_station_name(stn, stations)

    # download data
    if fill_missing_obs:
        knmi_df, variables, station_meta = \
            fill_missing_measurements(stn, stn_name, meteo_var, start, end,
                                      interval, use_api=use_api,
                                      raise_exceptions=raise_exceptions)
    else:
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, stn_name, meteo_var, start, end,
                               interval, inseason,
                               use_api=use_api,
                               raise_exceptions=raise_exceptions)
    if not station_meta is None:
        meta = station_meta.to_dict()
    else:
        meta = {}
    meta.update(variables)

    # set metadata
    x = stations.loc[stn, 'x']
    y = stations.loc[stn, 'y']
    meta.update({'x': x, 'y': y, 'station': stn, 'name': stn_name})

    return knmi_df, meta


def get_knmi_obslist(locations=None, stns=None, xmid=None, ymid=None,
                     meteo_vars=["RD"], start=[None], end=[None],
                     ObsClass=None, fill_missing_obs=True,
                     normalize_index=True, interval='daily', inseason=False,
                     cache=False, raise_exceptions=False):
    """Get a list of observations of knmi stations. Either specify a list of
    knmi stations (stns) or a dataframe with x, y coordinates (locations).

    Parameters
    ----------
    locations : pd.DataFrame or None
        dataframe with x and y coordinates. The default is None
    stns : list of str or None
        list of knmi stations. The default is None
    xmid : np.array, optional
        x coördinates of the cell centers of your grid shape(ncol)
    ymid : np.array, optional
        y coördinates of the cell centers of your grid shape(nrow)
    meteo_vars : list or tuple of str
        meteo variables e.g. ["RD", "EV24"]. The default is ("RD")
    start : list of str, datetime or None]
        start date of observations per meteo variable. The start date is
        included in the time series. If start is None the start date will
        be January 1st of the previous year. The default is [None]
    end : list of str, datetime or None]
        end date of observations per meteo variable. The end date is not
        included in the time series. If end is None the last date with
        measurements is used. The default is [None]
    ObsClass : type or None
        class of the observations, only KnmiObs is supported for now. The
        default is None
    fill_missing_obs : bool, optional
        if True nan values in time series are filled with nearby time series.
        The default is True.
    normalize_index : bool, optional
        if True the index is normalized.
    interval : str, optional
        desired time interval for observations. The default is 'daily'.
    inseason : boolean, optional
        flag to obtain inseason data. The default is False
    cache : boolean, optional
        if True the observation data will be cached or read from cache.
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.

    Returns
    -------
    obs_list : list of obsevation objects
        collection of multiple point observations
    """
    obs_list = []
    for i, meteo_var in enumerate(meteo_vars):
        start[i], end[i] = _start_end_to_datetime(start[i], end[i])
        if stns is None:
            stations = get_stations(meteo_var=meteo_var)
            if (locations is None) and (xmid is not None):
                _stns = get_nearest_station_grid(xmid, ymid,
                                                 stations=stations,
                                                 meteo_var=meteo_var)
            elif locations is not None:
                _stns = get_nearest_station_df(locations,
                                               stations=stations,
                                               meteo_var=meteo_var)
            else:
                raise ValueError('stns, location and xmid are all None'
                                 'please specify one of these')
        else:
            _stns = stns

        for stn in _stns:
            if cache:
                cache_dir = os.path.join(tempfile.gettempdir(), 'knmi')
                if not os.path.isdir(cache_dir):
                    os.mkdir(cache_dir)

                fname = (f'{stn}-{meteo_var}-{start[i].strftime("%Y%m%d")}'
                         f'-{end[i].strftime("%Y%m%d")}'
                         f'-{fill_missing_obs}' + '.pklz')
                pklz_path = os.path.join(cache_dir, fname)

                if os.path.isfile(pklz_path):
                    logger.info(f'reading {fname} from cache')
                    o = pd.read_pickle(pklz_path)
                else:
                    o = ObsClass.from_knmi(stn, meteo_var, start[i], end[i],
                                           fill_missing_obs=fill_missing_obs,
                                           interval=interval,
                                           inseason=inseason,
                                           raise_exceptions=raise_exceptions)
                    o.to_pickle(pklz_path)
            else:
                o = ObsClass.from_knmi(stn, meteo_var, start[i], end[i],
                                       fill_missing_obs=fill_missing_obs,
                                       interval=interval,
                                       inseason=inseason,
                                       raise_exceptions=raise_exceptions)
            if normalize_index:
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
    knmi_df : pd.DataFrame
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
    knmi_df : pd.DataFrame
        data from one station from one type of observation
    """
    # check if given dates are more or less similar than measurement dates
    if (knmi_df.index[0] - start).days < 2:
        new_start = knmi_df.index[0]
    else:
        new_start = pd.Timestamp(year=start.year, month=start.month,
                                 day=start.day, hour=knmi_df.index[0].hour,
                                 minute=knmi_df.index[0].minute,
                                 second=knmi_df.index[0].second)
        logger.info(f'station {stn} has no measurements before {knmi_df.index[0]}')

    if (end - knmi_df.index[-1]).days < 2:
        new_end = knmi_df.index[-1]
    else:
        new_end = pd.Timestamp(year=end.year, month=end.month, day=end.day,
                               hour=knmi_df.index[-1].hour,
                               minute=knmi_df.index[-1].minute,
                               second=knmi_df.index[-1].second)
        logger.info(f'station {stn} has no measurements after {knmi_df.index[-1]}')

    # add missing indices
    new_index = pd.date_range(new_start, new_end, freq='D')
    knmi_df = knmi_df.reindex(new_index)

    return knmi_df


def fill_missing_measurements(stn, stn_name=None, meteo_var='RD',
                              start=None, end=None,
                              interval='daily',
                              use_api=True,
                              raise_exceptions=False):
    """fill missing measurements in knmi data.

    Parameters
    ----------
    stn : int or str
        measurement station.
    meteo_var : str, optional
        observation type. The default is 'RD'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    interval : str, optional
        desired time interval for observations. The default is 'daily'.
    use_api : bool, optional
        if True the api is used to obtain the data, API documentation is here:
            https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
        Default is True as the API since (July 2021).
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.

    Returns
    -------
    knmi_df : pd.DataFrame
        data from one station from one type of observation, with additional
        column to see which station is used to fill the value
    variables : dictionary
        information about the observerd variables
    stations : pd.DataFrame
        information about the measurement station.
    """

    if not isinstance(meteo_var, str):
        raise (TypeError(f'meteo var should be string not {type(meteo_var)}'))
    # get the location of the stations
    stations = get_stations(meteo_var=meteo_var)
    if stn_name is None:
        stn_name = get_station_name(stn, stations)

    # get start and end date
    start, end = _start_end_to_datetime(start, end)
    if (meteo_var == 'RD') and (end > (dt.datetime.now() - pd.Timedelta(90, unit='D'))):
        end = min(end, _check_latest_measurement_date_RD_debilt(use_api=use_api))
        logger.info(f'changing end_date to {end.strftime("%Y-%m-%d")}')

    # download data from station
    knmi_df, variables, station_meta = \
        download_knmi_data(stn, stn_name, meteo_var, start=start,
                           end=end, interval=interval,
                           inseason=False,
                           use_api=use_api,
                           raise_exceptions=raise_exceptions)

    # if the first station cannot be read, read another station as the first
    ignore = [stn]
    while knmi_df.empty:
        logger.info(f'station {stn} has no measurements between {start} and {end}')
        logger.info('trying to get measurements from nearest station')
        stn = get_nearest_station_df(
            stations.loc[[stn]], meteo_var=meteo_var, ignore=ignore)[0]
        stn_name = get_station_name(stn, stations)
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, stn_name, meteo_var, start=start,
                               end=end, interval=interval,
                               inseason=False,
                               use_api=use_api,
                               raise_exceptions=raise_exceptions)
        ignore.append(stn)

    # find missing values
    knmi_df = add_missing_indices(knmi_df, stn, start, end)

    missing = knmi_df[meteo_var].isna()
    logger.info(f'station {stn} has {missing.sum()} missing measurements')

    # fill missing values
    while np.any(missing) and not np.all(missing):
        stn_comp = get_nearest_station_df(
            stations.loc[[stn]], meteo_var=meteo_var, ignore=ignore)

        logger.info(f'trying to fill {missing.sum()} '
                     f'measurements with station {stn_comp}')

        if stn_comp is None:
            logger.info('could not fill all missing measurements there are '
                         'no stations left to check')

            missing[:] = False
            break
        else:
            stn_comp = stn_comp[0]
        stn_name_comp = get_station_name(stn_comp, stations)
        knmi_df_comp, _, __ = \
            download_knmi_data(stn_comp, stn_name_comp, meteo_var,
                               start=start, end=end,
                               interval=interval,
                               inseason=False,
                               use_api=use_api,
                               raise_exceptions=raise_exceptions)

        if knmi_df_comp.empty:
            logger.warning(f'station {stn_comp} cannot be downloaded')

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
                knmi_df.loc[ix_idx,
                            meteo_var] = knmi_df_comp.loc[ix_idx, meteo_var]
                # add source station number
                knmi_df.loc[ix_idx, 'station_opvulwaarde'] = str(stn_comp)

        missing = knmi_df[meteo_var].isna()
        ignore.append(stn_comp)

    return knmi_df, variables, station_meta
