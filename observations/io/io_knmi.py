import os
import re
from io import StringIO

import numpy as np
import pandas as pd
import requests


def get_stations(variable='RD'):
    """get knmi stations from json files according to variable

    Parameters
    ----------
    variable : str, optional
        [description], by default 'RD'

    Returns
    -------
    pandas DataFrame with stations, names and coordinates (Lat/Lon & RD)
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))

    if variable == "RD":
        fname = "../../data/knmi_neerslagstation.json"
    else:
        fname = "../../data/knmi_meteostation.json"

    stations = pd.read_json(os.path.join(dir_path, fname))

    if variable == 'PG':
        # in Ell wordt geen luchtdruk gemeten
        stations.drop(377, inplace=True)
    if variable == 'EV24':
        # in Woensdrecht wordt geen verdamping gemeten
        stations.drop(340, inplace=True)

    return stations


def get_nearest_stations_xy(x, y, variable, n=1, stations=None, ignore=None):
    """find the KNMI stations that measure 'variable' closest to the
    x, y coordinates


    Parameters
    ----------
    x : int or float
        x coordinate in RD
    y : int or float
        x coordinate in RD
    variable : str
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
        stations = get_stations(variable=variable)
    if ignore is not None:
        stations.drop(ignore, inplace=True)
        if stations.empty:
            return None

    d = np.sqrt((stations.x - x)**2 + (stations.y - y)**2)

    return d.nsmallest(n).index.to_list()


def get_nearest_station_df(locations, stations=None, variable="RD", ignore=None):
    """find the KNMI stations that measure 'variable' closest to the
    stations in 'locations'.


    Parameters
    ----------
    locations : pd.DataFrame
        station number, x and y coordinates
    stations : pd.DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    variable : str
        measurement variable e.g. 'RD' or 'EV24'
    ignore : list, optional
        list of stations to ignore. The default is None.

    Returns
    -------
    stns : list
        station numbers.

    """
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


def _start_end_to_datetime(start, end):
    """convert start and endtime to datetime


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
        start = pd.Timestamp(pd.Timestamp.today().year, 1, 1)
    else:
        start = pd.to_datetime(start)
    # start date one day before because later the datetime index is modified
    start = start - pd.Timedelta(1, 'D')
    if end is None:
        end = pd.Timestamp.today()
    else:
        end = pd.to_datetime(end)

    return start, end


def download_knmi_data(stn, var='RD', start=None, end=None, interval='daily',
                       inseason=False, raise_exceptions=True, verbose=False):
    """download knmi data of a measurements station for certain observation
    type


    Parameters
    ----------
    stn : int or str
        number of measurements station
    var : str, optional
        measurement type 'RD' or 'EV24'. The default is 'RD'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    interval : str, optional
        time interval of observations. The default is 'daily'.
    inseason : bool, optional
        passed to the knmi api. The default is False.
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.
    verbose : boolean, optional
            Print additional information to the screen (default is False).

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
    if interval.startswith('hour') and var == 'RD':
        message = 'Interval can not be hourly for rainfall-stations'
        raise (ValueError(message))

    if interval != 'daily':
        raise NotImplementedError('only daily intervals are working now')

    if inseason:
        raise NotImplementedError('season stuff not implemented')

    start, end = _start_end_to_datetime(start, end)

    # convert possible integer to string
    stn = str(stn)

    # define variables
    knmi_df = pd.DataFrame()
    variables = {}
    stations = pd.DataFrame()

    # download and read data
    try:
        if interval.startswith('hour'):
            # hourly data from meteorological stations
            url = 'http://projects.knmi.nl/klimatologie/uurgegevens/getdata_uur.cgi'
            knmi_df = get_knmi_hourly(url, stn, var, start, end)

        elif var == 'RD':
            # daily data from rainfall-stations
            url = 'http://projects.knmi.nl/klimatologie/monv/reeksen/getdata_rr.cgi'
            knmi_df, variables = get_knmi_daily_rainfall(
                url, stn, var, start, end, inseason, verbose)
        else:
            # daily data from meteorological stations
            url = 'http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi'
            knmi_df, variables, stations = get_knmi_daily_meteo(
                url, stn, var, start, end, inseason, verbose)
    except ValueError as e:
        if verbose:
            print(e)
        if raise_exceptions:
            raise ValueError(e)

    return knmi_df, variables, stations


def get_knmi_daily_rainfall(url, stn, var, start, end, inseason, verbose=False):
    """download and read knmi daily rainfall


    Parameters
    ----------
    url : str
        download url.
    stn : str
        station number.
    var : str
        must be 'RD'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    inseason : boolean
        flag to obtain inseason data.
    verbose : boolean, optional
        Print additional information to the screen (default is False).

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
        'start': start.strftime('%Y%m%d'),
        'end': end.strftime('%Y%m%d'),
        'inseason': str(int(inseason)),
        'vars': var,
        'stns': stn,
    }

    result = requests.get(url, params=data).text

    f = StringIO(result)
    knmi_df, variables = read_knmi_daily_rainfall(f, var, verbose=verbose)
    if int(stn) not in knmi_df.STN.unique():
        raise ValueError(f'KNMI station {stn} not recognized please provide '
                         'valid rainfall station number')

    return knmi_df[[var]], variables


def _read_knmi_header(f, verbose=False):

    variables = {}
    line = f.readline()
    if 'DOCTYPE HTML PUBLIC' in line:
        if verbose:
            print(f.read())
        raise ValueError('Internal Server Error')
    for iline in range(500):
        if ' = ' in line:
            line = line.lstrip(' #').strip('\n')
            varDes = line.split(' = ')
            variables[varDes[0].strip()] = varDes[1].strip()

        if 'STN,YY' in line:
            header = line.split(',')
            header = [item.lstrip().rstrip() for item in header]
            break

        line = f.readline()

    if iline > 498:
        raise ValueError('cannot read measurements from file')

    return f, variables, header


def _transform_variables(df, variables, verbose=False):

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
            if verbose:
                print(f'transform {key}, {value} from 0.1 to 1')
            # transform 0.1 to 1
            df[key] = df[key] * 0.1
            value = value.replace('0.1 ', '')
        if ' tiende ' in value:
            if verbose:
                print(f'transform {key}, {value} from 0.1 to 1')
            # transform 0.1 to 1
            df[key] = df[key] * 0.1
            value = value.replace(' tiende ', ' ')
        if ' mm' in value:
            if verbose:
                print(f'transform {key}, {value} from mm to m')
            # transform mm to m
            df[key] = df[key] * 0.001
            value = value.replace(' mm', ' m')
        if ' millimeters' in value:
            if verbose:
                print(f'transform {key}, {value} from mm to m')
            # transform mm to m
            df[key] = df[key] * 0.001
            value = value.replace(' millimeters', ' m')
        # Store new variable
        variables[key] = value

    return df, variables


def read_knmi_daily_rainfall(f, var, verbose=False):

    f, variables, header = _read_knmi_header(f)

    df = pd.read_csv(f, header=None, names=header, na_values='     ')
    f.close()

    df.set_index(pd.to_datetime(df.YYYYMMDD, format='%Y%m%d'),
                 inplace=True)
    df = df.drop('YYYYMMDD', axis=1)

    # sometimes the last row is messed up, check for that and remove it
    if df.iloc[-1].isna().any():
        if verbose:
            print('last row contains no data, remove last row')
        df = df.drop(index=df.index[-1])
        df.loc[:, var] = df[var].astype(float)

    # daily precipitation amount in 0.1 mm over the period 08.00
    # preceding day - 08.00 UTC present day
    df.index = df.index + pd.to_timedelta(8, unit='h')

    # from UT to UT+1 (standard-time in the Netherlands)
    df.index = df.index + pd.to_timedelta(1, unit='h')

    df, variables = _transform_variables(df, variables, verbose)

    return df, variables


def _read_station_location(f, verbose=False):

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

            if ':' in f.readline():
                raise ValueError(
                    'KNMI station number not recognized please provide '
                    'valid meteo station number')
            break

        line = f.readline()

    if stations is None:
        if verbose:
            print('could not find stations')

    return f, stations


def get_knmi_daily_meteo(url, stn, var, start, end, inseason, verbose=False):
    """download and read knmi daily meteo data


    Parameters
    ----------
    url : str
        download url.
    stn : str
        station number.
    var : str
        e.g. 'EV24'.
    start : pd.TimeStamp
        start time of observations.
    end : pd.TimeStamp
        end time of observations.
    inseason : boolean
        flag to obtain inseason data.
    verbose : boolean, optional
        Print additional information to the screen (default is False).

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
        'vars': var,
        'stns': stn,
    }

    result = requests.get(url, params=data).text

    f = StringIO(result)
    knmi_df, variables, stations = read_knmi_daily_meteo(f)

    return knmi_df[[var]], variables, stations


def read_knmi_daily_meteo(f, verbose=False):

    f, stations = _read_station_location(f, verbose)
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

    df, variables = _transform_variables(df, variables, verbose)

    return df, variables, stations

   # return knmi_series


def get_knmi_hourly(url, stn, var, start, end):

    data = {
        'start': start.strftime('%Y%m%d') + '01',
        'end': end.strftime('%Y%m%d') + '24',
        'vars': var,
        'stns': stn,
    }

    result = requests.get(url, params=data).text
    f = StringIO(result)
    knmi_series = read_knmi_hourly(f)

    return knmi_series


def read_knmi_hourly(f):
    raise NotImplementedError('work in progress')
    knmi_df = pd.DataFrame()
    return knmi_df


def get_knmi_timeseries_xy(x, y, var, start, end, fill_missing_obs=True,
                           interval='daily', inseason=False,
                           raise_exceptions=False,
                           verbose=False):

    # get station
    stations = get_stations(variable=var)
    stn = get_nearest_stations_xy(x, y, var, stations=stations)[0]

    # download data
    if fill_missing_obs:
        knmi_df, variables, station_meta = \
            fill_missing_measurements(stn, var, start, end,
                                      interval, raise_exceptions,
                                      verbose=verbose)
    else:
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, var, start, end,
                               interval, inseason, raise_exceptions,
                               verbose=verbose)

    meta = station_meta.to_dict()
    meta.update(variables)

    # set metadata
    name = var + ' ' + stations.loc[stn, 'naam']
    x = stations.loc[stn, 'x']
    y = stations.loc[stn, 'y']
    meta.update({'x': x, 'y': y, 'station': stn, 'name': name})

    return knmi_df, meta


def get_knmi_timeseries_stn(stn, var, start, end,
                            fill_missing_obs=True, interval='daily',
                            inseason=False, raise_exceptions=False,
                            verbose=False):

    # get station
    stations = get_stations(variable=var)

    # download data
    if fill_missing_obs:
        knmi_df, variables, station_meta = \
            fill_missing_measurements(stn, var, start, end,
                                      interval, raise_exceptions,
                                      verbose=verbose)
    else:
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, var, start, end,
                               interval, inseason, raise_exceptions,
                               verbose=verbose)

    meta = station_meta.to_dict()
    meta.update(variables)

    # set metadata
    name = var + ' ' + stations.loc[stn, 'naam']
    x = stations.loc[stn, 'x']
    y = stations.loc[stn, 'y']
    meta.update({'x': x, 'y': y, 'station': stn, 'name': name})

    return knmi_df, meta


def add_missing_indices(knmi_df, stn, start, end, verbose=False):
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
    verbose : boolean, optional
        Print additional information to the screen (default is False).

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
        if verbose:
            print(
                f'station {stn} has no measurements before {knmi_df.index[0]}')
    if (knmi_df.index[-1] - end).days < 2:
        new_end = knmi_df.index[-1]
    else:
        new_end = pd.Timestamp(year=end.year, month=end.month, day=end.day,
                               hour=knmi_df.index[-1].hour,
                               minute=knmi_df.index[-1].minute,
                               second=knmi_df.index[-1].second)
        if verbose:
            print(
                f'station {stn} has no measurements after {knmi_df.index[-1]}')

    # add missing indices
    new_index = pd.date_range(new_start, new_end, freq='D')
    knmi_df = knmi_df.reindex(new_index)

    return knmi_df


def fill_missing_measurements(stn, var='RD', start=None, end=None,
                              interval='daily',
                              raise_exceptions=False, verbose=False):
    """fill missing measurements in knmi data


    Parameters
    ----------
    stn : int or str
        measurement station.
    interval : str, optional
        desired time interval for observations. The default is 'daily'.
    var : str, optional
        observation type. The default is 'RD'.
    start : str, datetime or None, optional
        start date of observations. The default is None.
    end : str, datetime or None, optional
        end date of observations. The default is None.
    raise_exceptions : bool, optional
        if True you get errors when no data is returned. The default is True.
    verbose : boolean, optional
        Print additional information to the screen (default is False).

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

    if type(var) is not str:
        raise (TypeError('Only one variable supported for now'))
    # get the location of the stations
    stations = get_stations(variable=var)

    if start is None:
        start = pd.Timestamp(pd.Timestamp.today().year, 1, 1)
    else:
        start = pd.to_datetime(start)
    if end is None:
        end = pd.Timestamp.today()
    else:
        end = pd.to_datetime(end)

    if verbose:
        print('Download ' + var + ' from ' +
              str(stn) + ' ' + stations.loc[stn, 'naam'])
    knmi_df, variables, station_meta = \
        download_knmi_data(stn, var, start=start,
                           end=end, interval=interval,
                           inseason=False,
                           raise_exceptions=raise_exceptions,
                           verbose=verbose)

    # if the first station cannot be read, read another station as the first
    ignore = [stn]
    while knmi_df.empty:
        stn = get_nearest_station_df(
            stations.loc[[stn]], variable=var, ignore=ignore)[0]
        if verbose:
            print('Download ' + var + ' from ' +
                  str(stn) + ' ' + stations.loc[stn, 'naam'])
        knmi_df, variables, station_meta = \
            download_knmi_data(stn, var, start=start,
                               end=end, interval=interval,
                               inseason=False,
                               raise_exceptions=raise_exceptions,
                               verbose=verbose)
        ignore.append(stn)

    # find missing values
    knmi_df = add_missing_indices(knmi_df, stn, start, end, verbose)

    missing = knmi_df[var].isna()
    if verbose:
        print(f'station {stn} has {missing.sum()} missing measurements')

    # fill missing values
    while np.any(missing) and not np.all(missing):
        stn_comp = get_nearest_station_df(
            stations.loc[[stn]], variable=var, ignore=ignore)
        if verbose:
            print(f'trying to fill {missing.sum()} '
                  f'measurements with station {stn_comp}')
        if stn_comp is None:
            if verbose:
                print('could not fill all missing measurements there are '
                      'no stations left to check')
            missing[:] = False
            break
        else:
            stn_comp = stn_comp[0]
        knmi_df_comp, _, __ = \
            download_knmi_data(stn_comp, var,
                               start=start, end=end,
                               interval=interval,
                               inseason=False,
                               raise_exceptions=raise_exceptions,
                               verbose=verbose)

        if knmi_df_comp.empty:
            if verbose:
                print(f'station {stn_comp} cannot be downloaded')
        else:
            # dropnans from new data
            knmi_df_comp = knmi_df_comp.loc[~knmi_df_comp[var].isna(), :]
            # get index of missing data in original timeseries
            missing_idx = missing.loc[missing].index
            # if any missing are in the new data, update
            if missing_idx.isin(knmi_df_comp.index).any():
                # index for missing but in newly downloaded data
                ix_idx = missing_idx.intersection(knmi_df_comp.index)
                # update missing data
                knmi_df.loc[ix_idx, var] = knmi_df_comp.loc[ix_idx, var]
                # add source station number
                knmi_df.loc[ix_idx, 'station_opvulwaarde'] = str(stn_comp)

        missing = knmi_df[var].isna()
        ignore.append(stn_comp)

    return knmi_df, variables, station_meta
