from functools import lru_cache
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

URL_WOW_KNMI = "https://wow.knmi.nl/sites"

meteo_vars_wow = (
    "rain_rate",
    "temperature",
    # "weather_code",
    "wind",
    "humidity",
    "pressure_msl",
)

meteo_vars_wow_translate = {
    "neerslagintensiteit (mm/uur)": "precipitation [mm]",
    "temperatuur (C)": "temperature [C]",
    "windsnelheid (m/s)": "windspeed [m/s]",
    "relatieve vochtigheid  (%)": "relative humidity [%]",
    "luchtdruk (hPa)": "air pressure [hPa]",
}

wow_filters = ("wow_observations", "official_observations")


def get_wow_stations(
    meteo_var: str = "rain_rate",
    date: Optional[pd.Timestamp] = None,
    bbox: Optional[List[float]] = None,
    obs_filter: Optional[str] = None,
) -> pd.DataFrame:
    """
    Get a DataFrame with the observation stations from WOW-KNMI.

    Parameters
    ----------
    meteo_var : str, optional
        The meteorological variable to query, one of 'rain_rate',
        'temperature', 'wind', 'humidity', 'pressure_msl' by default
        'rain_rate'.
    date : Optional[pd.Timestamp], optional
        The date and time of the observations, by default the latest available
        (10 minutes ago).
    bbox : Optional[List[float]], optional
        The bounding box of the spatial query in the format [min_lon, min_lat,
        max_lon, max_lat], by default the extent of the Netherlands.
    obs_filter : Optional[str], optional
        The filter to apply to the observations, one of 'wow_observations',
        'official_observations', by default None which selects all stations.

    Returns
    -------
    pd.DataFrame
        A dataframe with the station information and observations.

    Raises
    ------
    ValueError
        If meteo_var or obs_filter are not valid choices.
    requests.HTTPError
        If the request to the WOW-KNMI API fails.
    """
    if meteo_var not in meteo_vars_wow:
        raise ValueError(f"{meteo_var} must be one of {meteo_vars_wow}")

    if date is None:
        date = pd.Timestamp.today() - pd.Timedelta(days=0, hours=0, minutes=10)
    datestr = _wow_strftime(date)

    if bbox is None:
        bbox = [2.5, 50.0, 8.0, 54.0]  # Netherlands extent
    bboxstr = (
        str(bbox[0]) + "%20" + str(bbox[1]) + "," + str(bbox[2]) + "%20" + str(bbox[3])
    )  # lat lon,lat lon

    if obs_filter is None:
        obs_filter = ""
    else:
        if obs_filter not in wow_filters:
            raise ValueError(f"{obs_filter} must be one of {wow_filters}")

    try:
        url = (
            f"{URL_WOW_KNMI}?bbox={bboxstr}&layer={meteo_var}"
            f"&filter={obs_filter}&date={datestr}"
        )
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        sites = r.json()["sites"]
        lat_lon = np.array([x["geo"]["coordinates"] for x in sites])
        xy = pd.DataFrame(
            np.column_stack([lat_lon[:, 1], lat_lon[:, 0]]),
            columns=["x", "y"],
        )
        stations = pd.concat([pd.DataFrame(sites), xy], axis=1).set_index(
            "id", drop=False
        )
    except requests.HTTPError as ex:
        raise ex

    return stations


def _wow_strftime(timestamp: pd.Timestamp) -> str:
    """Convert pandast Timestamp to wow string time"""
    return (
        timestamp.strftime("%Y-%m-%dT%H")
        + "%3A"
        + timestamp.strftime("%M")[0]
        + "0%3A00"
    )


def get_wow(
    meteo_var: str,
    stn: str,
    xy: List[float] = None,
    start: pd.Timestamp = None,
    end: pd.Timestamp = None,
) -> Tuple[pd.DataFrame, dict]:
    """
    Get weather observations and metadata from a WOW-KNMI station based on the
    station name or lat, lon location

    Parameters
    ----------
    meteo_var : str
        The meteorological variable to query, one of 'rain_rate',
        'temperature', 'wind', 'humidity', 'pressure_msl'.
    stn : str
        The station ID of the WOW-KNMI station.
    xy : List[float], optional
        The coordinates of the location to query in the format [lon, lat], by
        default None. If specified, the nearest station within 1 degree will be
        used.
    start : pd.Timestamp, optional
        The start date and time of the observations, by default None. If None,
        the januari first of the year before will be used.
    end : pd.Timestamp, optional
        The end date and time of the observations, by default None. If None,
        the yesterday will be used.

    Returns
    -------
    Tuple[pd.DataFrame, dict]
        A tuple of a dataframe with the observations and a dictionary with the
        station metadata.

    Raises
    ------
    ValueError
        If both stn and xy are None or if meteo_var is not valid.
    """

    if stn is None and xy is None:
        raise ValueError("specify KNMI station (stn) or coordinates (xy)")

    if stn is None and xy is not None:
        bbox = [xy[0] - 1, xy[1] - 1, xy[0] + 1, xy[1] + 1]  # lat lon lat lon
        stations = get_wow_stations(
            meteo_var=meteo_var, date=start, bbox=bbox, obs_filter=None
        )
        nearest_stations = get_nearest_station_xy([xy], stations=stations, ignore=None)
        stn = nearest_stations[0]

    metadata = get_wow_metadata(stn=stn)
    measurements = get_wow_measurements(
        stn=stn, meteo_var=meteo_var, start=start, end=end
    )

    return measurements, metadata


def get_wow_metadata(stn: str) -> dict:
    """Get metadata from a WOW-KNMI station."""
    r = requests.get(f"{URL_WOW_KNMI}/{stn}")
    r.raise_for_status()
    metadata = r.json()
    return metadata


@lru_cache()
def get_wow_measurements(
    stn: str,
    meteo_var: str,
    start: pd.Timestamp = None,
    end: pd.Timestamp = None,
) -> pd.DataFrame:
    """
    Get measurements from a WOW-KNMI station.

    Parameters
    ----------
    stn : str
        The station ID of the WOW-KNMI station.
    meteo_var : str
        The meteorological variable to query, one of 'rain_rate',
        'temperature', 'wind', 'humidity', 'pressure_msl'.
    start : pd.Timestamp, optional
        The start date and time of the observations, by default None. If None,
        the januari first of the year before will be used.
    end : pd.Timestamp, optional
        The end date and time of the observations, by default None. If None,
        the yesterday will be used.

    Returns
    -------
    pd.DataFrame
        A dataframe with the measurements.

    Raises
    ------
    ValueError
        If meteo_var is not valid.
    requests.HTTPError
        If the request to the WOW-KNMI API fails.
    """
    if meteo_var not in meteo_vars_wow:
        raise ValueError(f"{meteo_var} must be one of {meteo_vars_wow}")

    start, end = _start_end_to_datetime(start, end)
    if (end - start) > pd.to_timedelta(365, unit="D"):
        t = list(pd.date_range(start, end, freq="365D"))
        if t[-1] != end:
            t.append(end)
        dfs = []
        for i in range(len(t) - 1):
            dfs.append(get_wow_measurements(stn, meteo_var, start=t[i], end=t[i + 1]))
        df = pd.concat(dfs)
        df = df[~df.index.duplicated()]
        return df

    startstr = _wow_strftime(start)
    endstr = _wow_strftime(end)
    # get station measurements
    try:
        url = (
            f"{URL_WOW_KNMI}/{stn}/export?start={startstr}"
            f"&end={endstr}&layer={meteo_var}"
        )
        meas = pd.read_csv(url, delimiter=";", index_col=["datum"])
        meas.index = pd.DatetimeIndex(pd.to_datetime(meas.index.values), name="date")
        measurements = meas.rename(
            columns={
                x: meteo_vars_wow_translate[x.replace(f" [{stn}]", "")]
                for x in meas.columns
            }
        )
    except requests.HTTPError as ex:
        raise ex

    if meteo_var == "rain_rate":
        weights = (
            (measurements.index[1:] - measurements.index[:-1])
            / pd.Timedelta(1, unit="H")
        ).values
        measurements = measurements.iloc[0:-1].multiply(weights, axis=0)

    return measurements


def get_nearest_station_xy(
    xy: List[List[float]], stations: pd.DataFrame, ignore: List[str] = None
) -> List[str]:
    """find the stations that measure closest to the given
    longitude and latitude coordinates.

    Parameters
    ----------
    xy : list or numpy array, optional
        longitude, latitude coordinates of the locations. e.g. [[10,25], [5,25]]
    stations : pandas DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
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

    locations = pd.DataFrame(data=xy, columns=["lon", "lat"])

    stns = get_nearest_station_df(
        locations=locations,
        stations=stations,
        xcol="lon",
        ycol="lat",
        ignore=ignore,
    )

    return stns


def get_nearest_station_df(
    locations: pd.DataFrame,
    stations: pd.DataFrame,
    xcol: str = "x",
    ycol: str = "y",
    ignore: List[str] = None,
) -> List[str]:
    """Find the nearest stations that measure 'meteo_var' closest to the
    coordinates in 'locations'.

    Parameters
    ----------
    locations : pandas.DataFrame
        DataFrame containing x and y coordinates
    stations : pandas DataFrame, optional
        if None stations will be obtained using the get_stations function.
        The default is None.
    xcol : str
        The name of the column in the locations dataframe with the x values, by
        default 'x'. If lon in xcol, longitude is assumed.
    ycol : str
        The name of the column in the locations dataframe with the y values, by
        default 'y'. If lat in ycol, latitude is assumed.
    ignore : list, optional
        A list of station IDs to ignore, by default None.

    Returns
    -------
    stns : list
        A list of station IDs that are closest to each location.

    """

    if ignore is not None:
        stations = stations.drop(ignore, inplace=False)
        if stations.empty:
            return None

    xo = pd.to_numeric(locations[xcol])
    xt = pd.to_numeric(stations.x)
    yo = pd.to_numeric(locations[ycol])
    yt = pd.to_numeric(stations.y)

    xh, xi = np.meshgrid(xt, xo)
    yh, yi = np.meshgrid(yt, yo)

    if "lon" in xcol.lower() or "lat" in ycol.lower():
        distances = pd.DataFrame(
            _latlon_distance(yh, xh, yi, xi),
            index=locations.index,
            columns=stations.index,
        )
    else:
        distances = pd.DataFrame(
            np.sqrt((xh - xi) ** 2 + (yh - yi) ** 2),
            index=locations.index,
            columns=stations.index,
        )

    stns = distances.idxmin(axis=1).to_list()

    return stns


def _latlon_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Haversine formula"""
    p = 0.017453292519943295
    hav = (
        0.5
        - np.cos((lat2 - lat1) * p) / 2
        + np.cos(lat1 * p) * np.cos(lat2 * p) * (1 - np.cos((lon2 - lon1) * p)) / 2
    )
    return 12742 * np.arcsin(np.sqrt(hav))


def _start_end_to_datetime(
    start: Optional[str], end: Optional[str]
) -> Tuple[pd.Timestamp]:
    """Convert start and endtime to datetime.

    Parameters
    ----------
    start : str, None
        start time
    end : str, None
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
