import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from shapely.geometry import box
from tqdm import tqdm

from ..util import get_transformer28992

logger = logging.getLogger(__name__)

GHCN_STATIONS_URL = "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
GHCN_DAILY_URL = "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/{station_id}.dly"

# GHCN daily depth-like elements are reported in 0.1 mm. Convert to m.
_DEPTH_ELEMENTS_TO_M = {
    "PRCP",
    "SNOW",
    "SNWD",
    "WESD",
    "WESF",
    "EVAP",
}


def get_stations(extent=None):
    """Get GHCN stations within a specific extent.

    Parameters
    ----------
    extent : list, tuple or numpy-array, optional
        get GHCN stations within this extent [xmin, xmax, ymin, ymax]. Coordinates should
        be in WGS84 (EPSG:4326). The default is None.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame containing the GHCN stations within the specified extent.
    """
    url = GHCN_STATIONS_URL
    colspecs = [
        (0, 11),  # ID
        (12, 20),  # LATITUDE
        (21, 30),  # LONGITUDE
        (31, 37),  # ELEVATION
        (38, 40),  # STATE
        (41, 71),  # NAME
        (72, 75),  # GSN FLAG
        (76, 79),  # HCN/CRN FLAG
        (80, 85),  # WMO ID
    ]

    colnames = [
        "id",
        "latitude",
        "longitude",
        "elevation",
        "state",
        "name",
        "gsn_flag",
        "hcn_crn_flag",
        "wmo_id",
    ]

    stations = pd.read_fwf(url, colspecs=colspecs, names=colnames)
    stations = stations.set_index("id")
    geometry = gpd.points_from_xy(stations["longitude"], stations["latitude"])
    stations_gdf = gpd.GeoDataFrame(stations, geometry=geometry, crs="EPSG:4326")
    if extent is not None:
        # Extent format is [xmin, ymin, xmax, ymax] in this function.
        polygon = box(extent[0], extent[2], extent[1], extent[3])
        stations_gdf = stations_gdf[stations_gdf.intersects(polygon)]

    return stations_gdf


def get_station_data(station_id, element=None, start_date=None, end_date=None):
    """Get daily GHCN data for one station.

    Parameters
    ----------
    station_id : str
        GHCN station ID for which to download the data.
    element : str, list of str, or None, optional
        GHCN element(s) to download (e.g. 'PRCP', 'TMAX', 'TMIN').
        If None all available elements per station are downloaded.
        Depth-like elements (e.g. PRCP, SNOW, SNWD, WESD, WESF, EVAP)
        are converted from 0.1 mm to m.
        The default is None.
    start_date : str or None, optional
        start date of observations (e.g. '2020-01-01'). The default is None.
    end_date : str or None, optional
        end date of observations (e.g. '2021-12-31'). The default is None.

    Returns
    -------
    DataFrame
        DataFrame containing the daily GHCN data for the specified station and element(s).

    Notes
    -----
    Timestamps are shifted by +1 day so the index represents the end of the
    daily period, consistent with KNMI daily indexing in hydropandas.
    """
    url = GHCN_DAILY_URL.format(station_id=station_id)
    colspecs = [
        (0, 11),  # ID
        (11, 15),  # YEAR
        (15, 17),  # MONTH
        (17, 21),  # ELEMENT
    ] + [(21 + i * 8, 26 + i * 8) for i in range(31)]  # VALUE1-31

    colnames = ["id", "year", "month", "element"] + [f"value{i}" for i in range(1, 32)]

    data = pd.read_fwf(url, colspecs=colspecs, names=colnames)
    data = data.melt(
        id_vars=["id", "year", "month", "element"],
        var_name="day_col",
        value_name="value",
    )
    data["day"] = data["day_col"].str.replace("value", "", regex=False).astype(int)
    data["date"] = pd.to_datetime(
        {
            "year": data["year"],
            "month": data["month"],
            "day": data["day"],
        },
        errors="coerce",
    )
    data = data.dropna(subset=["date"])
    data["date"] = data["date"] + pd.Timedelta(1, "d")
    data = data.drop(columns=["year", "month", "day_col", "day"])
    if start_date is not None:
        data = data[data["date"] >= pd.to_datetime(start_date)]
    if end_date is not None:
        data = data[data["date"] <= pd.to_datetime(end_date)]
    if element is not None:
        data = data[data["element"] == element]

    return data


def get_obs_list_from_extent(
    extent,
    ObsClass,
    elements=None,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=True,
    crs=4326,
):
    """Get GHCN observations within a specific extent.

    Parameters
    ----------
    extent : list, tuple or numpy-array
        get GHCN stations within this extent [xmin, xmax, ymin, ymax]
    ObsClass : type
        class of the observations, e.g. MeteoObs or PrecipitationObs
    elements : str, list of str, or None, optional
        GHCN element(s) to download (e.g. 'PRCP', 'TMAX', 'TMIN').
        If None all available elements per station are downloaded.
        Depth-like elements (e.g. PRCP, SNOW, SNWD, WESD, WESF, EVAP)
        are converted from 0.1 mm to m.
        The default is None.
    tmin : str or None, optional
        start date of observations (e.g. '2020-01-01'). The default is None.
    tmax : str or None, optional
        end date of observations (e.g. '2021-12-31'). The default is None.
    only_metadata : bool, optional
        if True download only station metadata, significantly faster.
        The default is False.
    keep_all_obs : bool, optional
        if False, only observations with measurements are kept.
        The default is True.
    crs : str, int or pyproj.CRS, optional
        The coordinate reference system of the extent, this crs is also
        used for the observations. The default is 4326 (WGS84).

    Returns
    -------
    list
        list with Obs objects
    """
    # transform extent corners to WGS84 for station selection
    crs = pyproj.CRS(crs)
    if crs != pyproj.CRS(4326):
        transformer = get_transformer28992(crs, pyproj.CRS(4326))

        lon_min, lat_min = transformer.transform(extent[0], extent[2])
        lon_max, lat_max = transformer.transform(extent[1], extent[3])

        transformer_from_wgs84 = get_transformer28992(pyproj.CRS(4326), crs)
    else:
        # standard hydropandas extent: [xmin, xmax, ymin, ymax]
        lon_min, lon_max = extent[0], extent[1]
        lat_min, lat_max = extent[2], extent[3]
        transformer_from_wgs84 = None

    stations_gdf = get_stations(extent=[lon_min, lon_max, lat_min, lat_max])

    if stations_gdf.empty:
        logger.warning(f"No GHCN stations found within extent {extent}")
        return []

    logger.info(f"downloading GHCN data from {len(stations_gdf)} stations")

    if isinstance(elements, str):
        elements = [elements]

    obs_list = []
    for station_id, row in tqdm(
        stations_gdf.iterrows(), total=len(stations_gdf), desc="station"
    ):
        lon = row.geometry.x
        lat = row.geometry.y
        if transformer_from_wgs84 is not None:
            x, y = transformer_from_wgs84.transform(lon, lat)
        else:
            x, y = lon, lat

        meta = {
            "station": station_id,
            "name": row.get("name", station_id),
            "elevation": row.get("elevation", np.nan),
            "latitude": lat,
            "longitude": lon,
            "x": x,
            "y": y,
            "crs": crs,
        }

        if only_metadata:
            o = ObsClass(
                name=station_id,
                x=x,
                y=y,
                station=station_id,
                source="GHCN",
                crs=crs,
                meta=meta,
            )
            obs_list.append(o)
            continue

        data = get_station_data(station_id, start_date=tmin, end_date=tmax)

        # replace GHCN missing-value flag with NaN
        data["value"] = data["value"].replace(-9999, np.nan)

        # filter by requested elements
        if elements is not None:
            data = data[data["element"].isin(elements)]

        if data.empty:
            if keep_all_obs:
                o = ObsClass(
                    name=station_id,
                    x=x,
                    y=y,
                    station=station_id,
                    source="GHCN",
                    crs=crs,
                    meta=meta,
                )
                obs_list.append(o)
            continue

        for element, element_data in data.groupby("element"):
            # GHCN depth-like elements are reported in 0.1 mm. Convert to meters.
            if element in _DEPTH_ELEMENTS_TO_M:
                element_data = element_data.copy()
                element_data["value"] = element_data["value"] * 1e-4

            ts = (
                element_data.set_index("date")[["value"]]
                .rename(columns={"value": element})
                .sort_index()
            )

            meta["meteo_var"] = element
            meta["unit"] = "m" if element in _DEPTH_ELEMENTS_TO_M else "unknown"
            o = ObsClass(
                ts,
                name=f"{station_id}_{element}",
                x=x,
                y=y,
                crs=crs,
                station=station_id,
                meteo_var=element,
                source="GHCN",
                unit=meta["unit"],
                meta=meta,
            )
            obs_list.append(o)

    return obs_list
