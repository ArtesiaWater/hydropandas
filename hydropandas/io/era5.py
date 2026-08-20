import logging

import numpy as np
import pandas as pd
import requests
from pyproj import Transformer
from tqdm import tqdm

logger = logging.getLogger(__name__)

ERA5_ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"

_SOURCE_TO_MODEL = {
    "era5": "era5",
    "era5_land": "era5_land",
    "era5_hourly": "era5",
    "era5_seamless": "era5_seamless",
}

# Open-Meteo ERA5 precipitation-like variables are provided in mm.
_MM_TO_M_VARS = {
    "precipitation_sum",
    "rain_sum",
    "showers_sum",
    "snowfall_sum",
}


def _extent_to_wgs84(extent, epsg):
    if epsg == 4326:
        return extent[0], extent[1], extent[2], extent[3]

    transformer = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    lon_min, lat_min = transformer.transform(extent[0], extent[2])
    lon_max, lat_max = transformer.transform(extent[1], extent[3])
    return lon_min, lon_max, lat_min, lat_max


def _xy_to_wgs84(xy, epsg):
    x, y = xy
    if epsg == 4326:
        return float(x), float(y)

    transformer = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(x, y)
    return float(lon), float(lat)


def _axis_values(vmin, vmax, step):
    if np.isclose(vmin, vmax):
        return np.array([vmin])

    start = min(vmin, vmax)
    end = max(vmin, vmax)
    values = np.arange(start, end + step * 0.5, step)
    return values


def _get_grid_points(extent, epsg=4326, grid_size=0.25):
    lon_min, lon_max, lat_min, lat_max = _extent_to_wgs84(extent, epsg)
    lons = _axis_values(lon_min, lon_max, grid_size)
    lats = _axis_values(lat_min, lat_max, grid_size)

    return [(float(lat), float(lon)) for lon in lons for lat in lats]


def _resolve_source_and_interval(source, interval):
    if source not in _SOURCE_TO_MODEL:
        raise ValueError(
            "source should be one of 'era5', 'era5_land', 'era5_hourly', 'era5_seamless'"
        )

    model = _SOURCE_TO_MODEL[source]
    if source == "era5_hourly":
        interval = "hourly"

    if interval not in ("daily", "hourly"):
        raise ValueError("interval should be 'daily' or 'hourly'")

    return model, interval


def request_era5_api(
    latitude,
    longitude,
    tmin,
    tmax,
    variables,
    model="era5",
    interval="daily",
    timeout=120,
):
    params = {
        "latitude": latitude,
        "longitude": longitude,
        "start_date": pd.Timestamp(tmin).strftime("%Y-%m-%d"),
        "end_date": pd.Timestamp(tmax).strftime("%Y-%m-%d"),
        "timezone": "UTC",
        "models": model,
    }

    if interval == "daily":
        params["daily"] = ",".join(variables)
    elif interval == "hourly":
        params["hourly"] = ",".join(variables)
    else:
        raise ValueError("interval should be 'daily' or 'hourly'")

    response = requests.get(ERA5_ARCHIVE_URL, params=params, timeout=timeout)
    response.raise_for_status()
    return response.json()


def get_obs_list_from_extent(
    ObsClass,
    extent=None,
    xy=None,
    variables=("precipitation_sum",),
    source="era5_seamless",
    tmin=None,
    tmax=None,
    interval="daily",
    only_metadata=False,
    keep_all_obs=False,
    epsg=4326,
    grid_size=0.25,
    timeout=120,
    max_points=200,
):
    """Get ERA5 observations from a regular grid within an extent or at one xy.

    Parameters
    ----------
    extent : list, tuple, numpy-array or None, optional
        get ERA5 grid points within this extent [xmin, xmax, ymin, ymax]
    xy : tuple, list or None, optional
        single point coordinates (x, y). If provided, extent is ignored and
        ERA5 data is downloaded for this point only.
    ObsClass : type
        class of the observations, e.g. MeteoObs
    variables : tuple, list or str, optional
        ERA5 variable(s) to download, default is ('precipitation_sum',)
    source : str, optional
        ERA5 product selection. Options are 'era5', 'era5_land',
        'era5_hourly', and 'era5_seamless'. The default is
        'era5_seamless'.
    tmin : str or None, optional
        start date of observations. If None, one month before today is used.
    tmax : str or None, optional
        end date of observations. If None, today is used.
    interval : str, optional
        one of 'daily' or 'hourly', by default 'daily'
        Returned timestamps are in UTC and shifted to the end of each
        aggregation period: +1 day for daily and +1 hour for hourly.
    only_metadata : bool, optional
        if True download only metadata, significantly faster.
        The default is False.
    keep_all_obs : bool, optional
        if False, only observations with measurements are kept.
        The default is False.
    epsg : int, optional
        epsg code of the supplied extent. Returned observation x/y
        coordinates are also in this CRS. The default is 4326 (WGS84).
    grid_size : float, optional
        ERA5 grid sampling size in degrees, default is 0.25
    timeout : int, optional
        request timeout in seconds, default is 120
    max_points : int, optional
        maximum number of grid points to download, default is 200

    Returns
    -------
    list
        list with Obs objects
    """
    if isinstance(variables, str):
        variables = (variables,)

    model, interval = _resolve_source_and_interval(source, interval)

    if tmin is None:
        tmin = pd.Timestamp.today().normalize() - pd.Timedelta(days=30)
    if tmax is None:
        tmax = pd.Timestamp.today().normalize()

    if xy is not None:
        lon, lat = _xy_to_wgs84(xy, epsg)
        points = [(lat, lon)]
    elif extent is not None:
        points = _get_grid_points(extent, epsg=epsg, grid_size=grid_size)
    else:
        raise ValueError("specify extent or xy for ERA5 data")

    if len(points) > max_points:
        raise ValueError(
            f"selected extent results in {len(points)} ERA5 points; "
            f"increase grid_size or max_points (current max_points={max_points})"
        )

    if epsg != 4326:
        transformer_from_wgs84 = Transformer.from_crs(
            "EPSG:4326", f"EPSG:{epsg}", always_xy=True
        )
    else:
        transformer_from_wgs84 = None

    obs_list = []
    for lat, lon in tqdm(points, total=len(points), desc="era5 point"):
        if transformer_from_wgs84 is not None:
            x, y = transformer_from_wgs84.transform(lon, lat)
        else:
            x, y = lon, lat

        base_meta = {
            "latitude": lat,
            "longitude": lon,
            "x": x,
            "y": y,
            "epsg": epsg,
            "source": "ERA5",
            "era5_source": source,
        }

        if only_metadata:
            for variable in variables:
                obs_list.append(
                    ObsClass(
                        name=f"ERA5_{lat:.4f}_{lon:.4f}_{variable}",
                        x=x,
                        y=y,
                        source="ERA5",
                        meteo_var=variable,
                        unit="",
                        meta=base_meta | {"meteo_var": variable},
                    )
                )
            continue

        data = request_era5_api(
            latitude=lat,
            longitude=lon,
            tmin=tmin,
            tmax=tmax,
            variables=variables,
            model=model,
            interval=interval,
            timeout=timeout,
        )

        values_key = "daily" if interval == "daily" else "hourly"
        units_key = "daily_units" if interval == "daily" else "hourly_units"

        if values_key not in data or "time" not in data[values_key]:
            logger.warning(f"no ERA5 data returned for lat={lat}, lon={lon}")
            continue

        index = pd.to_datetime(data[values_key]["time"])
        if interval == "daily":
            # Keep daily timestamps aligned with KNMI convention in hydropandas.
            index = index + pd.Timedelta(days=1)
        elif interval == "hourly":
            # Use end-of-hour timestamps for consistency with daily behavior.
            index = index + pd.Timedelta(hours=1)

        for variable in variables:
            if variable not in data[values_key]:
                if keep_all_obs:
                    obs_list.append(
                        ObsClass(
                            name=f"ERA5_{lat:.4f}_{lon:.4f}_{variable}",
                            x=x,
                            y=y,
                            source="ERA5",
                            meteo_var=variable,
                            unit="",
                            meta=base_meta | {"meteo_var": variable},
                        )
                    )
                continue

            values = pd.Series(data[values_key][variable], index=index, dtype=float)
            unit = data.get(units_key, {}).get(variable, "")

            if variable in _MM_TO_M_VARS and unit == "mm":
                values = values * 1e-3
                unit = "m"

            ts = pd.DataFrame({variable: values})
            ts = ts.dropna(how="all")

            if ts.empty and not keep_all_obs:
                continue

            meta = base_meta | {
                "meteo_var": variable,
                "unit": unit,
                "interval": interval,
            }
            obs_list.append(
                ObsClass(
                    ts,
                    name=f"ERA5_{lat:.4f}_{lon:.4f}_{variable}",
                    x=x,
                    y=y,
                    source="ERA5",
                    meteo_var=variable,
                    unit=unit,
                    meta=meta,
                )
            )

    return obs_list
