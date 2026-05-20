import html
import logging
import re

import pandas as pd
import requests
from pyproj import Transformer
from tqdm import tqdm

logger = logging.getLogger(__name__)

GGMN_WFS_URL = "https://ggis.un-igrac.org/geoserver/wfs"
GGMN_SITE_URL = "https://ggis.un-igrac.org"
GGMN_LEVEL_LAYER = "groundwater:GGMN_Levels_Data"


def _extent_to_wgs84(extent, epsg):
    if epsg == 4326:
        return extent[0], extent[1], extent[2], extent[3]

    transformer = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    lon_min, lat_min = transformer.transform(extent[0], extent[2])
    lon_max, lat_max = transformer.transform(extent[1], extent[3])
    return lon_min, lon_max, lat_min, lat_max


def get_locations_within_extent(extent, epsg=4326, max_locations=200, timeout=120):
    """Get GGMN monitoring locations within an extent."""
    lon_min, lon_max, lat_min, lat_max = _extent_to_wgs84(extent, epsg)

    params = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeNames": GGMN_LEVEL_LAYER,
        "srsName": "EPSG:4326",
        "outputFormat": "application/json",
        "count": max_locations,
        "bbox": f"{lon_min},{lat_min},{lon_max},{lat_max},EPSG:4326",
    }

    response = requests.get(GGMN_WFS_URL, params=params, timeout=timeout)
    response.raise_for_status()
    payload = response.json()

    return payload.get("features", [])


def _parse_measurement_row(row_html):
    time_match = re.search(r'name="time" value="([^"]*)"', row_html)
    value_match = re.search(r'name="value_value"\s+value="([^"]*)"', row_html)
    param_match = re.search(r"<option[^>]*selected[^>]*>([^<]*)</option>", row_html)
    unit_match = re.search(
        r'name="value_unit"[\s\S]*?<option[^>]*selected[^>]*>([^<]*)</option>',
        row_html,
    )

    if time_match is None or value_match is None:
        return None

    time_str = time_match.group(1).strip()
    value_str = value_match.group(1).strip()

    if (time_str == "") or (value_str == ""):
        return None

    try:
        value = float(value_str.replace(",", "."))
    except ValueError:
        return None

    parameter = html.unescape(param_match.group(1).strip()) if param_match else ""
    unit = html.unescape(unit_match.group(1).strip()) if unit_match else ""
    return {
        "time": pd.to_datetime(time_str),
        "value": value,
        "parameter": parameter,
        "unit": unit,
    }


def _parameter_matches(parameter_name, requested_parameter):
    if requested_parameter is None:
        return True

    if isinstance(requested_parameter, str):
        requested = [requested_parameter]
    else:
        requested = list(requested_parameter)

    pname = parameter_name.casefold()
    return any(req.casefold() in pname for req in requested)


def get_level_measurements(
    record_id,
    tmin=None,
    tmax=None,
    parameter=None,
    max_pages=20,
    timeout=120,
):
    """Get groundwater level measurements for a single GGMN record."""
    set_idx = 1
    rows = []

    for _ in range(max_pages):
        url = (
            f"{GGMN_SITE_URL}/groundwater/record/{record_id}/"
            f"WellLevelMeasurement/list?set={set_idx}"
        )
        response = requests.get(
            url,
            headers={"X-Requested-With": "XMLHttpRequest"},
            timeout=timeout,
        )
        response.raise_for_status()
        payload = response.json()

        for item in payload.get("data", []):
            parsed = _parse_measurement_row(item.get("html", ""))
            if parsed is not None and _parameter_matches(
                parsed["parameter"], parameter
            ):
                rows.append(parsed)

        if payload.get("end", True):
            break

        next_set = payload.get("set", set_idx + 1)
        if next_set == set_idx:
            break
        set_idx = next_set

    if not rows:
        return pd.DataFrame(columns=["groundwater_level"]), ""

    df = pd.DataFrame(rows)

    df = df.sort_values("time").drop_duplicates(subset="time", keep="first")
    df = df.set_index("time")[["value"]].rename(columns={"value": "groundwater_level"})

    if tmin is not None:
        df = df.loc[df.index >= pd.to_datetime(tmin)]
    if tmax is not None:
        df = df.loc[df.index <= pd.to_datetime(tmax)]

    parameter_name = rows[0].get("parameter", "")
    unit_match = re.search(r"\[(.*?)\]", parameter_name)
    unit = unit_match.group(1) if unit_match else rows[0].get("unit", "")

    return df, unit


def get_obs_list_from_extent(
    extent,
    ObsClass,
    tmin=None,
    tmax=None,
    parameter=None,
    only_metadata=False,
    keep_all_obs=True,
    epsg=4326,
    max_locations=200,
    max_pages=20,
    timeout=120,
):
    """Get GGMN observations within a specific extent.

    Parameters
    ----------
    extent : list, tuple, numpy-array
        get GGMN locations within this extent [xmin, xmax, ymin, ymax]
    ObsClass : type
        class of the observations, e.g. GroundwaterObs
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    parameter : str, iterable of str, or None, optional
        groundwater-level parameter name filter. Set to None (default)
        to include all available level parameters.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default is False.
    keep_all_obs : bool, optional
        if False, only observations with measurements are kept. The default is True.
    epsg : int, optional
        epsg code of the supplied extent. Returned observation x/y
        coordinates are also in this CRS. The default is 4326 (WGS84).
    max_locations : int, optional
        maximum number of locations to download, by default 200
    max_pages : int, optional
        maximum number of measurement pages per location, by default 20
    timeout : int, optional
        request timeout in seconds, by default 120

    Returns
    -------
    list
        list with Obs objects
    """
    features = get_locations_within_extent(
        extent, epsg=epsg, max_locations=max_locations, timeout=timeout
    )

    if not features:
        logger.warning(f"No GGMN locations found within extent {extent}")
        return []

    if epsg != 4326:
        transformer_from_wgs84 = Transformer.from_crs(
            "EPSG:4326", f"EPSG:{epsg}", always_xy=True
        )
    else:
        transformer_from_wgs84 = None

    obs_list = []
    for feature in tqdm(features, total=len(features), desc="ggmn location"):
        props = feature.get("properties", {})
        coords = feature.get("geometry", {}).get("coordinates", [None, None])
        lon, lat = coords[0], coords[1]

        if (lon is None) or (lat is None):
            continue

        if transformer_from_wgs84 is not None:
            x, y = transformer_from_wgs84.transform(lon, lat)
        else:
            x, y = lon, lat

        record_id = props.get("id")
        name = props.get("ggis_uid") or props.get("name") or f"GGMN_{record_id}"

        meta = {
            "record_id": record_id,
            "ggis_uid": props.get("ggis_uid"),
            "organisation": props.get("organisation"),
            "country": props.get("country"),
            "feature_type": props.get("feature_type"),
            "x": x,
            "y": y,
            "latitude": lat,
            "longitude": lon,
            "epsg": epsg,
            "source": "GGMN",
        }

        if only_metadata:
            obs_list.append(
                ObsClass(
                    name=name,
                    x=x,
                    y=y,
                    source="GGMN",
                    unit="m",
                    meta=meta,
                )
            )
            continue

        if record_id is None:
            if keep_all_obs:
                obs_list.append(
                    ObsClass(
                        name=name,
                        x=x,
                        y=y,
                        source="GGMN",
                        unit="m",
                        meta=meta,
                    )
                )
            continue

        ts, unit = get_level_measurements(
            record_id,
            tmin=tmin,
            tmax=tmax,
            parameter=parameter,
            max_pages=max_pages,
            timeout=timeout,
        )

        if ts.empty and not keep_all_obs:
            continue

        obs_list.append(
            ObsClass(
                ts,
                name=name,
                x=x,
                y=y,
                source="GGMN",
                unit=unit or "m",
                meta=meta,
            )
        )

    return obs_list
