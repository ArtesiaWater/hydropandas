import logging
import math
import pathlib
from concurrent.futures import ThreadPoolExecutor

import geopandas
import pandas as pd
import requests
from pyproj import Transformer
from shapely.geometry import Polygon
from tqdm import tqdm

logger = logging.getLogger(__name__)

# TODO:
# - check transformation from EPSG:28992 to WGS84 (elsewhere in hydropandas we use
#   another definition for EPSG:28992 that is provided in util.py)

# Generic Lizard API endpoint (with 'organisation' as placeholder, following the Lizard documentation
lizard_api_endpoint = "https://{organisation}.lizard.net/api/v4/"


def check_status_obs(metadata, timeseries):
    """Checks if a monitoring tube is still active.

    If there are no measurements in the last 180 days, the monitoring
    tube is considered inactive.

    Parameters
    ----------
    metadata : dict
        metadata of the monitoring tube
    timeseries : pandas DataFrame
        timeseries of the monitoring well

    Returns
    -------
    metadata DataFrame including the status of the monitoring well
    """
    if timeseries.empty:
        metadata["status"] = "no timeseries available"
        return metadata

    last_measurement_date = timeseries.last_valid_index()
    today = pd.to_datetime("today").normalize()

    if today - last_measurement_date < pd.Timedelta(days=180):
        metadata["status"] = "active"

    else:
        metadata["status"] = "inactive"

    return metadata


def extent_to_wgs84_polygon(extent):
    """Translates an extent (xmin, xmax, ymin, ymax) to a polygon with coordinate system
    WGS84.

    Parameters
    ----------
    extent : list or tuple
        extent in epsg 28992 within which the observations are collected.

    Returns
    -------
    polygon of the extent with coordinate system WGS84
    """
    transformer = Transformer.from_crs("EPSG:28992", "WGS84")

    lon_min, lat_min = transformer.transform(extent[0], extent[2])
    lon_max, lat_max = transformer.transform(extent[1], extent[3])

    poly_T = Polygon(
        [(lat_min, lon_min), (lat_max, lon_min), (lat_max, lon_max), (lat_min, lon_max)]
    )

    return poly_T


def translate_flag(timeseries):
    """Translates Vitens Lizard flags from integer to text.

    Parameters
    ----------
    timeseries : pandas.DataFrame
        timeseries of a monitoring well with flags

    Returns
    -------
    timeseries : pandas.DataFrame
        timeseries with translated quality flags
    """
    translate_dic = {
        0: "betrouwbaar",
        1: "betrouwbaar",
        3: "onbeslist",
        4: "onbeslist",
        6: "onbetrouwbaar",
        7: "onbetrouwbaar",
        99: "ongevalideerd",
        -99: "verwijderd",
    }
    timeseries["flag"] = timeseries["flag"].replace(translate_dic)

    return timeseries


def get_metadata_mw_from_code(code, organisation="vitens", auth=None):
    """Extracts the Groundwater Station parameters from a monitoring well based on the
    code of the monitoring well.

    Parameters
    ----------
    code : str
        code of the monitoring well
    organisation : str
        organisation indicating URL endpoint, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Raises
    ------
    ValueError
        if code of the monitoring well is not known

    Returns
    -------
    groundwaterstation_metadata : dict
        dictionary with all available metadata of the monitoring well and its filters
    """
    base_url = lizard_api_endpoint.format(organisation=organisation)
    lizard_GWS_endpoint = f"{base_url}groundwaterstations/"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}?code={code}"

    ValueError(url_groundwaterstation_code)
    r = requests.get(url_groundwaterstation_code, auth=auth, timeout=1200)
    r.raise_for_status()

    try:
        groundwaterstation_metadata = r.json()["results"][0]
    except IndexError:
        raise ValueError(r.json())
        raise ValueError(f"Code {code} is invalid")

    return groundwaterstation_metadata


def _prepare_API_input(nr_pages, url_groundwater):
    """Get API data pages within the defined extent.

    Parameters
    ----------
    nr_pages : int
        number of the pages on which the information is stored
    url_groundwater : str
        location of the used API to extract the data

    Returns
    -------
    urls : list
        list of the page number and the corresponding url
    """
    urls = []
    for page in range(nr_pages):
        true_page = page + 1  # The real page number is attached to the import thread
        urls += [url_groundwater + "&page={}".format(true_page)]
    return urls


def _download(url, timeout=1800, auth=None):
    """Function to download the data from the API using the ThreadPoolExecutor.

    Parameters
    ----------
    url : str
        url of an API page
    timeout : int, optional
        number of seconds to wait before terminating request
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    dictionary with timeseries data
    """
    try:
        r = requests.get(url=url, timeout=timeout, auth=auth)
        r.raise_for_status()
        data = r.json()["results"]
    except requests.exceptions.Timeout:
        logger.error(f"Timeout while requesting {url}. Please check your connection.")
        data = []

    return data


def _split_mw_tube_nr(code):
    """get the tube number from a code that consists of the name and the tube number.

    Parameters
    ----------
    code : str
        name + tube_nr. e.g. 'BUWP014-11' or 'BUWP014012'

    Returns
    -------
    monitoring well, tube_number (str, int)

    Notes
    -----
    The format of the name + tube_nr is not very consistent and this function may need
    further finetuning.
    """

    if code[-3:].isdigit():
        return code[:-3], int(code[-3:])
    else:
        # assume there is a '-' to split name and filter number
        tube_nr = code.split("-")[-1]
        return code.strip(f"-{tube_nr}"), int(tube_nr)


def _extract_timeseries_info_from_tube(mtd_tube, auth=None):
    """
    Extracts timeseries information (hand/diver UUIDs and types) from a tube/filter dict.

    Parameters
    ----------
    mtd_tube : dict
        Tube/filter metadata dictionary.
    auth : tuple, optional
        Authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    dict
        Dictionary with timeseries info and type.
    """
    info = {}
    if not mtd_tube["timeseries"]:
        info["timeseries_type"] = None
        return info
    for series in mtd_tube["timeseries"]:
        r = requests.get(series, auth=auth)
        r.raise_for_status()
        series_info = r.json()
        if series_info["code"] == "WNS9040.hand":
            info["uuid_hand"] = series_info["uuid"]
            info["start_hand"] = series_info["start"]
        elif series_info["code"] == "WNS9040":
            info["uuid_diver"] = series_info["uuid"]
            info["start_diver"] = series_info["start"]
        elif series_info["code"] == "WNS9040.val":
            info["uuid_diver_validated"] = series_info["uuid"]
            info["start_diver_validated"] = series_info["start"]
            info["end_diver_validated"] = series_info["end"]

    # Create string with all timeseries types
    ts_types = []

    if info.get("start_hand") is not None:
        ts_types.append("hand")
    if info.get("start_diver") is not None:
        ts_types.append("diver")
    if info.get("start_diver_validated") is not None:
        ts_types.append("diver_validated")
    info["timeseries_type"] = " + ".join(ts_types) if ts_types else None

    return info


def get_metadata_tube(metadata_mw, tube_nr, auth=None):
    """Extract the metadata for a specific tube from the monitoring well metadata.

    Parameters
    ----------
    metadata_mw : dict
        dictionary with all available metadata of the monitoring well and all its
        filters
    tube_nr : int or None
        select metadata from a specific tube number
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Raises
    ------
    ValueError
        if code of the monitoring well is invalid.

    Returns
    -------
    dictionary with metadata of a specific tube

    .. warning::
       This function assumes that there is only one 'hand' timeseries and one 'diver'
       timeseries for each tube. This seems to comply with the Vitens use of Lizard.
    """

    # Set default tube number if not provided
    if tube_nr is None:
        tube_nr = 1

    # Prepare a base metadata dict from the monitoring well
    metadata = {
        "location": metadata_mw["name"],
        "ground_level": metadata_mw["surface_level"],
        "source": "lizard",
        "organisation": metadata_mw["organisation"],
        "unit": "m NAP",
        "metadata_available": True,
        "status": None,
    }

    # Searches for filters matching the requested tube number
    metadata_tube_list = []
    for metadata_tube in metadata_mw["filters"]:
        # check if name+filternr ends with three digits
        code, tbnr = _split_mw_tube_nr(metadata_tube["code"])
        if tbnr == tube_nr:
            metadata_tube_list.append(metadata_tube)

    # Handles cases with no, one or multiple tubes with the same code and tube number
    if len(metadata_tube_list) == 0:
        raise ValueError(f"{metadata_mw['name']} doesn't have a tube number {tube_nr}")
    elif len(metadata_tube_list) == 1:
        mtd_tube = metadata_tube_list[0]
    elif len(metadata_tube_list) > 1:
        # tube has probably been replaced, multiple tubes with the same code and tube nr
        # merge metadata from all tubes
        logger.info(
            f"there are {len(metadata_tube_list)} instances of {code} and tube "
            f"{tube_nr}, trying to merge all in one observation object"
        )
        mtd_tube = metadata_tube_list[0].copy()
        relevant_keys = {
            "top_level",
            "filter_top_level",
            "filter_bottom_level",
            "timeseries",
        }
        for metadata_tube in metadata_tube_list:
            for key in set(metadata_tube.keys()) & relevant_keys:
                # check if properties are always the same for a tube number
                val = metadata_tube[key]
                if key in ["top_level", "filter_top_level", "filter_bottom_level"]:
                    if val != mtd_tube[key]:
                        logger.warning(
                            f"multiple {key} values found ({val} & {mtd_tube[key]})"
                            f" for {code} and tube {tube_nr}, using {mtd_tube[key]}"
                        )
                # merge time series from all tubes with the same code and tube number
                elif key == "timeseries":
                    mtd_tube[key] += val

        mtd_tube["code"] = f"{code}{tube_nr}"

    # Updates metadata with tube-specific information (top level, screen levels, coords)
    metadata.update(
        {
            "tube_nr": tube_nr,
            "name": mtd_tube["code"].replace("-", ""),
            "tube_top": mtd_tube["top_level"],
            "screen_top": mtd_tube["filter_top_level"],
            "screen_bottom": mtd_tube["filter_bottom_level"],
        }
    )

    lon, lat, _ = metadata_mw["geometry"]["coordinates"]
    transformer = Transformer.from_crs("WGS84", "EPSG:28992")
    metadata["x"], metadata["y"] = transformer.transform(lat, lon)

    # Extracts timeseries information (hand/diver UUIDs and types)
    metadata.update(_extract_timeseries_info_from_tube(mtd_tube, auth))

    return metadata


def get_timeseries_uuid(
    uuid, tmin, tmax, page_size=100000, organisation="vitens", auth=None
):
    """
    Get the time series (hand or diver) using the uuid.

    ----------
    uuid : str
        Universally Unique Identifier of the tube and type of time series.
    tmin : str YYYY-m-d
        start of the observations, by default the entire serie is returned
    tmax : str YYYY-m-d
        end of the observations, by default the entire serie is returned
    page_size : int, optional
        Query parameter which can extend the response size. The default is 100000.
    organisation : str, optional
        organisation as used by Lizard, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    pd.DataFrame
        pandas DataFrame with the timeseries of the monitoring well
    """
    base_url = lizard_api_endpoint.format(organisation=organisation)
    url_timeseries = f"{base_url}timeseries/{uuid}"

    if tmin is not None:
        tmin = pd.to_datetime(tmin).isoformat("T")

    if tmax is not None:
        tmax = pd.to_datetime(tmax).isoformat("T")

    params = {"start": tmin, "end": tmax, "page_size": page_size}
    url = url_timeseries + "/events/"

    r = requests.get(url=url, params=params, auth=auth)
    time_series_events = r.json()["results"]
    time_series_df = pd.DataFrame(time_series_events)

    if time_series_df.empty:
        return pd.DataFrame()

    else:
        time_series_df = translate_flag(time_series_df)

        timeseries_sel = time_series_df.loc[:, ["time", "value", "flag", "comment"]]
        timeseries_sel["time"] = pd.to_datetime(
            timeseries_sel["time"], format="%Y-%m-%dT%H:%M:%SZ", errors="coerce"
        ) + pd.DateOffset(hours=1)

        timeseries_sel = timeseries_sel[~timeseries_sel["time"].isnull()]

        timeseries_sel.set_index("time", inplace=True)
        timeseries_sel.index.rename("peil_datum_tijd", inplace=True)
        # timeseries_sel.dropna(inplace=True)

    return timeseries_sel


def _filter_timeseries(ts_dict, datafilters):
    """
    Generic filter function for multiple timeseries.

    Parameters
    ----------
    ts_dict : dict
        Dictionary of timeseries DataFrames, e.g. {'hand': df1, 'diver': df2, ...}
    datafilters : list of str
        List of datafilter names as strings, e.g. ["remove_unvalidated_diver_values_when_validated_available", "remove_hand_during_diver_period"]

    Returns
    -------
    dict
        Filtered ts_dict.
    """
    # Define standard filters by name. Note that the order may be relevant (uppermost filter is applied first).
    standard_datafilters = {
        "remove_unvalidated_diver_values_when_validated_available": {
            "target": "diver",
            "action": "remove_before",
            "reference": "diver_validated",
            "how": "max",
        },
        "remove_hand_during_diver_period": {
            "target": "hand",
            "action": "remove_between",
            "reference": ["diver", "diver_validated"],
            "how": "range",
        },
    }

    # Convert string datafilters to dicts using standard_datafilters
    datafilter_dicts = []
    for f in datafilters:
        if isinstance(f, str):
            if f in standard_datafilters:
                datafilter_dicts.append(standard_datafilters[f])
            else:
                raise ValueError(f"Unknown filter name: {f}")
        else:
            raise ValueError(
                "Each filter must be a string referring to a standard filter."
            )

    ts_dict = {k: v.copy() for k, v in ts_dict.items()}

    for f in datafilter_dicts:
        targets = f["target"] if isinstance(f["target"], list) else [f["target"]]
        refs = f["reference"] if isinstance(f["reference"], list) else [f["reference"]]
        how = f.get("how", "range")
        on = f.get("on", None)
        for target in targets:
            df = ts_dict.get(target)
            if df is None or df.empty:
                continue
            mask = pd.Series(True, index=df.index)
            for ref in refs:
                ref_df = ts_dict.get(ref)
                if ref_df is None or ref_df.empty:
                    continue
                idx = ref_df.index if on is None else ref_df[on]
                if f["action"] == "remove_before":
                    val = idx.max() if how == "max" else idx.min()
                    mask &= df.index > val
                elif f["action"] == "remove_after":
                    val = idx.min() if how == "min" else idx.max()
                    mask &= df.index < val
                elif f["action"] == "remove_between":
                    mask &= ~df.index.to_series().between(idx.min(), idx.max())
                elif f["action"] == "keep_only":
                    mask &= df.index.to_series().between(idx.min(), idx.max())
                # Add more actions as needed
            ts_dict[target] = df[mask]
    return ts_dict


def get_timeseries_tube(
    tube_metadata,
    tmin,
    tmax,
    type_timeseries=None,  # deprecated argument
    which_timeseries=["hand", "diver"],  # new preferred argument
    datafilters=None,
    combine_method="merge",
    organisation="vitens",
    auth=None,
):
    """
    Extracts specified timeseries for a tube and combines them as requested.

    Parameters
    ----------
    tube_metadata : dict
        metadata of a tube
    tmin : str YYYY-m-d, optional
        start of the observations
    tmax : str YYYY-m-d, optional
        end of the observations
    type_timeseries : str, optional (deprecated)
        Old keyword, use which_timeseries instead.
    which_timeseries : list of str, optional
        Which timeseries to retrieve. Options: "hand", "diver", "diver_validated".
        Defaults to ["hand", "diver"] (which should be correct for Vitens).
    datafilters : list of strings, optional
        Methods to filter the timeseries data.
        If None (default), all measurements will be shown.
        Currently implemented datafilter methods:
        "remove_unvalidated_diver_values_when_validated_available": Removes diver values before last date with validated diver.
        "remove_hand_during_diver_period": Removes hand measurements during periods where diver or diver_validated measurements are available.
    combine_method : str, optional
        "merge" (vertical stack with 'origin' column) or "combine" (side-by-side columns).
        If None, defaults to "merge".
    organisation : str, optional
        organisation as used by Lizard.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    measurements : pandas DataFrame
        timeseries of the monitoring well
    metadata_df : dict
        metadata of the monitoring well
    """
    # Deprecation warning for type_timeseries
    if type_timeseries is not None:
        logger.warning(
            "The 'type_timeseries' argument is deprecated. "
            "Please use 'which_timeseries' (a list, e.g. ['hand', 'diver']) and 'combine_method' instead."
        )
        # Map old type_timeseries to which_timeseries and combine_method
        if type_timeseries == "combine":
            combine_method = "combine"
        elif type_timeseries == "merge":
            combine_method = "merge"
        else:
            which_timeseries = [type_timeseries]
            combine_method = "merge"

    if tube_metadata["timeseries_type"] is None:
        return pd.DataFrame(), tube_metadata

    # Fetch requested timeseries
    ts_dict = {}
    for ts_type in which_timeseries:
        uuid_key = f"uuid_{ts_type}"
        start_key = f"start_{ts_type}"
        if tube_metadata.get(start_key) is not None:
            ts = get_timeseries_uuid(
                tube_metadata.get(uuid_key),
                tmin,
                tmax,
                organisation=organisation,
                auth=auth,
            )
        else:
            ts = pd.DataFrame()

        ts_dict[ts_type] = ts

    # Filter
    if datafilters is not None:
        ts_dict_filtered = _filter_timeseries(ts_dict, datafilters)
    else:
        ts_dict_filtered = ts_dict

    # Combine as requested
    if combine_method == "combine":
        # Side-by-side
        if not ts_dict_filtered.get("hand", pd.DataFrame()).empty:
            ts_dict_filtered["hand"] = ts_dict_filtered["hand"].rename(
                columns={"value": "value_hand", "flag": "flag_hand"}
            )
        if not ts_dict_filtered.get("diver", pd.DataFrame()).empty:
            ts_dict_filtered["diver"] = ts_dict_filtered["diver"].rename(
                columns={"value": "value_diver", "flag": "flag_diver"}
            )
        if not ts_dict_filtered.get("diver_validated", pd.DataFrame()).empty:
            ts_dict_filtered["diver_validated"] = ts_dict_filtered[
                "diver_validated"
            ].rename(
                columns={
                    "value": "value_diver_validated",
                    "flag": "flag_diver_validated",
                }
            )

        dfs = [
            df
            for key in ["hand", "diver_validated", "diver"]
            if key in ts_dict_filtered and not ts_dict_filtered[key].empty
            for df in [ts_dict_filtered[key]]
        ]
        measurements = pd.concat(dfs, axis=1) if dfs else pd.DataFrame()
        # Only keep present expected columns
        expected_cols = [
            "value_hand",
            "value_diver_validated",
            "value_diver",
            "flag_hand",
            "flag_diver_validated",
            "flag_diver",
        ]
        present_cols = [col for col in expected_cols if col in measurements.columns]
        if not measurements.empty:
            measurements = measurements.loc[:, present_cols]
    else:
        # Default: merge (vertical stack)
        dfs = []
        for key in which_timeseries:
            df = ts_dict_filtered.get(key)
            if df is not None and not df.empty:
                df = df.copy()
                df["origin"] = key
                dfs.append(df)
        measurements = pd.concat(dfs, axis=0).sort_index() if dfs else pd.DataFrame()

    return measurements, tube_metadata


def get_lizard_groundwater(
    code,
    tube_nr=None,
    tmin=None,
    tmax=None,
    type_timeseries=None,  # deprecated argument
    which_timeseries=["hand", "diver"],  # new preferred argument
    datafilters=None,
    combine_method="merge",
    only_metadata=False,
    organisation="vitens",
    auth=None,
):
    """Extracts the metadata and timeseries of an observation well from a LIZARD-API
    based on the code of a monitoring well.

    Parameters
    ----------
    code : str
        code of the measuring well, e.g. '27B-0444'
    tube_nr : int, optional
        select specific tube top
        Default selects tube_nr = 1
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : str YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional (deprecated)
        Old keyword, use which_timeseries instead.
    which_timeseries : list of str, optional
        Which timeseries to retrieve. Options: "hand", "diver", "diver_validated".
        Defaults to ["hand", "diver"] (which should be correct for Vitens).
    datafilters : list of strings, optional
        Methods to filter the timeseries data.
        If None (default), all measurements will be shown.
        Currently implemented datafilter methods:
        "remove_unvalidated_diver_values_when_validated_available": Removes diver values before last date with validated diver.
        "remove_hand_during_diver_period": Removes hand measurements during periods where diver or diver_validated measurements are available.
    combine_method : str, optional
        "merge" (vertical stack with 'origin' column) or "combine" (side-by-side columns).
        If None, defaults to "merge".
    only_metadata : bool, optional
        if True only metadata is returned and no time series data. The
        default is False.
    organisation : str, optional
        organisation as used by Lizard, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    measurements : pd.DataFrame
        returns a DataFrame with metadata and timeseries
    tube_metadata : dict
        dictionary containing metadata
    """

    # Deprecation warning for type_timeseries
    if type_timeseries is not None:
        logger.warning(
            "The 'type_timeseries' argument is deprecated. "
            "Please use 'which_timeseries' (a list, e.g. ['hand', 'diver']) and 'combine_method' instead."
        )

        # Map old type_timeseries to which_timeseries and combine_method
        if type_timeseries == "combine":
            combine_method = "combine"
        elif type_timeseries == "merge":
            combine_method = "merge"
        else:
            which_timeseries = [type_timeseries]
            combine_method = "merge"

    groundwaterstation_metadata = get_metadata_mw_from_code(
        code, organisation=organisation, auth=auth
    )

    tube_metadata = get_metadata_tube(groundwaterstation_metadata, tube_nr, auth=auth)

    if only_metadata:
        return pd.DataFrame(), tube_metadata

    measurements, tube_metadata = get_timeseries_tube(
        tube_metadata,
        tmin,
        tmax,
        which_timeseries=which_timeseries,
        datafilters=datafilters,
        combine_method=combine_method,
        organisation=organisation,
        auth=auth,
    )
    tube_metadata = check_status_obs(tube_metadata, measurements)

    return measurements, tube_metadata


def get_obs_list_from_codes(
    codes,
    ObsClass,
    tube_nr="all",
    tmin=None,
    tmax=None,
    type_timeseries=None,  # deprecated argument
    which_timeseries=["hand", "diver"],  # new preferred argument
    datafilters=None,
    combine_method="merge",
    only_metadata=False,
    organisation="vitens",
    auth=None,
):
    """Get all observations from a list of codes of the monitoring wells and a list of
    tube numbers.

    Parameters
    ----------
    codes : list of str or str
        codes of the monitoring wells
    ObsClass : type
        class of the observations, e.g. GroundwaterObs
    tube_nr : lst of str
        list of tube numbers of the monitoring wells that should be selected.
        By default 'all' available tubes are selected.
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : str YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional (deprecated)
        Old keyword, use which_timeseries instead.
    which_timeseries : list of str, optional
        Which timeseries to retrieve. Options: "hand", "diver", "diver_validated".
        Defaults to ["hand", "diver"] (which should be correct for Vitens).
    datafilters : list of strings, optional
        Methods to filter the timeseries data.
        If None (default), all measurements will be shown.
        Currently implemented datafilter methods:
        "remove_unvalidated_diver_values_when_validated_available": Removes diver values before last date with validated diver.
        "remove_hand_during_diver_period": Removes hand measurements during periods where diver or diver_validated measurements are available.
    combine_method : str, optional
        "merge" (vertical stack with 'origin' column) or "combine" (side-by-side columns).
        If None, defaults to "merge".
    only_metadata : bool, optional
        if True only metadata is returned and no time series data. The
        default is False.
    organisation : str, optional
        organisation as used by Lizard, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    obs_list
        list of observations
    """

    # Deprecation warning for type_timeseries
    if type_timeseries is not None:
        logger.warning(
            "The 'type_timeseries' argument is deprecated. "
            "Please use 'which_timeseries' (a list, e.g. ['hand', 'diver']) and 'combine_method' instead."
        )
        # Map old type_timeseries to which_timeseries and combine_method
        if type_timeseries == "combine":
            combine_method = "combine"
        elif type_timeseries == "merge":
            combine_method = "merge"
        else:
            which_timeseries = [type_timeseries]
            combine_method = "merge"

    if isinstance(codes, str):
        codes = [codes]

    if not hasattr(codes, "__iter__"):
        raise TypeError("argument 'codes' should be an iterable")

    obs_list = []
    for code in codes:
        groundwaterstation_metadata = get_metadata_mw_from_code(
            code, organisation=organisation, auth=auth
        )
        tubes = []
        if tube_nr == "all":
            for metadata_tube in groundwaterstation_metadata["filters"]:
                tnr = _split_mw_tube_nr(metadata_tube["code"])[-1]
                if tnr not in tubes:
                    logger.debug(f"get {code}-{tnr}")
                    o = ObsClass.from_lizard(
                        code,
                        tnr,
                        tmin,
                        tmax,
                        # type_timeseries,
                        which_timeseries=which_timeseries,
                        datafilters=datafilters,
                        combine_method=combine_method,
                        only_metadata=only_metadata,
                        organisation=organisation,
                        auth=auth,
                    )
                    obs_list.append(o)
                    tubes.append(tnr)

        else:
            o = ObsClass.from_lizard(
                code,
                tube_nr,
                tmin,
                tmax,
                which_timeseries=which_timeseries,
                datafilters=datafilters,
                combine_method=combine_method,
                only_metadata=only_metadata,
                organisation=organisation,
                auth=auth,
            )
            obs_list.append(o)

    return obs_list


def get_obs_list_from_extent(
    extent,
    ObsClass,
    tube_nr="all",
    tmin=None,
    tmax=None,
    type_timeseries=None,  # deprecated argument
    which_timeseries=["hand", "diver"],  # new preferred argument
    datafilters=None,
    combine_method="merge",
    only_metadata=False,
    page_size=100,
    nr_threads=10,
    organisation="vitens",
    auth=None,
):
    """Get all observations within a specified extent.

    Parameters
    ----------
    extent : list or shapefile
        get groundwater monitoring wells within this extent [xmin, xmax, ymin, ymax]
        or within a predefined Polygon from a shapefile
    ObsClass : type
        class of the observations, e.g. GroundwaterObs
    tube_nr : lst of str
        list of tube numbers of the monitoring wells that should be selected.
        By default 'all' available tubes are selected.
    tmin : str, optional
        start of the observations (format YYYY-m-d), by default the entire series
        is returned
    tmax : str, optional
        end of the observations (format YYYY-m-d), by default the entire series
        is returned
    type_timeseries : str, optional (deprecated)
        Old keyword, use which_timeseries instead.
    which_timeseries : list of str, optional
        Which timeseries to retrieve. Options: "hand", "diver", "diver_validated".
        Defaults to ["hand", "diver"] (which should be correct for Vitens).
    datafilters : list of strings, optional
        Methods to filter the timeseries data.
        If None (default), all measurements will be shown.
        Currently implemented datafilter methods:
        "remove_unvalidated_diver_values_when_validated_available": Removes diver values before last date with validated diver.
        "remove_hand_during_diver_period": Removes hand measurements during periods where diver or diver_validated measurements are available.
    combine_method : str, optional
        "merge" (vertical stack with 'origin' column) or "combine" (side-by-side columns).
        If None, defaults to "merge".
    only_metadata : bool, optional
        if True only metadata is returned and no time series data. The
        default is False.
    organisation : str
        organisation as used by Lizard, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)
    page_size : int, optional
        number of records to retrieve per page, default is 100
    nr_threads : int, optional
        number of threads to use for the API requests, default is 10


    Returns
    -------
    obs_col : TYPE
        ObsCollection DataFrame with the 'obs' column
    """

    # Deprecation warning for type_timeseries
    if type_timeseries is not None:
        logger.warning(
            "The 'type_timeseries' argument is deprecated. "
            "Please use 'which_timeseries' (a list, e.g. ['hand', 'diver']) and 'combine_method' instead."
        )
        # Map old type_timeseries to which_timeseries and combine_method
        if type_timeseries == "combine":
            combine_method = "combine"
        elif type_timeseries == "merge":
            combine_method = "merge"
        else:
            which_timeseries = [type_timeseries]
            combine_method = "merge"

    if isinstance(extent, (list, tuple)):
        polygon_T = extent_to_wgs84_polygon(extent)

    elif isinstance(extent, str) or isinstance(extent, pathlib.PurePath):
        polygon = geopandas.read_file(extent)
        # TODO: check this transformation
        polygon_T = polygon.to_crs("WGS84", "EPSG:28992").loc[0, "geometry"]
    else:
        raise TypeError("Extent should be a shapefile or a list of coordinates")

    base_url = lizard_api_endpoint.format(organisation=organisation)
    lizard_GWS_endpoint = f"{base_url}groundwaterstations/"
    url_groundwaterstation_extent = (
        f"{lizard_GWS_endpoint}?geometry__within={polygon_T}&page_size={page_size}"
    )
    r = requests.get(url_groundwaterstation_extent, auth=auth)
    r.raise_for_status()
    groundwaterstation_data = r.json()
    nr_results = groundwaterstation_data["count"]
    nr_pages = math.ceil(nr_results / page_size)

    logger.info("Number of monitoring wells: {}".format(nr_results))
    logger.info("Number of pages: {}".format(nr_pages))

    if nr_results == 0:
        ValueError(r.json())
        logger.warning(
            "No monitoring wells found in the specified extent. "
            "Please check the extent or the organisation."
        )
        return []

    if nr_threads > nr_pages:
        nr_threads = nr_pages

    urls = _prepare_API_input(nr_pages, url_groundwaterstation_extent)

    # Prepare arguments for get_obs_list_from_codes
    kwargs = {
        "ObsClass": ObsClass,
        "tube_nr": tube_nr,
        "tmin": tmin,
        "tmax": tmax,
        "which_timeseries": which_timeseries,
        "datafilters": datafilters,
        "combine_method": combine_method,
        "only_metadata": only_metadata,
        "organisation": organisation,
        "auth": auth,
    }

    codes = []
    with ThreadPoolExecutor(max_workers=nr_threads) as executor:
        for result in tqdm(
            executor.map(lambda url: _download(url, auth=auth), urls),
            total=nr_pages,
            desc="Page",
        ):
            codes += [d["code"] for d in result]

    obs_list = []
    with ThreadPoolExecutor() as executor:
        for obs_list_mw in tqdm(
            executor.map(lambda code: get_obs_list_from_codes(code, **kwargs), codes),
            total=len(codes),
            desc="monitoring well",
        ):
            obs_list += obs_list_mw

    return obs_list
