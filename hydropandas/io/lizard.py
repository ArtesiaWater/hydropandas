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

    try:
        groundwaterstation_metadata = requests.get(
            url_groundwaterstation_code, auth=auth
        ).json()["results"][0]

    except IndexError:
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
        data = requests.get(url=url, timeout=timeout, auth=auth)
        data = data.json()["results"]
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
        series_info = requests.get(series, auth=auth).json()
        if series_info["code"] == "WNS9040.hand":
            info["uuid_hand"] = series_info["uuid"]
            info["start_hand"] = series_info["start"]
        elif series_info["code"] == "WNS9040":
            info["uuid_diver"] = series_info["uuid"]
            info["start_diver"] = series_info["start"]
    if (info.get("start_hand") is None) and (info.get("start_diver") is None):
        info["timeseries_type"] = None
    elif (info.get("start_hand") is not None) and (info.get("start_diver") is not None):
        info["timeseries_type"] = "diver + hand"
    elif info.get("start_hand") is None:
        info["timeseries_type"] = "diver"
    elif info.get("start_diver") is None:
        info["timeseries_type"] = "hand"
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
    tmax : int YYYY-m-d
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

    time_series_events = requests.get(url=url, params=params, auth=auth).json()[
        "results"
    ]
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


def _merge_timeseries(hand_measurements, diver_measurements):
    """Merges the timeseries of the hand and diver measurements into one timeserie.

    Parameters
    ----------
    hand_measurements : DataFrame
        DataFrame containing the hand measurements of the monitoring well
    diver_measurements : DataFrame
        DataFrame containing the Diver measurements of the monitoring well

    Returns
    -------
    DataFrame where hand and diver measurements are merged in one timeseries
    """
    if hand_measurements.empty and diver_measurements.empty:
        measurements = pd.DataFrame()

    elif diver_measurements.first_valid_index() is None:
        measurements = hand_measurements
        logger.debug("no diver measurements available")

    else:
        hand_measurements_sel = hand_measurements.loc[
            hand_measurements.index < diver_measurements.first_valid_index()
        ]
        measurements = pd.concat([hand_measurements_sel, diver_measurements], axis=0)

    return measurements


def _combine_timeseries(hand_measurements, diver_measurements):
    """Combines the timeseries of the hand and diver measurements into one DataFrame.

    Parameters
    ----------
    hand_measurements : DataFrame
        DataFrame containing the hand measurements of the monitoring well
    diver_measurements : DataFrame
        DataFrame containing the Diver measurements of the monitoring well

    Returns
    -------
    a combined DataFrame with both hand, and diver measurements
        DESCRIPTION.
    """
    hand_measurements.rename(
        columns={"value": "value_hand", "flag": "flag_hand"}, inplace=True
    )
    diver_measurements.rename(
        columns={"value": "value_diver", "flag": "flag_diver"}, inplace=True
    )

    measurements = pd.concat([hand_measurements, diver_measurements], axis=1)
    measurements = measurements.loc[
        :, ["value_hand", "value_diver", "flag_hand", "flag_diver"]
    ]

    return measurements


def get_timeseries_tube(
    tube_metadata, tmin, tmax, type_timeseries, organisation="vitens", auth=None
):
    """Extracts multiple timeseries (hand and/or diver measurements) for a specific tube
    using the Lizard API.

    Parameters
    ----------
    tube_metadata : dict
        metadata of a tube
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : str YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements
        merge: the hand and diver measurements into one time series (default)
        combine: keeps hand and diver measurements separated
    organisation : str, optional
        organisation as used by Lizard, currently only "vitens" is officially supported.
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    measurements : pandas DataFrame
        timeseries of the monitoring well
    metadata_df : dict
        metadata of the monitoring well
    """

    if tube_metadata["timeseries_type"] is None:
        return pd.DataFrame(), tube_metadata

    if type_timeseries in ["hand", "merge", "combine"]:
        if "hand" in tube_metadata["timeseries_type"]:
            hand_measurements = get_timeseries_uuid(
                tube_metadata.pop("uuid_hand"),
                tmin,
                tmax,
                organisation=organisation,
                auth=auth,
            )
        else:
            hand_measurements = None

    if type_timeseries in ["diver", "merge", "combine"]:
        if "diver" in tube_metadata["timeseries_type"]:
            diver_measurements = get_timeseries_uuid(
                tube_metadata.pop("uuid_diver"),
                tmin,
                tmax,
                organisation=organisation,
                auth=auth,
            )
        else:
            diver_measurements = None

    if type_timeseries == "hand" and hand_measurements is not None:
        measurements = hand_measurements
    elif type_timeseries == "diver" and diver_measurements is not None:
        measurements = diver_measurements
    elif type_timeseries in ["merge", "combine"]:
        if (hand_measurements is not None) and (diver_measurements is not None):
            if type_timeseries == "merge":
                measurements = _merge_timeseries(hand_measurements, diver_measurements)
            elif type_timeseries == "combine":
                measurements = _combine_timeseries(
                    hand_measurements, diver_measurements
                )
        elif hand_measurements is not None:
            measurements = hand_measurements
        elif diver_measurements is not None:
            measurements = diver_measurements

    return measurements, tube_metadata


def get_lizard_groundwater(
    code,
    tube_nr=None,
    tmin=None,
    tmax=None,
    type_timeseries="merge",
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
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements
        merge: the hand and diver measurements into one time series (default)
        combine: keeps hand and diver measurements separated
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

    groundwaterstation_metadata = get_metadata_mw_from_code(
        code, organisation=organisation, auth=auth
    )

    tube_metadata = get_metadata_tube(groundwaterstation_metadata, tube_nr, auth=auth)

    if only_metadata:
        return pd.DataFrame(), tube_metadata

    measurements, tube_metadata = get_timeseries_tube(
        tube_metadata, tmin, tmax, type_timeseries, organisation=organisation, auth=auth
    )
    tube_metadata = check_status_obs(tube_metadata, measurements)

    return measurements, tube_metadata


def get_obs_list_from_codes(
    codes,
    ObsClass,
    tube_nr="all",
    tmin=None,
    tmax=None,
    type_timeseries="merge",
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
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements
        merge: the hand and diver measurements into one time series (default)
        combine: keeps hand and diver measurements separeted
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
                        type_timeseries,
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
                type_timeseries,
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
    type_timeseries="merge",
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
    type_timeseries : str, optional
        merge: the hand and diver measurements into one time series (merge; default) or
        combine: keeps hand and diver measurements separeted
        The default is merge.
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

    groundwaterstation_data = requests.get(
        url_groundwaterstation_extent, auth=auth
    ).json()
    nr_results = groundwaterstation_data["count"]
    nr_pages = math.ceil(nr_results / page_size)

    logger.info("Number of monitoring wells: {}".format(nr_results))
    logger.info("Number of pages: {}".format(nr_pages))

    if nr_results == 0:
        logger.warning(
            "No monitoring wells found in the specified extent. "
            "Please check the extent or the organisation."
        )
        return []

    if nr_threads > nr_pages:
        nr_threads = nr_pages

    urls = _prepare_API_input(nr_pages, url_groundwaterstation_extent)

    arg_tuple = (ObsClass, tube_nr, tmin, tmax, type_timeseries, only_metadata)
    codes = []
    with ThreadPoolExecutor(max_workers=nr_threads) as executor:
        for result in tqdm(
            executor.map(lambda url: _download(url, auth=auth), urls),
            total=nr_pages,
            desc="Page",
        ):
            codes += [(d["code"],) + arg_tuple for d in result]

    obs_list = []
    with ThreadPoolExecutor() as executor:
        for obs_list_mw in tqdm(
            executor.map(lambda args: get_obs_list_from_codes(*args, auth=auth), codes),
            total=len(codes),
            desc="monitoring well",
        ):
            obs_list += obs_list_mw

    return obs_list


def get_monitoring_networks(organisation="vitens", auth=None, page_size=1000):
    """
    Get all monitoring networks of the specified organisation.

    ----------
    organisation : str, optional
        organisation indicating URL endpoint, currently only "vitens" is officially supported
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)
    page_size : int, optional
        number of records to retrieve per page, default is 1000

    Returns:
    ----------
    pd.DataFrame
        pandas DataFrame containing the available monitoring networks
    """
    base_url = lizard_api_endpoint.format(organisation=organisation)
    url_monitoringnetworks = f"{base_url}monitoringnetworks/"
    params = {"page_size": page_size}

    monitoring_networks = requests.get(
        url_monitoringnetworks, params=params, auth=auth
    ).json()["results"]
    monitoring_networks_df = pd.DataFrame(monitoring_networks)

    return monitoring_networks_df


def get_locations_in_monitoring_network_uuid(
    monitoring_network_uuid, organisation="vitens", auth=None, page_size=1000
):
    """
    Get all locations within a given monitoring network (specified by UUID).

    ----------
    monitoring_network_uuid : str
        UUID of the monitoring network to query
    organisation : str, optional
        name of the organisation, default is "vitens"
    auth : tuple, optional
        authentication credentials for the API request, e.g.: ("__key__", your_api_key)
    page_size : int, optional
        number of records to retrieve per page, default is 1000

    Returns:
    ----------
    pd.DataFrame
        pandas DataFrame containing the available monitoring networks
    """

    base_url = lizard_api_endpoint.format(organisation=organisation)
    url_monitoring_network_locations = (
        f"{base_url}monitoringnetworks/{monitoring_network_uuid}/locations/"
    )
    params = {"page_size": page_size}

    locations_in_monitoring_network = requests.get(
        url=url_monitoring_network_locations, params=params, auth=auth
    ).json()["results"]
    locations_in_monitoring_network_df = pd.DataFrame(locations_in_monitoring_network)

    nr_results = len(locations_in_monitoring_network_df)
    logger.info(
        "Number of locations in current monitoring network: {}".format(nr_results)
    )

    return locations_in_monitoring_network_df


def get_locations_in_monitoring_networks(
    monitoring_networks, organisation="vitens", auth=None
):
    """
    Retrieve all locations in all given monitoring networks and return a concatenated DataFrame.

    Parameters
    ----------
    monitoring_networks : pd.DataFrame
        DataFrame containing monitoring network UUIDs in the 'uuid' column.
    organisation : str
        Organisation name for API requests, default is "vitens".
    auth : tuple
        Authentication credentials for the API request, e.g.: ("__key__", your_api_key)

    Returns
    -------
    pd.DataFrame
        DataFrame containing all locations from all monitoring networks.
    """
    locations_list = []
    for monitoring_network_uuid in monitoring_networks["uuid"]:
        locations_in_monitoring_network = get_locations_in_monitoring_network_uuid(
            monitoring_network_uuid=monitoring_network_uuid,
            organisation=organisation,
            auth=auth,
        )
        if not locations_in_monitoring_network.empty:
            locations_list.append(locations_in_monitoring_network)
        else:
            logger.info(
                f"No locations found in monitoring network {monitoring_network_uuid}"
            )
    if locations_list:
        locations_df = pd.concat(locations_list, ignore_index=True)
    else:
        logger.warning("No locations found in any monitoring network.")
        locations_df = pd.DataFrame()
    return locations_df
