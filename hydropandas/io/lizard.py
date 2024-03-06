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

# NOTE: currently only the vitens API is officially supported. If/when new endpoints
# are added we should check whether we want to add the URL as argument or add supported
# sources to this dictionary:
LIZARD_APIS = {"vitens": "https://vitens.lizard.net/api/v4/"}


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
        99: "onongevalideerd",
        -99: "verwijderd",
    }
    timeseries["flag"] = timeseries["flag"].replace(translate_dic)

    return timeseries


def get_metadata_mw_from_code(code, source="vitens"):
    """Extracts the Groundwater Station parameters from a monitoring well based on the
    code of the monitoring well.

    Parameters
    ----------
    code : str
        code of the monitoring well
    source : str
        source indicating URL endpoint, currently only "vitens" is officially supported.

    Raises
    ------
    ValueError
        if code of the monitoring well is not known

    Returns
    -------
    groundwaterstation_metadata : dict
        dictionary with all available metadata of the monitoring well and its filters
    """
    lizard_GWS_endpoint = f"{LIZARD_APIS[source]}groundwaterstations/"
    url_groundwaterstation_code = f"{lizard_GWS_endpoint}?code={code}"

    try:
        groundwaterstation_metadata = requests.get(url_groundwaterstation_code).json()[
            "results"
        ][0]

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


def _download(url, timeout=1800):
    """Function to download the data from the API using the ThreadPoolExecutor.

    Parameters
    ----------
    url : str
        url of an API page
    timeout : int, optional
        number of seconds to wait before terminating request

    Returns
    -------
    dictionary with timeseries data
    """
    data = requests.get(url=url, timeout=timeout)
    data = data.json()["results"]

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


def get_metadata_tube(metadata_mw, tube_nr):
    """Extract the metadata for a specific tube from the monitoring well metadata.

    Parameters
    ----------
    metadata_mw : dict
        dictionary with all available metadata of the monitoring well and all its
        filters
    tube_nr : int or None
        select metadata from a specific tube number

    Raises
    ------
    ValueError
        if code of the monitoring well is invalid.

    Returns
    -------
    dictionary with metadata of a specific tube
    """

    if tube_nr is None:
        tube_nr = 1

    metadata = {
        "monitoring_well": metadata_mw["name"],
        "ground_level": metadata_mw["surface_level"],
        "source": "lizard",
        "unit": "m NAP",
        "metadata_available": True,
        "status": None,
    }

    metadata_tube_list = []
    for metadata_tube in metadata_mw["filters"]:
        # check if name+filternr ends with three digits
        code, tbnr = _split_mw_tube_nr(metadata_tube["code"])
        if tbnr == tube_nr:
            metadata_tube_list.append(metadata_tube)

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

    if not mtd_tube["timeseries"]:
        metadata["timeseries_type"] = None
    else:
        for series in mtd_tube["timeseries"]:
            series_info = requests.get(series).json()
            if series_info["name"] == "WNS9040.hand":
                metadata["uuid_hand"] = series_info["uuid"]
                metadata["start_hand"] = series_info["start"]
            elif series_info["name"] == "WNS9040":
                metadata["uuid_diver"] = series_info["uuid"]
                metadata["start_diver"] = series_info["start"]

        if (metadata.get("start_hand") is None) and (
            metadata.get("start_diver") is None
        ):
            metadata["timeseries_type"] = None
        elif (metadata.get("start_hand") is not None) and (
            metadata.get("start_diver") is not None
        ):
            metadata["timeseries_type"] = "diver + hand"
        elif metadata.get("start_hand") is None:
            metadata["timeseries_type"] = "diver"
        elif metadata.get("start_diver") is None:
            metadata["timeseries_type"] = "hand"

    return metadata


def get_timeseries_uuid(uuid, tmin, tmax, page_size=100000, source="vitens"):
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
    source : str, optional
        source indicating URL endpoint, currently only "vitens" is officially supported

    Returns
    -------
    pd.DataFrame
        pandas DataFrame with the timeseries of the monitoring well
    """

    url_timeseries = LIZARD_APIS[source] + "timeseries/{}".format(uuid)

    if tmin is not None:
        tmin = pd.to_datetime(tmin).isoformat("T")

    if tmax is not None:
        tmax = pd.to_datetime(tmax).isoformat("T")

    params = {"start": tmin, "end": tmax, "page_size": page_size}
    url = url_timeseries + "/events/"

    time_series_events = requests.get(url=url, params=params).json()["results"]
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
        print(
            "no diver measuremets available for {}".format(
                hand_measurements.iloc[0]["name"]
            )
        )

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


def get_timeseries_tube(tube_metadata, tmin, tmax, type_timeseries):
    """Extracts multiple timeseries (hand and/or diver measurements) for a specific tube
    using the Lizard API.

    Parameters
    ----------
    tube_metadata : dict
        metadata of a tube
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements
        merge: the hand and diver measurements into one time series (default)
        combine: keeps hand and diver measurements separeted

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
            )
        else:
            hand_measurements = None

    if type_timeseries in ["diver", "merge", "combine"]:
        if "diver" in tube_metadata["timeseries_type"]:
            diver_measurements = get_timeseries_uuid(
                tube_metadata.pop("uuid_diver"),
                tmin,
                tmax,
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

    Returns
    -------
    measurements : pd.DataFrame
        returns a DataFrame with metadata and timeseries
    tube_metadata : dict
        dictionary containing metadata
    """

    groundwaterstation_metadata = get_metadata_mw_from_code(code)

    tube_metadata = get_metadata_tube(groundwaterstation_metadata, tube_nr)

    if only_metadata:
        return pd.DataFrame(), tube_metadata

    measurements, tube_metadata = get_timeseries_tube(
        tube_metadata, tmin, tmax, type_timeseries
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
        groundwaterstation_metadata = get_metadata_mw_from_code(code)
        tubes = []
        if tube_nr == "all":
            for metadata_tube in groundwaterstation_metadata["filters"]:
                tnr = _split_mw_tube_nr(metadata_tube["code"])[-1]
                if tnr not in tubes:
                    logger.info(f"get {code}{tnr}")
                    o = ObsClass.from_lizard(
                        code,
                        tnr,
                        tmin,
                        tmax,
                        type_timeseries,
                        only_metadata=only_metadata,
                    )
                    obs_list.append(o)
                    tubes.append(tnr)

        else:
            o = ObsClass.from_lizard(
                code, tube_nr, tmin, tmax, type_timeseries, only_metadata=only_metadata
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
    source="vitens",
):
    """Get all observations within a specified extent.

    Parameters
    ----------
    extent : list or shapefile
        get groundwater monitoring wells wihtin this extent [xmin, xmax, ymin, ymax]
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
    source : str
        source indicating URL endpoint, currently only "vitens" is officially supported.


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

    lizard_GWS_endpoint = f"{LIZARD_APIS[source]}groundwaterstations/"
    url_groundwaterstation_extent = (
        f"{lizard_GWS_endpoint}?geometry__within={polygon_T}&page_size={page_size}"
    )

    groundwaterstation_data = requests.get(url_groundwaterstation_extent).json()
    nr_results = groundwaterstation_data["count"]
    nr_pages = math.ceil(nr_results / page_size)

    print("Number of monitoring wells: {}".format(nr_results))
    print("Number of pages: {}".format(nr_pages))

    if nr_threads > nr_pages:
        nr_threads = nr_pages

    urls = _prepare_API_input(nr_pages, url_groundwaterstation_extent)

    arg_tuple = (ObsClass, tube_nr, tmin, tmax, type_timeseries, only_metadata)
    codes = []
    with ThreadPoolExecutor(max_workers=nr_threads) as executor:
        for result in tqdm(executor.map(_download, urls), total=nr_pages, desc="Page"):
            codes += [(d["code"],) + arg_tuple for d in result]

    obs_list = []
    with ThreadPoolExecutor() as executor:
        for obs_list_mw in tqdm(
            executor.map(lambda args: get_obs_list_from_codes(*args), codes),
            total=len(codes),
            desc="monitoring well",
        ):
            obs_list += obs_list_mw

    return obs_list
