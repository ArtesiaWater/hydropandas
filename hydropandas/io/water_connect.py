import json
import logging
import zipfile
from io import BytesIO, StringIO
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from platformdirs import user_data_dir
from shapely.geometry import box
from tqdm import tqdm

logger = logging.getLogger(__name__)

DH_zip_url = (
    "https://www.waterconnect.sa.gov.au/Content/Downloads/DEW/WATER_Drillholes_shp.zip"
)


def get_obs_list_from_extent(
    extent,
    ObsClass,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=False,
    location_gdf=None,
    update=False,
    **kwargs,
):
    """Get observations within a specific extent and optionally for a specific set of
    locations.

    Parameters
    ----------
    extent : list, tuple, numpy-array or None, optional
        get water connect measurements within this extent
        [xmin, xmax, ymin, ymax]
    ObsClass : type
        class of the observations, e.g. GroundwaterObs
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : bool, optional
        if False, only observations with measurements are kept. The default
        is True.
    location_gdf : GeoDataFrame, optional
        geodataframe with the locations of the water drill holes you want to include.
    update : bool, optional
        if True new locations are downloaded and stored locally (slow) otherwise a
        cached version of the locations is used. By default False
    **kwargs
        additional keyword arguments are passed to the ObsClass.from_waterconnect()
        method

    Returns
    -------
    obs_list : list
        list with Obs objects

    """

    if only_metadata and not keep_all_obs:
        logger.warning(
            "Option with only_metadata is True and keep_all_obs is False not supported"
            "setting keep_all_obs=True"
        )
        keep_all_obs = True

    if location_gdf is None:
        location_gdf = get_locations_gdf(update=update)

    if extent is not None:
        location_gdf = get_locations_within_extent(location_gdf, extent)

    if location_gdf.empty:
        msg = f"No water connect measurements found within extent {extent}"
        logger.warning(msg)
        return []

    nr_obs_wells = len(location_gdf)
    logger.info(f"downloading waterconnect data from {nr_obs_wells} observation points")

    obs_list = []
    for dh_no, meta_series in tqdm(
        location_gdf.iterrows(), total=nr_obs_wells, desc="obs well"
    ):
        o = ObsClass.from_waterconnect(
            dh_no,
            meta_series,
            tmin=tmin,
            tmax=tmax,
            only_metadata=only_metadata,
            **kwargs,
        )
        if o.empty and not keep_all_obs:
            continue

        obs_list.append(o)

    return obs_list


def get_locations_gdf(fdir=None, keep_cols="all", update=False):
    """get locations of water drillholes in Southern Australia from water connect.

    Parameters
    ----------
    fdir : Path or str
        directory to store geodataframe with locations. If None a directory is created.
        The default is None.
    keep_cols : list or str, optional
        if provided only these columns of the geodataframe are stored. If 'all' all the
        columns are returned, by default 'all'.
    update : bool, optional
        if True new locations are downloaded and stored locally (slow) otherwise a
        cached version of the locations is used. By default False

    Returns
    -------
    geopandas.GeoDataFrame
        locations of water drillholes in Southern Australia
    """

    # set directories
    if fdir is None:
        fdir = Path(user_data_dir("waterconnect", "hydropandas"))
    elif isinstance(str, fdir):
        fdir = Path(fdir)
    fdir.mkdir(parents=True, exist_ok=True)
    fname_feather = fdir / "WATER_Drillholes_GDA2020.feather"

    # read files from feather
    if fname_feather.exists() and not update:
        logger.debug(f"Reading Water Connect locations from {fname_feather}")
        gdf_sel = gpd.read_feather(fname_feather)
        return gdf_sel

    msg = (
        "Initializing Water Connect, this may take several minutes. This only has to "
        "be done the first time you use Water Connect via hydropandas."
    )
    logger.info(msg)

    # download data
    r = requests.get(DH_zip_url, timeout=3600)
    r.raise_for_status()

    # extract files from zip and write to disk
    zip_buffer = BytesIO(r.content)
    with zipfile.ZipFile(zip_buffer, "r") as z:
        target_files = [
            name for name in z.namelist() if name.startswith("WATER_Drillholes_GDA2020")
        ]
        # Extract only those files
        logger.debug(f"Writing Water Connect locations shape to {fdir}")
        for file in target_files:
            z.extract(file, fdir)

    # read shapefiles to geodataframe
    logger.debug(f"Reading Water Connect locations from {fdir}")
    gdf_gda2020 = gpd.read_file(fdir / "WATER_Drillholes_GDA2020.shp")

    # select columns
    if keep_cols is None:
        keep_cols = [
            "UNIT_NO",
            "DHNO",
            "NAME",
            "LAT",
            "LON",
            "REF_ELEV",
            "GRND_ELEV",
            "STATE",
            "STATUS",
            "STAT_DESC",
            "PURPOSE",
            "PURP_DESC",
            "DRILL_DATE",
            "ORIG_DEPTH",
            "MAX_DEPTH",
            "CASE_TO",
            "geometry",
        ]
    elif isinstance(keep_cols, str) and keep_cols == "all":
        keep_cols = gdf_gda2020.columns

    if "geometry" not in keep_cols:
        keep_cols = list(keep_cols)
        keep_cols.append("geometry")

    gdf_sel = gdf_gda2020[keep_cols]

    # set index column
    gdf_sel.set_index("DHNO", inplace=True)

    # write to feather
    logger.debug(f"Write Water Connect locations to feather {fname_feather}")
    gdf_sel.to_feather(fname_feather)

    return gdf_sel


def get_locations_within_extent(gdf, extent=None):
    """get drillhole locations within an extent

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        all drillhole locations in Southern Australia
    extent : list, tuple, np.array, optional
        the extent, by default None

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame with locations in the given extent
    """

    if extent is None:
        return gdf

    pol_extent = box(*tuple(np.asarray(extent)[[0, 2, 1, 3]]))
    gdf_extent = gdf.loc[gdf.within(pol_extent)]
    return gdf_extent


def get_waterconnect_obs(
    dh_no,
    meta_series=None,
    tmin=None,
    tmax=None,
    only_metadata=False,
    verify=True,
    pumping=True,
    anomalous=True,
    **kwargs,
):
    """get waterconnect observations using the API

    Parameters
    ----------
    dh_no : int or str
        drill hole number
    meta_series : pd.Series, optional
        series with metadata. Typically a row of the locations gdf.
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True only metadata and no measurements are returned. BY default False
    verify : bool, optional
        use verification to get a secure connection
    pumping : bool, optional
        return observations from pumping wells
    anomalous : bool, optional
        return anomalous observations
    **kwargs
        kwargs are passed to 'get_location_gdf'

    Returns
    -------
    pd.DataFrame: measurements

    dict: metadata

    """
    if meta_series is None:
        gdf = get_locations_gdf(**kwargs)
        meta_series = gdf.loc[dh_no]
    elif isinstance(meta_series, pd.DataFrame):
        meta_series = gdf.loc[dh_no]

    if only_metadata:
        measurements = pd.DataFrame()
        source = "water connect"
        unit = ""
    else:
        logger.debug(f"reading measurements for {dh_no}")
        measurements = get_measurements(
            dh_no,
            tmin=tmin,
            tmax=tmax,
            verify=verify,
            pumping=pumping,
            anomalous=anomalous,
        )

        if measurements.empty:
            logger.debug(f"no measurements for {dh_no}")
            source = "water connect"
            unit = ""
        else:
            if "data_source" in measurements.columns:
                source = (
                    ", ".join(measurements["data_source"].unique()) + " (water connect)"
                )
            else:
                source = "water connect"
            unit = "m AHD"

    meta = {
        "name": dh_no,
        "x": meta_series.pop("LON"),
        "y": meta_series.pop("LAT"),
        "unit": unit,
        "location": meta_series.pop("UNIT_NO"),
        "ground_level": meta_series.pop("GRND_ELEV"),
        "source": source,
        "meta": meta_series.to_dict(),
    }

    return measurements, meta


def get_measurements(
    dh_no, tmin=None, tmax=None, verify=True, pumping=True, anomalous=True
):
    """get measurement from a dh_no using the water connect api.

    Parameters
    ----------
    dh_no : int or str
        drill hole number
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    verify : bool, optional
        use verification to get a secure connection
    pumping : bool, optional
        return observations from pumping wells
    anomalous : bool, optional
        return anomalous observations

    Returns
    -------
    pd.DataFrame with measurements, an empty DataFrame is returned if there are no
    measurements
    """
    result_str = request_api(dh_no, verify=verify, pumping=pumping, anomalous=anomalous)
    df = pd.read_csv(result_str)

    if df.empty:
        return df

    df.index = pd.to_datetime(df.pop("obs_date"), format="%d/%m/%Y")
    df = df.loc[tmin:tmax]

    if "rswl" in df.columns:
        df.insert(0, "rswl", df.pop("rswl"))

    return df


def request_api(
    dh_no, verify=True, pumping=True, anomalous=True, fname=None
) -> StringIO:
    """Download water connect data from the API

    Parameters
    ----------
    dh_no : int or str
        drill hole number
    verify : bool, optional
        use verification to get a secure connection
    pumping : bool, optional
        return observations from pumping wells
    anomalous : bool, optional
        return anomalous observations
    fname : str or Path, optional
        write requested data to this filename. Only if fname is not None

    Returns
    -------
    StringIO

    """
    url = "https://www.waterconnect.sa.gov.au/_layouts/15/dfw.sharepoint.wdd/WDDDMS.ashx/GetWaterLevelDownload?bulkOutput=CSV"
    json_data = {"DHNOs": [dh_no], "Pumping": pumping, "Anomalous": anomalous}
    r = requests.post(
        url, verify=verify, data={"exportdata": json.dumps(json_data)}, timeout=60
    )
    r.raise_for_status()

    result_str = r.text

    if result_str.startswith("<!DOCTYPE html>"):
        raise RuntimeError("water connect API down\n" + result_str)

    if fname is not None:
        with open(fname, "w") as f:
            f.write(result_str)

    return StringIO(result_str)
