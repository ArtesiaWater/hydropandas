import datetime as dt
import logging
import os
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box
from tqdm import tqdm

logger = logging.getLogger(__name__)


def get_obs_list_from_extent(
    extent,
    ObsClass,
    grootheid_code=None,
    locatie=None,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=False,
    epsg=28992,
):
    """Get observations within a specific extent and optionally for a specific location
    and grootheid_code.

    Parameters
    ----------
    extent : list, tuple, numpy-array or None, optional
        get waterinfo measurements within this extent
        [xmin, xmax, ymin, ymax]
    ObsClass : type
        class of the observations, e.g. WaterlvlObs
    grootheid_code : str, optional
        select only measurement with this grootheid_code, e.g. 'WATHTE', default is None
    locatie : str, optional
        select only measurement with this location, e.g. 'SCHOONHVN', default is None
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
    epsg : int, optional
        epsg code of the extent. The default is 28992 (RD).

    Returns
    -------
    obs_list : list
        list with Obs objects

    """

    gdf = get_locations_gdf(extent=extent, epsg=epsg)
    if gdf.empty:
        msg = f"No waterinfo measurements found within extent {extent}"
        logger.warning(msg)
        return []

    if locatie is not None:
        gdf = gdf.loc[locatie]
    if grootheid_code is not None:
        gdf = gdf.loc[gdf["Grootheid.Code"] == grootheid_code]
    if gdf.empty:
        msg = f"No waterinfo measurements found for locatie {locatie} and grootheid_code {grootheid_code}"
        logger.warning(msg)
        return []

    if only_metadata:
        obs_list = []
        for _, row in gdf.iterrows():
            meta = _get_metadata_from_series(row)
            o = ObsClass(meta=meta, **meta)
            obs_list.append(o)
        return obs_list

    obs_list = []
    for _, row in gdf.iterrows():
        o = ObsClass.from_waterinfo(location_gdf=row, tmin=tmin, tmax=tmax)
        if not keep_all_obs and o.empty:
            continue
        obs_list.append(o)

    return obs_list


def get_waterinfo_obs(
    path=None,
    location_gdf=None,
    grootheid_code=None,
    locatie=None,
    tmin=None,
    tmax=None,
    **kwargs,
):
    """Get waterinfo observations from a file or ddlpy

    Parameters
    ----------
    path : str, optional
        path to waterinfo file (.zip or .csv), default is None
    location_gdf : geopandas.GeoDataFrame, optional
        geodataframe with locations, default is None
    grootheid_code : str, optional
        code of the grootheid, e.g. 'WATHTE', default is None
    locatie : str, optional
        name of the location, e.g. 'SCHOONHVN', default is None
    tmin : datetime, optional
        start date of the measurements, default is None
    tmax : datetime, optional
        end date of the measurements, default is None

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with measurements
    meta : dict
        dict with metadata
    """

    if path is not None:
        df, meta = read_waterinfo_file(path, **kwargs)
    elif location_gdf is not None or (
        grootheid_code is not None and locatie is not None
    ):
        df, meta = get_measurements_ddlpy(
            location_gdf, grootheid_code, locatie, tmin, tmax
        )
    else:
        raise ValueError("Provide path or grootheid_code and locatie!")

    return df, meta


def _get_metadata_from_series(selected):
    """Get metadata from a series with location information

    Parameters
    ----------
    selected : pandas.Series
        series with location information

    Returns
    -------
    meta : dict
        dict with metadata
    """
    d = selected.to_dict()
    p = d.pop("geometry")
    meta = {
        "name": d["Naam"] + " " + d.pop("Grootheid.Omschrijving"),
        "unit": d.pop("Eenheid.Code") + " " + d.pop("Hoedanigheid.Code"),
        "x": p.x,
        "y": p.y,
        "source": "waterinfo (ddlpy)",
        "location": d.pop("Naam"),
        "meta": d,
    }
    return meta


def get_measurements_ddlpy(
    location_gdf=None, grootheid_code=None, locatie=None, tmin=None, tmax=None
):
    """Get measurements from ddlpy for a specific location and grootheid_code

    Parameters
    ----------
    location_gdf : geopandas.GeoDataFrame, optional
        geodataframe with one or more locations, default is None
    grootheid_code : str
        code of the grootheid
    locatie : str
        name of the location
    tmin : datetime, optional
        start date of the measurements, default is 2025-01-01
    tmax : datetime, optional
        end date of the measurements, default is now

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with measurements
    meta : dict
        dict with metadata
    """

    import ddlpy

    if tmin is None:
        tmin = dt.datetime(2024, 1, 1)
    if tmax is None:
        tmax = dt.datetime.now()

    if location_gdf is None:
        location_gdf = get_locations_gdf()

    if isinstance(location_gdf, pd.Series):
        selected = location_gdf
        locatie = selected.name
        grootheid_code = selected["Grootheid.Code"]
    else:
        selected = location_gdf.loc[
            location_gdf["Grootheid.Code"] == grootheid_code
        ].loc[locatie]
        if selected.empty:
            raise ValueError(f"No location found for {locatie} and {grootheid_code}")
        elif isinstance(selected, pd.DataFrame):
            if len(selected) == 1:
                selected = selected.iloc[0]
            else:
                logger.warning("Multiple locations found, selecting last one")
                selected = selected.iloc[-1]

    df = ddlpy.measurements(selected, start_date=tmin, end_date=tmax)

    if df.empty:
        msg = f"no measurements for {locatie} en {grootheid_code} between {tmin} and {tmax}"
        logger.info(msg)
    else:
        if "Meetwaarde.Waarde_Numeriek" in df.columns:
            df = df[["Meetwaarde.Waarde_Numeriek"]]
            df.columns = ["value"]

    meta = _get_metadata_from_series(selected)
    return df, meta


def get_locations_gdf(extent=(482.06, 306602.42, 284182.97, 637049.52), epsg=28992):
    """Get locations from ddlpy and return as geodataframe

    Parameters
    ----------
    extent : tuple, optional
        extent of the locations. The default is the extent of the Netherlands (RD).

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        geodataframe with locations. This dataframe is needed to obtain measurements
        using ddlpy
    """

    import ddlpy

    locations = ddlpy.locations()
    geometries = gpd.points_from_xy(locations["X"], locations["Y"])
    gdf = gpd.GeoDataFrame(locations, geometry=geometries, crs=25831)
    gdf.to_crs(epsg, inplace=True)
    polygon_ext = box(*tuple(np.array(extent)[[0, 2, 1, 3]]))
    gdf = gdf.loc[gdf.within(polygon_ext)]

    return gdf


def read_waterinfo_file(
    path,
    index_cols=None,
    return_metadata=True,
    value_col=None,
    location_col=None,
    xcol=None,
    ycol=None,
    transform_coords=True,
):
    """Read waterinfo file (CSV or zip)

    Parameters
    ----------
    path : str
        path to waterinfo file (.zip or .csv)

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing file content
    metadata : dict, optional
        dict containing metadata, returned if return_metadata is
        True, default is False
    """
    from pyproj import Transformer

    name = os.path.splitext(os.path.basename(path))[0]

    if path.endswith(".csv"):
        f = path
    elif path.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        f = zf.open("{}.csv".format(name))
    else:
        raise NotImplementedError(
            "File type '{}' not supported!".format(os.path.splitext(path)[-1])
        )

    if value_col is None:
        value_col = "NUMERIEKEWAARDE"

    if location_col is None:
        location_col = "MEETPUNT_IDENTIFICATIE"

    if xcol is None:
        xcol = "X"

    if ycol is None:
        ycol = "Y"

    # read data
    df = pd.read_csv(
        f,
        sep=";",
        decimal=",",
        encoding="ISO-8859-1",
        dayfirst=True,
    )

    if index_cols is None:
        index_cols = ["WAARNEMINGDATUM"]
        if "WAARNEMINGTIJD (MET/CET)" in df.columns:
            index_cols += ["WAARNEMINGTIJD (MET/CET)"]
        elif "WAARNEMINGTIJD" in df.columns:
            index_cols += ["WAARNEMINGTIJD"]
        else:
            raise KeyError("expected column with WAARNEMINGTIJD but could not find one")

    df.index = pd.to_datetime(
        df[index_cols[0]] + " " + df[index_cols[1]], dayfirst=True
    )
    df.drop(columns=index_cols, inplace=True)

    # do some conversions
    df.loc[df[value_col] == 999999999, value_col] = np.nan
    df[value_col] = df[value_col] / 100.0

    # parse metadata into dict
    if return_metadata:
        if len(df[location_col].unique()) > 1:
            raise ValueError(
                "File contains data for more than one location!"
                " Use ObsCollection.from_waterinfo()!"
            )

        metadata = {}

        if transform_coords:
            transformer = Transformer.from_crs("epsg:25831", "epsg:28992")
            x, y = transformer.transform(df[xcol].iloc[-1], df[ycol].iloc[-1])
        else:
            x = df[xcol].iloc[-1] / 100.0
            y = df[ycol].iloc[-1] / 100.0
        metadata["name"] = df[location_col].iloc[-1]
        metadata["x"] = x
        metadata["y"] = y
        metadata["filename"] = f
        metadata["source"] = "waterinfo"

        return df, metadata
    else:
        return df


def read_waterinfo_obs(file_or_dir, ObsClass, progressbar=False, **kwargs):
    """Read waterinfo file or directory and extract locations and observations.

    Parameters
    ----------
    file_or_dir : str
        path to file or directory
    ObsClass: Obs type
        type of Obs to store data in
    progressbar : bool, optional
        show progressbar if True, default is False

    Returns
    -------
    obs_collection : list
        list of Obs objects
    """
    from pyproj import Transformer

    # Waterinfo file
    if os.path.isfile(file_or_dir):
        files = [file_or_dir]
    # directory with waterinfo files (zips or csvs)
    elif os.path.isdir(file_or_dir):
        files = [os.path.join(file_or_dir, f) for f in sorted(os.listdir(file_or_dir))]
    else:
        raise NotImplementedError("Provide path to file or directory!")

    location_col = kwargs.pop("location_col", "MEETPUNT_IDENTIFICATIE")

    # loop over files
    metadata = {}
    obs_collection = []

    transformer = Transformer.from_crs("epsg:25831", "epsg:28992")

    for filenm in tqdm(files) if progressbar else files:
        # read file or zip
        df = read_waterinfo_file(
            filenm, location_col=location_col, return_metadata=False, **kwargs
        )

        # get location and convert to m RD
        for stn in df[location_col].unique():
            mask = df[location_col] == stn
            x, y = transformer.transform(
                df.loc[mask, "X"].iloc[-1], df.loc[mask, "Y"].iloc[-1]
            )
            metadata = {
                "name": stn,
                "x": x,
                "y": y,
                "filename": filenm,
                "source": "waterinfo",
            }

            # add to list
            o = ObsClass(df.loc[mask, :], meta=metadata, **metadata)
            obs_collection.append(o)

    return obs_collection
