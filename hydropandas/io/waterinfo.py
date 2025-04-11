import logging
import os
import zipfile
from functools import lru_cache

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box
from tqdm import tqdm

logger = logging.getLogger(__name__)


def get_obs_list_from_extent(
    extent,
    ObsClass,
    locatie=None,
    grootheid_code=None,
    groepering_code=None,
    parameter_code=None,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=False,
    epsg=28992,
    location_gdf=None,
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
    locatie : str or list of str, optional
        select only measurement with this location(s), e.g. 'SCHOONHVN', default is None
    grootheid_code : str or list of str, optional
        select only measurement with this grootheid_code, e.g. 'WATHTE', default is None
    groepering_code : str or list of str, optional
        select only measurement with this groepering_code, e.g. 'GETETBRKD2', default is None
    parameter_code :  str or list of str, optional
        select only measurement with this parameter_code, e.g. 'Cl', default is None
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
    location_gdf : GeoDataFrame, optional
        geodataframe with the locations of the measurements you want to include. If
        location_gdf is provided the provided extent and epgs will be ignored.

    Returns
    -------
    obs_list : list
        list with Obs objects

    """
    if location_gdf is None:
        gdf = get_locations_gdf(epsg=epsg)
        gdf = get_locations_within_extent(gdf, extent=extent)
    else:
        gdf = location_gdf
        if gdf.empty:
            msg = f"No waterinfo measurements found within extent {extent}"
            logger.warning(msg)
            return []

    gdf = _select_location(
        gdf, locatie, grootheid_code, groepering_code, parameter_code
    )

    if gdf.empty:
        return []

    logger.info(
        f"downloading waterinfo measurements from {len(gdf)} observation points"
    )

    obs_list = []
    onames = []
    for _, row in gdf.iterrows():
        if only_metadata:
            meta = _get_metadata_from_series(row)
            o = ObsClass(meta=meta, **meta)
        else:
            o = ObsClass.from_waterinfo(location_gdf=row, tmin=tmin, tmax=tmax)
            if not keep_all_obs and o.empty:
                continue
        # rename if observation name already exists
        if o.name in onames:
            counter = 1
            new_name = f"{o.name} ({counter})"
            while new_name in onames:
                counter += 1
                new_name = f"{o.name} ({counter})"
            setattr(o, "name", new_name)
        obs_list.append(o)
        onames.append(o.name)

    return obs_list


def get_waterinfo_obs(
    path=None,
    location_gdf=None,
    locatie=None,
    grootheid_code=None,
    groepering_code=None,
    parameter_code=None,
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
    locatie : str or list of str, optional
        name of the location, e.g. 'SCHOONHVN', default is None
    grootheid_code : str or list of str, optional
        code(s) of the grootheid, e.g. 'WATHTE', default is None
    groepering_code : str or list of str, optional
        code(s) of the groepering, e.g. 'GETETBRKD2', default is None
    parameter_code : str or list of str, optional
        code(s) of the parameter, e.g. 'Cl', default is None
    tmin : pd.Timestamp, optional
        start date of the measurements, default is None
    tmax : pd.Timestamp, optional
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
    else:
        df, meta = get_measurements_ddlpy(
            location_gdf,
            locatie,
            grootheid_code,
            groepering_code,
            parameter_code,
            tmin,
            tmax,
        )

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
    unit = d.pop("Eenheid.Code")
    if d["Hoedanigheid.Code"] not in ["", "NVT"]:
        unit += " " + d.pop("Hoedanigheid.Code")
    meta = {
        "name": d["Naam"] + " " + d.pop("Grootheid.Omschrijving"),
        "unit": unit,
        "x": p.x,
        "y": p.y,
        "source": "waterinfo (ddlpy)",
        "location": d.pop("Naam"),
        "meta": d,
    }
    return meta


def _select_location(
    location_gdf, locatie, grootheid_code, groepering_code, parameter_code
):
    """Select location from a geodataframe with locations

    Parameters
    ----------
    location_gdf : geopandas.GeoDataFrame
        geodataframe with locations
    locatie : str or list of str
        name(s) of the location
    grootheid_code : str or list of str
        code(s) of the grootheid
    groepering_code : str or list of str
        code(s) of the groepering, e.g. 'GETETBRKD2', default is None
    parameter_code : str or list of str
        code(s) of the parameter, e.g. 'Cl', default is None

    Returns
    -------
    selected : pandas.Series
        series with location information
    """

    if locatie is not None:
        if isinstance(locatie, str):
            locatie = [locatie]
        location_gdf = location_gdf.loc[locatie]
    if grootheid_code is not None:
        if isinstance(grootheid_code, str):
            grootheid_code = [grootheid_code]
        location_gdf = location_gdf.loc[
            location_gdf["Grootheid.Code"].isin(grootheid_code)
        ]
    if groepering_code is not None:
        if isinstance(groepering_code, str):
            groepering_code = [groepering_code]
        location_gdf = location_gdf.loc[
            location_gdf["Groepering.Code"].isin(groepering_code)
        ]

    if parameter_code is not None:
        if isinstance(parameter_code, str):
            parameter_code = [parameter_code]
        location_gdf = location_gdf.loc[
            location_gdf["Parameter.Code"].isin(parameter_code)
        ]

    if location_gdf.empty:
        raise ValueError(
            f"No location found for {locatie=}, {grootheid_code=} and {groepering_code=}"
        )

    return location_gdf


def get_measurements_ddlpy(
    location_gdf=None,
    locatie=None,
    grootheid_code=None,
    groepering_code=None,
    parameter_code=None,
    tmin=None,
    tmax=None,
):
    """Get measurements from ddlpy for a specific location and grootheid_code

    Parameters
    ----------
    location_gdf : geopandas.GeoDataFrame, optional
        geodataframe with one or more locations, default is None
    locatie : str or list of str, optional
        name(s) of the location
    grootheid_code : str or list of str, optional
        code(s) of the grootheid
    groepering_code : str or list of str, optional
        code(s) of the groepering
    parameter_code : str or list of str, optional
        code(s) of the parameter
    tmin : pd.Timestamp, optional
        start date of the measurements, default is 2025-01-01
    tmax : pd.Timestamp, optional
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
        tmin = pd.Timestamp("2024")
    else:
        tmin = pd.to_datetime(tmin)
    if tmax is None:
        tmax = pd.Timestamp.now()
    else:
        tmax = pd.to_datetime(tmax)

    if location_gdf is None:
        location_gdf = get_locations_gdf()

    if isinstance(location_gdf, pd.Series):
        selected = location_gdf
        locatie = selected.name
        grootheid_code = selected["Grootheid.Code"]
        groepering_code = selected["Groepering.Code"]
        parameter_code = selected["Parameter.Code"]
        df = ddlpy.measurements(selected, start_date=tmin, end_date=tmax)
    else:
        selected = _select_location(
            location_gdf, locatie, grootheid_code, groepering_code, parameter_code
        )
        if selected.empty:
            msg = (
                f"No location found for {locatie=}, {grootheid_code=}, "
                f"{groepering_code=} and {parameter_code=}"
            )
            raise ValueError(msg)
        elif isinstance(selected, pd.DataFrame):
            if len(selected) == 1:
                selected = selected.iloc[0]
                df = ddlpy.measurements(selected, start_date=tmin, end_date=tmax)
            else:
                logger.info(
                    "Multiple observation points match critera, select first one with measurements"
                )
                for _, row in selected.iterrows():
                    df = ddlpy.measurements(row, start_date=tmin, end_date=tmax)
                    if not df.empty:
                        break
                selected = row
        else:
            df = ddlpy.measurements(selected, start_date=tmin, end_date=tmax)

    if df.empty:
        msg = (
            f"No measurements for {locatie=}, {grootheid_code=}, "
            f"{groepering_code=} and {parameter_code=} between {tmin} and {tmax}"
        )
        logger.info(msg)
    else:
        if "Meetwaarde.Waarde_Numeriek" in df.columns:
            df = df[["Meetwaarde.Waarde_Numeriek"]]
            df.columns = ["value"]

        # remove time zone information by transforming to dutch winter time
        df.index = pd.to_datetime(df.index, utc=True).tz_localize(None) + pd.Timedelta(
            1, unit="h"
        )

    meta = _get_metadata_from_series(selected)
    return df, meta


@lru_cache()
def get_locations_gdf(epsg=28992):
    """Get locations from ddlpy and return as geodataframe

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

    return gdf


def get_locations_within_extent(gdf, extent=(482.06, 306602.42, 284182.97, 637049.52)):
    """Get locations from ddlpy and return as geodataframe

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        geodataframe with locations. This dataframe is needed to obtain measurements
        using ddlpy
    extent : tuple, optional
        extent of the locations. The default is the extent of the Netherlands (RD).

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        geodataframe with locations. This dataframe is needed to obtain measurements
        using ddlpy
    """

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
