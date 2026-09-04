"""Created on Wed Sep 12 12:15:42 2018.

@author: Artesia
"""

import logging
import os
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import pandas as pd
import pyproj
from scipy.interpolate import RBFInterpolator

logger = logging.getLogger(__name__)

EPSG_28992 = (
    "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 "
    "+x_0=155000 +y_0=463000 +ellps=bessel "
    "+towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m "
    "+no_defs"
)


def get_transformer28992(crs_from, crs_to, always_xy=True, **kwargs):
    """This is simply a wrapper around pyproj.Transformer.from_crs in order
    to handle the special case of EPSG:28992.

    Parameters
    ----------
    crs_from : pyproj.CRS
        source coordinate reference system.
    crs_to : pyproj.CRS
        target coordinate reference system.
    **kwargs are passed to pyproj.Transformer.from_crs.

    Returns
    -------
    pyproj.Transformer
        transformer object to convert coordinates from crs_from to crs_to.
    """
    if crs_to == pyproj.CRS(28992):
        crs_to = pyproj.CRS(EPSG_28992)
    if crs_from == pyproj.CRS(28992):
        crs_from = pyproj.CRS(EPSG_28992)

    transformer = pyproj.Transformer.from_crs(
        crs_from, crs_to, always_xy=always_xy, **kwargs
    )

    return transformer


def _obslist_to_frame(obs_list):
    """Convert a list of observations to a pandas DataFrame.

    Parameters
    ----------
    obs_list : list of hydropandas.*Obs
        list containing *Obs objects that will be stored in DataFrame. All
        observations must have the same crs value (or an empty string).

    Returns
    -------
    obs_df : pandas.DataFrame
        DataFrame containing all data
    crs : pyproj.CRS
        coordinate reference system of the observations, if available
    """
    if len(obs_list) > 0:
        obs_df = pd.DataFrame(
            [o.to_collection_dict() for o in obs_list],
            columns=obs_list[0].to_collection_dict().keys(),
        )

        # infer crs from the observations
        crs = {c for c in obs_df.pop("crs") if c != ""}
        if len(crs) > 1:
            raise ValueError(
                "multiple crs values in observations, an ObsCollection can only have one crs value"
            )
        crs = next(iter(crs), "")

        obs_df.set_index("name", inplace=True)
        if obs_df.index.duplicated().any():
            logger.warning("multiple observations with the same name")
    else:
        obs_df = pd.DataFrame()
        crs = ""

    return obs_df, crs


def unzip_file(src, dst, force=False, preserve_datetime=False):
    """Unzip file.

    Parameters
    ----------
    src : str or pathlib.Path
        source zip file
    dst : str or pathlib.Path
        destination directory
    force : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        use date of the zipfile for the destination file

    Returns
    -------
    int
        1 of True
    """
    dst = Path(dst)
    if dst.exists() and not force:
        print(
            "File not unzipped. Destination already exists. Use'force=True' to unzip."
        )
        return
    if preserve_datetime:
        zipf = zipfile.ZipFile(src, "r")
        for f in zipf.infolist():
            zipf.extract(f, path=dst)
            date_time = time.mktime(f.date_time + (0, 0, -1))
            os.utime(dst / f.filename, (date_time, date_time))
        zipf.close()
    else:
        zipf = zipfile.ZipFile(src, "r")
        zipf.extractall(dst)
        zipf.close()
    return 1


def get_files(
    file_or_dir, ext, unpackdir=None, force_unpack=False, preserve_datetime=False
):
    """Internal method to get list of files with specific extension from dirname.

    Parameters
    ----------
    file_or_dir : str or pathlib.Path
        file or path to data.
    ext : str
        extension of filenames to store in list.
    unpackdir : str or pathlib.Path
        directory to story unpacked zip file, only used in case of a zipfile.
    force_unpack : bool, optional
        force unzip, by default False.
    preserve_datetime : bool, optional
        preserve datetime of unzipped files, by default False. Used for
        checking whether data has changed.
    """
    file_or_dir = Path(file_or_dir)
    if unpackdir is not None:
        unpackdir = Path(unpackdir)

    # check if unpackdir is same as file_or_dir, if same, this can cause
    # problems when the unpackdir still contains zips that will be unpacked
    # again.
    if (unpackdir is not None) and (unpackdir.resolve() == file_or_dir.resolve()):
        raise ValueError("Please specify a different folder to unpack files!")

    # identify whether file_or_dir started as zip
    iszip = file_or_dir.suffix == ".zip"

    # unzip dir
    if iszip:
        zipf = file_or_dir
        if unpackdir is None:
            file_or_dir = Path(tempfile.TemporaryDirectory().name)
        else:
            file_or_dir = unpackdir
        unzip_file(
            zipf, file_or_dir, force=force_unpack, preserve_datetime=preserve_datetime
        )

    # file_or_dir is directory
    if file_or_dir.is_dir():
        # check for zips in dir
        zip_fnames = [i for i in file_or_dir.iterdir() if i.suffix == ".zip"]
        if len(zip_fnames) > 0:
            # unzip zips
            if unpackdir is None:
                dirname = Path(tempfile.TemporaryDirectory().name)
            else:
                dirname = unpackdir
            for zipf in zip_fnames:
                unzip_file(
                    zipf,
                    dirname,
                    force=True,
                    preserve_datetime=preserve_datetime,
                )
                # remove intermediate zipfiles if initial file_or_dir was zip
                if iszip:
                    zipf.unlink()
        else:
            dirname = file_or_dir
        # get all files with extension ext
        unzip_fnames = [i.name for i in dirname.iterdir() if i.name.endswith(ext)]
    elif file_or_dir.is_file():
        # file_or_dir is actually an xml
        unzip_fnames = [file_or_dir.name]  # get file name
        dirname = file_or_dir.parent  # get directory path
    else:
        raise NotImplementedError(f"Cannot parse 'file_or_dir': {file_or_dir}!")

    return dirname, unzip_fnames


def df2gdf(df, xcol="x", ycol="y", crs=28992, custom_crs_28992=False):
    """Create a GeoDataFrame from a DataFrame with xy points.

    Parameters
    ----------
    df : pd.DataFrame
        input dataframe
    xcol : str, optional
        column name with x values. The default is  'x'.
    ycol : str, optional
        column name with y values. The default is  'x'.
    crs : int, optional
        coordinate reference system, by default 28992 (RD new).
    custom_crs_28992 : bool, optional
        if True, use a custom definition for EPSG:28992 instead of the default one.
        In some cases the default EPSG:28992 definition gives incorrect results
        when converting to another crs, so a custom definition may be necessary.
        The default is False.

    Returns
    -------
    geopandas GeoDataFrame
        geodataframe
    """
    from geopandas import GeoDataFrame
    from shapely.geometry import Point

    if pyproj.CRS(crs) == pyproj.CRS(28992) and custom_crs_28992:
        crs = pyproj.CRS(EPSG_28992)

    gdf = GeoDataFrame(
        df.copy(),
        geometry=[Point((s[xcol], s[ycol])) for _, s in df.iterrows()],
        crs=crs,
    )
    return gdf


class ColoredFormatter(logging.Formatter):
    """Colored log formatter.

    Taken from
    https://gist.github.com/joshbode/58fac7ababc700f51e2a9ecdebe563ad
    """

    def __init__(self, *args, colors: dict[str, str] | None = None, **kwargs) -> None:
        """Initialize the formatter with specified format strings."""

        super().__init__(*args, **kwargs)

        self.colors = colors if colors else {}

    def format(self, record) -> str:
        """Format the specified record as text."""

        record.color = self.colors.get(record.levelname, "")
        record.reset = "\x1b[0m"

        return super().format(record)


def get_color_logger(level="INFO", logger_name=None):
    """Get a logger with colored output.

    Parameters
    ----------
    level : str, optional
        The logging level to set for the logger. Default is "INFO".

    Returns
    -------
    logger : logging.Logger
        The configured logger object.
    """
    if level == "DEBUG":
        FORMAT = "{color}{levelname}:{name}.{funcName}:{lineno}:{message}{reset}"
    else:
        FORMAT = "{color}{levelname}:{name}.{funcName}:{message}{reset}"
    formatter = ColoredFormatter(
        FORMAT,
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        colors={
            "DEBUG": "\x1b[36m",
            "INFO": "\x1b[32m",
            "WARNING": "\x1b[33m",
            "ERROR": "\x1b[31m",
            "CRITICAL": "\x1b[31m" + "\x1b[47m" + "\x1b[1m",
        },
    )

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    clogger = logging.getLogger(logger_name)
    clogger.handlers[:] = []
    clogger.addHandler(handler)
    clogger.setLevel(getattr(logging, level))

    logging.captureWarnings(True)
    return clogger


def oc_to_df(oc, col: str | None = None) -> pd.DataFrame:
    """Convert an observation collection to a DataFrame where every column has one
    observation.

    Parameters
    ----------
    oc : hydropandas ObsCollection
        observation collection
    col : Optional[str], optional
        Name of a column in hte observation collection, by default None

    Returns
    -------
    DataFrame
        _description_
    """
    df_list = []
    for o in oc.obs.values:
        if not o.empty:
            if col is None:
                vals = o.loc[:, o._get_first_numeric_col_name()]
            else:
                vals = o.loc[:, col]
            vals.name = o.name
            df_list.append(vals)
    return pd.concat(df_list, axis=1)


def interpolate(
    xy: list[list[float]],
    obsdf: pd.DataFrame,
    obsloc: pd.DataFrame,
    kernel: str = "thin_plate_spline",
    kernel2: str = "linear",
    epsilon: int | None = None,
) -> pd.DataFrame:
    """Interpolation method using the Scipy radial basis function (RBF)

    Parameters
    ----------
    xy : List[List[float]]
        xy coordinates of locations of interest e.g. [[10,25], [5,25]]
    obsdf : DataFrame
        Dataframe containing the observation locations as columns and
        the observations at a measurement time in each row.
    obsloc : DataFrame
        Dataframe containing the observation locations coordinates
        with observation locations as index and columns ["x", "y"]
    kernel : str, optional
        Type of radial basis funtion, by default thin_plate_spline.
        Other options are linear, gaussian, inverse_quadratic,
        multiquadric, inverse_multiquadric, cubic or quintic.
    kernel2 : str, optional
        Kernel in case there are not enough observations (3 or 6) for
        time step, by default linear. Other options are gaussian,
        inverse_quadratic, multiquadric, or inverse_multiquadric.
    epsilon : Optional[int], optional
        Shape parameter that scales the input to the RBF. If kernel is
        linear, thin_plate_spline, cubic, or quintic, this defaults to 1.
        Otherwise this must be specified.

    Returns
    -------
    DataFrame
        DataFrame with locations of interest as columns and interpolated values
        at a measurement time in each row.
    """

    if kernel in {"thin_plate_spline", "cubic"}:
        min_val = 3
    elif kernel == "quintic":
        min_val = 6
    else:
        min_val = len(obsdf.index)

    fill_df = (
        pd.DataFrame(
            index=obsdf.index,
            columns=[f"{int(_xy[0])}_{int(_xy[1])}" for _xy in xy],
        )
        .sort_index()
        .astype(float)
    )

    for idx in obsdf.index:
        # get all stations with values for this date
        val = obsdf.loc[idx].dropna()
        # get stations for this date
        coor = obsloc.loc[val.index]

        if len(val) < min_val:
            kernel = kernel2

        # create an scipy interpolator
        rbf = RBFInterpolator(
            coor.to_numpy(), val.to_numpy(), epsilon=epsilon, kernel=kernel
        )

        val_rbf = rbf(xy)
        fill_df.loc[idx] = val_rbf

    return fill_df
