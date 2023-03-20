# -*- coding: utf-8 -*-
"""Created on Wed Sep 12 12:15:42 2018.

@author: Artesia
"""

import logging
import os
import sys
import tempfile
import time
import zipfile
from typing import Dict, List, Optional

from colorama import Back, Fore, Style
from numpy import unique
from pandas import DataFrame, DatetimeIndex, Timedelta, Timestamp
from scipy.interpolate import RBFInterpolator

logger = logging.getLogger(__name__)


def _obslist_to_frame(obs_list):
    """convert a list of observations to a pandas DataFrame.

    Parameters
    ----------
    obs_list : list of hydropandas.*Obs
        list containing *Obs objects that will be stored in DataFrame.

    Returns
    -------
    obs_df : pandas.DataFrame
        DataFrame containing all data
    """
    if len(obs_list) > 0:
        obs_df = DataFrame(
            [o.to_collection_dict() for o in obs_list],
            columns=obs_list[0].to_collection_dict().keys(),
        )
        obs_df.set_index("name", inplace=True)
        if obs_df.index.duplicated().any():
            logger.warning("multiple observations with the same name")
    else:
        obs_df = DataFrame()

    return obs_df


def unzip_file(src, dst, force=False, preserve_datetime=False):
    """Unzip file.

    Parameters
    ----------
    src : str
        source zip file
    dst : str
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
    if os.path.exists(dst):
        if not force:
            print(
                "File not unzipped. Destination already exists. Use"
                "'force=True' to unzip."
            )
            return
    if preserve_datetime:
        zipf = zipfile.ZipFile(src, "r")
        for f in zipf.infolist():
            zipf.extract(f, path=dst)
            date_time = time.mktime(f.date_time + (0, 0, -1))
            os.utime(os.path.join(dst, f.filename), (date_time, date_time))
        zipf.close()
    else:
        zipf = zipfile.ZipFile(src, "r")
        zipf.extractall(dst)
        zipf.close()
    return 1


def unzip_changed_files(zipname, pathname, check_time=True, check_size=False):
    # Extract each file in a zip-file only when the properties are different
    # With the default arguments this method only checks the modification time
    with zipfile.ZipFile(zipname) as zf:
        infolist = zf.infolist()
        for info in infolist:
            fname = os.path.join(pathname, info.filename)
            extract = False
            if os.path.exists(fname):
                if check_time:
                    tz = time.mktime(info.date_time + (0, 0, -1))
                    tf = os.path.getmtime(fname)
                    if tz != tf:
                        extract = True
                if check_size:
                    sz = info.file_size
                    sf = os.path.getsize(fname)
                    if sz != sf:
                        extract = True
            else:
                extract = True
            if extract:
                logging.info("extracting {}".format(info.filename))
                zf.extract(info.filename, pathname)
                # set the correct modification time
                # (which is the time of extraction by default)
                tz = time.mktime(info.date_time + (0, 0, -1))
                os.utime(os.path.join(pathname, info.filename), (tz, tz))


def matlab2datetime(tindex):
    """Transform a matlab time to a datetime, rounded to seconds."""
    day = Timestamp.fromordinal(int(tindex))
    dayfrac = Timedelta(days=float(tindex) % 1) - Timedelta(days=366)
    return day + dayfrac


def get_files(
    file_or_dir, ext, unpackdir=None, force_unpack=False, preserve_datetime=False
):
    """internal method to get list of files with specific extension from
    dirname.

    Parameters
    ----------
    file_or_dir : str
        file or path to data
    ext : str
        extension of filenames to store in list
    force_unpack : bool, optional
        force unzip, by default False
    preserve_datetime : bool, optional
        preserve datetime of unzipped files, by default False
        (useful for checking whether data has changed)
    """
    # check if unpackdir is same as file_or_dir, if same, this can cause
    # problems when the unpackdir still contains zips that will be unpacked
    # again.
    if unpackdir is not None:
        if os.path.normcase(unpackdir) == os.path.normcase(file_or_dir):
            raise ValueError("Please specify a different folder to unpack" " files!")
    # identify whether file_or_dir started as zip
    if file_or_dir.endswith(".zip"):
        iszip = True
    else:
        iszip = False

    # unzip dir
    if file_or_dir.endswith(".zip"):
        zipf = file_or_dir
        if unpackdir is None:
            file_or_dir = tempfile.TemporaryDirectory().name
        else:
            file_or_dir = unpackdir
        unzip_file(
            zipf, file_or_dir, force=force_unpack, preserve_datetime=preserve_datetime
        )

    # file_or_dir is directory
    if os.path.isdir(file_or_dir):
        # check for zips in dir
        zip_fnames = [i for i in os.listdir(file_or_dir) if i.endswith(".zip")]
        if len(zip_fnames) > 0:
            # unzip zips
            if unpackdir is None:
                dirname = tempfile.TemporaryDirectory().name
            else:
                dirname = unpackdir
            for zipf in zip_fnames:
                unzip_file(
                    os.path.join(file_or_dir, zipf),
                    dirname,
                    force=True,
                    preserve_datetime=preserve_datetime,
                )
                # remove intermediate zipfiles if initial file_or_dir was zip
                if iszip:
                    os.remove(os.path.join(file_or_dir, zipf))
        else:
            dirname = file_or_dir
        # get all files with extension ext
        unzip_fnames = [i for i in os.listdir(dirname) if i.endswith(ext)]
    elif os.path.isfile(file_or_dir):
        # file_or_dir is actually an xml
        unzip_fnames = [os.path.basename(file_or_dir)]  # get file name
        dirname = os.path.dirname(file_or_dir)  # get directory path
    else:
        raise NotImplementedError("Cannot parse 'file_or_dir': " f"{file_or_dir}!")

    return dirname, unzip_fnames


def df2gdf(df, xcol="x", ycol="y"):
    """Make a GeoDataFrame from a DataFrame, assuming the geometry are
    points."""
    from geopandas import GeoDataFrame
    from shapely.geometry import Point

    gdf = GeoDataFrame(
        df.copy(), geometry=[Point((s[xcol], s[ycol])) for i, s in df.iterrows()]
    )
    return gdf


def show_versions():
    """Method to print the version of dependencies."""
    from sys import version as os_version

    from matplotlib import __version__ as mpl_version
    from numpy import __version__ as np_version
    from pandas import __version__ as pd_version
    from scipy import __version__ as sc_version

    msg = (
        f"Python version    : {os_version}\n"
        f"Numpy version     : {np_version}\n"
        f"Scipy version     : {sc_version}\n"
        f"Pandas version    : {pd_version}\n"
        f"Matplotlib version: {mpl_version}"
    )

    return print(msg)


class ColoredFormatter(logging.Formatter):
    """Colored log formatter.

    Taken from
    https://gist.github.com/joshbode/58fac7ababc700f51e2a9ecdebe563ad
    """

    def __init__(
        self, *args, colors: Optional[Dict[str, str]] = None, **kwargs
    ) -> None:
        """Initialize the formatter with specified format strings."""

        super().__init__(*args, **kwargs)

        self.colors = colors if colors else {}

    def format(self, record) -> str:
        """Format the specified record as text."""

        record.color = self.colors.get(record.levelname, "")
        record.reset = Style.RESET_ALL

        return super().format(record)


def get_color_logger(level="INFO"):
    formatter = ColoredFormatter(
        "{color}{levelname}:{name}:{message}{reset}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        colors={
            "DEBUG": Fore.CYAN,
            "INFO": Fore.GREEN,
            "WARNING": Fore.YELLOW,
            "ERROR": Fore.RED,
            "CRITICAL": Fore.RED + Back.WHITE + Style.BRIGHT,
        },
    )

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.handlers[:] = []
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, level))

    logging.captureWarnings(True)
    return logger


def oc_to_df(oc, col: Optional[str] = None) -> DataFrame:
    obs_times = DatetimeIndex(unique([obs.index for obs in oc.obs]))
    df = DataFrame(index=obs_times, columns=oc.index, dtype="float").sort_index()
    for obs in oc.obs:
        if not obs.empty:
            if col is None:
                vals = obs.loc[:, obs._get_first_numeric_col_name()]
            else:
                vals = obs.loc[:, col]
            df.loc[obs.index, obs.name] = vals.values

    return df


def interpolate(
    xy: List[List[float]],
    obsdf: DataFrame,
    obsloc: DataFrame,
    kernel: str = "thin_plate_spline",
    kernel2: str = "linear",
    epsilon: Optional[int] = None,
) -> DataFrame:
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

    if (kernel == "thin_plate_spline") or (kernel == "cubic"):
        min_val = 3
    elif kernel == "quintic":
        min_val = 6
    else:
        min_val = len(obsdf.index)

    fill_df = (
        DataFrame(
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

        if len(val) >= min_val:
            kernel = kernel
        else:
            kernel = kernel2

        # create an scipy interpolator
        rbf = RBFInterpolator(
            coor.to_numpy(), val.to_numpy(), epsilon=epsilon, kernel=kernel
        )

        val_rbf = rbf(xy)
        fill_df.loc[idx] = val_rbf

    return fill_df
