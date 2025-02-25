import logging
import os

import numpy as np
import pandas as pd

from ..util import get_files

logger = logging.getLogger(__name__)


def _read_wiski_header(f, header_sep=":", header_identifier="#", end_header_str=None):
    line = f.readline()
    header = dict()
    while header_identifier in line:
        prop, val = line.split(header_sep)
        prop = prop.strip().replace(header_identifier, "")
        val = val.strip()
        try:
            val = float(val)
        except ValueError:
            pass
        header[prop] = val
        line = f.readline()
        if end_header_str is not None:
            if end_header_str in line:
                break

    return line, header


def read_wiski_file(
    path,
    sep=";",
    header_sep=None,
    header_identifier="#",
    read_series=True,
    translate_dic=None,
    tz_localize=True,
    unit="",
    **kwargs,
):
    """
    Read data from a WISKI file.

    Parameters:
    -----------
    path : str
        The path of the file to be read.
    sep : str, optional (default=";")
        The delimiter used to separate fields in the file.
    header_sep : str, optional (default=None)
        The delimiter used to separate fields in the header. If None, the
        function will try to automatically detect the separator.
    header_identifier : str, optional (default="#")
        The character used to identify header lines.
    read_series : bool, optional (default=True)
        Whether to read the time series data from the file.
    translate_dic : dict, optional (default=None)
        A dictionary mapping header field names to the desired output names.
    tz_localize : bool, optional (default=True)
        Whether to localize the datetime index to the machine's timezone.
    unit : str, optional (default="")
        The unit of measurement of the data.
    **kwargs : keyword arguments
        Additional arguments to pass to the pandas `read_csv` function.

    Returns:
    --------
    data : pandas.DataFrame or None
        A dataframe containing the time series data from the file. Returns None
        if `read_series` is False.
    metadata : dict
        A dictionary containing metadata about the data in the file.
    """
    logger.info("reading -> {}".format(os.path.split(path)[-1]))

    if translate_dic is None:
        translate_dic = {}

    # manually break header parse at certain point
    if "end_header_str" in kwargs.keys():
        end_header_str = kwargs.pop("end_header_str")
    else:
        end_header_str = None

    # read header
    with open(path, "r") as f:
        if header_sep is None:
            line, header = _read_wiski_header(
                f, end_header_str=end_header_str, header_identifier=header_identifier
            )
        else:
            line, header = _read_wiski_header(
                f,
                header_sep=header_sep,
                header_identifier=header_identifier,
                end_header_str=end_header_str,
            )

        # specify names
        for key, item in translate_dic.items():
            header[key] = header[item]

        # get column names of data
        # columns = split(r'\s{2,}', line)
        if sep == r"\s+":
            columns = line.split("\t")
        else:
            columns = line.split(sep)
        columns = [icol.strip("\n") for icol in columns]
        columns = [icol.replace("[", "_[") for icol in columns]

        validator = np.lib._iotools.NameValidator()
        columns = validator(columns)

        if read_series:
            # read data
            data = pd.read_csv(
                f,
                sep=sep,
                header=None,
                names=columns,
                **kwargs,
            )

            if tz_localize:
                data.index = data.index.tz_localize(None)

            # convert Value to float
            col = [icol for icol in data.columns if icol.lower().startswith("value")][0]

            data[col] = pd.to_numeric(data[col], errors="coerce")
        else:
            data = None

        # translate some header keys
        metadata = {"source": "wiski", "unit": unit}
        for key, val in header.items():
            if key == "Station Site":
                metadata["location"] = val
            elif key == "x":
                metadata["x"] = val
            elif key == "y":
                metadata["y"] = val
            elif key == "name":
                metadata["name"] = val
            # this adds the other header keys to metadata
            # but this will not work if keys contain spaces etc.
            # because it tries adding them as attributes.
            # else:
            #     metadata[key] = val

    return data, metadata


def read_wiski_dir(
    dirname,
    ObsClass=None,
    suffix=".csv",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    **kwargs,
):
    """
    Reads WISKI CSV files from a directory and returns a list of observation
    objects.

    Parameters
    ----------
    dirname : str
        The path of the directory containing the WISKI CSV files.
    ObsClass : object, optional
        The observation class to use for creating observation objects. Default
        is None.
    suffix : str, optional
        The file extension of the WISKI CSV files. Default is ".csv".
    unpackdir : str, optional
        The directory to which the files should be unpacked. Default is None.
    force_unpack : bool, optional
        If True, forces the files to be unpacked even if they are already in the
        target directory. Default is False.
    preserve_datetime : bool, optional
        If True, preserves the original modification times of the files when
        unpacking them. Default is False.
    keep_all_obs : bool, optional
        If True, keeps all observation objects even if they have no metadata
        available. Default is True.
    **kwargs
        Additional keyword arguments to pass to the `from_wiski` method of the
        `ObsClass` object.

    Returns
    -------
    list
        A list of observation objects created from the WISKI CSV files in the
        directory.

    Raises
    ------
    FileNotFoundError
        If no WISKI CSV files are found in the directory.

    """
    # get files
    dirname, unzip_fnames = get_files(
        dirname,
        ext=suffix,
        force_unpack=force_unpack,
        preserve_datetime=preserve_datetime,
    )

    if not unzip_fnames:
        raise FileNotFoundError(
            "no files were found in '{}' that end with '{}'".format(
                os.path.join(dirname), suffix
            )
        )

    # gather all obs in list
    obs_list = []
    for i, csv in enumerate(unzip_fnames):
        logger.info("reading {0}/{1} -> {2}".format(i + 1, len(unzip_fnames), csv))
        obs = ObsClass.from_wiski(os.path.join(dirname, csv), **kwargs)

        if obs.metadata_available:
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            logger.info(f"not added to collection -> {csv}")

    return obs_list
