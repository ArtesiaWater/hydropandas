import logging
import os

import numpy as np
from pandas import DataFrame, Series, Timedelta, Timestamp
from scipy.io import loadmat

from ..observation import GroundwaterObs, WaterlvlObs

logger = logging.getLogger(__name__)


def matlab2datetime(tindex):
    """
    Transform a MATLAB serial date number to a Python datetime object, rounded
    to seconds.

    Parameters
    ----------
    tindex : float
        The MATLAB serial date number to convert.

    Returns
    -------
    datetime : datetime.datetime
        The equivalent datetime object in Python.

    Notes
    -----
    MATLAB serial date numbers represent the number of days elapsed since
    January 1, 0000 (the proleptic Gregorian calendar), with January 1, 0000 as
    day 1. Fractions of a day can be represented as a decimal.

    The returned datetime object is rounded to the nearest second.

    Examples
    --------
    >>> matlab2datetime(719529.496527778)
    datetime.datetime(2019, 1, 1, 11, 55, 2)

    """
    day = Timestamp.fromordinal(int(tindex))
    dayfrac = Timedelta(days=float(tindex) % 1) - Timedelta(days=366)
    return day + dayfrac


def read_file(path, ObsClass, load_oseries=True, load_stresses=True):
    """
    Read data from a Menyanthes file and create observation objects.

    Parameters
    ----------
    path : str
        Full path of the Menyanthes file (.men) to read.
    ObsClass : GroundwaterObs or WaterlvlObs
        Class of observation object to create.
    load_oseries : bool, optional
        Flag indicating whether to load observation series or not, by default
        True.
    load_stresses : bool, optional
        Flag indicating whether to load stresses or not, by default True.

    Returns
    -------
    obs_list : list
        List of observation objects created from the Menyanthes file.
    """

    logger.info(f"reading menyanthes file {path}")

    if ObsClass == GroundwaterObs:
        _rename_dic = {
            "xcoord": "x",
            "ycoord": "y",
            "upfiltlev": "screen_top",
            "lowfiltlev": "screen_bottom",
            "surflev": "ground_level",
            "filtnr": "tube_nr",
            "measpointlev": "tube_top",
        }

        _keys_o = [
            "name",
            "x",
            "y",
            "source",
            "unit",
            "monitoring_well",
            "tube_nr",
            "metadata_available",
            "ground_level",
            "tube_top",
            "screen_top",
            "screen_bottom",
        ]
        unit = "m NAP"
    elif ObsClass == WaterlvlObs:
        _rename_dic = {"xcoord": "x", "ycoord": "y", "measpointlev": "tube_top"}
        _keys_o = ["name", "x", "y", "source", "unit", "monitoring_well"]
        unit = "m NAP"
    else:
        _rename_dic = {
            "xcoord": "x",
            "ycoord": "y",
        }
        _keys_o = ["name", "x", "y", "source", "unit"]
        unit = ""

    # Check if file is present
    if not (os.path.isfile(path)):
        print("Could not find file ", path)

    mat = loadmat(path, struct_as_record=False, squeeze_me=True, chars_as_strings=True)

    obs_list = []
    if load_oseries and ("H" in mat.keys()):
        d_h = read_oseries(mat)

        locations = d_h.keys()
        for location in locations:
            metadata = d_h[location]
            metadata["projection"] = "epsg:28992"
            metadata["metadata_available"] = True
            metadata["source"] = "Menyanthes"
            metadata["unit"] = unit

            df = DataFrame(metadata.pop("values"), columns=["values"])
            for key in _rename_dic.keys():
                if key in metadata.keys():
                    metadata[_rename_dic[key]] = metadata.pop(key)

            meta_o = {k: metadata[k] for k in _keys_o if k in metadata}

            o = ObsClass(df, meta=metadata, **meta_o, filename=path)
            obs_list.append(o)

    if load_stresses and ("IN" in mat.keys()):
        d_in = read_stresses(mat)
        stresses = d_in.keys()
        for stress in stresses:
            metadata = d_in[stress]
            metadata["projection"] = "epsg:28992"
            metadata["metadata_available"] = True
            metadata["source"] = "Menyanthes"
            metadata["unit"] = unit
            s = metadata.pop("values")
            df = DataFrame(s, columns=["values"])
            for key in _rename_dic.keys():
                if key in metadata.keys():
                    metadata[_rename_dic[key]] = metadata.pop(key)
            o = ObsClass(
                df,
                meta=metadata,
                name=metadata["name"],
                x=metadata["x"],
                y=metadata["y"],
                source=metadata["source"],
                unit=metadata["unit"],
                filename=path,
            )
            obs_list.append(o)

    return obs_list


def read_oseries(mat):
    """Read the oseries from a mat file from menyanthes.

    Parameters
    ----------
    mat : dict
        A dictionary object containing the Menyanthes file data.

    Returns
    -------
    dict
        A dictionary containing oseries data, with oseries names as keys and
        their corresponding metadata and values as values.

    Notes
    -----
    This function reads the oseries data from a Menyanthes file in .mat format
    and returns it in a dictionary format. The oseries data contains the
    following metadata:
        - name: The name of the oseries.
        - x: The x-coordinate of the oseries location.
        - y: The y-coordinate of the oseries location.
        - source: The data source.
        - unit: The unit of measurement.

    In addition to the metadata, the oseries data also contains a pandas Series
    object named 'values', which contains the time series data for the oseries.

    Examples
    --------
    >>> mat = loadmat('menyanthes_file.mat')
    >>> d_h = read_oseries(mat)
    """
    d_h = {}

    # Check if more then one time series model is present
    if not isinstance(mat["H"], np.ndarray):
        mat["H"] = [mat["H"]]

    # Read all the time series models
    for i, H in enumerate(mat["H"]):
        if not hasattr(H, "Name") and not hasattr(H, "name"):
            H.Name = "H" + str(i)  # Give it the index name
        if hasattr(H, "name"):
            H.Name = H.name
        if len(H.Name) == 0:
            H.Name = H.tnocode
        logger.info(f"reading oseries -> {H.Name}")

        data = {}

        for name in H._fieldnames:
            val = getattr(H, name)
            if name != "values":
                # if value is an empty numpy array set value to nan
                if isinstance(val, np.ndarray) and val.size == 0:
                    val = np.nan
                data[name.lower()] = val
            else:
                if H.values.size == 0:
                    # when diver-files are used, values will be empty
                    series = Series()
                else:
                    tindex = map(matlab2datetime, H.values[:, 0])
                    # measurement is used as is
                    series = Series(H.values[:, 1], index=tindex)
                    # round on seconds, to get rid of conversion milliseconds
                    series.index = series.index.round("s")
                data["values"] = series

        # add to self.H
        d_h[H.Name] = data

    return d_h


def read_stresses(mat):
    """Reads the stresses from a mat file from menyanthes.

    Parameters
    ----------
    mat : dict
        A dictionary object containing the mat file.

    Returns
    -------
    dict
        A dictionary object containing the stresses data.
    """
    d_in = {}

    # Read all the time series
    for i, IN in enumerate(mat["IN"]):
        if not hasattr(IN, "Name") and not hasattr(IN, "name"):
            IN.Name = "IN" + str(i)  # Give it the index name
        if hasattr(IN, "name"):
            IN.Name = IN.name
        if len(IN.Name) == 0:
            IN.Name = IN.tnocode

        logger.info(f"reading stress -> {IN.Name}")

        data = {}
        for name in IN._fieldnames:
            val = getattr(IN, name)
            if name != "values":
                # if value is an empty numpy array set value to nan
                if isinstance(val, np.ndarray) and val.size == 0:
                    val = np.nan
                data[name.lower()] = val
            else:
                if IN.values.size == 0:
                    # when diver-files are used, values will be empty
                    series = Series()
                else:
                    tindex = map(matlab2datetime, IN.values[:, 0])
                    # measurement is used as is
                    series = Series(IN.values[:, 1], index=tindex)
                    # round on seconds, to get rid of conversion milliseconds
                    series.index = series.index.round("s")
                data["values"] = series

        # add to self.H
        d_in[IN.Name] = data

    return d_in
