# -*- coding: utf-8 -*-
"""Created on Thu Oct 10 11:01:22 2019.

@author: oebbe
"""

import logging
import os

import numpy as np
from pandas import DataFrame, Series
from scipy.io import loadmat

from ..observation import GroundwaterObs, WaterlvlObs
from ..util import matlab2datetime

logger = logging.getLogger(__name__)


def read_file(fname, ObsClass, load_oseries=True, load_stresses=True):
    """This method is used to read the file."""

    logger.info(f"reading menyanthes file {fname}")

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
    if not (os.path.isfile(fname)):
        print("Could not find file ", fname)

    mat = loadmat(fname, struct_as_record=False, squeeze_me=True, chars_as_strings=True)

    obs_list = []
    if load_oseries:
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

            o = ObsClass(df, meta=metadata, **meta_o)
            obs_list.append(o)

    if load_stresses:
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
            )
            obs_list.append(o)

    return obs_list


def read_oseries(mat):
    """Read the oseries from a mat file from menyanthes."""
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
            if name != "values":
                data[name.lower()] = getattr(H, name)
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
    d_in = {}

    # Check if more then one time series is present
    # if not isinstance(mat["IN"], np.ndarray):
    #     mat["IN"] = [mat["IN"]]

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
            if name != "values":
                data[name.lower()] = getattr(IN, name)
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
