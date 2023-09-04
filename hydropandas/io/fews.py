import datetime
import logging
import os
import xml.etree.ElementTree as etree
from io import StringIO
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .. import observation
from ..observation import Obs

logger = logging.getLogger(__name__)


def read_xml_fname(
    fname: str,
    ObsClass: Union[Obs, Dict[str, Obs]],
    translate_dic: Optional[Dict[str, str]] = None,
    low_memory: bool = True,
    locationIds: Optional[List[str]] = None,
    filterdict: Optional[Dict[str, List[str]]] = None,
    return_events: bool = True,
    keep_flags: Tuple[int] = (0, 1),
    return_df: bool = False,
    tags: Tuple[str] = ("series", "header", "event"),
    remove_nan: bool = False,
    **kwargs: dict,  # unused
):
    """Read an xml filename into a list of observations objects.

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass: Union[Obs, Dict[str, Obs]]
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'monitoring_well'}.
    low_memory : bool, optional
        whether to use xml-parsing method with lower memory footprint,
        default is True
    locationIds : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    filterdict : dict, optional
        dictionary with tag name to apply filter to as keys, and list of
        accepted names as dictionary values to keep in final result,
        i.e. {"locationId": ["B001", "B002"]}
    return_events : bool, optional
        return all event-information in a DataFrame per location, instead of
        just a Series (defaults to False). Overrules keep_flags kwarg.
    keep_flags : list of ints, optional
        keep the values with these flags (defaults to 0 and 1). Only used
        when return_events is False.
    tags : list of strings, optional
        Select the tags to be parsed. Defaults to series, header and event
    return_df : bool, optional
        return a DataFame with the data, instead of two lists (default is
        False)
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost, only used if low_memory=False

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects
    """
    if translate_dic is None:
        translate_dic = {"locationId": "monitoring_well"}

    if low_memory is True:
        obs_list = iterparse_pi_xml(
            fname,
            ObsClass,
            translate_dic=translate_dic,
            locationIds=locationIds,
            filterdict=filterdict,
            return_events=return_events,
            keep_flags=keep_flags,
            return_df=return_df,
            tags=tags,
        )
    else:
        tree = etree.parse(fname)
        root = tree.getroot()
        obs_list = read_xml_root(
            root,
            ObsClass,
            translate_dic=translate_dic,
            locationIds=locationIds,
            remove_nan=remove_nan,
        )

    return obs_list


def iterparse_pi_xml(
    fname: str,
    ObsClass: Union[Obs, Dict[str, Obs]],
    translate_dic: Optional[Dict[str, str]] = None,
    filterdict: Optional[Dict[str, List[str]]] = None,
    locationIds: Optional[List[str]] = None,
    return_events: bool = True,
    keep_flags: Tuple[int] = (0, 1),
    return_df: bool = False,
    tags: Tuple[str] = ("series", "header", "event"),
):
    """Read a FEWS XML-file with measurements, memory efficient.

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass: Union[Obs, Dict[str, Obs]],
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'monitoring_well'}.
    locationIds : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    filterdict : dict, optional
        dictionary with tag name to apply filter to as keys, and list of
        accepted names as dictionary values to keep in final result,
        i.e. {"locationId": ["B001", "B002"]}
    return_events : bool, optional
        return all event-information in a DataFrame per location, instead of
        just a Series (defaults to False). Overrules keep_flags kwarg.
    keep_flags : list of ints, optional
        keep the values with these flags (defaults to 0 and 1). Only used
        when return_events is False.
    tags : list of strings, optional
        Select the tags to be parsed. Defaults to series, header and event
    return_df : bool, optional
        return a DataFame with the data, instead of two lists (default is
        False)

    Returns
    -------
    df : pandas.DataFrame
        a DataFrame containing the metadata and the series if 'return_df'
        is True
    obs_list : list of pandas Series
        list of timeseries if 'return_df' is False
    """
    from lxml.etree import iterparse

    if translate_dic is None:
        translate_dic = {"locationId": "monitoring_well"}

    tags = ["{{http://www.wldelft.nl/fews/PI}}{}".format(tag) for tag in tags]

    context = iterparse(fname, tag=tags)
    # _, root = next(context)

    header_list = []
    obs_list = []

    keep_flags = [str(flag) for flag in keep_flags]

    for _, element in context:
        if element.tag.endswith("header"):
            header = {}
            for h_attr in element:
                tag = h_attr.tag.replace(
                    "{{{0}}}".format("http://www.wldelft.nl/fews/PI"), ""
                )

                if tag.startswith("locationId"):
                    logger.info(f"reading {h_attr.text}")

                # if specific locations are provided only read those
                if locationIds is not None and tag.startswith("locationId"):
                    loc = h_attr.text
                    if loc not in locationIds:
                        element.clear()
                        logger.info(f" ... skipping '{loc}', not in locationIds")
                        continue

                if filterdict is not None:
                    for k, v in filterdict.items():
                        if tag.startswith(k):
                            attr = h_attr.text
                            if attr not in v:
                                element.clear()
                                logger.info(
                                    f" ... skipping '{attr}' not "
                                    f"in accepted values for '{k}'"
                                )
                                continue

                if h_attr.text is not None:
                    header[tag] = h_attr.text
                elif len(h_attr.attrib) != 0:
                    header[tag] = {**h_attr.attrib}
                else:
                    header[tag] = None
            events = []

        elif element.tag.endswith("event"):
            # if specific locations are provided only read those
            if locationIds is not None:
                if loc not in locationIds:
                    element.clear()
                    continue

            if filterdict is not None:
                skip = False
                for (
                    k,
                    v,
                ) in filterdict.items():
                    if header.get(k, None) not in v:
                        skip = True
                if skip:
                    element.clear()
                    continue
            events.append({**element.attrib})

        elif element.tag.endswith("series"):
            # if specific locations are provided only read those
            if locationIds is not None:
                if loc not in locationIds:
                    element.clear()
                    continue

            if filterdict is not None:
                skip = False
                for (
                    k,
                    v,
                ) in filterdict.items():
                    if header.get(k, None) not in v:
                        skip = True
                if skip:
                    element.clear()
                    continue

            if len(events) == 0:
                if return_events:
                    ts = pd.DataFrame()
                else:
                    ts = pd.Series()
            else:
                df = pd.DataFrame(events)
                df.index = pd.to_datetime(
                    [d + " " + t for d, t in zip(df["date"], df["time"])],
                    errors="coerce",
                )
                df.drop(columns=["date", "time"], inplace=True)
                if return_events:
                    df["value"] = pd.to_numeric(df["value"], errors="coerce")
                    df["flag"] = pd.to_numeric(df["flag"])
                    ts = df
                else:
                    mask = df["flag"].isin(keep_flags)
                    ts = pd.to_numeric(df.loc[mask, "value"], errors="coerce")

            o, header = _obs_from_meta(ts, header, translate_dic, ObsClass)
            header_list.append(header)
            obs_list.append(o)

        # Free memory.
        element.clear()

    if return_df:
        for h, s in zip(header_list, obs_list):
            h["series"] = s
        return pd.DataFrame(header_list)
    else:
        return obs_list


def read_xmlstring(
    xmlstring: str,
    ObsClass: Union[Obs, Dict[str, Obs]],
    translate_dic: Optional[Dict[str, str]] = None,
    filterdict: Optional[Dict[str, List[str]]] = None,
    locationIds: Optional[List[str]] = None,
    low_memory: bool = True,
    remove_nan: bool = False,
):
    """Read xmlstring into an list of Obs objects. Xmlstrings are usually
    obtained using a fews api.

    Parameters
    ----------
    xmlstring : str
        xml string to be parsed. Typically from a fews api.
    ObsClass: Union[Obs, Dict[str, Obs]]
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'monitoring_well'}.
    locationIds : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    low_memory : bool, optional
        whether to use xml-parsing method with lower memory footprint,
        default is True
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost, only used if low_memory=False

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects
    """
    if translate_dic is None:
        translate_dic = {"locationId": "monitoring_well"}

    if low_memory:
        obs_list = iterparse_pi_xml(
            StringIO(xmlstring),
            ObsClass,
            translate_dic=translate_dic,
            filterdict=filterdict,
            locationIds=locationIds,
        )
    else:
        root = etree.fromstring(xmlstring)
        obs_list = read_xml_root(
            root,
            ObsClass,
            translate_dic=translate_dic,
            locationIds=locationIds,
            remove_nan=remove_nan,
        )

    return obs_list


def read_xml_root(
    root: etree.Element,
    ObsClass: Union[Obs, Dict[str, Obs]],
    translate_dic: Dict[str, str] = None,
    locationIds: List[str] = None,
    remove_nan: bool = False,
):
    """Read a FEWS XML-file with measurements, return list of ObsClass objects.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        root element of a fews xml
    ObsClass: Union[Obs, Dict[str, Obs]],
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'monitoring_well'}.
    locationIds : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects
    """
    if translate_dic is None:
        translate_dic = {"locationId": "monitoring_well"}

    obs_list = []
    for item in root:
        if item.tag.endswith("series"):
            header = {}
            date = []
            time = []
            events = []
            for subitem in item:
                if subitem.tag.endswith("header"):
                    for subsubitem in subitem:
                        prop = subsubitem.tag.split("}")[-1]
                        val = subsubitem.text
                        if prop == "x" or prop == "y" or prop == "lat" or prop == "lon":
                            val = float(val)
                        header[prop] = val
                        if prop == "locationId":
                            logger.info(f"read {val}")
                elif subitem.tag.endswith("event"):
                    date.append(subitem.attrib.pop("date"))
                    time.append(subitem.attrib.pop("time"))
                    events.append({**subitem.attrib})

            # combine events in a dataframe
            index = pd.to_datetime(
                [d + " " + t for d, t in zip(date, time)], errors="coerce"
            )
            ts = pd.DataFrame(events, index=index)
            ts.loc[:, "value"] = ts.loc[:, "value"].astype(float)

            if remove_nan and (not ts.empty):
                ts.dropna(subset=["value"], inplace=True)
                header["unit"] = "m NAP"

            o, header = _obs_from_meta(ts, header, translate_dic, ObsClass)
            if locationIds is not None:
                if header["monitoring_well"] in locationIds:
                    obs_list.append(o)
            else:
                obs_list.append(o)

    return obs_list


def _obs_from_meta(
    ts: pd.DataFrame,
    header: Dict[str, str],
    translate_dic: Dict[str, str],
    ObsClass: Union[Obs, Dict[str, Obs]],
):
    """Internal function to convert timeseries and header into Obs objects.

    Parameters
    ----------
    ts : pd.DataFrame
        timeseries data.
    header : dictionary
        metadata.
    translate_dic : dictionary
        translate dictionary.
    ObsClass: Union[Obs, Dict[str, Obs]],
        class of the observations, e.g. GroundwaterObs or WaterlvlObs

    Returns
    -------
    o : GroundwaterObs or WaterlvlObs
        hyrdopandas observation object.
    header : dictionary
        metadata.
    """
    for key, item in translate_dic.items():
        header[item] = header.pop(key)

    if "x" in header.keys():
        x = float(header["x"])
    else:
        x = np.nan
    if "y" in header.keys():
        y = float(header["y"])
    else:
        y = np.nan

    if "units" in header.keys():
        unit = str(header["units"])
    else:
        unit = np.nan

    if np.isnan(x) or np.isnan(y):
        metadata_available = False
    else:
        metadata_available = True

    if "parameterId" in header:
        parid = header["parameterId"]
        name = header["monitoring_well"] + "_" + parid
    else:
        name = header["monitoring_well"]

    if isinstance(ObsClass, dict):
        if parid in ObsClass.keys():
            ObsC = ObsClass[parid]
        else:
            ObsC = Obs
    else:
        ObsC = ObsClass

    if ObsC in (observation.WaterlvlObs,):
        o = ObsC(
            ts,
            x=x,
            y=y,
            unit=unit,
            meta=header,
            name=name,
            monitoring_well=header["monitoring_well"],
            metadata_available=metadata_available,
            source="FEWS",
        )
    elif ObsC in (observation.GroundwaterObs,):
        if "z" in header.keys():
            z = float(header["z"])
        else:
            z = np.nan
        o = ObsC(
            ts,
            x=x,
            y=y,
            ground_level=z,
            unit=unit,
            meta=header,
            name=name,
            monitoring_well=header["monitoring_well"],
            metadata_available=metadata_available,
            source="FEWS",
        )
    elif ObsC in (
        observation.MeteoObs,
        observation.PrecipitationObs,
        observation.EvaporationObs,
    ):
        o = ObsC(
            ts,
            x=x,
            y=y,
            unit=unit,
            meta=header,
            name=name,
            meteo_var=parid,
            source="FEWS",
        )
    else:
        o = ObsC(ts, x=x, y=y, unit=unit, meta=header, name=name, source="FEWS")

    return o, header


def write_pi_xml(obs_coll, fname: str, timezone: float = 1.0, version: str = "1.24"):
    """Write TimeSeries object to PI-XML file.

    Parameters
    ----------
    fname: path
        path to XML file
    """

    assert fname.endswith(".xml"), "Output file should have '.xml' extension!"

    # first line of XML file
    line0 = '<?xml version="1.0" encoding="UTF-8"?>\n'

    # some definitions for timeseries XML file
    NS = r"http://www.wldelft.nl/fews/PI"
    FS = r"http://www.wldelft.nl/fews/fs"
    XSI = r"http://www.w3.org/2001/XMLSchema-instance"
    schemaLocation = (
        r"http://fews.wldelft.nl/schemas/version1.0" r"/Pi-schemas/pi_timeseries.xsd"
    )
    timeseriesline = (
        '<TimeSeries xmlns="{NS}" xmlns:xsi="{XSI}" '
        'xsi:schemaLocation="{NS} {schema}" version="{version}" '
        'xmlns:fs="{FS}">\n'
    )

    # line templates
    paramline = "<{tag}>{param}</{tag}>\n"

    # write file
    with open(fname, "w") as f:
        f.write(line0)
        f.write(
            timeseriesline.format(
                NS=NS, FS=FS, XSI=XSI, schema=schemaLocation, version=version
            )
        )
        tzline = "\t" + paramline.format(tag="timeZone", param=timezone)
        f.write(tzline)

        for o in obs_coll.obs:
            # start series
            start = "\t" + "<series>\n"
            f.write(start)
            # write header
            hlines = []
            hstart = 2 * "\t" + "<header>\n"
            hlines.append(hstart)
            for htag, hval in o.meta.items():
                if htag.endswith("Date"):
                    try:
                        hdate = hval.strftime("%Y-%m-%d")
                        htime = hval.strftime("%H:%M:%S")
                    except AttributeError as e:
                        if htag.startswith("start"):
                            hdate = o.index[0].strftime("%Y-%m-%d")
                            htime = o.index[0].strftime("%H:%M:%S")
                        elif htag.startswith("end"):
                            hdate = o.index[-1].strftime("%Y-%m-%d")
                            htime = o.index[-1].strftime("%H:%M:%S")
                        else:
                            raise (e)
                    hline = '<{tag} date="{date}" time="{time}"/>\n'.format(
                        tag=htag, date=hdate, time=htime
                    )
                elif htag.endswith("timeStep"):
                    hline = '<{tag} unit="{unit}"/>\n'.format(tag=htag, unit=hval)
                else:
                    hline = paramline.format(tag=htag, param=hval)
                hlines.append(3 * "\t" + hline)
            hlines.append(2 * "\t" + "</header>\n")
            f.writelines(hlines)

            # write timeseries
            dates = o.reset_index()["index"].apply(
                lambda s: datetime.datetime.strftime(s, "%Y-%m-%d")
            )
            times = o.reset_index()["index"].apply(
                lambda s: datetime.datetime.strftime(s, "%H:%M:%S")
            )
            # set date and time attributes
            events = (
                2 * "\t" + '<event date="' + dates.values + '" time="' + times.values
            )
            # loop through columns and add to event
            for icol in o.columns:
                val = o[icol].astype(str)
                events += '" {}="'.format(icol) + val.values
            # close event
            events += '"/>\n'
            # write to file
            f.writelines(events)
            # end series
            f.write("\t" + "</series>\n")
        # end Timeseries
        f.write("</TimeSeries>\n")


def read_xml_filelist(
    fnames: List[str],
    ObsClass: Union[Obs, Dict[str, Obs]],
    directory: Optional[str] = None,
    locations: Optional[List[str]] = None,
    translate_dic: Optional[Dict[str, str]] = None,
    filterdict: Optional[Dict[str, List[str]]] = None,
    remove_nan: bool = False,
    low_memory: bool = True,
    **kwargs: dict,
):
    """Read a list of xml files into a list of observation objects.

    Parameters
    ----------
    fnames : TYPE
        DESCRIPTION.
    ObsClass: Union[Obs, Dict[str, Obs]]
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    directory : TYPE, optional
        DESCRIPTION. The default is None.
    locations : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'monitoring_well'}.
    filterdict : dict, optional
        dictionary with tag name to apply filter to as keys, and list of
        accepted names as dictionary values to keep in final result,
        i.e. {"locationId": ["B001", "B002"]}
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost, only used if low_memory=False
    low_memory : bool, optional
        whether to use xml-parsing method with lower memory footprint,
        default is True

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects
    """
    if translate_dic is None:
        translate_dic = {"locationId": "monitoring_well"}

    obs_list = []
    nfiles = len(fnames)
    for j, ixml in enumerate(fnames):
        # print message
        logger.info(f"{j+1}/{nfiles} read {ixml}")

        # join directory to filename if provided
        if directory is None:
            fullpath = ixml
        else:
            fullpath = os.path.join(directory, ixml)

        # read xml fname
        obs_list += read_xml_fname(
            fullpath,
            ObsClass,
            translate_dic=translate_dic,
            filterdict=filterdict,
            low_memory=low_memory,
            locationIds=locations,
            remove_nan=remove_nan,
            **kwargs,
        )

    return obs_list


def get_fews_pid(name: str) -> Dict[str, Obs]:
    """Get matching ParameterId's and HydroPandas Observation Classes

    Parameters
    ----------
    name : str
        Waterboard name

    Returns
    -------
    Dict[str, Obs]
        Dictonary with ParameterId and the resulting Observation Class
    """
    from ..data.fews_parameterid import pid

    pid_sel = pid[name.lower()]

    return {key: getattr(observation, value) for (key, value) in pid_sel.items()}
