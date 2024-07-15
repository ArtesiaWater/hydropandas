"""
https://www.bro-productomgeving.nl/bpo/latest/informatie-voor-softwareleveranciers/url-s-publieke-rest-services
"""

import json
import logging
import os
import xml.etree.ElementTree
from functools import lru_cache

import numpy as np
import pandas as pd
import requests
import tempfile
from pyproj import Proj, Transformer
from tqdm import tqdm

from ..rcparams import rcParams
from ..util import EPSG_28992

logger = logging.getLogger(__name__)

ns_gld = {
        "ns11": "http://www.broservices.nl/xsd/dsgld/1.0",
        "gldcommon": "http://www.broservices.nl/xsd/gldcommon/1.0",
        "waterml": "http://www.opengis.net/waterml/2.0",
        "swe": "http://www.opengis.net/swe/2.0",
        "om": "http://www.opengis.net/om/2.0",
    }

ns_gmw = {
        "dsgmw": "http://www.broservices.nl/xsd/dsgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "gml": "http://www.opengis.net/gml/3.2",
        "brocom":  "http://www.broservices.nl/xsd/brocommon/3.0"
    }

ns_gmn = {"xmlns": "http://www.broservices.nl/xsd/dsgmn/1.0"}


def get_obs_list_from_gmn(bro_id_gmn, ObsClass, only_metadata=False, keep_all_obs=True):
    """get a list of observation from a groundwater monitoring network (GMN).

    Parameters
    ----------
    bro_id_gmn : str
        starts with 'GMN' e.g. 'GMN000000000163'.
    ObsClass : type
        class of the observations, so far only GroundwaterObs is supported
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without
        measurements

    Returns
    -------
    obs_list : list
        list with observation objects.
    meta : dict
        metadata of the groundwater monitoring net.

    """
    gmn = get_gmn_elem_api(bro_id_gmn)
    return get_obs_list_from_gmn_api(gmn, bro_id_gmn, ObsClass, only_metadata=only_metadata, 
                                     keep_all_obs=keep_all_obs)


def get_gmn_elem_api(bro_id_gmn):
    """get a list of observation from a groundwater monitoring network.

    Parameters
    ----------
    bro_id_gmn : str
        starts with 'GMN' e.g. 'GMN000000000163'.

    Raises
    ------
    ValueError
        if bro_id_gmn is invalid.

    Returns
    -------
    gmn : Element
        xml reference to gmn object
    """

    if not bro_id_gmn.startswith("GMN"):
        raise ValueError("bro id should start with GMN")

    url = f"https://publiek.broservices.nl/gm/gmn/v1/objects/{bro_id_gmn}"
    req = requests.get(url)

    if req.status_code > 200:
        print(req.json()["errors"][0]["message"])

    tree = xml.etree.ElementTree.fromstring(req.text)
    gmn = tree.find(".//xmlns:GMN_PO", ns_gmn)

    return gmn

def get_obs_list_from_gmn_api(gmn, bro_id_gmn, ObsClass, only_metadata=False, keep_all_obs=True):
    """get a list of observations from a groundwater monitoring network (GMN) xml.

    Parameters
    ----------
    gmn : Element
        xml reference to gmn object
    bro_id_gmn : str
        starts with 'GMN' e.g. 'GMN000000000163'.
    ObsClass : type
        class of the observations, so far only GroundwaterObs is supported
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without
        measurements

    Returns
    -------
    obs_list : list
        list with observation objects.
    meta : dict
        metadata of the groundwater monitoring net.
    """
    gmws = gmn.findall("xmlns:measuringPoint", ns_gmn)

    logger.info(f"{len(gmws)} groundwater monitoring wells within groundwater meetnet")

    obs_list = []
    for gmw in tqdm(gmws):
        tags = ["MeasuringPoint", "monitoringTube", "GroundwaterMonitoringTube"]
        tube = gmw.find(f".//xmlns:{'//xmlns:'.join(tags)}", ns_gmn)
        bro_id_gmw = tube.find("xmlns:broId", ns_gmn).text
        tube_nr = int(tube.find("xmlns:tubeNumber", ns_gmn).text)

        o = ObsClass.from_bro(
            bro_id=bro_id_gmw, tube_nr=tube_nr, only_metadata=only_metadata
        )
        if o.empty:
            logger.debug(
                f"no measurements found for gmw_id {bro_id_gmw} and tube number {tube_nr}"
            )
            if keep_all_obs:
                obs_list.append(o)
        else:
            obs_list.append(o)
        obs_list.append(o)

    meta = {}
    meta["name"] = gmn.find("xmlns:name", ns_gmn).text
    meta["doel"] = gmn.find("xmlns:monitoringPurpose", ns_gmn).text
    meta["bro_id"] = bro_id_gmn

    return obs_list, meta


def get_bro_groundwater(bro_id=None, fname=None, tube_nr=None, 
                        only_metadata=False, tmin=None, tmax=None, **kwargs):
    """get bro groundwater measurements and/or metadata.

    Parameters
    ----------
    bro_id : str, optional
        starts with 'GLD' or 'GMW' e.g. 'GLD000000012893'.
    fname : str, optional
        fname xml file e.g. 'GLD000000012893.xml'
    tube_nr : int or None, optional
        tube number, required if bro_id starts with 'GMW'. The default is
        None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    tmin : str or None, optional
        start date in format YYYY-MM-DD
    tmax : str or None, optional
        end date in format YYYY-MM-DD
    **kwargs :
        passes to measurements_from_gld_elem.

    Returns
    -------
    dataframe
        measurements.
    dictionary
        metadata.

    Raises
    ------
    ValueError
        if fname and bro_id are not specified
    ValueError
        if both fname and bro_id are specified
    """
    if fname is None and bro_id is None:
        raise ValueError('please specify either fname or bro_id')
    elif fname is not None and bro_id is not None:
        raise ValueError('please specify fname or bro_id not both')
    elif bro_id is not None:
        return get_bro_groundwater_api(bro_id, tube_nr=tube_nr, only_metadata=only_metadata, 
                                       tmin=tmin, tmax=tmax, **kwargs)
    elif fname is not None:
        logger.warning('reading data from a single xml file, this either gives only'
                       'the measurements (for GLD files) or only the metadata (for '
                       'GMW files)')
        return get_bro_groundwater_fname(fname, tube_nr=tube_nr, only_metadata=only_metadata, 
                                         tmin=tmin, tmax=tmax, **kwargs)


def get_bro_groundwater_api(bro_id, tube_nr=None, only_metadata=False, 
                            tmin=None, tmax=None, **kwargs):
    """get bro groundwater measurement from a GLD id or a GMW id with a
    filter number.

    Parameters
    ----------
    bro_id : str, optional
        starts with 'GLD' or 'GMW' e.g. 'GLD000000012893'.
    tube_nr : int or None, optional
        tube number, required if bro_id starts with 'GMW'. The default is
        None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. Should be False if
        bro_id is a GLD id. The default is False.
    tmin : str or None, optional
        start date in format YYYY-MM-DD
    tmax : str or None, optional
        end date in format YYYY-MM-DD
    **kwargs :
        passes to measurements_from_gld_elem.

    Raises
    ------
    ValueError
        if bro_id is a GLD id and only_metadata is True

    Returns
    -------
    dataframe
        measurements.
    dictionary
        metadata.

    """
    logger.debug(f"reading bro_id {bro_id}")

    if bro_id.startswith("GLD"):
        if only_metadata:
            raise ValueError("cannot get metadata from gld id")
        gld = get_gld_elem_api(bro_id, tmin=tmin, tmax=tmax)
        df, meta = measurements_from_gld_elem(gld, bro_id, tmin=tmin, tmax=tmax, **kwargs)
        
        # add metadata from gmw
        gmw = get_gmw_elem_api(meta["monitoring_well"])
        meta.update(get_metadata_from_gmw_elem(gmw, meta["monitoring_well"], meta["tube_nr"]))
        return df, meta
    elif bro_id.startswith("GMW"):
        if tube_nr is None:
            raise ValueError("if bro_id is GMW a tube_nr should be specified")

        gmw = get_gmw_elem_api(bro_id)
        meta = get_metadata_from_gmw_elem(gmw, bro_id, tube_nr)
        bro_id_glds = get_gld_ids_from_gmw(bro_id, tube_nr)

        if bro_id_glds is None:
            meta["name"] = f"{bro_id}_{tube_nr}"
            only_metadata = True  # cannot get time series without gld id
        else:
            meta["name"] = f"{bro_id}_{tube_nr}"
            meta["gld_ids"] = bro_id_glds

        if only_metadata:
            empty_df = pd.DataFrame()
            return empty_df, meta

        dfl = []
        for bro_id_gld in bro_id_glds:
            gld = get_gld_elem_api(bro_id_gld, tmin=tmin, tmax=tmax)
            df, meta_new = measurements_from_gld_elem(gld, bro_id_gld, 
                                                      tmin=tmin,
                                                        tmax=tmax, **kwargs)
            meta.update(meta_new)
            dfl.append(df)
        df = pd.concat(dfl, axis=0)

        return df, meta

def get_bro_groundwater_fname(fname, tube_nr=None, only_metadata=False, 
                              tmin=None, tmax=None, **kwargs):
    """get bro groundwater measurement from an xml GLD or GMW file.

    Parameters
    ----------
    fname : str, optional
        fname xml file e.g. 'GLD000000012893.xml'
    tube_nr : int or None, optional
        tube number, required if file is an gmw file. The default is
        None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    tmin : str or None, optional
        start date in format YYYY-MM-DD
    tmax : str or None, optional
        end date in format YYYY-MM-DD
    **kwargs :
        passes to measurements_from_gld_elem.

    Raises
    ------
    ValueError
        for invalid filename.

    Returns
    -------
    dataframe
        measurements.
    dictionary
        metadata.

    """
    bro_id = os.path.basename(fname).split(".")[0]
    if bro_id.startswith("GLD"):
        logger.info(f'reading GLD data from file {fname}, do not read metadata')
        gld = get_gld_elem_file(fname)
        df, meta = measurements_from_gld_elem(gld, bro_id, tmin=tmin, tmax=tmax, **kwargs)
        meta.update({'x':np.nan, 'y':np.nan, 'unit': '','screen_bottom':np.nan, 
                     'screen_top':np.nan, 
                     'ground_level':np.nan, 'metadata_available':False, 
                     'tube_top':np.nan})
    elif bro_id.startswith("GMW"):
        if tube_nr is None:
            raise ValueError("if fname is GMW file a tube_nr should be specified")
        logger.info(f'reading GMW data from file {fname}, do not read measurements')
        
        gmw = get_gmw_elem_file(fname)
        meta = get_metadata_from_gmw_elem(gmw, bro_id, tube_nr)

        meta["name"] = f"{bro_id}_{tube_nr}"
        df = pd.DataFrame()
    meta['filename'] = fname
        
    return df, meta


def get_gld_ids_from_gmw(bro_id_gmw, tube_nr):
    """get bro_ids of multiple grondwaterstandendossier (gld) from a bro_id of a
    grondwatermonitoringsput (gmw).

    Parameters
    ----------
    bro_id_gmw : str
        starts with 'GMW' e.g. 'GMW000000036287'.
    tube_nr : int
        tube number.

    Raises
    ------
    ValueError
        if bro_id_gld is invalid.

    Returns
    -------
    list of str or None
        bro_ids of a grondwaterstandonderzoek (gld).

    """
    if not bro_id_gmw.startswith("GMW"):
        raise ValueError("bro id should start with GMW")

    url = f"https://publiek.broservices.nl/gm/v1/gmw-relations/{bro_id_gmw}"
    req = requests.get(url)

    if req.status_code > 200:
        print(req.json()["errors"][0]["message"])

    d = json.loads(req.text)

    if len(d["monitoringTubeReferences"]) == 0:
        logger.debug(
            f"no groundwater level dossier for {bro_id_gmw} and tube number {tube_nr}"
        )
        return None

    for tube in d["monitoringTubeReferences"]:
        if tube["tubeNumber"] == tube_nr:
            if len(tube["gldReferences"]) == 0:
                logger.debug(
                    f"no groundwater level dossier for {bro_id_gmw} and tube number"
                    f"{tube_nr}"
                )
                return None
            else:
                return [gldref["broId"] for gldref in tube["gldReferences"]]

def get_gld_elem_api(bro_id_gld, tmin=None, tmax=None):
    """ get a gld element from a bro_gld_id

    Parameters
    ----------
    bro_id_gld : str
        e.g. 'GLD000000012893'.
    tmin : str or None, optional
        start date in format YYYY-MM-DD
    tmax : str or None, optional
        end date in format YYYY-MM-DD

    Returns
    -------
    gld : Element
        xml reference to gld object

    Raises
    ------
    ValueError
        if bro_id_gld is invalid.

    """

    if not bro_id_gld.startswith("GLD"):
        raise ValueError("can only get observations if bro id starts with GLD")

    url = "https://publiek.broservices.nl/gm/gld/v1/objects/{}"
    params = {}
    if tmin is not None:
        tmin = pd.to_datetime(tmin)
        params["observationPeriodBeginDate"] = tmin.strftime("%Y-%m-%d")
    if tmax is not None:
        tmax = pd.to_datetime(tmax)
        params["observationPeriodEndDate"] = tmax.strftime("%Y-%m-%d")
    req = requests.get(url.format(bro_id_gld), params=params)

    if req.status_code > 200:
        print(req.json()["errors"][0]["message"])

    tree = xml.etree.ElementTree.fromstring(req.text)

    glds = tree.findall(".//ns11:GLD_O", ns_gld)
    if len(glds) != 1:
        raise (Exception("Only one gld supported"))
    gld = glds[0]

    return gld


def get_gld_elem_file(fname):
    """ get a gld element from an xml file

    Parameters
    ----------
    fname : str
        filename of gld xml file e.g. 'GLD000000012893.xml'.


    Returns
    -------
    gld : Element
        xml reference to gld object
    """

    tree = xml.etree.ElementTree.parse(fname)

    glds = tree.findall(".//ns11:GLD_O", ns_gld)
    if len(glds) != 1:
        raise (Exception("Only one gld supported"))
    gld = glds[0]

    return gld


def measurements_from_gld_elem(
    gld, bro_id_gld, tmin=None, tmax=None, to_wintertime=True, drop_duplicate_times=True
):
    """get measurements and metadata from a grondwaterstandonderzoek (gld)
    xml element

    Parameters
    ----------
    gld : Element
        xml reference to gld object
    bro_id_gld : str
        e.g. 'GLD000000012893'.
    tmin : str or None, optional
        start date in format YYYY-MM-DD
    tmax : str or None, optional
        end date in format YYYY-MM-DD
    to_wintertime : bool, optional
        if True the time index is converted to Dutch winter time. The default
        is True.
    drop_duplicate_times : bool, optional
        if True rows with a duplicate time stamp are removed keeping only the
        first row. The default is True.
    add_registration_history : bool, optional
        if True the registration history is added to the metadata. The defualt
        is True.

    Returns
    -------
    df : pd.DataFrame
        measurements.
    meta : dict
        metadata.

    """
    

    meta = {"source": "BRO"}
    meta["monitoring_well"] = gld.find("ns11:monitoringPoint//gldcommon:broId", ns_gld).text
    meta["tube_nr"] = int(
        gld.find("ns11:monitoringPoint//gldcommon:tubeNumber", ns_gld).text
    )
    meta["name"] = f"{meta['monitoring_well']}_{meta['tube_nr']}"
    gmn = gld.find("ns11:groundwaterMonitoringNet//gldcommon:broId", ns_gld)
    if gmn is None:
        meta["monitoringsnet"] = None
    else:
        meta["monitoringsnet"] = gmn.text

    # get observations
    msts = "ns11:observation//om:result//waterml:MeasurementTimeseries"
    times = [time.text for time in gld.findall(f"{msts}//waterml:time", ns_gld)]
    values = [
        np.nan if value.text is None else float(value.text)
        for value in gld.findall(f"{msts}//waterml:value", ns_gld)
    ]
    qualifiers = [q.text for q in gld.findall(f"{msts}//swe:Category//swe:value", ns_gld)]

    # to dataframe
    df = pd.DataFrame(
        index=pd.to_datetime(times, utc=True).tz_convert("CET"),
        data={"values": values, "qualifier": qualifiers},
    )

    # wintertime
    if to_wintertime:
        # remove time zone information by transforming to dutch winter time
        df.index = pd.to_datetime(df.index, utc=True).tz_localize(None) + pd.Timedelta(
            1, unit="h"
        )

    # duplicates
    if df.index.has_duplicates and drop_duplicate_times:
        duplicates = df.index.duplicated(keep="first")
        message = "{} contains {} duplicates (of {}). Keeping only first values."
        logger.debug(message.format(bro_id_gld, duplicates.sum(), len(df)))
        df = df[~duplicates]

    df = df.sort_index()

    # slice to tmin and tmax
    df = df.loc[tmin:tmax]

    return df, meta


def get_full_metadata_from_gmw(bro_id_gmw, tube_nr):
    """get metadata for a groundwater monitoring well.


    Parameters
    ----------
    bro_id_gmw : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.
    tube_nr : int
        filter number you want metadata for.

    Raises
    ------
    ValueError
        if bro_id_gmw is invalid.

    Returns
    -------
    meta : dict
        dictionary with metadata.

    """

    if not bro_id_gmw.startswith("GMW"):
        raise ValueError("can only get metadata if bro id starts with GMW")

    url = f"https://publiek.broservices.nl/gm/gmw/v1/objects/{bro_id_gmw}"
    req = requests.get(url)

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)
    ns = "{http://www.broservices.nl/xsd/dsgmw/1.1}"

    gmws = tree.findall(f".//{ns}GMW_PO")
    if len(gmws) != 1:
        raise (Exception("Only one gmw supported"))
    gmw = gmws[0]
    meta = {"monitoring_well": bro_id_gmw, "tube_nr": tube_nr, "source": "BRO"}
    for child in gmw:
        key = child.tag.split("}", 1)[1]
        if len(child) == 0:
            meta[key] = child.text
        elif key == "deliveredLocation":
            ns = "{http://www.broservices.nl/xsd/gmwcommon/1.1}"
            point = child.find(f"{ns}location")
            ns = "{http://www.opengis.net/gml/3.2}"
            xy = [float(x) for x in point.find(f"{ns}pos").text.split()]
            meta["x"], meta["y"] = xy
        elif key in ["wellConstructionDate"]:
            meta[key] = child[0].text
        elif key == "wellHistory":
            for grandchild in child:
                key = grandchild.tag.split("}", 1)[1]
                meta[key] = grandchild[0].text
        elif key in ["deliveredVerticalPosition", "registrationHistory"]:
            for grandchild in child:
                key = grandchild.tag.split("}", 1)[1]
                meta[key] = grandchild.text
        elif key in ["monitoringTube"]:
            tube = False
            tube_dic = {}
            for grandchild in child:
                if len(grandchild) == 0:
                    tube_key = grandchild.tag.split("}", 1)[1]
                    if tube_key == "tubeNumber":
                        if int(grandchild.text) == tube_nr:
                            tube = True
                    else:
                        tube_dic[tube_key] = grandchild.text

                else:
                    for greatgrandchild in grandchild:
                        tube_key = greatgrandchild.tag.split("}", 1)[1]
                        tube_dic[tube_key] = greatgrandchild.text
            if tube:
                meta.update(tube_dic)

    rename_dic = {
        "groundLevelPosition": "ground_level",
        "screenTopPosition": "screen_top",
        "screenBottomPosition": "screen_bottom",
        "tubeTopPosition": "tube_top",
    }

    for key, val in rename_dic.items():
        meta[val] = meta.pop(key)

    return meta


@lru_cache()
def get_gmw_elem_api(bro_id_gmw, retries=0):
    """get a gmw object from a bro_id

    Parameters
    ----------
    bro_id_gmw : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.
    retries : int, optional
        internal counter used to count number of api retries

    Returns
    -------
    Element
        xml reference to gmw object

    Raises
    ------
    ValueError
        if bro_id_gmw is invalid.
    """

    if not bro_id_gmw.startswith("GMW"):
        raise ValueError("can only get metadata if bro id starts with GMW")

    url = f"https://publiek.broservices.nl/gm/gmw/v1/objects/{bro_id_gmw}"
    req = requests.get(url)

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)

    gmws = tree.findall(".//dsgmw:GMW_PO", ns_gmw)
    if len(gmws) != 1:
        max_retries = rcParams["bro"]["max_retries"]
        val_ind = req.text.find("valid")
        valid = req.text[(val_ind + 9) : (val_ind + 14)]
        if valid == "false" and retries < max_retries:
            logger.debug(
                f"got invalid response for {bro_id_gmw}, trying again {retries+1}/{max_retries}"
            )
            return get_gmw_elem_api(bro_id_gmw, retries=retries + 1)
        elif valid == "false":
            raise Exception(
                f"got invalid response for {bro_id_gmw} after trying {retries} times"
            )
        raise (Exception("Only one gmw supported"))
    gmw = gmws[0]

    return gmw


def get_gmw_elem_file(fname):
    """get a gmw object from an xml file

    Parameters
    ----------
    fname : str
        filename of gmw xml file e.g. 'GMW000000036287.xml'.

    Returns
    -------
    Element
        xml reference to gmw object
    """

    # read results
    tree = xml.etree.ElementTree.parse(fname)

    if tree.find("brocom:requestReference", ns_gmw).text == 'uitgifteloket':
        gmwtag = "GMW_PPO"
    else:
        gmwtag = "GMW_PO"

    gmws = tree.findall(f".//dsgmw:{gmwtag}", ns_gmw)
    if len(gmws) != 1:
        raise (Exception("Only one gmw supported"))
    gmw = gmws[0]

    return gmw

def get_tube_nrs_from_gmw(gmw):
    """returns all tube numbers from a groundwater monitoring well (gmw)

    Parameters
    ----------
    gmw : element
        xml reference to gmw object

    Returns
    -------
    list of int
        tube numbers
    """

    # get tube nrs
    tube_numbers = [
        int(tube.text)
        for tube in gmw.findall("dsgmw:monitoringTube/dsgmw:tubeNumber", ns_gmw)
    ]

    return tube_numbers


def get_metadata_from_gmw_elem(gmw, bro_id_gmw, tube_nr):
    """get selection of metadata for a groundwater monitoring well.
    coordinates, ground_level, tube_top and tube screen

    Parameters
    ----------
    gmw : element
        xml reference to gmw object
    bro_id_gmw : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.
    tube_nr : int
        tube number you want metadata for.

    Raises
    ------
    TypeError
        if tube_nr is not an int

    Returns
    -------
    meta : dict
        dictionary with metadata.

    """
    if not isinstance(tube_nr, int):
        raise TypeError(f"expected integer got {type(tube_nr)}")

    meta = {"monitoring_well": bro_id_gmw, "tube_nr": tube_nr, "source": "BRO",
            "filename":''}

    # x and y
    xy_elem = gmw.find("dsgmw:deliveredLocation//gmwcommon:location//gml:pos", ns_gmw)
    xy = [float(val) for val in xy_elem.text.split()]

    # convert crs
    srsname = gmw.find("dsgmw:deliveredLocation//gmwcommon:location", ns_gmw).attrib[
        "srsName"
    ]
    epsg_gwm = int(srsname.split(":")[-1])
    proj_from = Proj(f"EPSG:{epsg_gwm}")
    proj_to = Proj(EPSG_28992)
    transformer = Transformer.from_proj(proj_from, proj_to)
    xy = transformer.transform(xy[0], xy[1])

    meta["x"], meta["y"] = xy

    # ground_level
    vert_pos = gmw.find("dsgmw:deliveredVerticalPosition", ns_gmw)
    mv = vert_pos.find("gmwcommon:groundLevelPosition", ns_gmw)
    datum = vert_pos.find("gmwcommon:verticalDatum", ns_gmw)
    if datum.text == "NAP" and mv.attrib["uom"] == "m":
        meta["unit"] = "m NAP"
        if mv.text is None:
            meta["ground_level"] = np.nan
        else:
            meta["ground_level"] = float(mv.text)
    else:
        raise ValueError("invalid ground_level datum or unit")

    # buis eigenschappen
    tubes = gmw.findall("dsgmw:monitoringTube", ns_gmw)
    tube_nrs = [int(tube.find("dsgmw:tubeNumber", ns_gmw).text) for tube in tubes]
    if tube_nr not in tube_nrs:
        raise ValueError(
            f"gmw {bro_id_gmw} has no tube_nr {tube_nr} please choose a tube_nr from"
            f"{tube_nrs}"
        )
    tube = tubes[tube_nrs.index(tube_nr)]

    # tube_top
    mp = tube.find("dsgmw:tubeTopPosition", ns_gmw)
    if mp.attrib["uom"] == "m":
        meta["tube_top"] = float(mp.text)

    # bovenkant filter
    bkf = tube.find("dsgmw:screen//dsgmw:screenTopPosition", ns_gmw)
    if bkf.attrib["uom"] == "m":
        meta["screen_top"] = float(bkf.text)

    # onderkant filter
    okf = tube.find("dsgmw:screen//dsgmw:screenBottomPosition", ns_gmw)
    if okf.attrib["uom"] == "m":
        meta["screen_bottom"] = float(okf.text)

    meta["metadata_available"] = True

    return meta


def get_obs_list_from_extent(
    extent,
    ObsClass,
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=True,
    epsg=28992,
    ignore_max_obs=False,
):
    """get a list of gmw observations within an extent.

    Parameters
    ----------
    extent : list, tuple, numpy-array or None, optional
        get groundwater monitoring wells within this extent
        [xmin, xmax, ymin, ymax]
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the list, even without measurements
    epsg : int, optional
        epsg code of the extent. The default is 28992 (RD).
    ignore_max_obs : bool, optional
        by default you get a prompt if you want to download over a 1000
        observations at once. if ignore_max_obs is True you won't get the
        prompt. The default is False

    Returns
    -------
    obs_list : list
        list of observations.

    """

    if only_metadata and not keep_all_obs:
        logger.error(
            "you will get an empty ObsCollection with only_metadata is True and"
            "keep_all_obs is False"
        )

    url = "https://publiek.broservices.nl/gm/gmw/v1/characteristics/searches?"

    data = {}
    if tmin is None or tmax is None:
        data["registrationPeriod"] = {}
        if tmin is not None:
            beginDate = pd.to_datetime(tmin).strftime("%Y-%m-%d")
            data["registrationPeriod"]["beginDate"] = beginDate
        if tmax is not None:
            endDate = pd.to_datetime(tmax).strftime("%Y-%m-%d")
            data["registrationPeriod"]["endDate"] = endDate

    transformer = Transformer.from_crs(epsg, 4326)
    data["area"] = {}
    if extent is not None:
        lat1, lon1 = transformer.transform(extent[0], extent[2])
        lat2, lon2 = transformer.transform(extent[1], extent[3])
        data["area"]["boundingBox"] = {
            "lowerCorner": {"lat": lat1, "lon": lon1},
            "upperCorner": {"lat": lat2, "lon": lon2},
        }
    req = requests.post(url, json=data)
    if req.status_code > 200:
        logger.error(
            "could not get monitoring wells, your extent is probably too big."
            "Try a smaller extent"
        )
        req.raise_for_status()
        # print(req.json()["errors"][0]["message"])

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)

    if tree.find(".//brocom:responseType", ns_gmw).text == "rejection":
        raise RuntimeError(tree.find(".//brocom:rejectionReason", ns_gmw).text)

    gmws_ids = np.unique(
        [gmw.text for gmw in tree.findall(".//dsgmw:GMW_C//brocom:broId", ns_gmw)]
    )

    if len(gmws_ids) > 1000 and not ignore_max_obs:
        ans = input(
            f"You requested to download {len(gmws_ids)} observations, this can"
            "take a while. Are you sure you want to continue [Y/n]? "
        )
        if ans not in ["Y", "y", "yes", "Yes", "YES"]:
            return []

    obs_list = []
    for bro_id_gmw in tqdm(gmws_ids):
        gmw = get_gmw_elem_api(bro_id_gmw)
        tube_nrs = get_tube_nrs_from_gmw(gmw)
        for tube_nr in tube_nrs:
            o = ObsClass.from_bro(
                bro_id_gmw,
                tube_nr=tube_nr,
                tmin=tmin,
                tmax=tmax,
                only_metadata=only_metadata,
            )
            if o.empty:
                logger.debug(
                    f"no measurements found for gmw_id {bro_id_gmw} and tube number"
                    f"{tube_nr}"
                )
                if keep_all_obs:
                    obs_list.append(o)
            else:
                obs_list.append(o)

    return obs_list



def get_obs_list_from_directory(dirname, ObsClass,
                                gmwsdir='BRO_Grondwatermonitoringput',
                                gldsdir='BRO_Grondwaterstandonderzoek',
                                unpackdir=None, force_unpack=False,
                                tmin=None, tmax=None, 
                                only_metadata=False, keep_all_obs=False,
                                **kwargs):
    """ get a list of observations from a directory. Tries to read all the xml files
    that start with GMW and GLD.

    Parameters
    ----------
    dirname : str
        directory name, can be a .zip file or the parent directory of gmwsdir en gldsdir
    ObsClass : type
        class of the observations, e.g. hpd.GroundwaterObs
    gmwsdir : str, optional
        subdirectory with xml files starting with GMW, by default 'BRO_Grondwatermonitoringput'
    gldsdir : str, optional
        subdirectory with xml files starting with GLD, by default 'BRO_Grondwaterstandonderzoek'
    unpackdir : str
        destination directory of the unzipped file
    force_unpack : boolean, optional
        force unpack if dst already exists
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the list, even without measurements
    only_metadata : 

    Returns
    -------
    obs_list : list
        list of observations.
    """

    if dirname.endswith(".zip"):
        from ..util import unzip_file

        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(
            zipf, dirname, force=force_unpack
        )

    # gmnfiles = os.listdir(os.path.join(dirname, gmnsdir))
    gmwfiles = os.listdir(os.path.join(dirname, gmwsdir))
    gldfiles = os.listdir(os.path.join(dirname, gldsdir))

    obs_list = []
    obs_names = []

    if not only_metadata:
        for fname in gldfiles:
            # read gld file
            fpath = os.path.join(dirname, gldsdir, fname)
            bro_id_gld = fname.strip('.xml')
            gld = get_gld_elem_file(fpath)
            measurements, meta = measurements_from_gld_elem(gld, 
                                                            bro_id_gld, tmin=tmin, tmax=tmax, **kwargs)
            
            # get metadata from gmw
            bro_id_gmw = meta["monitoring_well"]
            fpath_gmw = os.path.join(dirname, gmwsdir, f'{bro_id_gmw}.xml')
            gmw = get_gmw_elem_file(fpath_gmw)
            meta.update(get_metadata_from_gmw_elem(gmw, bro_id_gmw, 
                                                meta["tube_nr"]))
            obs_names.append(meta['name'])

            o = ObsClass(
                measurements,
                meta=meta,
                name=meta.pop("name"),
                filename=fpath,
                x=meta.pop("x"),
                y=meta.pop("y"),
                source=meta.pop("source"),
                unit=meta.pop("unit"),
                screen_bottom=meta.pop("screen_bottom"),
                screen_top=meta.pop("screen_top"),
                ground_level=meta.pop("ground_level"),
                metadata_available=meta.pop("metadata_available"),
                monitoring_well=meta.pop("monitoring_well"),
                tube_nr=meta.pop("tube_nr"),
                tube_top=meta.pop("tube_top"),)
            
            obs_list.append(o)
    
    if not keep_all_obs:
        return obs_list
    
    # check for tubes without measurements (GLD files)
    for fname in gmwfiles:
        fpath_gmw = os.path.join(dirname, gmwsdir, fname)
        bro_id_gmw = fname.strip('.xml')
        gmw = get_gmw_elem_file(fpath_gmw)

        tubes = gmw.findall("dsgmw:monitoringTube", ns_gmw)
        tube_nrs = [int(tube.find("dsgmw:tubeNumber", ns_gmw).text) for tube in tubes]
        for tube_nr in tube_nrs:
            name = f"{bro_id_gmw}_{tube_nr}"
            if name not in obs_names:
                meta = get_metadata_from_gmw_elem(gmw, bro_id_gmw, tube_nr)
                meta["name"] = name
                empty_df = pd.DataFrame()
                
                o = ObsClass(
                        empty_df,
                        meta=meta,
                        name=meta.pop("name"),
                        filename=fpath_gmw,
                        x=meta.pop("x"),
                        y=meta.pop("y"),
                        source=meta.pop("source"),
                        unit=meta.pop("unit"),
                        screen_bottom=meta.pop("screen_bottom"),
                        screen_top=meta.pop("screen_top"),
                        ground_level=meta.pop("ground_level"),
                        metadata_available=meta.pop("metadata_available"),
                        monitoring_well=meta.pop("monitoring_well"),
                        tube_nr=meta.pop("tube_nr"),
                        tube_top=meta.pop("tube_top"),)
                obs_list.append(o)
    
    return obs_list