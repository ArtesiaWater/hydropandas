"""
https://www.bro-productomgeving.nl/bpo/latest/informatie-voor-softwareleveranciers/url-s-publieke-rest-services

function levels:
1. get_obs_list_from_gmn: list of observations from grondwatermonitoring net
    2. get_bro_groundwater: single observation
        3. measurements_from_gld:
        4. get_metadata_from_gmw:
        5. get_gld_ids_from_gmw:

"""

import json
import logging
import xml.etree.ElementTree
from functools import lru_cache

import numpy as np
import pandas as pd
import requests
from pyproj import Proj, Transformer
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

from ..rcparams import rcParams
from ..util import EPSG_28992

logger = logging.getLogger(__name__)


def get_obs_list_from_gmn(bro_id, ObsClass, only_metadata=False, keep_all_obs=True):
    """get a list of observation from a groundwater monitoring network.

    Parameters
    ----------
    bro_id : str
        starts with 'GMN' e.g. 'GMN000000000163'.
    ObsClass : type
        class of the observations, so far only GroundwaterObs is supported
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without
        measurements

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    obs_list : list
        list with observation objects.
    meta : dict
        metadata of the groundwater monitoring net.

    """

    if not bro_id.startswith("GMN"):
        raise ValueError("bro id should start with GMN")

    url = f"https://publiek.broservices.nl/gm/gmn/v1/objects/{bro_id}"
    req = requests.get(url)

    if req.status_code > 200:
        logger.error(
            "could not get monitoring wells, this groundwater monitoring net is"
            "probably too big. Try a smaller extent"
        )
        req.raise_for_status()

    ns = {"xmlns": "http://www.broservices.nl/xsd/dsgmn/1.0"}

    tree = xml.etree.ElementTree.fromstring(req.text)
    gmn = tree.find(".//xmlns:GMN_PO", ns)
    gmws = gmn.findall("xmlns:measuringPoint", ns)

    logger.info(f"{len(gmws)} groundwater monitoring wells within groundwater meetnet")

    obs_list = []
    for gmw in tqdm(gmws):
        tags = ["MeasuringPoint", "monitoringTube", "GroundwaterMonitoringTube"]
        tube = gmw.find(f".//xmlns:{'//xmlns:'.join(tags)}", ns)
        gmw_id = tube.find("xmlns:broId", ns).text
        tube_nr = int(tube.find("xmlns:tubeNumber", ns).text)

        o = ObsClass.from_bro(
            bro_id=gmw_id, tube_nr=tube_nr, only_metadata=only_metadata
        )
        if o.empty:
            logger.debug(
                f"no measurements found for gmw_id {gmw_id} and tube number {tube_nr}"
            )
            if keep_all_obs:
                obs_list.append(o)
        else:
            obs_list.append(o)
        obs_list.append(o)

    meta = {}
    meta["name"] = gmn.find("xmlns:name", ns).text
    meta["doel"] = gmn.find("xmlns:monitoringPurpose", ns).text
    meta["bro_id"] = bro_id

    return obs_list, meta


def get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False, **kwargs):
    """get bro groundwater measurement from a GLD id or a GMW id with a
    filter number.


    Parameters
    ----------
    bro_id : str
        starts with 'GLD' or 'GMW' e.g. 'GLD000000012893'.
    tube_nr : int or None, optional
        tube number, required if bro_id starts with 'GMW'. The default is
        None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    **kwargs :
        passes to measurements_from_gld.

    Raises
    ------
    ValueError
        DESCRIPTION.

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
        return measurements_from_gld(bro_id, **kwargs)

    elif bro_id.startswith("GMW"):
        if tube_nr is None:
            raise ValueError("if bro_id is GMW a tube_nr should be specified")

        meta = get_metadata_from_gmw(bro_id, tube_nr)
        gld_ids = get_gld_ids_from_gmw(bro_id, tube_nr)

        if gld_ids is None:
            meta["name"] = f"{bro_id}_{tube_nr}"
            only_metadata = True  # cannot get time series without gld id
        else:
            meta["name"] = f"{bro_id}_{tube_nr}"
            meta["gld_ids"] = gld_ids

        if only_metadata:
            empty_df = pd.DataFrame()
            return empty_df, meta

        dfl = []
        for i, gld_id in enumerate(gld_ids):
            df, meta_new = measurements_from_gld(gld_id, **kwargs)
            meta.update(meta_new)
            dfl.append(df)
        df = pd.concat(dfl, axis=0).sort_index()

        return df, meta


def get_gld_ids_from_gmw(bro_id, tube_nr):
    """get bro_ids of multiple grondwaterstandendossier (gld) from a bro_id of a
    grondwatermonitoringsput (gmw).

    Parameters
    ----------
    bro_id : str
        starts with 'GLD' or 'GMW' e.g. 'GMW000000036287'.
    tube_nr : int
        tube number.

    Raises
    ------
    ValueError
        DESCRIPTION.
    RuntimeError
        DESCRIPTION.

    Returns
    -------
    list of str or None
        bro_ids of a grondwaterstandonderzoek (gld).

    """
    if not bro_id.startswith("GMW"):
        raise ValueError("bro id should start with GMW")

    url = f"https://publiek.broservices.nl/gm/v1/gmw-relations/{bro_id}"
    req = requests.get(url)

    if req.status_code > 200:
        logger.error(f"could not get gld ids {req.text}")
        req.raise_for_status()

    d = json.loads(req.text)

    if len(d["monitoringTubeReferences"]) == 0:
        logger.debug(
            f"no groundwater level dossier for {bro_id} and tube number {tube_nr}"
        )
        return None

    for tube in d["monitoringTubeReferences"]:
        if tube["tubeNumber"] == tube_nr:
            if len(tube["gldReferences"]) == 0:
                logger.debug(
                    f"no groundwater level dossier for {bro_id} and tube number"
                    f"{tube_nr}"
                )
                return None
            else:
                return [gldref["broId"] for gldref in tube["gldReferences"]]


def measurements_from_gld(
    bro_id, tmin=None, tmax=None, to_wintertime=True, drop_duplicate_times=True
):
    """get measurements and metadata from a grondwaterstandonderzoek (gld)
    bro_id


    Parameters
    ----------
    bro_id : str
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

    Raises
    ------
    ValueError
        if bro_id is invalid.

    Returns
    -------
    df : pd.DataFrame
        measurements.
    meta : dict
        metadata.

    """
    if not bro_id.startswith("GLD"):
        raise ValueError("can only get observations if bro id starts with GLD")

    url = "https://publiek.broservices.nl/gm/gld/v1/objects/{}"
    params = {}
    if tmin is not None:
        tmin = pd.to_datetime(tmin)
        params["observationPeriodBeginDate"] = tmin.strftime("%Y-%m-%d")
    if tmax is not None:
        tmax = pd.to_datetime(tmax)
        params["observationPeriodEndDate"] = tmax.strftime("%Y-%m-%d")

    # add some logic to retry in case of a 429 response
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5, status_forcelist=[429])
    s.mount("https://", HTTPAdapter(max_retries=retries))
    req = s.get(url.format(bro_id), params=params)

    if req.status_code > 200:
        req.raise_for_status()

    ns = {
        "ns11": "http://www.broservices.nl/xsd/dsgld/1.0",
        "gldcommon": "http://www.broservices.nl/xsd/gldcommon/1.0",
        "waterml": "http://www.opengis.net/waterml/2.0",
        "swe": "http://www.opengis.net/swe/2.0",
        "om": "http://www.opengis.net/om/2.0",
    }

    tree = xml.etree.ElementTree.fromstring(req.text)

    glds = tree.findall(".//ns11:GLD_O", ns)
    if len(glds) != 1:
        raise (Exception("Only one gld supported"))
    gld = glds[0]

    meta = {"source": "BRO"}
    meta["location"] = gld.find("ns11:monitoringPoint//gldcommon:broId", ns).text
    meta["tube_nr"] = int(
        gld.find("ns11:monitoringPoint//gldcommon:tubeNumber", ns).text
    )
    meta["name"] = f"{meta['location']}_{meta['tube_nr']}"
    gmn = gld.find("ns11:groundwaterMonitoringNet//gldcommon:broId", ns)
    if gmn is None:
        meta["monitoringsnet"] = None
    else:
        meta["monitoringsnet"] = gmn.text

    # get observations
    msts = "ns11:observation//om:result//waterml:MeasurementTimeseries"
    times = [time.text for time in gld.findall(f"{msts}//waterml:time", ns)]
    values = [
        np.nan if value.text is None else float(value.text)
        for value in gld.findall(f"{msts}//waterml:value", ns)
    ]
    qualifiers = [q.text for q in gld.findall(f"{msts}//swe:Category//swe:value", ns)]

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
        logger.debug(message.format(bro_id, duplicates.sum(), len(df)))
        df = df[~duplicates]

    df = df.sort_index()

    # slice to tmin and tmax
    df = df.loc[tmin:tmax]

    # add metadata from gmw
    meta.update(get_metadata_from_gmw(meta["location"], meta["tube_nr"]))

    return df, meta


def get_full_metadata_from_gmw(bro_id, tube_nr):
    """get metadata for a groundwater monitoring well.


    Parameters
    ----------
    bro_id : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.
    tube_nr : int
        filter number you want metadata for.

    Raises
    ------
    ValueError
        if bro_id is invalid.

    Returns
    -------
    meta : dict
        dictionary with metadata.

    """

    if not bro_id.startswith("GMW"):
        raise ValueError("can only get metadata if bro id starts with GMW")

    url = f"https://publiek.broservices.nl/gm/gmw/v1/objects/{bro_id}"
    req = requests.get(url)

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)
    ns = "{http://www.broservices.nl/xsd/dsgmw/1.1}"

    gmws = tree.findall(f".//{ns}GMW_PO")
    if len(gmws) != 1:
        raise (Exception("Only one gmw supported"))
    gmw = gmws[0]
    meta = {"location": bro_id, "tube_nr": tube_nr, "source": "BRO"}
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
def _get_gmw_from_bro_id(bro_id, retries=0):
    """get a gmw object from a bro_id

    Parameters
    ----------
    bro_id : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.

    Returns
    -------
    Element
        xml reference to gmw object

    Raises
    ------
    ValueError
        if bro_id is invalid.
    """

    if not bro_id.startswith("GMW"):
        raise ValueError("can only get metadata if bro id starts with GMW")

    ns = {
        "dsgmw": "http://www.broservices.nl/xsd/dsgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "gml": "http://www.opengis.net/gml/3.2",
    }

    url = f"https://publiek.broservices.nl/gm/gmw/v1/objects/{bro_id}"
    req = requests.get(url)

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)

    gmws = tree.findall(".//dsgmw:GMW_PO", ns)
    if len(gmws) != 1:
        max_retries = rcParams["bro"]["max_retries"]
        val_ind = req.text.find("valid")
        valid = req.text[(val_ind + 9) : (val_ind + 14)]
        if valid == "false" and retries < max_retries:
            logger.debug(
                f"got invalid response for {bro_id}, trying again {retries + 1}/{max_retries}"
            )
            return _get_gmw_from_bro_id(bro_id, retries=retries + 1)
        elif valid == "false":
            raise Exception(
                f"got invalid response for {bro_id} after trying {retries} times"
            )
        raise (Exception("Only one gmw supported"))
    gmw = gmws[0]

    return gmw


def get_tube_nrs_from_gmw(bro_id):
    """returns all tube numbers from a groundwater monitoring well (gmw)

    Parameters
    ----------
    bro_id : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.

    Returns
    -------
    list of int
        tube numbers
    """
    ns = {
        "dsgmw": "http://www.broservices.nl/xsd/dsgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "gml": "http://www.opengis.net/gml/3.2",
    }

    gmw = _get_gmw_from_bro_id(bro_id)

    # get tube nrs
    tube_numbers = [
        int(tube.text)
        for tube in gmw.findall("dsgmw:monitoringTube/dsgmw:tubeNumber", ns)
    ]

    return tube_numbers


def get_metadata_from_gmw(bro_id, tube_nr):
    """get selection of metadata for a groundwater monitoring well.
    coordinates, ground_level, tube_top and tube screen

    Parameters
    ----------
    bro_id : str
        bro id of groundwater monitoring well e.g. 'GMW000000036287'.
    tube_nr : int
        tube number you want metadata for.

    Raises
    ------
    ValueError
        if bro_id is invalid.
    TypeError
        if tube_nr is not an int

    Returns
    -------
    meta : dict
        dictionary with metadata.

    """
    if not isinstance(tube_nr, int):
        raise TypeError(f"expected integer got {type(tube_nr)}")

    ns = {
        "dsgmw": "http://www.broservices.nl/xsd/dsgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "gml": "http://www.opengis.net/gml/3.2",
    }

    gmw = _get_gmw_from_bro_id(bro_id)

    meta = {"location": bro_id, "tube_nr": tube_nr, "source": "BRO"}

    # x and y
    xy_elem = gmw.find("dsgmw:deliveredLocation//gmwcommon:location//gml:pos", ns)
    xy = [float(val) for val in xy_elem.text.split()]

    # convert crs
    srsname = gmw.find("dsgmw:deliveredLocation//gmwcommon:location", ns).attrib[
        "srsName"
    ]
    epsg_gwm = int(srsname.split(":")[-1])
    proj_from = Proj(f"EPSG:{epsg_gwm}")
    proj_to = Proj(EPSG_28992)
    transformer = Transformer.from_proj(proj_from, proj_to)
    xy = transformer.transform(xy[0], xy[1])

    meta["x"], meta["y"] = xy

    # ground_level
    vert_pos = gmw.find("dsgmw:deliveredVerticalPosition", ns)
    mv = vert_pos.find("gmwcommon:groundLevelPosition", ns)
    datum = vert_pos.find("gmwcommon:verticalDatum", ns)
    if datum.text == "NAP" and mv.attrib["uom"] == "m":
        meta["unit"] = "m NAP"
        if mv.text is None:
            meta["ground_level"] = np.nan
        else:
            meta["ground_level"] = float(mv.text)
    else:
        raise ValueError("invalid ground_level datum or unit")

    # buis eigenschappen
    tubes = gmw.findall("dsgmw:monitoringTube", ns)
    tube_nrs = [int(tube.find("dsgmw:tubeNumber", ns).text) for tube in tubes]
    if tube_nr not in tube_nrs:
        raise ValueError(
            f"gmw {bro_id} has no tube_nr {tube_nr} please choose a tube_nr from"
            f"{tube_nrs}"
        )
    tube = tubes[tube_nrs.index(tube_nr)]

    # tube_top
    mp = tube.find("dsgmw:tubeTopPosition", ns)
    if mp.attrib["uom"] == "m":
        meta["tube_top"] = float(mp.text)

    # bovenkant filter
    bkf = tube.find("dsgmw:screen//dsgmw:screenTopPosition", ns)
    if bkf.attrib["uom"] == "m":
        meta["screen_top"] = float(bkf.text)

    # onderkant filter
    okf = tube.find("dsgmw:screen//dsgmw:screenBottomPosition", ns)
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
    epsg : int, optional
        epsg code of the extent. The default is 28992 (RD).
    ignore_max_obs : bool, optional
        by default you get a prompt if you want to download over a 1000
        observations at once. if ignore_max_obs is True you won't get the
        prompt. The default is False

    Raises
    ------

        DESCRIPTION.

    Returns
    -------
    obs_list : TYPE
        DESCRIPTION.

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

    data["area"] = {}
    if epsg == 4326:
        data["area"]["boundingBox"] = {
            "lowerCorner": {"lat": extent[2], "lon": extent[0]},
            "upperCorner": {"lat": extent[3], "lon": extent[1]},
        }
    else:
        transformer = Transformer.from_crs(epsg, 4326)
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

    # read results
    tree = xml.etree.ElementTree.fromstring(req.text)

    ns = {
        "dsgmw": "http://www.broservices.nl/xsd/dsgmw/1.1",
        "gml": "http://www.opengis.net/gml/3.2",
        "brocom": "http://www.broservices.nl/xsd/brocommon/3.0",
    }

    if tree.find(".//brocom:responseType", ns).text == "rejection":
        raise RuntimeError(tree.find(".//brocom:rejectionReason", ns).text)

    gmws_ids = np.unique(
        [gmw.text for gmw in tree.findall(".//dsgmw:GMW_C//brocom:broId", ns)]
    )

    if len(gmws_ids) > 1000 and not ignore_max_obs:
        ans = input(
            f"You requested to download {len(gmws_ids)} observations, this can"
            "take a while. Are you sure you want to continue [Y/n]? "
        )
        if ans not in ["Y", "y", "yes", "Yes", "YES"]:
            return []

    obs_list = []
    for gmw_id in tqdm(gmws_ids):
        gmws = tree.findall(f'.//*[brocom:broId="{gmw_id}"]', ns)
        if len(gmws) < 1:
            raise RuntimeError("unexpected")

        tube_nrs = get_tube_nrs_from_gmw(gmw_id)
        for tube_nr in tube_nrs:
            o = ObsClass.from_bro(
                gmw_id,
                tube_nr=tube_nr,
                tmin=tmin,
                tmax=tmax,
                only_metadata=only_metadata,
            )
            if o.empty:
                logger.debug(
                    f"no measurements found for gmw_id {gmw_id} and tube number"
                    f"{tube_nr}"
                )
                if keep_all_obs:
                    obs_list.append(o)
            else:
                obs_list.append(o)

    return obs_list
