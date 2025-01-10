"""
This module is all about data is to be submitted to
https://www.bronhouderportaal-bro.nl, i.e. XMLs of GMWs that are newly
installed by a fieldwork company
"""

import logging
import pathlib
import xml.etree.ElementTree

import pandas as pd

logger = logging.getLogger(__name__)


def get_tube_nrs_from_xml(tree, ns):
    """
    Get tube numbers from xml data

    Parameters
    ----------
    tree : xml.etree.ElementTree.ElementTree
        The XML input.
    ns : dictonary
        Namespaces in XML input

    Returns
    -------
    all_tube_nrs : list
        The numbers of all filters in the xml data.

    """
    # get numbers of individual filters from XML file
    all_tube_nrs = []
    tubes = tree.findall(
        "isgmw:sourceDocument//isgmw:GMW_Construction//isgmw:monitoringTube", ns
    )
    for tube in tubes:
        all_tube_nrs.append(int(tube.find("isgmw:tubeNumber", ns).text))

    return all_tube_nrs


def get_obs_list_from_dir(
    dirname,
    ObsClass,
    full_meta=False,
):
    """get a list of gmw observations within a dirname.

    Parameters
    ----------
    dirname : str
        name of directory with XML files
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    full_meta : bool, optional
        process not only the standard metadata to ObsCollection


    Returns
    -------
    obs_list : list
        all observation wells

    """
    ns = {
        "brocommon": "http://www.broservices.nl/xsd/brocommon/3.0",
        "isgmw": "http://www.broservices.nl/xsd/isgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "d5p1": "http://www.opengis.net/gml/3.2",
    }

    obs_list = []

    # loop over files
    for path_xml in pathlib.Path(dirname).glob("*.xml"):
        # open xml and get number of monitoring tubes
        tree = xml.etree.ElementTree.parse(path_xml)
        # get number of tubes from XML file
        ntubes_from_tree = int(
            tree.find(
                "isgmw:sourceDocument//"
                "isgmw:GMW_Construction//"
                "isgmw:numberOfMonitoringTubes",
                ns,
            ).text
        )
        all_tube_nrs = get_tube_nrs_from_xml(tree, ns)
        if len(all_tube_nrs) != ntubes_from_tree:
            logger.warning(
                f"{path_xml} numberOfMonitoringTubes suggests that "
                f"{ntubes_from_tree} filters are present, "
                f"however meta data for {len(all_tube_nrs)} filters is "
                f"in the xml. Meta data for {len(all_tube_nrs)} filters "
                "is processed."
            )

        for tube_nr in all_tube_nrs:
            # append each tube to the obs_list
            o = ObsClass.from_bronhouderportaal_bro(
                path_xml,
                tube_nr=tube_nr,
                full_meta=full_meta,
            )
            obs_list.append(o)

    return obs_list


def get_metadata_from_gmw(path_xml, tube_nr, full_meta=False):
    """get metadata for a groundwater monitoring well.
    standard: coordinates, ground_level, tube_top and tube screen

    Function is based upon get_metadata_from_gmw in bro.py, changes are made because:
        - bro_id is not available, because data is not yet imported to bro
        - XML layout of input is slightly different (example of Fugro is used)


    Parameters
    ----------
    path_xml : pathlib.WindowsPath
        path of groundwater monitoring well XML file.
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
    if not isinstance(path_xml, pathlib.WindowsPath):
        try:
            path_xml = pathlib.Path(path_xml)
        except ValueError:
            print(f"invalid path_xml {path_xml}")

    # open meta dictonary
    meta = {
        "filename": path_xml.name,
        "tube_nr": tube_nr,
        "source": "bronhouderportaal-bro",
    }

    # open xml file
    tree = xml.etree.ElementTree.parse(path_xml)
    ns = {
        "brocommon": "http://www.broservices.nl/xsd/brocommon/3.0",
        "isgmw": "http://www.broservices.nl/xsd/isgmw/1.1",
        "gmwcommon": "http://www.broservices.nl/xsd/gmwcommon/1.1",
        "d5p1": "http://www.opengis.net/gml/3.2",
    }

    # is tube_nr valid
    all_tube_nrs = get_tube_nrs_from_xml(tree, ns)
    if tube_nr not in all_tube_nrs:
        raise ValueError(
            f"Request is to process tube_nr {tube_nr}, that number is not "
            f"in the XML data. Valid tube_nrs are {all_tube_nrs} "
        )

    # name of well via requestReference
    # bro_id is not used, because well is not yet imported in Broloket
    if tree.find("brocommon:requestReference", ns) is not None:
        name = tree.find("brocommon:requestReference", ns).text
    else:
        name = path_xml.stem  # filname without extention
    meta["name"] = f"{name}-{tube_nr}"
    meta["monitoring_well"] = name

    if full_meta:
        if tree.find("brocommon:deliveryAccountableParty", ns) is not None:
            # deliveryAccountableParty is KvK nummer
            meta["deliveryAccountableParty"] = tree.find(
                "brocommon:deliveryAccountableParty", ns
            ).text
        if tree.find("brocommon:qualityRegime", ns) is not None:
            meta["qualityRegime"] = tree.find("brocommon:qualityRegime", ns).text

    # open sourceDocument part
    sourceDocument = tree.find("isgmw:sourceDocument", ns)

    # GMW_Construction element
    GMW_c = sourceDocument.find("isgmw:GMW_Construction", ns)

    if full_meta:
        if GMW_c.find("isgmw:objectIdAccountableParty", ns) is not None:
            meta["objectIdAccountableParty"] = GMW_c.find(
                "isgmw:objectIdAccountableParty", ns
            ).text
        if GMW_c.find("isgmw:deliveryContext", ns) is not None:
            meta["deliveryContext"] = GMW_c.find("isgmw:deliveryContext", ns).text
        if GMW_c.find("isgmw:constructionStandard", ns) is not None:
            meta["constructionStandard"] = GMW_c.find(
                "isgmw:constructionStandard", ns
            ).text
        if GMW_c.find("isgmw:initialFunction", ns) is not None:
            meta["initialFunction"] = GMW_c.find("isgmw:initialFunction", ns).text
        if GMW_c.find("isgmw:groundLevelStable", ns) is not None:
            meta["groundLevelStable"] = GMW_c.find("isgmw:groundLevelStable", ns).text
        if GMW_c.find("isgmw:owner", ns) is not None:
            meta["owner"] = GMW_c.find("isgmw:owner", ns).text
        if GMW_c.find("isgmw:wellHeadProtector", ns) is not None:
            meta["wellHeadProtector"] = GMW_c.find("isgmw:wellHeadProtector", ns).text
        if GMW_c.find("isgmw:wellConstructionDate", ns) is not None:
            xml_text = GMW_c.find("isgmw:wellConstructionDate//brocommon:date", ns).text
            meta["wellConstructionDate"] = pd.Timestamp(xml_text)

    # x and y
    xy = GMW_c.find("isgmw:deliveredLocation//gmwcommon:location//d5p1:pos", ns)
    meta["x"], meta["y"] = [float(val) for val in xy.text.split()]
    if full_meta:
        if (
            GMW_c.find(
                "isgmw:deliveredLocation//gmwcommon:horizontalPositioningMethod", ns
            )
            is not None
        ):
            meta["horizontalPositioningMethod"] = GMW_c.find(
                "isgmw:deliveredLocation//gmwcommon:horizontalPositioningMethod", ns
            ).text
        if (
            GMW_c.find(
                "isgmw:deliveredVerticalPosition//"
                "gmwcommon:localVerticalReferencePoint",
                ns,
            )
            is not None
        ):
            meta["localVerticalReferencePoint"] = GMW_c.find(
                "isgmw:deliveredVerticalPosition//"
                "gmwcommon:localVerticalReferencePoint",
                ns,
            ).text
        if (
            GMW_c.find("isgmw:deliveredVerticalPosition//gmwcommon:offset", ns)
            is not None
        ):
            meta["deliveredVerticalPosition_offset"] = float(
                GMW_c.find("isgmw:deliveredVerticalPosition//gmwcommon:offset", ns).text
            )
        if (
            GMW_c.find(
                "isgmw:deliveredVerticalPosition//"
                "gmwcommon:groundLevelPositioningMethod",
                ns,
            )
            is not None
        ):
            meta["groundLevelPositioningMethod"] = GMW_c.find(
                "isgmw:deliveredVerticalPosition//"
                "gmwcommon:groundLevelPositioningMethod",
                ns,
            ).text

    # ground_level
    glp_xml = GMW_c.find(
        "isgmw:deliveredVerticalPosition//gmwcommon:groundLevelPosition", ns
    )
    vert_datum = GMW_c.find(
        "isgmw:deliveredVerticalPosition//gmwcommon:verticalDatum", ns
    ).text
    meta["unit"] = glp_xml.attrib["uom"] + " " + vert_datum
    if glp_xml.attrib["uom"].lower() != "m":
        logger.info(
            f"groundlevel unit is unexpected {glp_xml.attrib['uom']}, m is expected"
        )
    if vert_datum.lower() != "nap":
        logger.info(f"datum has unexpected value {vert_datum}, NAP is expected")
    meta["ground_level"] = float(glp_xml.text)

    # tube properties
    tubes = GMW_c.findall("isgmw:monitoringTube", ns)
    for tube in tubes:
        # open only requested tube
        if int(tube.find("isgmw:tubeNumber", ns).text) == tube_nr:
            break

    if full_meta:
        if tube.find("isgmw:tubeType", ns) is not None:
            meta["tubeType"] = tube.find("isgmw:tubeType", ns).text
        if tube.find("isgmw:artesianWellCapPresent", ns) is not None:
            meta["artesianWellCapPresent"] = tube.find(
                "isgmw:artesianWellCapPresent", ns
            ).text
        if tube.find("isgmw:sedimentSumpPresent", ns) is not None:
            meta["sedimentSumpPresent"] = tube.find(
                "isgmw:sedimentSumpPresent", ns
            ).text
        if tube.find("isgmw:numberOfGeoOhmCables", ns) is not None:
            meta["numberOfGeoOhmCables"] = int(
                tube.find("isgmw:numberOfGeoOhmCables", ns).text
            )
        if tube.find("isgmw:tubeTopDiameter", ns) is not None:
            meta["tubeTopDiameter"] = float(tube.find("isgmw:tubeTopDiameter", ns).text)
            meta["tubeTopDiameter_unit"] = tube.find(
                "isgmw:tubeTopDiameter", ns
            ).attrib["uom"]
        if tube.find("isgmw:variableDiameter", ns) is not None:
            meta["variableDiameter"] = tube.find("isgmw:variableDiameter", ns).text
        if tube.find("isgmw:tubeStatus", ns) is not None:
            meta["tubeStatus"] = tube.find("isgmw:tubeStatus", ns).text
        if tube.find("isgmw:tubeTopPositioningMethod", ns) is not None:
            meta["tubeTopPositioningMethod"] = tube.find(
                "isgmw:tubeTopPositioningMethod", ns
            ).text
        if (
            tube.find("isgmw:materialUsed//gmwcommon:tubePackingMaterial", ns)
            is not None
        ):
            meta["tubePackingMaterial"] = tube.find(
                "isgmw:materialUsed//gmwcommon:tubePackingMaterial", ns
            ).text
        if tube.find("isgmw:materialUsed//gmwcommon:tubeMaterial", ns) is not None:
            meta["tubeMaterial"] = tube.find(
                "isgmw:materialUsed//gmwcommon:tubeMaterial", ns
            ).text
        if tube.find("isgmw:materialUsed//gmwcommon:glue", ns) is not None:
            meta["tubeGlue"] = tube.find("isgmw:materialUsed//gmwcommon:glue", ns).text
        if tube.find("isgmw:screen//isgmw:sockMaterial", ns) is not None:
            meta["sockMaterial"] = tube.find(
                "isgmw:screen//isgmw:sockMaterial", ns
            ).text

    # tube_top
    tube_top_xml = tube.find("isgmw:tubeTopPosition", ns)
    meta["tube_top"] = float(tube_top_xml.text)
    tube_top_unit = tube_top_xml.attrib["uom"]
    if tube_top_unit != "m":
        logger.info(f"tube_top unit is unexpected {tube_top_unit}, m expected")

    # some lenghts to calculate screen levels
    screenLength_xml = tube.find("isgmw:screen//isgmw:screenLength", ns)
    screenLength = float(screenLength_xml.text)
    if full_meta:
        meta["screenLength"] = screenLength
    screenLength_unit = screenLength_xml.attrib["uom"]
    if screenLength_unit != "m":
        logger.info(
            f"screenLength unit is unexpected {screenLength.attrib['uom']},m expected"
        )

    plainTubePartLength_xml = tube.find(
        "isgmw:plainTubePart//gmwcommon:plainTubePartLength", ns
    )
    plainTubePartLength = float(plainTubePartLength_xml.text)
    if full_meta:
        meta["plainTubePartLength"] = plainTubePartLength
    plainTubePartLength_unit = plainTubePartLength_xml.attrib["uom"]
    if plainTubePartLength_unit != "m":
        logger.info(
            "plainTubePartLength unit is unexpected"
            f"{plainTubePartLength_unit}, m expected"
        )

    if tube.find("isgmw:sedimentSumpPresent", ns).text.lower() in ["ja", "yes"]:
        sedimentSumpLength_xml = tube.find(
            "isgmw:sedimentSump//gmwcommon:sedimentSumpLength", ns
        )
        sedimentSumpLength = float(sedimentSumpLength_xml.text)
        sedimentSumpLength_unit = sedimentSumpLength_xml.attrib["uom"]
        if sedimentSumpLength_unit != "m":
            logger.info(
                "sedimentSumpLength unit is unexpected"
                f"{sedimentSumpLength_unit}, m expected"
            )

        if full_meta:
            meta["sedimentSumpLength"] = sedimentSumpLength

    if (
        ("tube_top" in meta)
        and ("plainTubePartLength" in locals())
        and (tube_top_unit == plainTubePartLength_unit)
    ):
        meta["screen_top"] = round(meta["tube_top"] - plainTubePartLength, 2)
    if (
        ("screen_top" in meta)
        and ("screenLength" in locals())
        and (plainTubePartLength_unit == screenLength_unit)
    ):
        meta["screen_bottom"] = round(meta["screen_top"] - screenLength, 2)
    if (
        ("screen_bottom" in meta)
        and ("sedimentSumpLength" in locals())
        and (screenLength_unit == sedimentSumpLength_unit)
    ):
        meta["tube_bottom"] = round(meta["screen_bottom"] - sedimentSumpLength, 2)

    meta["metadata_available"] = True

    return meta
