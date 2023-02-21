import datetime as dt
import logging
import os
import re
import tempfile
import warnings
from timeit import default_timer

import numpy as np
import pandas as pd
import requests
import zeep
from requests.exceptions import HTTPError
from zeep import Plugin
from zeep.plugins import HistoryPlugin
from zeep.wsa import WsAddressingPlugin

from ..util import unzip_file

logger = logging.getLogger(__name__)


def _read_dino_groundwater_header(f):
    line = f.readline()
    header = dict()
    while line not in ["\n", "", "\r\n"]:
        if "," in line:  # for csv from dinoloket
            propval = line.split(",")
            prop = propval[0]
            prop = prop.replace(":", "")
            prop = prop.strip()
            val = propval[1]
            if propval[2] != "":
                val = val + " " + propval[2].replace(":", "") + " " + propval[3]
            header[prop] = val
        else:  # for artdiver dino-csv
            propval = line.split(":")
            prop = propval[0]
            val = ":".join(propval[1:])
            header[prop] = val
        line = f.readline()

    return line, header


def _read_empty(f, line):
    while (line == "\n") or (line == "\r\n"):
        line = f.readline()
    return line


def _read_dino_groundwater_referencelvl(f, line):
    ref = {}
    while line not in ["\n", "", "\r\n"]:
        propval = line.split(",")
        prop = propval[0]
        prop = prop.replace(":", "")
        prop = prop.strip()
        if len(propval) > 1:
            val = propval[1]
            ref[prop] = val
        line = f.readline()
    return line, ref


def _read_dino_groundwater_metadata(f, line):
    _translate_dic_float = {
        "x-coordinaat": "x",
        "y-coordinaat": "y",
        "filternummer": "tube_nr",
    }
    _translate_dic_div_100 = {
        "meetpunt (cm t.o.v. nap)": "tube_top",
        "bovenkant filter (cm t.o.v. nap)": "screen_top",
        "onderkant filter (cm t.o.v. nap)": "screen_bottom",
        "maaiveld (cm t.o.v. nap)": "ground_level",
    }
    metalist = list()
    line = line.strip()
    properties = line.split(",")
    line = f.readline()

    while line not in ["\n", "", "\r\n"]:
        meta = dict()
        line = line.strip()
        values = line.split(",")
        for i in range(0, len(values)):
            meta[properties[i].lower()] = values[i]
        metalist.append(meta)
        line = f.readline()

    meta_ts = {}
    if metalist:
        # add time dependent metadata to meta_ts
        for i, meta in enumerate(metalist):
            meta_tsi = {}
            start_date = pd.to_datetime(meta.pop("startdatum"), dayfirst=True)
            # end_date = pd.to_datetime(meta.pop('einddatum'), dayfirst=True)
            meta_tsi["monitoring_well"] = meta["locatie"]
            for key in _translate_dic_float.keys():
                if meta[key] == "":
                    meta_tsi[_translate_dic_float[key]] = np.nan
                else:
                    meta_tsi[_translate_dic_float[key]] = float(meta[key])

            for key in _translate_dic_div_100.keys():
                if meta[key] == "":
                    meta_tsi[_translate_dic_div_100[key]] = np.nan
                else:
                    meta_tsi[_translate_dic_div_100[key]] = float(meta[key]) / 100.0
            if i == 0:
                for key in meta_tsi.keys():
                    meta_ts[key] = pd.Series(name=key, dtype=type(meta_tsi[key]))

            for key in meta_tsi.keys():
                meta_ts[key].loc[start_date] = meta_tsi[key]

        # remove series with non time variant metadata from meta_ts
        ts_keys = (
            ["monitoring_well"]
            + list(_translate_dic_float.values())
            + list(_translate_dic_div_100.values())
        )
        for key in ts_keys:
            unique_values = meta_ts[key].unique()
            if len(unique_values) == 1:
                meta_ts.pop(key)

        obs_att = meta_tsi.copy()
        obs_att["name"] = f'{obs_att["monitoring_well"]}-{int(obs_att["tube_nr"]):03d}'
        obs_att["metadata_available"] = True
    else:
        # no metadata
        obs_att = {}
        obs_att["monitoring_well"] = ""
        obs_att["tube_nr"] = np.nan
        obs_att["name"] = "unknown"
        obs_att["x"] = np.nan
        obs_att["y"] = np.nan
        obs_att["tube_top"] = np.nan
        obs_att["ground_level"] = np.nan
        obs_att["screen_top"] = np.nan
        obs_att["screen_bottom"] = np.nan
        obs_att["metadata_available"] = False

    return line, obs_att, meta_ts


def _read_dino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(",")
    while "" in titel:
        titel.remove("")

    if line != "":
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(
                f,
                header=None,
                names=titel,
                parse_dates=["peildatum"],
                index_col="peildatum",
                dayfirst=True,
                usecols=usecols,
            )
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_dino_groundwater_quality_txt(fname):
    """Read dino groundwater quality (grondwatersamenstelling) from a dinoloket
    txt file.

    Notes
    -----
    this function has not been tested thoroughly

    Parameters
    ----------
    fname : str
        path to txt file

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    logger.info("reading -> {}".format(os.path.split(fname)[-1]))
    with open(fname, "r") as f:
        # LOCATIE gegevens
        line = f.readline().rstrip("\n")
        assert line == "LOCATIE gegevens"

        strt_locatie = f.tell()
        # determine the number of rows
        nrows = -1  # the header does not count
        line = f.readline()
        while line not in ["\n", ""]:
            nrows += 1
            line = f.readline()
        eind_locatie = f.tell()
        # go back to where we were before
        f.seek(strt_locatie)
        # read the location-information
        locatie = pd.read_csv(f, sep="\t", nrows=nrows)
        # there is always only one location (change if this is not the case)
        assert nrows == 1

        locatie = locatie.squeeze()

        # KWALITEIT gegevens VLOEIBAAR
        f.seek(eind_locatie)
        line = f.readline().rstrip("\n")
        assert line == "KWALITEIT gegevens VLOEIBAAR"

        measurements = pd.read_csv(
            f,
            sep="\t",
            parse_dates=["Monster datum", "Analyse datum"],
            dayfirst=True,
            index_col="Monster datum",
        )

    meta = {
        "filename": fname,
        "source": "dino",
        "monitoring_well": locatie["NITG-nr"],
        "name": locatie["NITG-nr"],
        "x": locatie["X-coord"],
        "y": locatie["Y-coord"],
    }
    try:
        meta["ground_level"] = locatie["Maaiveldhoogte (m tov NAP)"]
    except KeyError:
        meta["ground_level"] = np.nan

    return measurements, meta


def read_dino_groundwater_csv(
    fname, to_mnap=True, read_series=True, remove_duplicates=False, keep_dup="last"
):
    """Read dino groundwater quantity data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    remove_duplicates : boolean, optional
        if True duplicate indices are removed. Default is False.
    keep_dup : str, optional
        indicate which duplicate indices should be kept, only used when
        remove_duplicates is True. Default is 'last'

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    logger.info(f"reading -> {os.path.split(fname)[-1]}")

    with open(fname, "r") as f:
        # read header
        line, header = _read_dino_groundwater_header(f)
        line = _read_empty(f, line)
        if not header:
            logger.warning(f"could not read header -> {fname}")

        # read reference level
        line, ref = _read_dino_groundwater_referencelvl(f, line)
        line = _read_empty(f, line)
        if not ref:
            logger.warning(f"could not read reference level -> {fname}")

        # read metadata
        line, meta, meta_ts = _read_dino_groundwater_metadata(f, line)
        line = _read_empty(f, line)
        if not meta["metadata_available"]:
            logger.warning(f"could not read metadata -> {fname}")
        meta["filename"] = fname
        meta["source"] = "dino"

        # read measurements
        if read_series:
            line, measurements = _read_dino_groundwater_measurements(f, line)
            if measurements is None:
                logger.warning(f"could not read measurements -> {fname}")
            elif measurements[~measurements.stand_cm_tov_nap.isna()].empty:
                logger.warning(f"no NAP measurements available -> {fname}")
            if to_mnap and measurements is not None:
                measurements.insert(
                    0, "stand_m_tov_nap", measurements["stand_cm_tov_nap"] / 100.0
                )
                meta["unit"] = "m NAP"
            elif not to_mnap:
                meta["unit"] = "cm NAP"
            if remove_duplicates:
                measurements = measurements[
                    ~measurements.index.duplicated(keep=keep_dup)
                ]

            # add time variant metadata to measurements
            for s in meta_ts.values():
                if measurements is None:
                    measurements = pd.DataFrame(data=s.copy(), columns=[s.name])
                else:
                    measurements = measurements.join(s, how="outer")
                measurements.loc[:, s.name] = measurements.loc[:, s.name].ffill()

        else:
            measurements = None

    return measurements, meta


def _read_artdino_groundwater_metadata(f, line):
    metalist = list()
    line = line.strip()
    properties = line.split(",")
    line = f.readline()
    while line not in ["\n", "", "\r\n"]:
        meta = dict()
        line = line.strip()
        values = line.split(",")
        for i in range(0, len(values)):
            meta[properties[i].lower()] = values[i]
        metalist.append(meta)
        line = f.readline()

    meta = {}
    if metalist:
        meta["monitoring_well"] = metalist[-1]["locatie"]
        meta["tube_nr"] = int(float(metalist[-1]["filternummer"]))
        meta["name"] = "-".join([meta["monitoring_well"], metalist[-1]["filternummer"]])
        meta["x"] = float(metalist[-1]["x-coordinaat"])
        meta["y"] = float(metalist[-1]["y-coordinaat"])
        meetpunt = metalist[-1]["meetpunt nap"]
        if meetpunt == "":
            meta["tube_top"] = np.nan
        else:
            meta["tube_top"] = float(meetpunt) / 100.0
        maaiveld = metalist[-1]["maaiveld nap"]
        if maaiveld == "":
            meta["ground_level"] = np.nan
        else:
            meta["ground_level"] = float(maaiveld) / 100
        bovenkant_filter = metalist[-1]["bovenkant filter"]
        if bovenkant_filter == "":
            meta["screen_top"] = np.nan
        else:
            meta["screen_top"] = float(bovenkant_filter) / 100
        onderkant_filter = metalist[-1]["onderkant filter"]
        if onderkant_filter == "":
            meta["screen_bottom"] = np.nan
        else:
            meta["screen_bottom"] = float(onderkant_filter) / 100
        meta["metadata_available"] = True
    else:
        # no metadata
        meta["monitoring_well"] = ""
        meta["tube_nr"] = np.nan
        meta["name"] = "unknown"
        meta["x"] = np.nan
        meta["y"] = np.nan
        meta["tube_top"] = np.nan
        meta["ground_level"] = np.nan
        meta["screen_top"] = np.nan
        meta["screen_bottom"] = np.nan
        meta["metadata_available"] = False

    return line, meta


def _read_artdino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(",")
    while "" in titel:
        titel.remove("")

    if line != "":
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(
                f,
                header=None,
                names=titel,
                parse_dates=["peil_datum_tijd"],
                index_col="peil_datum_tijd",
                dayfirst=True,
                usecols=usecols,
            )
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_artdino_groundwater_csv(fname, to_mnap=True, read_series=True):
    """Read dino groundwater quantity data from a CSV file as exported by
    ArtDiver.

    Parameters
    ----------
    fname : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    logger.info(f"reading -> {os.path.split(fname)[-1]}")

    with open(fname, "r") as f:
        # read header
        line, header = _read_dino_groundwater_header(f)
        line = _read_empty(f, line)
        if not header:
            logger.warning(f"could not read header -> {fname}")

        # read metadata
        line, meta = _read_artdino_groundwater_metadata(f, line)
        line = _read_empty(f, line)

        if not meta["metadata_available"]:
            logger.warning(f"could not read metadata -> {fname}")

        meta["filename"] = fname
        meta["source"] = "dino"

        # read measurements
        if read_series:
            line, measurements = _read_artdino_groundwater_measurements(f, line)
            if measurements is None:
                logger.warning(f"could not read measurements -> {fname}")
            elif measurements[~measurements["stand_cm_nap"].isna()].empty:
                logger.warning(f"no NAP measurements available -> {fname}")
            if to_mnap and measurements is not None:
                measurements["stand_m_tov_nap"] = measurements["stand_cm_nap"] / 100.0
                meta["unit"] = "m NAP"
            elif not to_mnap:
                meta["unit"] = "cm NAP"
        else:
            measurements = None

    return measurements, meta


def read_artdino_dir(
    dirname,
    ObsClass=None,
    subdir="csv",
    suffix=".csv",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    **kwargs,
):
    """Read Dino directory with point observations.

    to do:
    - Evt. nog verbeteren door meteen Dataframe te vullen op het moment dat een observatie
    wordt ingelezen. Nu wordt eerst alles ingelezen in een lijst en daar een dataframe van gemaakt.
    - aparte unzip functie maken en toch de juiste tijdelijke directory krijgen.

    Parameters
    ----------
    dirname : str
        directory name, can be a .zip file or the parent directory of subdir
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    subdir : str
        subdirectory of dirname with data files
    suffix : str
        suffix of files in subdir that will be read
    unpackdir : str
        destination directory of the unzipped file
    force_unpack : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        use date of the zipfile for the destination file
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_dino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations
    """

    # unzip dir
    if dirname.endswith(".zip"):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(
            zipf, dirname, force=force_unpack, preserve_datetime=preserve_datetime
        )

    # read filenames
    files = os.listdir(os.path.join(dirname, subdir))
    if suffix:
        files = [file for file in files if file.endswith(suffix)]

    if not files:
        raise FileNotFoundError(
            "no files were found in {} that end with {}".format(
                os.path.join(dirname, subdir), suffix
            )
        )

    # read individual files
    obs_list = []
    for _, file in enumerate(files):
        fname = os.path.join(dirname, subdir, file)
        obs = ObsClass.from_artdino_file(fname=fname, **kwargs)
        if obs.metadata_available and (not obs.empty):
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            logging.info(f"not added to collection -> {fname}")

    return obs_list


# %% DINO download methods


class RemoveWSA(Plugin):
    """Helper class to remove wsa tags from header in zeep XML post.

    As described by:
    https://github.com/mvantellingen/python-zeep/issues/330#issuecomment-321781637
    """

    def egress(self, envelope, http_headers, operation, binding_options):
        envelope.find("{http://schemas.xmlsoap.org/soap/envelope/}Header").clear()
        return envelope, http_headers


class DinoREST:
    """This object is a first draft for retrieving data from the Dinoloket REST
    API.

    Currently only tested for piezometers.
    """

    aliases = {
        "grondwatermonitoring": "lks_gwo_rd",
        "sondering": "lks_gso_rd",
        "boring": "lks_gbo_rd",
        "peilschaal": "lks_owo_rd",
    }

    def __init__(self, query_url=None, gwo_url=None):
        """DinoREST object for getting metadata from dinoloket.

        Parameters
        ----------
        query_url : str, optional
            url for passing queries to, by default None
        gwo_url : str, optional
            url for getting details for dino piezometers, by default None
        """
        if query_url is None:
            self.query_url = (
                "https://www.dinoloket.nl/arcgis/rest/services/dinoloket/"
                "{layer}/MapServer/0/query"
            )
        else:
            self.query_url = query_url

        if gwo_url is None:
            self.gwo_url = (
                "https://www.dinoloket.nl/javascriptmapviewer-web/"
                "rest/gws/gwo/details"
            )
        else:
            self.gwo_url = gwo_url

    def get(self, url, query):
        """GET method.

        Parameters
        ----------
        url : str
            url
        query : dict
            dictionary containing query parameters

        Returns
        -------
        requests.models.Response
            response from dinoloket REST API
        """
        try:
            response = requests.get(url, params=query)
            response.raise_for_status()
        except HTTPError as http_err:
            raise http_err
        except Exception as err:
            raise err
        else:
            return response

    def post(self, url, json_d):
        """POST method.

        Parameters
        ----------
        url : str
            url
        json_d : dict
            dictionary containing payload as JSON

        Returns
        -------
        requests.models.Response
            response from dinoloket REST API
        """
        try:
            response = requests.post(url, json=json_d)
            response.raise_for_status()
        except HTTPError as http_err:
            raise http_err
        except Exception as err:
            raise err
        else:
            return response

    def query_locations_by_extent(self, extent, layer="grondwatermonitoring"):
        """Get unique locations by extent.

        Parameters
        ----------
        extent : list of floats
            xmin, xmax, ymin, ymax
        layer : str
            dino layer to get data from, see DinoREST.aliases for more
            information on which layer codes to use

        Returns
        -------
        response : requests.models.Response
            query response
        """
        xmin, xmax, ymin, ymax = extent

        query = {
            "f": "pjson",
            "outFields": "*",
            "spatialRel": "esriSpatialRelContains",
            "geometryType": "esriGeometryEnvelope",
            "geometry": f"{xmin},{ymin},{xmax},{ymax}",
            "inSR": 28992,
            "outSR": 28992,
        }

        url = self.query_url.format(layer=self.aliases[layer])
        response = self.get(url, query)
        return response

    def get_gwo_details(self, locations):
        """Get metadata details for locations.

        Parameters
        ----------
        locations : list of str
            list of dinoId strings

        Returns
        -------
        requests.models.Response
            response from the Dino REST API
        """
        response = self.post(self.gwo_url, locations)
        return response

    def get_gwo_metadata(self, location, tube_nr, raw_response=False):
        """Get metadata details for locations.

        Parameters
        ----------
        locations : list of str
            list of dinoId strings

        Returns
        -------
        meta dictionary or requests.models.Response
        """
        response = self.get_gwo_details([location])
        # response = self.post(self.gwo_url, [location])

        if raw_response:
            return response
        else:
            data = response.json()
            if data == []:
                meta = {"metadata_available": False}
            else:
                meta = self._parse_json_single_gwo_filter(data[0], location, tube_nr)
            return meta

    def get_locations_by_extent(self, extent, layer="grondwatermonitoring", raw=False):
        """get locations within an extent.

        Parameters
        ----------
        extent : list of floats
            xmin, xmax, ymin, ymax
        layer : str, optional
            which dino layer to query, by default "grondwatermonitoring",
            see DinoREST.aliases for more information on which layer codes
            are available
        raw : bool, optional
            return raw JSON response, by default False

        Returns
        -------
        pandas.DataFrame
            DataFrame containing names and locations
        """
        r = self.query_locations_by_extent(extent, layer=layer)

        if raw:
            return r.json()
        else:
            return self._parse_json_features(r.json())

    def get_metadata_by_extent(self, extent, layer="grondwatermonitoring", raw=False):
        """Get metadata details for all objects in extent.

        Parameters
        ----------
        extent : list of floats
            xmin, xmax, ymin, ymax
        layer : str, optional
            which dino layer to query, by default "grondwatermonitoring",
            see DinoREST.aliases for more information on which layer codes
            are available
        raw : bool, optional
            return raw json response, by default False

        Returns
        -------
        pandas.DataFrame
            DataFrame containing metadata per object
        """
        # get locations
        r = self.query_locations_by_extent(extent, layer=layer)
        # get features as json
        loc_json = r.json()["features"]
        # loop through names to create list of locations
        nloc = len(loc_json)
        locations = [loc_json[i]["attributes"]["DINO_NR"] for i in range(nloc)]

        if not locations:
            return pd.DataFrame()
        # get details as json
        rdetails = self.get_gwo_details(locations)
        json_details = rdetails.json()

        if not json_details:
            return pd.DataFrame()

        # return raw or parsed as dataframe
        if raw:
            return json_details
        else:
            if layer == "grondwatermonitoring":
                return self._parse_json_gwo_details(json_details)
            else:
                raise NotImplementedError("Not yet implemented for '{}'!".format(layer))

    @staticmethod
    def _parse_json_features(js, field="features"):
        data = js[field]
        if len(data) == 0:
            return pd.DataFrame()
        df = pd.DataFrame()
        for idx, idat in enumerate(data):
            for name, value in idat["attributes"].items():
                df.loc[idx, name] = value
            x = idat["geometry"]["x"]
            y = idat["geometry"]["y"]
            df.loc[idx, "x"] = x
            df.loc[idx, "y"] = y
        df.set_index("DINO_NR", inplace=True)
        return df

    @staticmethod
    def _parse_json_single_gwo_filter(data, location, tube_nr):
        """parse metadata obtained with get_dino_piezometer_metadata.

        Parameters
        ----------
        data : dictionary
            response from dinoloket
        location : str
            location of piezometer, e.g. B57F0077
        tube_nr : str
            tube number should be a string of 3 characters, e.g. 002

        Returns
        -------
        meta : dictionary
            parsed dictionary
        """
        _translate_dic_location = {
            "xcoord": "x",
            "ycoord": "y",
            "surfaceElevation": "ground_level",
        }

        _translate_dic_filter = {
            "bottomHeightNap": "screen_bottom",
            "topHeightNap": "screen_top",
        }

        meta = {
            "monitoring_well": location,
            "name": "-".join([location, tube_nr]),
            "tube_nr": int(tube_nr),
            "metadata_available": True,
        }

        piezometers = data.pop("piezoMeters")
        levels = data.pop("levels")
        samples = data.pop("samples")

        if piezometers:
            try:
                meta.update(
                    next(
                        item for item in piezometers if item["piezometerNr"] == tube_nr
                    )
                )
            except StopIteration:
                logging.warning(
                    "no piezometer metadata in metadata json "
                    f"{location} with tube_nr {tube_nr}"
                )
        if levels:
            try:
                meta.update(
                    next(item for item in levels if item["piezometerNr"] == tube_nr)
                )
            except StopIteration:
                logger.warning(
                    "no level metadata in metadata json "
                    f"{location} with tube_nr {tube_nr}"
                )
        if samples:
            try:
                meta.update(
                    next(item for item in samples if item["piezometerNr"] == tube_nr)
                )
            except StopIteration:
                logger.warning(
                    "no sample metadata in metadata json "
                    f"{location} with tube_nr {tube_nr}"
                )

        meta.update(data)
        meta.pop("dinoId")
        for key, item in _translate_dic_location.items():
            meta[item] = meta.pop(key)

        if "clusterList" in meta.keys():
            if isinstance(meta["clusterList"], list):
                meta["clusterList"] = ";".join(meta["clusterList"])

        if "piezometerNr" not in meta.keys():
            return meta

        meta.pop("piezometerNr")

        for key, item in _translate_dic_filter.items():
            meta[item] = meta.pop(key)
            if meta[item] is None:
                meta[item] = np.nan

        return meta

    def _parse_json_gwo_details(self, json_details, field="levels"):
        """convert json response to dataframe.

        TODO:
        - figure out how to deal with this dictionary to convert to a
        dataframe. There are 3 attributes in the JSON that contain lists:
          - levels
          - samples
          - piezoMeters
          I have no idea what the difference between these is... but it makes
          converting the result to a DataFrame annoying. Currently, I choose one
          of these fields and then skip the rest.
        - currently clusterList field is joined with ';' and added as a string
        instead of the list in JSON.

        Parameters
        ----------
        data : JSON
            json response
        field : str, optional
            default field to use to populate the dataframe, by default "levels"

        Returns
        -------
        pandas.DataFrame
            DataFrame containing metadata for all dino piezometers
        """
        dflist = []
        location_list = []
        for i in range(len(json_details)):
            data = json_details[i]
            location = data["dinoId"]

            if location in location_list:
                raise ValueError
            else:
                location_list.append(location)

            piezometers = data["piezoMeters"]
            levels = data["levels"]
            samples = data["samples"]

            # check number of filters
            if piezometers is None:
                piezometers = []
            if levels is None:
                levels = []
            if samples is None:
                samples = []

            # merge piezometers, levels and samples
            filter_dic = {}
            for l in [piezometers, levels, samples]:
                for d in l:
                    if "piezometerNr" in d.keys():
                        if d["piezometerNr"] not in filter_dic.keys():
                            filter_dic[d["piezometerNr"]] = d
                        else:
                            filter_dic[d["piezometerNr"]].update(d)

            # get metadata for every filter
            for tube_nr in filter_dic.keys():
                meta = self._parse_json_single_gwo_filter(
                    data.copy(), location, str(tube_nr).zfill(3)
                )
                idf = pd.DataFrame(meta, index=[meta.pop("name")])
                dflist.append(idf)

        df = pd.concat(dflist, axis=0, sort=True)

        return df


def _read_dino_waterlvl_metadata(f, line):
    """read dino waterlevel metadata

    Parameters
    ----------
    f : text wrapper
    line : str
        line with meta dictionary keys
    meta_dic : dict (optional)
        dictionary with metadata

    Returns
    -------
    meta : dict
        dictionary with metadata

    """
    meta_keys = line.strip().split(",")
    meta_values = f.readline().strip().split(",")
    meta = {}
    for key, value in zip(meta_keys, meta_values):
        key = key.strip()
        if key in ["X-coordinaat", "Y-coordinaat"]:
            if key == "X-coordinaat":
                meta["x"] = float(value)
            elif key == "Y-coordinaat":
                meta["y"] = float(value)
        elif key == "Locatie":
            meta["monitoring_well"] = value
            meta["name"] = value

    return meta


def _read_dino_waterlvl_measurements(f, line):
    """

    Parameters
    ----------
    f : text wrapper
    line: str
        header of dataframe

    Returns
    -------
    measurements : pd.DataFrame
    """

    titel = line.strip().split(",")
    while "" in titel:
        titel.remove("")

    validator = np.lib._iotools.NameValidator()
    titel = [i.lower() for i in validator(titel)]
    usecols = range(0, len(titel))

    measurements = pd.read_csv(
        f,
        header=None,
        names=titel,
        parse_dates=["peildatum"],
        index_col="peildatum",
        dayfirst=True,
        usecols=usecols,
    )

    # measurements['Stand (m t.o.v. NAP)'] = measurements['Stand (cm t.o.v. NAP)'] /100.

    # measurements.set_index('Peildatum', inplace=True)

    return measurements


def read_dino_waterlvl_csv(fname, to_mnap=True, read_series=True):
    """Read dino waterlevel data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    """
    logging.info(f"reading -> {os.path.split(fname)[-1]}")

    p_meta = re.compile(
        "Locatie,Externe aanduiding,X-coordinaat,Y-coordinaat, Startdatum, Einddatum"
    )
    p_data = re.compile(r"Locatie,Peildatum,Stand \(cm t.o.v. NAP\),Bijzonderheid")

    with open(fname, "r") as f:
        line = f.readline()
        while line != "":
            line = f.readline()
            if p_meta.match(line):
                meta = _read_dino_waterlvl_metadata(f, line)
                if meta:
                    meta["metadata_available"] = True
                else:
                    meta["metadata_available"] = False
                meta["filename"] = fname
                meta["source"] = "dino"
            elif p_data.match(line):
                if read_series:
                    measurements = _read_dino_waterlvl_measurements(f, line)
                    if to_mnap and measurements is not None:
                        measurements["stand_m_tov_nap"] = (
                            measurements["stand_cm_tov_nap"] / 100.0
                        )
                        meta["unit"] = "m NAP"
                    elif not to_mnap:
                        meta["unit"] = "cm NAP"
                else:
                    measurements = None

                return measurements, meta


def read_dino_dir(
    dirname,
    ObsClass=None,
    subdir="Boormonsterprofiel_Geologisch booronderzoek",
    suffix=".txt",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    **kwargs,
):
    """Read Dino directory with point observations.

    to do:
    - Evt. nog verbeteren door meteen Dataframe te vullen op het moment dat een observatie
    wordt ingelezen. Nu wordt eerst alles ingelezen in een lijst en daar een dataframe van gemaakt.
    - aparte unzip functie maken en toch de juiste tijdelijke directory krijgen.

    Parameters
    ----------
    dirname : str
        directory name, can be a .zip file or the parent directory of subdir
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    subdir : str
        subdirectory of dirname with data files
    suffix : str
        suffix of files in subdir that will be read
    unpackdir : str
        destination directory of the unzipped file
    force_unpack : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        use date of the zipfile for the destination file
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_dino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations
    """

    # unzip dir
    if dirname.endswith(".zip"):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(
            zipf, dirname, force=force_unpack, preserve_datetime=preserve_datetime
        )

    # read filenames
    files = os.listdir(os.path.join(dirname, subdir))
    if suffix:
        files = [file for file in files if file.endswith(suffix)]

    if not files:
        raise FileNotFoundError(
            "no files were found in {} that end with {}".format(
                os.path.join(dirname, subdir), suffix
            )
        )

    # read individual files
    obs_list = []
    for _, file in enumerate(files):
        fname = os.path.join(dirname, subdir, file)
        obs = ObsClass.from_dino(fname=fname, **kwargs)
        if obs.metadata_available and (not obs.empty):
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            logging.info(f"not added to collection -> {fname}")

    return obs_list
