import os
import warnings
import re
import tempfile

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

import requests
from requests.exceptions import HTTPError

import json

import zeep
from zeep import Plugin
from zeep.wsa import WsAddressingPlugin
from zeep.plugins import HistoryPlugin

from ..util import unzip_file


# %% DINO groundwater CSV methods
def _read_dino_groundwater_header(f):
    line = f.readline()
    header = dict()
    while line not in ['\n', '', '\r\n']:
        if "," in line:  # for csv from dinoloket
            propval = line.split(',')
            prop = propval[0]
            prop = prop.replace(':', '')
            prop = prop.strip()
            val = propval[1]
            if propval[2] != '':
                val = val + ' ' + propval[2].replace(':', '') + ' ' + \
                    propval[3]
            header[prop] = val
        else:  # for artdiver dino-csv
            propval = line.split(":")
            prop = propval[0]
            val = ":".join(propval[1:])
            header[prop] = val
        line = f.readline()

    return line, header


def _read_empty(f, line):
    while (line == '\n') or (line == '\r\n'):
        line = f.readline()
    return line


def _read_dino_groundwater_referencelvl(f, line):
    ref = {}
    while line not in ['\n', '', '\r\n']:
        propval = line.split(',')
        prop = propval[0]
        prop = prop.replace(':', '')
        prop = prop.strip()
        if len(propval) > 1:
            val = propval[1]
            ref[prop] = val
        line = f.readline()
    return line, ref


def _read_dino_groundwater_metadata(f, line):
    metalist = list()
    line = line.strip()
    properties = line.split(',')
    line = f.readline()
    while line not in ['\n', '', '\r\n']:
        meta = dict()
        line = line.strip()
        values = line.split(',')
        for i in range(0, len(values)):
            meta[properties[i].lower()] = values[i]
        metalist.append(meta)
        line = f.readline()

    meta = {}
    if metalist:
        meta["locatie"] = metalist[-1]['locatie']
        meta["filternr"] = int(float(metalist[-1]['filternummer']))
        meta["name"] = '-'.join([meta["locatie"],
                                 metalist[-1]['filternummer']])
        meta["x"] = float(metalist[-1]['x-coordinaat'])
        meta["y"] = float(metalist[-1]['y-coordinaat'])
        meetpunt = metalist[-1]['meetpunt (cm t.o.v. nap)']
        if meetpunt == '':
            meta["meetpunt"] = np.nan
        else:
            meta["meetpunt"] = float(meetpunt) / 100.
        maaiveld = metalist[-1]['maaiveld (cm t.o.v. nap)']
        if maaiveld == '':
            meta["maaiveld"] = np.nan
        else:
            meta["maaiveld"] = float(maaiveld) / 100
        bovenkant_filter = metalist[-1]['bovenkant filter (cm t.o.v. nap)']
        if bovenkant_filter == '':
            meta["bovenkant_filter"] = np.nan
        else:
            meta["bovenkant_filter"] = float(bovenkant_filter) / 100
        onderkant_filter = metalist[-1][
            'onderkant filter (cm t.o.v. nap)']
        if onderkant_filter == '':
            meta["onderkant_filter"] = np.nan
        else:
            meta["onderkant_filter"] = float(onderkant_filter) / 100
        meta["metadata_available"] = True
    else:
        # de metadata is leeg
        meta["locatie"] = ''
        meta["filternr"] = np.nan
        meta["name"] = 'unknown'
        meta["x"] = np.nan
        meta["y"] = np.nan
        meta["meetpunt"] = np.nan
        meta["maaiveld"] = np.nan
        meta["bovenkant_filter"] = np.nan
        meta["onderkant_filter"] = np.nan
        meta["metadata_available"] = False

    return line, meta


def _read_dino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(',')
    while '' in titel:
        titel.remove('')

    if line != '':
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(f, header=None, names=titel,
                                       parse_dates=['peildatum'],
                                       index_col='peildatum',
                                       dayfirst=True,
                                       usecols=usecols)
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_dino_groundwater_quality_txt(fname, verbose=False):
    """Read dino groundwater quality (grondwatersamenstelling) from a dinoloket
    txt file

    Notes
    -----
    this function has not been tested thoroughly

    Parameters
    ----------
    fname : str
        path to txt file
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """

    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))
    with open(fname, 'r') as f:

        # LOCATIE gegevens
        line = f.readline().rstrip('\n')
        assert line == 'LOCATIE gegevens'

        strt_locatie = f.tell()
        # determine the number of rows
        nrows = -1  # the header does not count
        line = f.readline()
        while line not in ['\n', '']:
            nrows += 1
            line = f.readline()
        eind_locatie = f.tell()
        # go back to where we were before
        f.seek(strt_locatie)
        # read the location-information
        locatie = pd.read_csv(f, sep='\t', nrows=nrows)
        # there is always only one location (change if this is not the case)
        assert nrows == 1

        locatie = locatie.squeeze()

        # KWALITEIT gegevens VLOEIBAAR
        f.seek(eind_locatie)
        line = f.readline().rstrip('\n')
        assert line == 'KWALITEIT gegevens VLOEIBAAR'

        measurements = pd.read_csv(f, sep='\t', parse_dates=['Monster datum', 'Analyse datum'],
                                   dayfirst=True, index_col='Monster datum')

    meta = {}
    meta['filename'] = fname
    meta['locatie'] = locatie['NITG-nr']
    meta['name'] = locatie['NITG-nr']
    meta['x'] = locatie['X-coord']
    meta['y'] = locatie['Y-coord']
    try:
        meta['maaiveld'] = locatie['Maaiveldhoogte (m tov NAP)']
    except KeyError:
        meta['maaiveld'] = np.nan

    return measurements, meta


def read_dino_groundwater_csv(fname, to_mnap=True,
                              read_series=True, verbose=False,):
    """Read dino groundwater quantity data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))

    with open(fname, 'r') as f:
        # read header
        line, header = _read_dino_groundwater_header(f)
        line = _read_empty(f, line)
        if verbose and (header == dict()):
            print('could not read header -> {}'.format(fname))

        # read reference level
        line, ref = _read_dino_groundwater_referencelvl(f, line)
        line = _read_empty(f, line)
        if verbose and (ref == dict()):
            print('could not read reference level -> {}'.format(fname))

        # read metadata
        line, meta = _read_dino_groundwater_metadata(f, line)
        line = _read_empty(f, line)
        if verbose and (meta['metadata_available'] == False):
            print('could not read metadata -> {}'.format(fname))
        meta['filename'] = fname

        # read measurements
        if read_series:
            line, measurements = _read_dino_groundwater_measurements(f, line)
            if verbose and (measurements is None):
                print('could not read measurements -> {}'.format(fname))
            elif verbose and measurements[~measurements.stand_cm_tov_nap.isna()].empty:
                print('no NAP measurements available -> {}'.format(fname))
            if to_mnap and measurements is not None:
                measurements['stand_m_tov_nap'] = measurements['stand_cm_tov_nap'] / 100.
        else:
            measurements = None

    return measurements, meta


# %% DINO CSV as exported by ArtDiver methods

def _read_artdino_groundwater_metadata(f, line):
    metalist = list()
    line = line.strip()
    properties = line.split(',')
    line = f.readline()
    while line not in ['\n', '', '\r\n']:
        meta = dict()
        line = line.strip()
        values = line.split(',')
        for i in range(0, len(values)):
            meta[properties[i].lower()] = values[i]
        metalist.append(meta)
        line = f.readline()

    meta = {}
    if metalist:
        meta["locatie"] = metalist[-1]['locatie']
        meta["filternr"] = int(float(metalist[-1]['filternummer']))
        meta["name"] = '-'.join([meta["locatie"],
                                 metalist[-1]['filternummer']])
        meta["x"] = float(metalist[-1]['x-coordinaat'])
        meta["y"] = float(metalist[-1]['y-coordinaat'])
        meetpunt = metalist[-1]['meetpunt nap']
        if meetpunt == '':
            meta["meetpunt"] = np.nan
        else:
            meta["meetpunt"] = float(meetpunt) / 100.
        maaiveld = metalist[-1]['maaiveld nap']
        if maaiveld == '':
            meta["maaiveld"] = np.nan
        else:
            meta["maaiveld"] = float(maaiveld) / 100
        bovenkant_filter = metalist[-1]['bovenkant filter']
        if bovenkant_filter == '':
            meta["bovenkant_filter"] = np.nan
        else:
            meta["bovenkant_filter"] = float(bovenkant_filter) / 100
        onderkant_filter = metalist[-1][
            'onderkant filter']
        if onderkant_filter == '':
            meta["onderkant_filter"] = np.nan
        else:
            meta["onderkant_filter"] = float(onderkant_filter) / 100
        meta["metadata_available"] = True
    else:
        # de metadata is leeg
        meta["locatie"] = ''
        meta["filternr"] = np.nan
        meta["name"] = 'unknown'
        meta["x"] = np.nan
        meta["y"] = np.nan
        meta["meetpunt"] = np.nan
        meta["maaiveld"] = np.nan
        meta["bovenkant_filter"] = np.nan
        meta["onderkant_filter"] = np.nan
        meta["metadata_available"] = False

    return line, meta


def _read_artdino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(',')
    while '' in titel:
        titel.remove('')

    if line != '':
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(f, header=None, names=titel,
                                       parse_dates=['peil_datum_tijd'],
                                       index_col='peil_datum_tijd',
                                       dayfirst=True,
                                       usecols=usecols)
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_artdino_groundwater_csv(fname, to_mnap=True,
                                 read_series=True, verbose=False,):
    """Read dino groundwater quantity data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))

    with open(fname, 'r') as f:
        # read header
        line, header = _read_dino_groundwater_header(f)
        line = _read_empty(f, line)
        if verbose and (header == dict()):
            print('could not read header -> {}'.format(fname))

        # read metadata
        line, meta = _read_artdino_groundwater_metadata(f, line)
        line = _read_empty(f, line)
        if verbose and (meta['metadata_available'] == False):
            print('could not read metadata -> {}'.format(fname))
        meta['filename'] = fname

        # read measurements
        if read_series:
            line, measurements = _read_artdino_groundwater_measurements(
                f, line)
            if verbose and (measurements is None):
                print('could not read measurements -> {}'.format(fname))
            elif verbose and measurements[~measurements["stand_cm_nap"].isna()].empty:
                print('no NAP measurements available -> {}'.format(fname))
            if to_mnap and measurements is not None:
                measurements['stand_m_tov_nap'] = measurements['stand_cm_nap'] / 100.
        else:
            measurements = None

    return measurements, meta


def read_artdino_dir(dirname, ObsClass=None,
                     subdir='csv', suffix='.csv',
                     unpackdir=None, force_unpack=False, preserve_datetime=False,
                     verbose=False, keep_all_obs=True, **kwargs
                     ):
    '''Read Dino directory with point observations

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
    verbose : boolean, optional
        Print additional information to the screen (default is False)
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_dino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations

    '''

    # unzip dir
    if dirname.endswith('.zip'):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(zipf, dirname, force=force_unpack,
                   preserve_datetime=preserve_datetime)

    # read filenames
    files = os.listdir(os.path.join(dirname, subdir))
    if suffix:
        files = [file for file in files if file.endswith(suffix)]

    if not files:
        raise FileNotFoundError('no files were found in {} that end with {}'.format(
            os.path.join(dirname, subdir), suffix))

    # read individual files
    obs_list = []
    for _, file in enumerate(files):
        fname = os.path.join(dirname, subdir, file)
        obs = ObsClass.from_artdino_file(
            fname=fname, verbose=verbose, **kwargs)
        if obs.metadata_available and (not obs.empty):
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            if verbose:
                print('not added to collection -> {}'.format(fname))

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df


# %% DINO download methods  

class RemoveWSA(Plugin):
    """Helper class to remove wsa tags
    from header in zeep XML post

    as described by: https://github.com/mvantellingen/python-zeep/issues/330#issuecomment-321781637

    """

    def egress(self, envelope, http_headers, operation, binding_options):
        envelope.find(
            '{http://schemas.xmlsoap.org/soap/envelope/}Header').clear()
        return envelope, http_headers


class DinoREST:
    """This object is a first draft for retrieving data from the
    Dinoloket REST API. Currently only tested for piezometers.

    """

    aliases = {
        "grondwatermonitoring": "lks_gwo_rd",
        "sondering": "lks_gso_rd",
        "boring": "lks_gbo_rd",
        "peilschaal": "lks_owo_rd"
    }

    def __init__(self, query_url=None, gwo_url=None):
        """DinoREST object for getting metadata from dinoloket

        Parameters
        ----------
        query_url : str, optional
            url for passing queries to, by default None
        gwo_url : str, optional
            url for getting details for dino piezometers, by default None
        """
        if query_url is None:
            self.query_url = ("https://www.dinoloket.nl/arcgis/rest/services/dinoloket/"
                              "{layer}/MapServer/0/query")
        else:
            self.query_url = query_url

        if gwo_url is None:
            self.gwo_url = ("https://www.dinoloket.nl/javascriptmapviewer-web/"
                            "rest/gws/gwo/details")
        else:
            self.gwo_url = gwo_url

    def get(self, url, query):
        """GET method

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

    def post(self, url, json):
        """POST method

        Parameters
        ----------
        url : str
            url
        json : dict
            dictionary containing payload as JSON

        Returns
        -------
        requests.models.Response
            response from dinoloket REST API

        """
        try:
            response = requests.post(url, json=json)
            response.raise_for_status()
        except HTTPError as http_err:
            raise http_err
        except Exception as err:
            raise err
        else:
            return response

    def query_locations_by_extent(self, extent, layer="grondwatermonitoring"):
        """Get unique locations by extent

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
            "outSR": 28992
        }

        url = self.query_url.format(layer=self.aliases[layer])
        response = self.get(url, query)
        return response
    
    def get_gwo_details(self, locations):
        """Get metadata details for locations

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

    def get_gwo_metadata(self, location, filternr, raw_response=False,
                         verbose=False):
        """Get metadata details for locations

        Parameters
        ----------
        locations : list of str
            list of dinoId strings

        Returns
        -------
        meta dictionary or requests.models.Response
            

        """
        response = self.get_gwo_details([location])
        #response = self.post(self.gwo_url, [location])
        
        if raw_response:
            return response
        else:
            data = json.loads(response.content)
            if data == []:
                meta = {'metadata_available':False}
            else:
                meta = self._parse_json_single_gwo_filter(data[0], location, 
                                                          filternr, 
                                                          verbose=verbose)
            return meta
    

    def get_locations_by_extent(self, extent, layer="grondwatermonitoring",
                                raw=False):
        """get locations within an extent

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
        
    def get_metadata_by_extent(self, extent, layer="grondwatermonitoring",
                               raw=False):
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
                raise NotImplementedError(
                    "Not yet implemented for '{}'!".format(layer))


    def _parse_json_features(self, js, field="features"):
        data = js[field]
        if len(data) == 0:
            return pd.DataFrame()
        gdf = gpd.GeoDataFrame()
        geometry = []
        for idx, idat in enumerate(data):
            for name, value in idat["attributes"].items():
                gdf.loc[idx, name] = value
            x = idat["geometry"]["x"]
            y = idat["geometry"]["y"]
            geometry.append(Point(x, y))
            gdf.loc[idx, "x"] = x
            gdf.loc[idx, "y"] = y
        gdf["geometry"] = geometry
        gdf.set_index("DINO_NR", inplace=True)
        return gdf

    def _parse_json_single_gwo_filter(self, data, location, filternr, verbose=False):
        """ parse metadata obtained with get_dino_piezometer_metadata
        
        Parameters
        ----------
        data : dictionary
            response from dinoloket
        location : str
            location of piezometer, e.g. B57F0077
        filternr : str
            filter number should be a string of 3 characters, e.g. 002
    
        Returns
        -------
        meta : dictionary
            parsed dictionary
    
        """

        meta = {
            "locatie": location,
            "name": '-'.join([location, filternr]),
            "filternr": int(filternr),
            "metadata_available": True
        }
        
        piezometers = data.pop('piezoMeters')
        levels = data.pop('levels')
        samples = data.pop('samples')
        
        
        if piezometers:
            try:
                meta.update(next(item for item in piezometers if item["piezometerNr"] == filternr))
            except StopIteration:
                if verbose:
                    print(f'cannot find piezometer metadata for location {location} and filternr {filternr}')
        if levels:
            try:
                meta.update(next(item for item in levels if item["piezometerNr"] == filternr))
            except StopIteration:
                if verbose:
                    print(f'cannot find level metadata for location {location} and filternr {filternr}')
        if samples:
            try:
                meta.update(next(item for item in samples if item["piezometerNr"] == filternr))
            except StopIteration:
                if verbose:
                    print(f'cannot find sample metadata for location {location} and filternr {filternr}')
        
        meta.update(data)
        meta.pop('piezometerNr')
        meta.pop('dinoId')
        
        translate_dic = {'bottomHeightNap':'onderkant_filter',
                         'topHeightNap':'bovenkant_filter',
                         'xcoord':'x',
                         'ycoord':'y',
                         'surfaceElevation':'maaiveld'}
        
        for key, item in translate_dic.items():
            meta[item] = meta.pop(key)
            
        if 'clusterList' in meta.keys():
            if isinstance(meta['clusterList'], list):
                meta['clusterList'] = ";".join(meta['clusterList'])
            
        return meta
    
    def _parse_json_gwo_details(self, json_details, field="levels",
                                verbose=False):
        """convert json response to dataframe

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
        for i in range(len(json_details)):
            data = json_details[i]
            location = data['dinoId']
            
            piezometers = data['piezoMeters']
            levels = data['levels']
            samples = data['samples']
            
            #check number of filters
            if piezometers is None:
                piezometers = {}
            if levels is None:
                levels = {}
            if samples is None:
                samples = {}
            
            no_filters = max(len(piezometers), len(levels), len(samples))
            
            
            for filternr in range(1,no_filters+1):
                meta = self._parse_json_single_gwo_filter(data.copy(), location, 
                                                          str(filternr).zfill(3), 
                                                          verbose=verbose)
                idf = pd.DataFrame(meta, index=[meta.pop('name')])
                dflist.append(idf)
                
        df = pd.concat(dflist, axis=0, sort=True)
       
        return df

class DinoWSDL:
    """Class for querying DINOLoket. Currently available methods are:
        - findMeetreeks: get timeseries for piezometer
        - findWaarnemingGegevens: get some simple all time statistics
          for a piezometer
        - findTechnischeGegevens: get metadata for a piezometer
        - findGrondwaterStatistiek: get stats for 14th and 28th day of
          each month for a given period

    Returns
    -------
    DinoWSDL :
        returns an object that can be used to query the DINO Webservice.

    """

    def __init__(self, wsdl="http://www.dinoservices.nl/gwservices/gws-v11?wsdl"):
        """Initialize DinoWSDL object that can be used to query DINO Webservice

        Parameters
        ----------
        wsdl : str, optional
            wsdl url, by default "http://www.dinoservices.nl/gwservices/gws-v11?wsdl"
        """

        # Create some plugins, some currently unused but left here as a reminder.
        history = HistoryPlugin()
        wsa = WsAddressingPlugin()
        rwsa = RemoveWSA()

        # Configure the client
        self.client = zeep.Client(
            wsdl=wsdl,
            plugins=[history, wsa, rwsa])

        # Set prefix
        self.client.set_ns_prefix("v11", "http://v11.ws.gws.dino.tno.nl/")    
        
    def findMeetreeks(self, location, filternr, tmin, tmax, unit="NAP",
                      raw_response=False):
        """Get a timeseries for a piezometer.

        Parameters
        ----------
        location : str
            location str of the piezometer, i.e. B57F0077
        filternr : str, int, float
            filter number, is converted to str, i.e. 004
        tmin : str or pandas.Timestamp
            start date in format YYYY-MM-DD (will be converted if Timestamp)
        tmax : str or pandas.Timestamp
            end date in format YYYY-MM-DD (will be converted if Timestamp)
        unit : str, optional
            unit in which timeseries is returned, by default "NAP"
        raw_response : bool, optional
            if True, return the zeep parsed XML response, by default
            False which will return a DataFrame

        Returns
        -------
        df or response
            returns DataFrame or zeep parsed XML response

        Note
        ----
        - method sometimes returns data

        """
        if isinstance(tmin,pd.Timestamp):
            tmin = tmin.strftime('%Y-%m-%d')

        if isinstance(tmax,pd.Timestamp):
            tmax = tmax.strftime('%Y-%m-%d')

        assert unit in ["NAP", "SFL",
                        "MP"], "'unit' must be one of 'NAP', 'SFL', 'MP'!"

        data = {"WELL_NITG_NR": location,
                "WELL_TUBE_NR": filternr,
                "START_DATE": tmin,
                "END_DATE": tmax,
                "UNIT": unit}

        try:
            r = self.client.service.findMeetreeks(**data)
        except Exception as e:
            raise e

        if raw_response:
            return r
        else:
            
            df = self._parse_grondwaterstand(r, f'stand_m_tov_{unit.lower()}')
            
            return df
    
    def findTechnischeGegevens(self, location, filter_nr, raw_response=False):
        """Get metadata for a piezometer.

        Parameters
        ----------
        location : str
            location str of the piezometer, i.e. B57F0077
        filter_nr : str, int, float
            filter number, is converted to str, i.e. 004
        raw_response : bool, optional
            if True, return the zeep parsed XML response, by default
            False which will return a DataFrame

        Returns
        -------
        df or response
            returns DataFrame or zeep parsed XML response

        """

        if isinstance(filter_nr, float) or isinstance(filter_nr, int):
            filter_nr = "{0:03g}".format(filter_nr)

        data = {"WELL_NITG_NR": location,
                "WELL_TUBE_NR": filter_nr}

        try:
            r = self.client.service.findTechnischeGegevens(**data)
        except Exception as e:
            raise e

        if raw_response:
            return r
        else:
            meta = self._parse_technische_gegevens(r)
            return meta
    
    @staticmethod
    def _parse_grondwaterstand(r, column_name='stand_m_tov_nap'):
        """Parse the response from findMeetreeks

        Parameters
        ----------
        r : zeep.objects.GrondwaterStanden
            the response of findMeetreeks()

        Returns
        -------
        pandas.DataFrame
            pandas.DataFrame containing the data

        """
        if isinstance(r, list):
            gws = r[0]

        index = [i["DATE"] for i in gws.LEVELS]
        values = [i["LEVEL"] for i in gws.LEVELS]
        remarks = [i["REMARK"] for i in gws.LEVELS]
        df = pd.DataFrame(index=index, data={column_name:values,
                                             'remarks':remarks})
        # convert from cm to m
        df.loc[:,column_name] = df.loc[:,column_name]/100
                
        return df
    
    @staticmethod
    def _parse_technische_gegevens(response):
        """Parse the response from findTechnischeGegevens

        Parameters
        ----------
        response : zeep.objects.TechnischeGegevens
            the response of findTechnischeGegevens()

        Returns
        -------
        pandas.DataFrame
            pandas.DataFrame containing the data

        """
        if isinstance(response, list):
            response = response[0]
          
        locatie = response.WELL_NITG_NR
        filternr = response.WELL_TUBE_NR
        meta = {
            "locatie": locatie,
            "name": '-'.join([locatie, filternr]),
            "filternr": int(filternr)
        }

        if response.SFL_HEIGHT is not None:
            meta["maaiveld"] = response.SFL_HEIGHT / 100.
        else:
            meta["maaiveld"] = np.nan

        if response.BOTTOM_FILTER is not None:
            meta["onderkant_filter"] = response.BOTTOM_FILTER / 100.
        else:
            meta["onderkant_filter"] = np.nan

        if response.FILTER_LENGTH is not None:
            meta["filterlengte"] = response.FILTER_LENGTH / 100.
        else:
            meta["filterlengte"] = np.nan

        meta["bovenkant_filter"] = meta["onderkant_filter"] + meta["filterlengte"]

        return meta


  

def download_dino_groundwater(location, filternr, tmin, tmax, 
                              verbose=False,
                             **kwargs):
    """ download measurements and metadata from a dino groundwater 
    observation well

    Parameters
    ----------
    location : str
        location str of the piezometer, i.e. B57F0077
    filternr : str, int, float
        filter number, is converted to str, i.e. 004
    tmin : str or pandas.Timestamp
        start date in format YYYY-MM-DD (will be converted if Timestamp)
    tmax : str or pandas.Timestamp
        end date in format YYYY-MM-DD (will be converted if Timestamp)
    verbose : 
    kwargs : key-word arguments
            these arguments are passed to dino.findMeetreeks functie

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    
    if isinstance(filternr, float) or isinstance(filternr, int):
        filternr = "{0:03g}".format(filternr)
    

    # download data from dino
    dino = DinoWSDL()

    # measurements
    measurements = dino.findMeetreeks(location, filternr, tmin, tmax,
                                      **kwargs)

    
    # old metadata method
    #meta = dino.findTechnischeGegevens(location, filternr)
    
    # new metadata method
    dinorest = DinoREST()
    meta = dinorest.get_gwo_metadata(location, filternr, verbose=verbose)
    
    return measurements, meta


def get_dino_locations(extent=None, bbox=None, layer='grondwatermonitoring',
                       response_fname=None, split_dx=None):
    """Get a GeoDataFrame with dinoloket-locations from a wfs

    Parameters
    ----------
    extent : list, tuple or numpy-array (user must specify extent or bbox)
        The extent, in RD-coordinates, for which you want to retreive locations
        [xmin, xmax, ymin, ymax]
    bbox : list, tuple or numpy-array (user must specify extent or bbox)
        The bounding box, in RD-coordinates, for which you want to retreive locations
        [xmin, ymin, xmax, ymax]
    layer : str, optional
        The type of timeseries (grondwatermonitoring or boring)

    Returns
    -------
    df : pandas DataFrame or a tuple of GeoDataFrame

    if kind=None
        a geodataframe with the properties of the locations

    """

    # get extent
    if extent is None and bbox is None:
        raise (Exception('Either bbox or extent needs to be specified'))
    if extent is not None and bbox is not None:
        raise (Exception('bbox and extent cannot both be specified'))
    if bbox is not None:
        extent = [bbox[0], bbox[2], bbox[1], bbox[3]]
        if split_dx is not None:
            gdf = split_extent(get_dino_locations, extent,
                               split_dx, layer=layer)
            return gdf
    if isinstance(extent, np.ndarray):
        extent = extent.tolist()

    # use dino REST api
    if layer == 'grondwatermonitoring' or layer == 'boring':
        pass
    else:
        raise(Exception('Unknown layer: {}'.format(layer)))
    dinorest = DinoREST()
    if layer == 'boring':
        df = dinorest.get_locations_by_extent(extent, layer=layer)
    elif layer == 'grondwatermonitoring':
        df = dinorest.get_metadata_by_extent(extent, layer=layer)
    else:
        raise NotImplementedError(
            "Not yet implemented for '{}'!".format(layer))

    # check if there are any measurements
    if df.empty:
        warnings.warn('no observation points were found')
        return df

    return df

def split_extent(method, extent0, dx=1000, dy=None, **kwargs):
    if dy is None:
        dy = dx
    ncol = int(np.ceil((extent0[1] - extent0[0]) / dx))
    nrow = int(np.ceil((extent0[3] - extent0[2]) / dy))
    x = np.linspace(extent0[0], extent0[1], ncol)
    y = np.linspace(extent0[2], extent0[3], nrow)
    for ix in range(len(x) - 1):
        for iy in range(len(y) - 1):
            print('Downloading part {} of {}'.format(
                ix * iy + iy + 1, len(x) * len(y)))
            extent = [x[ix], x[ix + 1], y[iy], y[iy + 1]]
            gdft = method(extent=extent, **kwargs)
            if ix == 0 and iy == 0:
                gdf = gdft
            else:
                gdf = pd.concat((gdf, gdft))


def download_dino_within_extent(extent=None, bbox=None, ObsClass=None,
                                layer='grondwatermonitoring',
                                tmin="1900-01-01", tmax="2040-01-01",
                                zmin=None, zmax=None, unit="NAP",
                                keep_all_obs=True,
                                verbose=False):
    """Download DINO data within a certain extent (or a bounding box)


    Parameters
    ----------
    extent : list, tuple or numpy-array (user must specify extent or bbox)
        The extent, in RD-coordinates, for which you want to retreive locations
        [xmin, xmax, ymin, ymax]
    bbox : list, tuple or numpy-array (user must specify extent or bbox)
        The bounding box, in RD-coordinates, for which you want to retreive locations
        [xmin, ymin, xmax, ymax]
    ObsClass : type
        class of the observations, so far only GroundwaterObs is supported
    layer : str
        The layer of timeseries ('grondwatermonitoring' of 'boring')
    tmin : str or pandas.Timestamp
        start date in format YYYY-MM-DD (will be converted if Timestamp)
    tmax : str or pandas.Timestamp
        end date in format YYYY-MM-DD (will be converted if Timestamp)
    zmin : float, optional
        minimum filter height in m NAP for which piezometers are returned
    zmax : float, optional
        maximum filter height in m NAP for which piezometers are returned
    unit : str, optional
        unit in which the timeseries are returned, by default "NAP"
    keep_all_obs : boolean, optional
            add all observation points to the collection, even without data or
            metadata
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations

    """

    # read locations
    gdf_loc = get_dino_locations(extent=extent, bbox=bbox, layer=layer)

    if verbose:
        print('\ndownload {} data from dino within:\n'
              '- extent: {} or\n'
              '- bbox: {}'.format(layer, extent, bbox))

    if gdf_loc.empty:
        return pd.DataFrame()

    gdf_loc["startDate"] = gdf_loc.startDate.astype('datetime64[ns]')
    gdf_loc["endDate"] = gdf_loc.endDate.astype('datetime64[ns]')
    
    # slice by properties
    if tmin is not None:
        tmin = pd.to_datetime(tmin)
        mask = gdf_loc.endDate > tmin
        gdf_loc = gdf_loc.loc[mask]
    if tmax is not None:
        tmax = pd.to_datetime(tmax)
        mask = gdf_loc.startDate < tmax
        gdf_loc = gdf_loc.loc[mask]
    if zmin is not None:
        mask = gdf_loc.topHeightNap >= zmin
        gdf_loc = gdf_loc.loc[mask]
    if zmax is not None:
        mask = gdf_loc.bottomHeightNap <= zmax
        gdf_loc = gdf_loc.loc[mask]

    if gdf_loc.empty:
        return pd.DataFrame()

    # read measurements
    obs_list = []
    for index, loc in gdf_loc.iterrows():
        if tmin is None:
            tmin_t = loc.startDate
        else:
            tmin_t = tmin
        if tmax is None:
            tmax_t = loc.endDate
        else:
            tmax_t = tmax

        if verbose:
            print('reading -> {}'.format(index))

        o = ObsClass.from_dino_server(location=loc.locatie,
                                      filternr=float(loc.filternr),
                                      tmin=tmin_t,
                                      tmax=tmax_t,
                                      unit=unit)
        
        if o.metadata_available and (not o.empty):
            obs_list.append(o)
        elif keep_all_obs:
            obs_list.append(o)
        else:
            if verbose:
                print('not added to collection -> {}'.format(o.name))

    # create dataframe
    if len(obs_list)>0:
        obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                              columns=obs_list[0].to_collection_dict().keys())
        obs_df.set_index('name', inplace=True)
        
        return obs_df
    
    else:
        return pd.DataFrame()

# %% DINO waterlevel observations
def _read_dino_waterlvl_metadata(f, line):
    '''

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

    '''
    meta_keys = line.strip().split(',')
    meta_values = f.readline().strip().split(',')
    meta = {}
    for key, value in zip(meta_keys, meta_values):
        key = key.strip()
        if key in ['X-coordinaat', 'Y-coordinaat']:
            if key == 'X-coordinaat':
                meta['x'] = float(value)
            elif key == 'Y-coordinaat':
                meta['y'] = float(value)
        elif key == 'Locatie':
            meta['locatie'] = value
            meta['name'] = value

    return meta


def _read_dino_waterlvl_measurements(f, line):
    '''

    Parameters
    ----------
    f : text wrapper
    line: str
        header of dataframe

    Returns
    -------
    measurements : pd.DataFrame
    '''

    titel = line.strip().split(',')
    while '' in titel:
        titel.remove('')

    validator = np.lib._iotools.NameValidator()
    titel = [i.lower() for i in validator(titel)]
    usecols = range(0, len(titel))

    measurements = pd.read_csv(f, header=None, names=titel,
                               parse_dates=['peildatum'],
                               index_col='peildatum',
                               dayfirst=True,
                               usecols=usecols)

    #measurements['Stand (m t.o.v. NAP)'] = measurements['Stand (cm t.o.v. NAP)'] /100.

    #measurements.set_index('Peildatum', inplace=True)

    return measurements


def read_dino_waterlvl_csv(fname, to_mnap=True, read_series=True, verbose=False):
    '''Read dino waterlevel data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    verbose : boolean, optional
        print additional information to the screen (default is False).
    '''
    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))

    p_meta = re.compile(
        'Locatie,Externe aanduiding,X-coordinaat,Y-coordinaat, Startdatum, Einddatum')
    p_data = re.compile(
        'Locatie,Peildatum,Stand \(cm t.o.v. NAP\),Bijzonderheid')

    with open(fname, 'r') as f:
        line = f.readline()
        while line != '':
            line = f.readline()
            if p_meta.match(line):
                meta = _read_dino_waterlvl_metadata(f, line)
                meta['filename'] = fname
            elif p_data.match(line):
                if read_series:
                    measurements = _read_dino_waterlvl_measurements(f, line)
                    if to_mnap and measurements is not None:
                        measurements['stand_m_tov_nap'] = measurements['stand_cm_tov_nap'] / 100.
                else:
                    measurements = None

                return measurements, meta


# %% General DINO methods
def read_dino_dir(dirname, ObsClass=None,
                  subdir='Boormonsterprofiel_Geologisch booronderzoek', suffix='.txt',
                  unpackdir=None, force_unpack=False, preserve_datetime=False, verbose=False,
                  keep_all_obs=True, **kwargs
                  ):
    '''Read Dino directory with point observations

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
    verbose : boolean, optional
        Print additional information to the screen (default is False)
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_dino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations

    '''

    # unzip dir
    if dirname.endswith('.zip'):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(zipf, dirname, force=force_unpack,
                   preserve_datetime=preserve_datetime)

    # read filenames
    files = os.listdir(os.path.join(dirname, subdir))
    if suffix:
        files = [file for file in files if file.endswith(suffix)]

    if not files:
        raise FileNotFoundError('no files were found in {} that end with {}'.format(
            os.path.join(dirname, subdir), suffix))

    # read individual files
    obs_list = []
    for _, file in enumerate(files):
        fname = os.path.join(dirname, subdir, file)
        obs = ObsClass.from_dino_file(fname=fname, verbose=verbose, **kwargs)
        if obs.metadata_available and (not obs.empty):
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            if verbose:
                print('not added to collection -> {}'.format(fname))

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df


# %% Testing
if __name__ == '__main__':

    fname = r'..\data\2019-Dino-test\Grondwatersamenstellingen_Put\B52C0057.txt'

    me1, meta = read_dino_groundwater_quality_txt(fname, verbose=True)

    fname = r'..\data\2019-Dino-test\Peilschaal\P58A0001.csv'
    measurements, meta = read_dino_waterlvl_csv(fname, verbose=True)

    fname = r'art_tools\data\2019-Dino-test\Grondwaterstanden_Put\B33F0080001_1.csv'
    obs_df = read_dino_groundwater_csv(fname, verbose=True)

    fname = r'..\data\2019-Dino-test\Peilschaal\P58A0005.csv'
    obs_df = read_dino_waterlvl_csv(fname, read_series=False, verbose=True)
