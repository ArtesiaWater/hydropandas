import os
import re
import tempfile

import numpy as np
import pandas as pd

import requests
import json

import zeep
from zeep import Plugin
from zeep.wsa import WsAddressingPlugin
from zeep.plugins import HistoryPlugin

from ..util import unzip_file, _import_art_tools


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
    
def get_dino_piezometer_metadata(location, filternr, raw_response=False,
                                 verbose=False):
    """ get metadata of a dino piezometer
    

    Parameters
    ----------
    location : str
        location of piezometer, e.g. B57F0077
    filternr : str
        filter number should be a string of 3 characters, e.g. 002
    raw_response : bool, optional
        if True, return the requests response

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    url = "https://www.dinoloket.nl/javascriptmapviewer-web/rest/gws/gwo/details"
    response = requests.post(url, data=json.dumps([location]),
                             headers={'content-type': 'application/json'})
    if raw_response:
        return response
    else:
        data = json.loads(response.content)
        if data == []:
            meta = {'metadata_available':False}
        else:
            meta = _parse_metadata(data[0], location, filternr, verbose=verbose)
        return meta
    
def _parse_metadata(data, location, filternr, verbose=False):
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
        
    return meta



class RemoveWSA(Plugin):
    """Helper class to remove wsa tags
    from header in zeep XML post

    as described by: https://github.com/mvantellingen/python-zeep/issues/330#issuecomment-321781637

    """

    def egress(self, envelope, http_headers, operation, binding_options):
        envelope.find(
            '{http://schemas.xmlsoap.org/soap/envelope/}Header').clear()
        return envelope, http_headers


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
    meta = get_dino_piezometer_metadata(location, filternr)

    return measurements, meta


def download_dino_within_extent(extent=None, bbox=None, ObsClass=None,
                                layer='grondwatermonitoring',
                                tmin="1900-01-01", tmax="2040-01-01",
                                zmin=None, zmax=None, unit="NAP",
                                get_metadata=True,
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
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations

    """
    # attempt art_tools import
    art = _import_art_tools()

    # read locations
    gdf_loc = art.dino_wfs.get_dino_locations(extent=extent, bbox=bbox,
                                              layer=layer)

    if verbose:
        print('\ndownload {} data from dino within:\n'
              '- extent: {} or\n'
              '- bbox: {}'.format(layer, extent, bbox))

    if gdf_loc.empty:
        return pd.DataFrame()

    gdf_loc["startDate"] = gdf_loc.startDate.astype('datetime64[ns]')
    gdf_loc["endDate"] = gdf_loc.endDate.astype('datetime64[ns]')
    gdf_loc["topHeightNap"] = gdf_loc.topHeightNap.astype(float)
    gdf_loc["bottomHeightNap"] = gdf_loc.bottomHeightNap.astype(float)

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

        o = ObsClass.from_dino_server(name=index.split('-')[0],
                                      filternr=float(loc.piezometerNr),
                                      tmin=tmin_t,
                                      tmax=tmax_t,
                                      unit=unit,
                                      get_metadata=get_metadata)

        obs_list.append(o)

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df


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
