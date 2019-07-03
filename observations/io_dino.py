import os
import re
import tempfile

import numpy as np
import pandas as pd

from .util import unzip_file, _import_art_tools


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


def download_dino_groundwater(name, filternr, tmin, tmax,
                              x=np.nan, y=np.nan, **kwargs):
    """

    Parameters
    ----------
    name : str
        location str of the piezometer, i.e. B57F0077
    filter_nr : str, int, float
        filter number, is converted to str, i.e. 004
    tmin : str or pandas.Timestamp
        start date in format YYYY-MM-DD (will be converted if Timestamp)
    tmax : str or pandas.Timestamp
        end date in format YYYY-MM-DD (will be converted if Timestamp)
    x : int, float, optional
        the x coördinate of the measurement point (not read from server)
    y : int, float, optional
        the y coördinate of the measurement point (not read from server)
    kwargs : key-word arguments
            these arguments are passed to dino.findMeetreeks functie

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    # attempt art_tools import
    art = _import_art_tools()

    # download data from dino
    dino = art.dino_wsdl.DinoWSDL()

    if isinstance(filternr, float) or isinstance(filternr, int):
        filternr = "{0:03g}".format(filternr)

    # measurements
    measurements = dino.findMeetreeks(name, filternr, tmin, tmax,
                                      **kwargs)

    measurements.rename(columns={'_'.join([name, filternr]): 'Stand_m_tov_NAP'},
                        inplace=True)

    meta = dino.findTechnischeGegevens(name, filternr)

    meta['onderkant_filter'] = meta.pop('BOTTOM_FILTER')
    meta['bovenkant_filter'] = meta.pop('TOP_FILTER')
    meta['meetpunt'] = meta.pop('SFL_HEIGHT')

    filternr = meta.pop('WELL_TUBE_NR')
    name = meta.pop('WELL_NITG_NR')
    name = '-'.join([name, filternr])
    meta['filternr'] = int(filternr)
    meta['name'] = name
    meta['x'] = x
    meta['y'] = y
    meta['locatie'] = name

    return measurements, meta


def read_dino_groundwater_csv(fname, to_mnap=True,
                              read_series=True, verbose=False,):
    """Read dino groundwater quantity data from a dinoloket csv file.

    Parameters
    ----------
    fname : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'Stand_m_tov_NAP' is added to the dataframe
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
            elif verbose and measurements[~measurements.Stand_cm_tov_NAP.isna()].empty:
                print('no NAP measurements available -> {}'.format(fname))
            if to_mnap and measurements is not None:
                measurements['Stand_m_tov_NAP'] = measurements['Stand_cm_tov_NAP']/100.
        else:
            measurements = None

    return measurements, meta


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
    titel = validator(titel)
    usecols = range(0, len(titel))

    measurements = pd.read_csv(f, header=None, names=titel,
                               parse_dates=['Peildatum'],
                               index_col='Peildatum',
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
        if True a column with 'Stand_m_tov_NAP' is added to the dataframe
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
                        measurements['Stand_m_tov_NAP'] = measurements['Stand_cm_tov_NAP']/100.
                else:
                    measurements = None

                return measurements, meta


def download_dino_within_extent(extent=None, bbox=None, ObsClass=None,
                                kind='Grondwaterstand',
                                tmin=None, tmax=None,
                                zmin=None, zmax=None, unit="NAP",
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
    kind : str
        The type of timeseries (Grondwaterstand or Grondwatersamenstelling)
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
                                              kind=kind)

    if verbose:
        print('\ndownload {} data from dino within:\n\
               - extent: {} or\n\
               - bbox: {}'.format(kind, extent, bbox))

    # slice by properties
    if tmin is not None:
        tmin = pd.to_datetime(tmin)
        mask = gdf_loc.end_date > tmin
        gdf_loc = gdf_loc.loc[mask]
    if tmax is not None:
        tmax = pd.to_datetime(tmax)
        mask = gdf_loc.start_date < tmax
        gdf_loc = gdf_loc.loc[mask]
    if zmin is not None:
        mask = gdf_loc.top_height_nap >= zmin
        gdf_loc = gdf_loc.loc[mask]
    if zmax is not None:
        mask = gdf_loc.bottom_height_nap <= zmax
        gdf_loc = gdf_loc.loc[mask]
        
    if gdf_loc.empty:
        return pd.DataFrame()

    # read measurements
    obs_list = []
    for index, loc in gdf_loc.iterrows():
        if tmin is None:
            tmin_t = loc.start_date
        else:
            tmin_t = tmin
        if tmax is None:
            tmax_t = loc.end_date
        else:
            tmax_t = tmax

        if verbose:
            print('reading -> {}-{}'.format(loc.dino_nr, loc.piezometer_nr))

        o = ObsClass.from_dino_server(name=loc.dino_nr,
                                      filternr=float(loc.piezometer_nr),
                                      tmin=tmin_t,
                                      tmax=tmax_t, x=loc['x_rd_crd'],
                                      y=loc['y_rd_crd'],
                                      unit=unit)

        obs_list.append(o)

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df


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
        add all observation points to the collection, even without metadata
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
        if obs.metadata_available:
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


def _read_dino_groundwater_header(f):
    line = f.readline()
    header = dict()
    while line not in ['\n', '', '\r\n']:
        propval = line.split(',')
        prop = propval[0]
        prop = prop.replace(':', '')
        prop = prop.strip()
        val = propval[1]
        if propval[2] != '':
            val = val + ' ' + propval[2].replace(':', '') + ' ' + \
                propval[3]
        header[prop] = val
        line = f.readline()

    return line, header


def _read_empty(f, line):
    while (line == '\n') or (line == '\r\n'):
        line = f.readline()

    return line


def _read_dino_groundwater_referencelvl(f, line):
    ref = dict()
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
            meta[properties[i]] = values[i]
        metalist.append(meta)
        line = f.readline()

    meta = {}
    if metalist:
        meta["locatie"] = metalist[-1]['Locatie']
        meta["filternr"] = int(float(metalist[-1]['Filternummer']))
        meta["name"] = '-'.join([meta["locatie"],
                                 metalist[-1]['Filternummer']])
        meta["x"] = float(metalist[-1]['X-coordinaat'])
        meta["y"] = float(metalist[-1]['Y-coordinaat'])
        meetpunt = metalist[-1]['Meetpunt (cm t.o.v. NAP)']
        if meetpunt == '':
            meta["meetpunt"] = np.nan
        else:
            meta["meetpunt"] = float(meetpunt) / 100.
        maaiveld = metalist[-1]['Maaiveld (cm t.o.v. NAP)']
        if maaiveld == '':
            meta["maaiveld"] = np.nan
        else:
            meta["maaiveld"] = float(maaiveld) / 100
        bovenkant_filter = metalist[-1]['Bovenkant filter (cm t.o.v. NAP)']
        if bovenkant_filter == '':
            meta["bovenkant_filter"] = np.nan
        else:
            meta["bovenkant_filter"] = float(bovenkant_filter) / 100
        onderkant_filter = metalist[-1][
            'Onderkant filter (cm t.o.v. NAP)']
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
        titel = validator(titel)
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(f, header=None, names=titel,
                                       parse_dates=['Peildatum'],
                                       index_col='Peildatum',
                                       dayfirst=True,
                                       usecols=usecols)
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


if __name__ == '__main__':

    fname = r'..\data\2019-Dino-test\Grondwatersamenstellingen_Put\B52C0057.txt'

    me1, meta = read_dino_groundwater_quality_txt(fname, verbose=True)

    fname = r'..\data\2019-Dino-test\Peilschaal\P58A0001.csv'
    measurements, meta = read_dino_waterlvl_csv(fname, verbose=True)

    fname = r'art_tools\data\2019-Dino-test\Grondwaterstanden_Put\B33F0080001_1.csv'
    obs_df = read_dino_groundwater_csv(fname, verbose=True)

    fname = r'..\data\2019-Dino-test\Peilschaal\P58A0005.csv'
    obs_df = read_dino_waterlvl_csv(fname, read_series=False, verbose=True)
