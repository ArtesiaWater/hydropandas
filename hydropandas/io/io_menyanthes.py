# -*- coding: utf-8 -*-
"""Created on Thu Oct 10 11:01:22 2019.

@author: oebbe
"""

import os
import logging

import numpy as np
from hydropandas import observation
from pandas import DataFrame, Series
from scipy.io import loadmat

from ..util import matlab2datetime
logger = logging.getLogger(__name__)

def read_file(fname, ObsClass):
    """This method is used to read the file."""

    logger.info(f'reading menyanthes file {fname}')

    if ObsClass == observation.GroundwaterObs:
        _rename_dic = {'xcoord': 'x',
                       'ycoord': 'y',
                       'upfiltlev': 'bovenkant_filter',
                       'lowfiltlev': 'onderkant_filter',
                       'surflev': 'maaiveld',
                       'filtnr': 'filternr',
                       'meetpunt': 'measpointlev'
                       }

        _keys_o = ['name', 'x', 'y', 'locatie', 'filternr',
                   'metadata_available', 'maaiveld', 'meetpunt',
                   'bovenkant_filter', 'onderkant_filter']

    elif ObsClass == observation.WaterlvlObs:
        _rename_dic = {'xcoord': 'x',
                       'ycoord': 'y',
                       'meetpunt': 'measpointlev'
                       }

        _keys_o = ['name', 'x', 'y', 'locatie']

    # Check if file is present
    if not (os.path.isfile(fname)):
        print('Could not find file ', fname)

    mat = loadmat(fname, struct_as_record=False, squeeze_me=True,
                  chars_as_strings=True)

    d_h = read_oseries(mat)

    locations = d_h.keys()
    obs_list = []
    for location in locations:
        logger.info(f'reading location -> {location}')
        metadata = d_h[location]
        metadata['projection'] = 'epsg:28992'
        metadata['metadata_available'] = True

        s = metadata.pop('values')
        df = DataFrame(s, columns=['stand_m_tov_nap'])
        for key in _rename_dic.keys():
            if key in metadata.keys():
                metadata[_rename_dic[key]] = metadata.pop(key)

        meta_o = {k: metadata[k] for k in _keys_o if k in metadata}

        o = ObsClass(df, meta=metadata, **meta_o)
        obs_list.append(o)

    return obs_list


def read_oseries(mat):
    """Read the oseries from a mat file from menyanthes."""
    d_h = {}
    # Check if more then one time series model is present
    if not isinstance(mat['H'], np.ndarray):
        mat['H'] = [mat['H']]

    # Read all the time series models
    for i, H in enumerate(mat['H']):
        data = {}

        for name in H._fieldnames:
            if name != 'values':
                data[name] = getattr(H, name)
            else:
                if H.values.size == 0:
                    # when diver-files are used, values will be empty
                    series = Series()
                else:
                    tindex = map(matlab2datetime, H.values[:, 0])
                    # measurement is used as is
                    series = Series(H.values[:, 1], index=tindex)
                    # round on seconds, to get rid of conversion milliseconds
                    series.index = series.index.round('s')
                data['values'] = series

        # add to self.H
        if not hasattr(H, 'Name') and not hasattr(H, 'name'):
            H.Name = 'H' + str(i)  # Give it the index name
        if hasattr(H, 'name'):
            H.Name = H.name
        if len(H.Name) == 0:
            H.Name = H.tnocode

        d_h[H.Name] = data

    return d_h
