# -*- coding: utf-8 -*-
"""Created on Wed Sep 12 12:15:42 2018.

@author: Artesia
"""
import logging
import numbers

import pandas as pd
import pastastore as pst

logger = logging.getLogger(__name__)


def _get_metadata_from_obs(o):
    """internal method to get metadata in the right format for a pastas series.
    A pastas timeseries cannot handle the same metadata format as an
    observation object.

    Parameters
    ----------
    o : observations.Obs
        observation time series with metadata

    Returns
    -------
    meta : dictionary
        meta dictionary.
    """
    meta = dict()
    for attr_key in o._metadata:
        val = getattr(o, attr_key)
        if isinstance(val, (int, float, str, bool)):
            meta[attr_key] = val
        elif isinstance(val, numbers.Number):
            meta[attr_key] = float(val)
        elif isinstance(val, dict):
            for k, v in val.items():
                if isinstance(v, (int, float, str, bool)):
                    meta[k] = v
                elif isinstance(v, numbers.Number):
                    meta[k] = float(v)
                else:
                    logger.info(
                        f'did not add {k} to metadata because datatype is {type(v)}')
        else:
            logger.info(
                f'did not add {attr_key} to metadata because datatype is {type(val)}')

    return meta


def create_pastastore(oc,
                      pstore,
                      pstore_name='',
                      conn=None,
                      add_metadata=True,
                      obs_column='stand_m_tov_nap',
                      kind='oseries',
                      overwrite=False):
    """add observations to a new or existing pastastore.

    Parameters
    ----------
    oc : observation.ObsCollection
        collection of observations
    pstore : pastastore.PastaStore, optional
            Existing pastastore, if None a new pastastore is created
    pstore_name : str, optional
        Name of the pastastore only used if pstore is None
    conn : pastastore.connectors
        connector for database
    obs_column : str, optional
        Name of the column in the Obs dataframe to be used
    kind : str, optional
        The kind of series that is added to the pastastore
    add_metadata : boolean, optional
        If True metadata from the observations added to the pastastore
    overwrite : boolean, optional
        if True, overwrite existing series in pastastore, default is False

    Returns
    -------
    pstore : pastastore.PastaStore
        the pastastore with the series from the ObsCollection

    """
    if pstore is None:
        if conn is None:
            conn = pst.DictConnector("my_conn")
        pstore = pst.PastaStore(pstore_name, connector=conn)

    for o in oc.obs.values:
        logger.info('add to pastastore -> {}'.format(o.name))

        if add_metadata:
            meta = _get_metadata_from_obs(o)
        else:
            meta = dict()

        if kind == 'oseries':
            pstore.conn.add_oseries(o[obs_column], o.name,
                                    metadata=meta,
                                    overwrite=overwrite)
        else:
            pstore.conn.add_stress(o[obs_column], o.name,
                                   kind, metadata=meta,
                                   overwrite=overwrite)

    return pstore
