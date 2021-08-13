import os
import warnings

import numpy as np

from .. import util
from ..observation import ModelObs

import logging
logger = logging.getLogger(__name__)


def read_imod_results(obs_collection, ml, runfile, mtime, model_ws,
                      modelname='', nlay=None, exclude_layers=0):
    """Read imod model results at point locations.

    Parameters
    ----------
    obs_collection : ObsCollection
        collection of observations at which points imod results will be read
    ml : flopy.modflow.mf.model
        modflow model
    runfile : Runfile
        imod runfile object
    mtime : list of datetimes
        datetimes corresponding to the model periods
    model_ws : str
        model workspace with imod model
    nlay : int, optional
        number of layers if None the number of layers from ml is used.
    modelname : str
        modelname
    exclude_layers : int
        exclude modellayers from being read from imod
    """
    import imod
    if ml.modelgrid.xoffset == 0 or ml.modelgrid.yoffset == 0:
        warnings.warn('you probably want to set the xll and/or yll '
                      'attributes of ml.modelgrid')

    if nlay is None:
        nlay = ml.modelgrid.nlay

    xmid, ymid, _ = ml.modelgrid.xyzcellcenters

    xy = np.array([xmid.ravel(), ymid.ravel()]).T
    uv = obs_collection.loc[:, ("x", "y")].dropna(how="any", axis=0).values
    vtx, wts = util.interp_weights(xy, uv)

    hm_ts = np.zeros((obs_collection.shape[0], len(mtime)))

    # loop over layers
    for m in range(nlay):
        if m < exclude_layers:
            continue
        mask = obs_collection.modellayer.values == m
        # loop over timesteps
        for t, date in enumerate(mtime):
            head_idf = 'head_{}_l{}.idf'.format(
                date.strftime('%Y%m%d'), m + 1)
            fname = os.path.join(
                model_ws,
                runfile.data['OUTPUTDIRECTORY'],
                'head',
                head_idf)

            logger.info(f'read {fname}')
            ihds, _attrs = imod.idf.read(fname)
            hm = util.interpolate(ihds, vtx, wts)
            hm_ts[mask, t] = hm[mask]

    mo_list = []
    for i, name in enumerate(obs_collection.index):
        mo = ModelObs(index=mtime,
                      data=hm_ts[i],
                      name=name,
                      model=modelname,
                      x=obs_collection.loc[name, 'x'],
                      y=obs_collection.loc[name, 'y'],
                      meta=obs_collection.loc[name, 'obs'].meta)
        mo_list.append(mo)

    return mo_list


def read_modflow_results(obs_collection, ml, hds_arr, mtime,
                         modelname='', nlay=None, exclude_layers=None):
    """Read modflow groundwater heads at points in obs_collection.

    Parameters
    ----------
    obs_collection : ObsCollection
        locations of model observation
    ml : flopy.modflow.mf.model
        modflow model
    hds_arr : numpy array
        heads with shape (ntimesteps, nlayers, nrow, ncol)
    mtime : list of datetimes
        dates for each model timestep
    modelname : str, optional
        modelname
    nlay : int, optional
        number of layers if None the number of layers from ml is used.
    exclude_layers : list of int, optional
        exclude the observations in these modellayers

    """
    if ml.modelgrid.grid_type == 'structured':
        if ml.modelgrid.xoffset == 0 or ml.modelgrid.yoffset == 0:
            warnings.warn('you probably want to set the xll '
                          'and/or yll attributes in DIS!')

    if nlay is None:
        nlay = ml.modelgrid.nlay

    if modelname == '':
        modelname = ml.name

    xmid, ymid, _ = ml.modelgrid.xyzcellcenters

    if ml.modelgrid.grid_type == 'structured':
        xy = np.array([xmid.ravel(), ymid.ravel()]).T
    elif ml.modelgrid.grid_type == 'vertex':
        xy = np.array([xmid, ymid]).T

    uv = obs_collection.loc[:, ("x", "y")].dropna(how="any", axis=0).values
    vtx, wts = util.interp_weights(xy, uv)

    # get interpolated timeseries from hds_arr
    hm_ts = np.nan * np.ones((obs_collection.shape[0], hds_arr.shape[0]))

    if exclude_layers is None:
        exclude_layers = []

    # loop over layers
    for m in range(nlay):
        if m in exclude_layers:
            continue
        mask = obs_collection["modellayer"].values == m
        # loop over timesteps
        for t in range(hds_arr.shape[0]):
            ihds = hds_arr[t, m]
            ihds[ihds <= -999.] = np.nan
            hm = util.interpolate(ihds, vtx, wts)
            hm_ts[mask, t] = hm[mask]

    mo_list = []
    for i, name in enumerate(obs_collection.index):
        mo = ModelObs(index=mtime,
                      data=hm_ts[i],
                      name=name,
                      model=modelname,
                      x=obs_collection.loc[name, 'x'],
                      y=obs_collection.loc[name, 'y'],
                      meta=obs_collection.loc[name, 'obs'].meta)
        mo_list.append(mo)

    return mo_list
