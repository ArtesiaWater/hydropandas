import numpy as np
import pandas as pd
import scipy.interpolate as intp

from . import accessor


def findrowcolumn(xpts, ypts, ml):
    """
    find row and column number for points in an irregular grid. Adapted from Flopy.

    Parameters
    ----------
    xpts : 1d array
        x-coordinates of points
    ypts : 1d array
        y-coordinates of points
    ml : flopy.modflow.mf.Modflow
        modflow model

    Returns
    -------
    r, c: arrays with dtype: int
        array containing row and column indices of the points

    """
    xedge, yedge = ml.modelgrid.xyedges

    if not isinstance(xpts, np.ndarray):
        xpts = np.asarray(xpts)
    if not isinstance(ypts, np.ndarray):
        ypts = np.asarray(ypts)
    if not isinstance(xedge, np.ndarray):
        xedge = np.asarray(xedge)
    if not isinstance(yedge, np.ndarray):
        yedge = np.asarray(yedge)

    fx = intp.interp1d(xedge, range(len(xedge)), bounds_error=False)
    fy = intp.interp1d(yedge, range(len(yedge)), bounds_error=False)

    c = np.asarray(np.floor(fx(xpts)), dtype=int)
    r = np.asarray(np.floor(fy(ypts)), dtype=int)

    return r, c


def get_model_layer(x, y, z, ml, dis, zgr=None, left=-999, right=999):
    """get index of model layer based on elevation
    at a specific location.

    Parameters
    ----------
    x : np.array
        x-coordinate of locations
    y : np.array
        y-coordinate of locations
    z : np.array
        elevations corresponding to x, y locations
    ml : flopy.modflow.mf.Modflow
        modflow model
    zgr : np.3darray, optional
        3D array containing layer elevations for if dis
        does not contain this information. (the default is
        None, which uses dis to determine layer elevations)

    Returns
    -------
    ilay: np.array
        array containing 0-indexed model layer for each point

    TODO:
        - avoid loop if model layer elevations are equal everywhere
        -
    """

    r, c = findrowcolumn(x, y, ml)

    if zgr is None:
        zgr = np.concatenate(
            [dis.top.array[np.newaxis], dis.botm.array], axis=0)

    ilay = np.nan * np.ones(x.shape, dtype=np.int)

    for i in range(x.shape[0]):
        ir = r[i]
        ic = c[i]

        if np.isnan(z[i]):
            continue
        elif (ir >= 0) & (ic >= 0):
            zvec = zgr[:, ir, ic]
            ilay[i] = int(np.interp(z[i], zvec[::-1], np.arange(len(zvec))[::-1],
                                    left=left, right=right))

    return ilay


def get_pb_modellayer(x, y, ftop, fbot, ml, zgr=None, left=-999, right=999, 
                      verbose=False):
    """get index of model layer based on filterdepth of
    piezometer.

    Parameters
    ----------
    x : np.array
        x-coordinate of piezometers
    y : np.array
        y-coordinate of piezometers
    ftop : np.array
        top elevation of filter
    fbot : np.array
        bottom elevation of filter
    ml : flopy.modflow.mf.Modflow
        modflow model
    zgr : np.3darray, optional
        array containing model layer elevation
        information (the default is None, which
        gets this information from the dis object)
    verbose : bool, optional
        print information about how layer is
        selected when piezometer is screened in
        multiple layers (the default is False,
        which prints no info)

    Returns
    -------
    ilay: np.array
        array containing index of model layer per piezometer

    TODO:
        - some things are done double in this function and in get_model_layer, necessary?
        - speed up if model layer elevation is constant everywhere?

    """
    if zgr is None:
        dis = ml.get_package('DIS')
    
    ilay_ftop = get_model_layer(x, y, ftop, ml, dis, zgr=zgr, 
                                left=left, right=right)
    ilay_fbot = get_model_layer(x, y, fbot, ml, dis, zgr=zgr, 
                                left=left, right=right)
    
    

    if zgr is None:
        zgr = np.concatenate(
            [dis.top.array[np.newaxis], dis.botm.array], axis=0)

    r, c = findrowcolumn(x, y, ml)

    ilay = np.nan * np.ones(x.shape, dtype=np.int)

    for i in range(x.shape[0]):
        ir = r[i]
        ic = c[i]

        if (ir >= 0) & (ic >= 0):
            zvec = zgr[:, ir, ic]

            fti = ilay_ftop[i]
            fbi = ilay_fbot[i]

            if np.isnan(fti) and np.isnan(fbi):
                continue
            elif np.isnan(fti):
                ilay[i] = fbi
            elif np.isnan(fbi):
                ilay[i] = fti
            # Warn if piezometer is on boundary. Piezometer in layer in which majority of screen is located
            elif fbi != fti:
                fbi = int(fbi)
                fti = int(fti)
                if verbose:
                    print("filter on layer boundary:")
                    print("-layers: {0}, {1}".format(fbi, fti))
                    print("-layer elev in between: {0:.2f}".format(zvec[fbi]))
                    print("-fb, ft: {0:.2f}, {1:.2f}".format(fbot[i], ftop[i]))
                    print("-length in layer: {0:.2f}, {1:.2f}".format(
                          zvec[fbi] - fbot[i], ftop[i] - zvec[fbi]))
                if fti==right and fbi==left:
                    if verbose:
                        print('filter spans all layers. '
                              'not sure which layer to select')
                elif fti==right:
                    if verbose:
                        print('filter top higher than top layer. '
                              'selected layer {}'.format(fbi))
                    ilay[i] = fbi
                elif fbi==left:
                    if verbose:
                        print('filter bot lower than bottom layer. '
                              'selected layer {}'.format(fti))
                    ilay[i] = fti
                elif fti - fbi > 2:
                    ilay[i] = np.nan  # set to unknown
                    if verbose:
                        print("Piezometer filter spans {} layers. "
                              "Not sure which layer to select".format(fti - fbi+1))
                elif fti - fbi == 2:  # if filter spans 3 layers
                    ilay[i] = fbi + 1  # use middle layer
                elif np.abs(fbot[i] - zvec[fbi]) > np.abs(ftop[i] - zvec[fbi]):
                    ilay[i] = fbi
                    if verbose:
                        print("selected layer:", fbi)
                elif np.abs(fbot[i] - zvec[fbi]) < np.abs(ftop[i] - zvec[fbi]):
                    ilay[i] = fti
                    if verbose:
                        print("selected layer:", fti)
            else:
                ilay[i] = fbi

    return ilay


@accessor.register_obscollection_accessor("gwobs")
class GwObsAccessor:

    def __init__(self, oc_obj):
        self._obj = oc_obj

    def set_filter_num(self, radius=1, xcol='x', ycol='y', if_exists='error',
                       add_to_meta=False):
        """This method computes the filternumbers based on the location of the
        observations. Then it sets the value of the filternumber:
            - in the ObsCollection dataframe
            - as the attribute of an Obs object
            - in the meta dictionary of the Obs object (only if add_to_meta is
            True)

        This method is useful for groundwater observations. If two or more
        observation points are close to each other they will be seen as one
        location with multiple filters. The filternr is based on the
        'onderkant_filter' attribute of the observations in such a way that
        the deepest filter has the highest filter number.

        Parameters
        ----------
        radius : int, optional
            max distance between two observations to be seen as one location, by default 1
        xcol : str, optional
            column name with x coordinates, by default 'x'
        ycol : str, optional
            column name with y coordinates, by default 'y'
        if_exists : str, optional
            what to do if an observation point already has a filternr, options:
            'error', 'replace' or 'keep', by default 'error'
        add_to_meta : bool, optional
            if True the filternr is added to the meta dictionary
            of an observation. The default is False.

        Raises
        ------
        RuntimeError
            if the column filternr exists and if_exists='error' an error is raised
        """

        if self._obj['filternr'].dtype != np.number:
            self._obj['filternr'] = pd.to_numeric(
                self._obj['filternr'], errors='coerce')

        # check if column exists in obscollection
        if 'filternr' in self._obj.columns:
            if if_exists == 'error':
                raise RuntimeError(
                    "the column 'filternr' already exist, set if_exists='replace' to replace the current values")
            elif if_exists == 'replace':
                self._obj['filternr'] = np.nan
        else:
            self._obj['filternr'] = np.nan

        # ken filternummers toe aan peilbuizen die dicht bij elkaar staan
        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, 'filternr']):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_filters = np.sqrt(
                    (self._obj[xcol] - x)**2 + (self._obj[ycol] - y)**2)
                dup_x = self._obj.loc[distance_to_other_filters < radius]
                if dup_x.shape[0] == 1:
                    self._obj._set_metadata_value(name, 'filternr', 1,
                                                  add_to_meta=add_to_meta)
                else:
                    dup_x2 = dup_x.sort_values(
                        'onderkant_filter', ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        self._obj._set_metadata_value(pb_dub, 'filternr', i + 1,
                                                      add_to_meta=add_to_meta)

    def set_filter_num_location(self, loc_col, radius=1, xcol='x', ycol='y',
                                if_exists='error', add_to_meta=False):
        """This method sets the filternr and locatie name of an observation
        point based on the locatie of the observations. When two or more
        filters are located close, defined by a radius, to eachother they get
        the same `locatie` and a different `filternr`.

        The value of the filternumber and the location are set:
            - in the ObsCollection dataframe
            - as the attribute of an Obs object
            - in the meta dictionary of the Obs object (only if add_to_meta is
            True)

        This method is useful for groundwater observations. If two or more
        observation points are close to each other they will be seen as one
        location with multiple filters. The filternr is based on the
        'onderkant_filter' attribute of the observations in such a way that
        the deepest filter has the highest filter number. The locatie is
        based on the named of the loc_col of the filter with the lowest
        filternr.

        Parameters
        ----------
        loc_col : str
            the column name with the names to use for the locatie
        radius : int, optional
            max distance between two observations to be seen as one location, by default 1
        xcol : str, optional
            column name with x coordinates, by default 'x'
        ycol : str, optional
            column name with y coordinates, by default 'y'
        if_exists : str, optional
            what to do if an observation point already has a filternr, options:
            'error', 'replace' or 'keep', by default 'error'
        add_to_meta : bool, optional
            if True the filternr and location are added to the meta dictionary
            of an observation. The default is False.

        Raises
        ------
        RuntimeError
            if the column filternr exists and if_exists='error' an error is raised
        """

        # check if columns exists in obscollection
        if 'filternr' in self._obj.columns or 'locatie' in self._obj.columns:
            if if_exists == 'error':
                raise RuntimeError(
                    "the column 'filternr or locatie' already exist, set if_exists='replace' to replace the current values")
            elif if_exists == 'replace':
                self._obj['filternr'] = np.nan
                self._obj['locatie'] = np.nan
        else:
            self._obj['filternr'] = np.nan
            self._obj['locatie'] = np.nan

        # ken filternummers toe aan peilbuizen die dicht bij elkaar staan
        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, 'filternr']):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_filters = np.sqrt(
                    (self._obj[xcol] - x)**2 + (self._obj[ycol] - y)**2)
                dup_x = self._obj.loc[distance_to_other_filters < radius]
                if dup_x.shape[0] == 1:
                    self._obj._set_metadata_value(name, 'filternr', 1,
                                                  add_to_meta=add_to_meta)
                    locatie = self._obj.loc[name, loc_col]
                    self._obj._set_metadata_value(name, 'locatie', locatie,
                                                  add_to_meta=add_to_meta)
                else:
                    dup_x2 = dup_x.sort_values(
                        'onderkant_filter', ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        if i == 0:
                            locatie = self._obj.loc[pb_dub, loc_col]
                        self._obj._set_metadata_value(pb_dub, 'filternr', i + 1,
                                                      add_to_meta=add_to_meta)
                        self._obj._set_metadata_value(pb_dub, 'locatie', locatie,
                                                      add_to_meta=add_to_meta)

    def get_modellayers(self, ml, zgr=None, verbose=False):
        """Get the modellayer per observation. The layers can be obtained
        from the modflow model or can be defined in zgr.

        Parameters
        ----------
        ml : flopy.modflow.mf.Modflow
            modflow model
        zgr : np.3darray, optional
            array containing model layer elevation
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        pd.Series with the modellayers of each observation
        """
        modellayers = []
        for o in self._obj.obs.values:
            modellayers.append(o.gwobs.get_modellayer(ml, zgr, verbose))
            if verbose:
                print(o.name)

        modellayers = pd.Series(index=self._obj.index, data=modellayers,
                                name='modellayer')

        return modellayers
    
    


@accessor.register_obs_accessor("gwobs")
class GeoAccessorObs:
    def __init__(self, obs):
        self._obj = obs
        
    def get_modellayer(self, ml, zgr=None, verbose=False):
        """Add modellayer to meta dictionary

        Parameters
        ----------
        ml : flopy.modflow.mf.Modflow
            modflow model
        zgr : np.3darray, optional
            array containing model layer elevation
            information (the default is None, which
            gets this information from the dis object)
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """
        modellayer = get_pb_modellayer(np.array([self._obj.x]) - ml.modelgrid.xoffset,
                                       np.array([self._obj.y]) -
                                       ml.modelgrid.yoffset,
                                       np.array([self._obj.bovenkant_filter]),
                                       np.array([self._obj.onderkant_filter]),
                                       ml, zgr, verbose=verbose)[0]

        return modellayer
        
    