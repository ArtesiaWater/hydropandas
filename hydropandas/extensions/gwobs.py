import numpy as np
import pandas as pd
import scipy.interpolate as intp

from . import accessor


def get_model_layer_z(z, zvec, left=-999, right=999,
                      verbose=False):
    """get index of model layer based on elevation. Assumptions:
        - the highest value in zvec is the top of model layer 0.
        - if z is equal to the bottom of a layer, the model layer above that
        layer is assigned.

    Parameters
    ----------
    z : int, float
        elevation.
    zvec : list, np.array
        elevations of model layers. shape is nlay + 1
    left : int, optional
        if z is below the lowest value in zvec, this value is returned. 
        The default is -999.
    right : TYPE, optional
        if z is above the highest value in zvec, this value is returned. 
        The default is 999.

    Returns
    -------
    int
        model layer.
        
    Examples
    --------
    >>> zvec = [0, -10, -20, -30]
    >>> get_model_layer_z(-5, zvec)
    0
    
    >>> get_model_layer_z(-25, zvec)
    2
    
    >>> get_model_layer_z(-50, zvec)
    -999
    
    >>> get_model_layer_z(100, zvec)
    999
    
    >>> get_model_layer_z(-20, zvec)
    1
    """
    # make sure zvec is descending
    zvec = np.sort(zvec)[::-1]
    if z in zvec:
        z+= 1e-10
    lay = int(np.interp(z, zvec[::-1], np.arange(len(zvec))[::-1],
              left=left, right=right))
    
    return lay

def check_if_var_is_invalid(var):
    
    if var is None:
        return True
    elif np.isnan(var):
        return True
    
    return False

def get_modellayer_from_filter(ftop, fbot, zvec, left=-999, right=999,
                               verbose=False):
    """
    

    Parameters
    ----------
    ftop : int or float
        top of filter.
    fbot : int or float
        bottom of filter, has to be lower than ftop.
    zvec : list or np.array
        elevations of the modellayers at the location of the filter.
    left : int, optional
        value to return if filter is below the modellayers. 
        The default is -999.
    right : int, optional
        value to return if filter is above the modellayers. 
        The default is-999.
    verbose : bool, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        raised if something unexpected happens.

    Returns
    -------
    int or np.nan
        modellayer.
        
    examples
    --------
    >>> zvec = [0, -10, -20, -30, -40]
    >>> get_modellayer_from_filter(-5, -7, zvec)
    0
    
    >>> get_modellayer_from_filter(-25, -27, zvec)
    2
    
    >>> get_modellayer_from_filter(-15, -27, zvec, verbose=True)
    2
    
    >>> get_modellayer_from_filter(-5, -27, zvec, verbose=True)
    1
    
    >>> get_modellayer_from_filter(-5, -37, zvec, verbose=True)
    nan
    
    >>> get_modellayer_from_filter(15, -5, zvec, verbose=True)
    0
    
    >>> get_modellayer_from_filter(15, 5, zvec, verbose=True)
    999
    
    >>> get_modellayer_from_filter(-55, -65, zvec, verbose=True)
    -999
    
    >>> get_modellayer_from_filter(15, -65, zvec, verbose=True)
    nan
    
    >>> get_modellayer_from_filter(None, -7, zvec, verbose=True)
    0
    
    >>> get_modellayer_from_filter(None, None, zvec, verbose=True)
    nan
    
    """

    zvec = np.sort(zvec)[::-1]
    ftop_invalid = check_if_var_is_invalid(ftop)
    fbot_invalid = check_if_var_is_invalid(fbot)
    if ftop_invalid and fbot_invalid:
        return np.nan
    elif fbot_invalid:
        lay_ftop = get_model_layer_z(ftop, zvec,
                                     left=left, right=right)
        return lay_ftop
    elif ftop_invalid:
        lay_fbot = get_model_layer_z(fbot, zvec,
                                     left=left, right=right)
        return lay_fbot
    
    if ftop<fbot:
        if verbose:
            print('filter top below filter bot, switch top and bot')
        fbot, ftop = ftop, fbot
        
    lay_ftop = get_model_layer_z(ftop, zvec,
                                 left=left, right=right)
    lay_fbot = get_model_layer_z(fbot, zvec,
                                 left=left, right=right)
    
    # Piezometer in layer in which majority of screen is located
    if lay_fbot == lay_ftop:
        return lay_fbot
    
    else:
        if lay_fbot==left and lay_ftop==right:
            if verbose:
                print('filter spans all layers. '
                      'return nan')
            return np.nan
        elif lay_ftop==right:
            if verbose:
                print('filter top higher than top layer. '
                      'selected layer {}'.format(lay_fbot))
            return lay_fbot
        
        elif lay_fbot==left:
            if verbose:
                print('filter bot lower than bottom layer. '
                      'selected layer {}'.format(lay_ftop))
            return lay_ftop
        
        if verbose:
            print("filter crosses layer boundary:")
            print("-layers: {0}, {1}".format(lay_ftop, lay_fbot))
            print("-layer elev in between: {0:.2f} - {1:.2f}".format(zvec[lay_fbot+1], zvec[lay_ftop]))
            print("-filter elev in between: {0:.2f} - {1:.2f}".format(fbot, ftop)) 
        
        if lay_fbot - lay_ftop > 2:
            if verbose:
                print(f"Piezometer filter spans {lay_fbot - lay_ftop +1} layers."
                      " return nan")
            return np.nan
        
        elif lay_fbot - lay_ftop == 2:  # if filter spans 3 layers
            if verbose:
                print(f"Piezometer filter spans {lay_fbot - lay_ftop +1} layers."
                      f" Use middle layer {lay_ftop+1}")
            return lay_ftop + 1  # use middle layer
        
        # check which layer has the biggest length of the filter 
        length_lay_bot = zvec[lay_fbot] - fbot
        length_lay_top = ftop - zvec[lay_fbot]
        if verbose:
            print(f"length per layer:\n- lay {lay_ftop}: {length_lay_top:.2f}"
                  f"\n- lay {lay_fbot}: {length_lay_bot:.2f}")
            
        if length_lay_top <= length_lay_bot:
            if verbose:
                print("selected layer:", lay_fbot)
            return lay_fbot
        
        elif length_lay_top > length_lay_bot:
            if verbose:
                print("selected layer:", lay_fbot)
            return lay_ftop
    
    raise ValueError('something is wrong with the input please contact Artesia')

    return np.nan

def get_zvec(x, y, gwf=None):
    """ get a list with the vertical layer boundaries at a point in the model
    

    Parameters
    ----------
    x : int or float
        x coördinate.
    y : int or float
        y coordinate.
    gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
        modflow model with top and bottoms

    Raises
    ------
    NotImplementedError
        not all grid types are supported yet.

    Returns
    -------
    zvec : list
        list of vertical layer boundaries. length is nlay + 1.

    """
    
    if gwf is not None:
        if gwf.modelgrid.grid_type == 'structured':
            r,c = gwf.modelgrid.intersect(x,y)
            zvec = [gwf.modelgrid.top[r,c]] + [gwf.modelgrid.botm[i,r,c] for i in range(gwf.modelgrid.nlay)]
        elif gwf.modelgrid.grid_type == 'vertex':
            idx = gwf.modelgrid.intersect(x,y)
            zvec = [gwf.modelgrid.top[idx]] +[gwf.modelgrid.botm[i, idx] for i in range(gwf.modelgrid.nlay)]
        else:
            raise NotImplementedError(f'gridtype {gwf.modelgrid.grid_type} not (yet) implemented' )
        
    return zvec


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

    def get_modellayers(self, gwf, verbose=False):
        """Get the modellayer per observation. The layers can be obtained
        from the modflow model or can be defined in zgr.

        Parameters
        ----------
        gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
            modflow model
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        pd.Series with the modellayers of each observation
        """
        modellayers = []
        for o in self._obj.obs.values:
            modellayers.append(o.gwobs.get_modellayer_modflow(gwf, 
                                                              verbose=verbose))
            if verbose:
                print(o.name)

        modellayers = pd.Series(index=self._obj.index, data=modellayers,
                                name='modellayer')

        return modellayers
    
    


@accessor.register_obs_accessor("gwobs")
class GeoAccessorObs:
    def __init__(self, obs):
        self._obj = obs
        
    def get_modellayer_modflow(self, gwf, left=-999, right=999, verbose=False):
        """Add modellayer to meta dictionary

        Parameters
        ----------
        gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
            modflow model
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        
        Returns
        -------
        int
            modellayer
        """
    
        zvec = get_zvec(self._obj.x, self._obj.y, gwf)
        
        modellayer = get_modellayer_from_filter(self._obj.bovenkant_filter, 
                                                self._obj.onderkant_filter, 
                                                zvec, 
                                                left=left, right=right,
                                                verbose=verbose)
        return modellayer
        
    