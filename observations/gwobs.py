import numpy as np
import pandas as pd

from . import accessor


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
        """his method computes the filternumbers and location name based on 
        the location of the observations. Then it sets the value of the 
        filternumber and the location:
            - in the ObsCollection dataframe
            - as the attribute of an Obs object
            - in the meta dictionary of the Obs object (only if add_to_meta is
            True)
        
        This method is useful for groundwater observations. If two or more
        observation points are close to each other they will be seen as one
        location with multiple filters. The filternr is based on the
        'onderkant_filter' attribute of the observations in such a way that 
        the deepest filter has the highest filter number. The location is
        based on the 

        Parameters
        ----------
        loc_col : str
            the column name with the names to use for the location
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
        modellayers=[]
        for o in self._obj.obs.values:
            modellayers.append(o.get_modellayer(ml, zgr, verbose))
            if verbose:
                print(o.name)
        
        modellayers = pd.Series(index=self._obj.index, data=modellayers, 
                                name='modellayer')
        
        return modellayers
