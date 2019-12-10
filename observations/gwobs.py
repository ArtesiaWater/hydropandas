import numpy as np

from . import accessor


@accessor.register_obscollection_accessor("gwobs")
class GwObsAccessor:

    def __init__(self, oc_obj):
        self._obj = oc_obj

    def get_filter_num(self, radius=1, xcol='x', ycol='y', if_exists='error'):
        """This method will add a column to the ObsCollection with the
        filternr. This is useful for groundwater observations. If two
        observation points are close to each other they will be seen as one
        location with multiple filters. The filternr will be added based on the
        'onderkant_filter' attribute of the observations.

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

        Raises
        ------
        RuntimeError
            if the column filternr exists and if_exists='error' an error is raised
        """
        # ken filternummers toe aan peilbuizen die dicht bij elkaar staan
        if 'filternr' in self._obj.columns:
            if if_exists == 'error':
                raise RuntimeError(
                    "the column 'filternr' already exist, set if_exists='replace' to replace the current values")
            elif if_exists == 'replace':
                self._obj['filternr'] = np.nan
        else:
            self._obj['filternr'] = np.nan

        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, 'filternr']):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_filters = np.sqrt(
                    (self._obj[xcol] - x)**2 + (self._obj[ycol] - y)**2)
                dup_x = self._obj.loc[distance_to_other_filters < radius]
                if dup_x.shape[0] == 1:
                    self._obj.loc[name, 'filternr'] = 1
                else:
                    dup_x2 = dup_x.sort_values(
                        'onderkant_filter', ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        self._obj.loc[pb_dub, 'filternr'] = i + 1
                        self._obj.loc[pb_dub, 'obs'].filternr = i + 1
                        self._obj.loc[pb_dub, 'obs'].meta['filternr'] = i + 1

    def get_filter_num_location(self, loc_col, radius=1, xcol='x', ycol='y',
                                if_exists='error'):
        """This method will add two columns to the ObsCollection with the
        filternr and the location. This is useful for groundwater observations.
        If two observation points are close to each other they will be seen as
        one location with multiple filters. The filternr will be added based on
        the 'onderkant_filter' attribute of the observations.

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

        Raises
        ------
        RuntimeError
            if the column filternr exists and if_exists='error' an error is raised
        """
        # ken filternummers toe aan peilbuizen die dicht bij elkaar staan
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

        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, 'filternr']):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_filters = np.sqrt(
                    (self._obj[xcol] - x)**2 + (self._obj[ycol] - y)**2)
                dup_x = self._obj.loc[distance_to_other_filters < radius]
                if dup_x.shape[0] == 1:
                    self._obj.loc[name, 'filternr'] = 1
                    self._obj.loc[name, 'obs'].filternr = 1
                    self._obj.loc[name, 'obs'].meta['filternr'] = 1
                    locatie = self._obj.loc[name, loc_col]
                    self._obj.loc[name, 'locatie'] = locatie
                    self._obj.loc[name, 'obs'].locatie = locatie
                    self._obj.loc[name, 'obs'].meta['locatie'] = locatie
                else:
                    dup_x2 = dup_x.sort_values(
                        'onderkant_filter', ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        if i == 0:
                            locatie = self._obj.loc[pb_dub, loc_col]
                        self._obj.loc[pb_dub, 'filternr'] = i + 1
                        self._obj.loc[pb_dub, 'obs'].filternr = i + 1
                        self._obj.loc[pb_dub, 'obs'].meta['filternr'] = i + 1

                        self._obj.loc[pb_dub, 'locatie'] = locatie
                        self._obj.loc[pb_dub, 'obs'].locatie = locatie
                        self._obj.loc[pb_dub, 'obs'].meta['locatie'] = locatie

    def get_pb_modellayers(self, ml, zgr=None, verbose=False):
        """Get the modellayers from the dis file

        Parameters
        ----------
        ml : flopy.modflow.mf.Modflow
            modflow model
        zgr : np.3darray, optional
            array containing model layer elevation
            information (if None , this information is obtained
            from the dis object)
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """

        for o in self._obj.obs.values:
            o.get_pb_modellayer(ml, zgr, verbose)
            if verbose:
                print(o.name)

        self._obj.add_meta_to_df('modellayer')

        # this does not work because there are nan values in the modellayer
        # column
#        if self._obj['modellayer'].dtype != int:
#            self._obj['modellayer'] == self._obj['modellayer'].astype(int)
