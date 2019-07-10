'''
module with ObsCollection class for a collection of observations.

The ObsCollection class is a subclass of a pandas DataFrame with
additional attributes and methods.

More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

'''
import os
import tempfile
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import observation as obs
from . import util


class ObsCollection(pd.DataFrame):
    """class for a collection of point observations

    An ObsCollection object is a subclass of a pandas.DataFrame and allows for
    additional attributes and methods. Additional attributes have to be
    defined in the '_metadata' attribute.

    Parameters
    ----------
    bbox : tuple
        bounding box of the obs collection
    name : str
        name of the observation collection
    meta : dic
        metadata of the observatino collection
    """
    # temporary properties
    _internal_names = pd.DataFrame._internal_names + ['none']
    _internal_names_set = set(_internal_names)

    # normal properties
    _metadata = ['bbox',
                 'name',
                 'meta',
                 ]

    def __init__(self, *args, **kwargs):
        """ constructor of the ObsCollection

        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.
        """
        self.bbox = kwargs.pop('bbox', ())
        self.name = kwargs.pop('name', '')
        self.meta = kwargs.pop('meta', {})

        super(ObsCollection, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return ObsCollection

    def _infer_otype(self, verbose=False):
        """Infer observation type from the obs column

        Parameters
        ----------
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        otype, otypes
            type or list of types of the observation objects

        """
        otypes = self.obs.apply(lambda x: type(x)).unique()
        if otypes.shape[0] == 1:
            if verbose:
                print('inferred observation type: {}'.format(otypes[0]))
            return otypes[0]
        elif otypes.shape[0] > 1:
            if verbose:
                print('inferred multiple otypes, types: {}'.format(otypes))
            return otypes
        else:
            raise TypeError('could not infer observation type')

    @classmethod
    def from_dino_server(cls, extent=None, bbox=None,
                         ObsClass=obs.GroundwaterObs,
                         name=None, verbose=False, **kwargs
                         ):
        """ Read dino data from a server

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
        name : str, optional
            the name of the observation collection
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs:
            kwargs are passed to the io_dino.download_dino_within_extent() function


        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """

        from . import io_dino

        if ObsClass == obs.GroundwaterObs:
            kind = 'Grondwaterstand'
        elif ObsClass == obs.GroundwaterQualityObs:
            kind = 'Grondwatersamenstelling'
        else:
            raise ValueError('cannot download {} from Dino'.format(ObsClass))

        obs_df = io_dino.download_dino_within_extent(extent, bbox,
                                                     ObsClass, kind=kind,
                                                     verbose=verbose, **kwargs)

        if bbox is None:
            bbox = [extent[0], extent[2], extent[1], extent[3]]

        if name is None:
            name = '{} from DINO'.format(kind)

        meta = kwargs
        meta.update({'verbose': verbose})

        return cls(obs_df, name=name, bbox=bbox, meta=meta)

    @classmethod
    def from_dino_dir(cls, dirname=None, ObsClass=obs.GroundwaterObs,
                      subdir='Grondwaterstanden_Put', suffix='1.csv',
                      unpackdir=None, force_unpack=False, preserve_datetime=False,
                      keep_all_obs=True, name=None, verbose=False,
                      **kwargs):
        """ Read a dino directory

        Parameters
        ----------
        extent : list, optional
            get dinodata online within this extent [xmin, xmax, ymin, ymax]
        dirname : str, optional
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
        keep_all_obs : boolean, optional
            add all observation points to the collection, even without metadata
        name : str, optional
            the name of the observation collection
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs:
            kwargs are passed to the io_dino.read_dino_dir() function

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """

        from . import io_dino

        if name is None:
            name = subdir

        meta = {'dirname': dirname,
                'type': ObsClass,
                'suffix': suffix,
                'unpackdir': unpackdir,
                'force_unpack': force_unpack,
                'preserver_datetime': preserve_datetime,
                'verbose': verbose,
                'keep_all_obs': keep_all_obs
                }

        obs_df = io_dino.read_dino_dir(dirname, ObsClass, subdir, suffix,
                                       unpackdir, force_unpack, preserve_datetime,
                                       verbose, keep_all_obs,
                                       **kwargs)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews(cls, dirname, ObsClass=obs.GroundwaterObs, name='fews',
                  to_mnap=True, remove_nan=True,
                  unpackdir=None, force_unpack=False,
                  preserve_datetime=False, verbose=False):
        """ read one or several XML-files with measurements from FEWS

        Parameters
        ----------
        dirname :  str
            De directory met xml bestanden die je wil inlezen
        ObsClass : type
            class of the observations, e.g. GroundwaterObs or WaterlvlObs
        name : str, optional
            name of the observation collection, 'fews' by default
        to_mnap : boolean, optional
            if True a column with 'Stand_m_tov_NAP' is added to the dataframe
        remove_nan : boolean, optional
            remove nan values from measurements, flag information about the
            nan values is also lost
        unpackdir : str
            destination directory to unzip file if fname is a .zip
        force_unpack : boolean, optional
            force unpack if dst already exists
        preserve_datetime : boolean, optional
            
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """
        from . import io_xml

        # unzip dir
        if dirname.endswith('.zip'):
            zipf = dirname
            if unpackdir is None:
                dirname = tempfile.TemporaryDirectory().name
            else:
                dirname = unpackdir
            util.unzip_file(zipf, dirname, force=force_unpack,
                            preserve_datetime=preserve_datetime)

        unzip_fnames = [i for i in os.listdir(dirname) if i.endswith(".xml")]

        meta = {'dirname': dirname,
                'type': ObsClass,
                'verbose': verbose
                }

        obs_list = []
        nfiles = len(unzip_fnames)
        for j, ixml in enumerate(unzip_fnames):
            if verbose:
                print("{0}/{1} read {2}".format(j+1, nfiles, ixml))
            fullpath = os.path.join(dirname, ixml)
            olist = io_xml.read_xml(fullpath, ObsClass=ObsClass, to_mnap=to_mnap,
                                    remove_nan=remove_nan, verbose=False)
            obs_list += olist

        obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                              columns=obs_list[0].to_collection_dict().keys())

        obs_df.set_index('name', inplace=True)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews2(cls, fname, ObsClass=obs.GroundwaterObs, name='fews',
                   to_mnap=True, remove_nan=True,
                   unpackdir=None, force_unpack=False,
                   preserve_datetime=False,
                   verbose=False, ):
        """ read a XML-file with measurements from FEWS

        Parameters
        ----------
        fname :  str
            De bestandsnaam van het bestand dat je wilt inlezen
        ObsClass : type
            class of the observations, e.g. GroundwaterObs or WaterlvlObs
        name : str, optional
            name of the observation collection, 'fews' by default
        to_mnap : boolean, optional
            if True a column with 'Stand_m_tov_NAP' is added to the dataframe
        remove_nan : boolean, optional
            remove nan values from measurements, flag information about the
            nan values is also lost
        unpackdir : str
            destination directory to unzip file if fname is a .zip
        force_unpack : boolean, optional
            force unpack if dst already exists
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """
        from . import io_xml

        # unzip dir
        if fname.endswith('.zip'):
            zipf = fname
            if unpackdir is None:
                dirname = tempfile.TemporaryDirectory().name
            else:
                dirname = unpackdir
            util.unzip_file(zipf, dirname, force=force_unpack,
                            preserve_datetime=preserve_datetime)

            unzip_fnames = os.listdir(dirname)
            if len(unzip_fnames) == 1:
                fname = os.path.join(dirname, unzip_fnames[0])
            elif len(unzip_fnames) > 1:
                raise ValueError(
                    'more than one file in zip, use from_fews_dir()')
            else:
                raise ValueError('empty zip')

        meta = {'fname': fname,
                'type': ObsClass,
                'verbose': verbose
                }

        obs_list = io_xml.read_xml_alternative(fname, ObsClass, to_mnap=to_mnap,
                                               remove_nan=remove_nan, verbose=verbose)

        obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                              columns=obs_list[0].to_collection_dict().keys())

        obs_df.set_index('name', inplace=True)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fieldlogger(cls, fname, name='', ObsClass=obs.GroundwaterObs):
        """Read a fieldlogger file into a list of observation objects

        Parameter
        ---------
        fname : str
            name of the fieldlogger location csv file
        name : str, optional
            name of the observation collection
        ObsClass : observation class
            type of fieldlogger observations


        """

        from . import fieldlogger

        obs_list, fieldlogger_meta = fieldlogger.fieldlogger_csv_to_obs_list(
            fname, ObsClass=obs.GroundwaterObs)

        return cls(cls.from_list(obs_list, name=name), meta=fieldlogger_meta)

    @classmethod
    def from_imod(cls, obs_collection, dis, runfile, mtime, model_ws,
                  modelname='', exclude_layers=0, verbose=False):
        """read 'observations' from an imod model

        Parameters
        ----------
        obs_collection : ObsCollection
            collection of observations at which points imod results will be read
        dis : modflow.Dis
            discretization from modflow
        runfile : Runfile
            imod runfile object
        mtime : list of datetimes
            datetimes corresponding to the model periods
        model_ws : str
            model workspace with imod model
        modelname : str
            modelname
        exclude_layers : int
            exclude modellayers from being read from imod
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        """

        import imod

        if dis.sr.xll == 0 or dis.sr.yll == 0:
            warnings.warn(
                'you probably want to set the xll and/or yll attributes of dis.sr')

        xmid = dis.sr.xcentergrid.ravel()
        ymid = dis.sr.ycentergrid.ravel()

        xy = np.c_[xmid, ymid]
        uv = obs_collection.loc[:, ("x", "y")].dropna(how="any", axis=0).values
        vtx, wts = util.interp_weights(xy, uv)

        hm_ts = np.zeros((obs_collection.shape[0], len(mtime)))

        # loop over layers
        for m in range(dis.nlay):
            if m < exclude_layers:
                continue
            mask = obs_collection.modellayer.values == m
            # loop over timesteps
            for t, date in enumerate(mtime):
                head_idf = 'head_{}_l{}.idf'.format(
                    date.strftime('%Y%m%d'), m+1)
                fname = os.path.join(
                    model_ws, runfile.data['OUTPUTDIRECTORY'], 'head', head_idf)
                if verbose:
                    print('read {}'.format(fname))
                ihds, _attrs = imod.idf.read(fname)
                hm = util.interpolate(ihds, vtx, wts)
                hm_ts[mask, t] = hm[mask]

        mo_list = []
        for i, name in enumerate(obs_collection.index):
            mo = obs.ModelObs(index=mtime, data=hm_ts[i], name=name, model=modelname,
                              x=obs_collection.loc[name,
                                                   'x'], y=obs_collection.loc[name, 'y'],
                              meta=obs_collection.loc[name, 'obs'].meta)
            mo_list.append(mo)

        mobs_df = pd.DataFrame([mo.to_collection_dict() for mo in mo_list],
                               columns=mo_list[0].to_collection_dict().keys())
        mobs_df.set_index('name', inplace=True)
        if verbose:
            print(mobs_df)

        return cls(mobs_df, name=modelname)

    @classmethod
    def from_modflow(cls, obs_collection, dis, hds_arr, mtime,
                     modelname='', exclude_layers=0, verbose=False):
        """ Read modflow data

        Parameters
        ----------
        obs_collection : ObsCollection
            locations of model observation
        dis : flopy.modflow.mfdis.ModflowDis
            discretisation package modflow
        hds_arr : numpy array
            heads with shape (ntimesteps, nlayers, nrow, ncol)
        mtime : list of datetimes
            dates for each model timestep
        modelname : str, optional
            modelname
        exclude_layers : int, optional
            exclude the observations up to these modellayers
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """

        if dis.sr.xll == 0 or dis.sr.yll == 0:
            warnings.warn(
                'you probably want to set the xll and/or yll attributes of dis.sr')

        xmid = dis.sr.xcentergrid.ravel()
        ymid = dis.sr.ycentergrid.ravel()

        xy = np.c_[xmid, ymid]
        uv = obs_collection.loc[:, ("x", "y")].dropna(how="any", axis=0).values
        vtx, wts = util.interp_weights(xy, uv)

        # get interpolated timeseries from hds_arr
        hm_ts = np.zeros((obs_collection.shape[0], hds_arr.shape[0]))

        # loop over layers
        for m in range(dis.nlay):
            if m < exclude_layers:
                continue
            mask = obs_collection.modellayer.values == m
            # loop over timesteps
            for t in range(hds_arr.shape[0]):
                ihds = hds_arr[t, m]
                ihds[ihds < -999.] = np.nan
                hm = util.interpolate(ihds, vtx, wts)
                hm_ts[mask, t] = hm[mask]

        mo_list = []
        for i, name in enumerate(obs_collection.index):
            mo = obs.ModelObs(index=mtime, data=hm_ts[i], name=name, model=modelname,
                              x=obs_collection.loc[name,
                                                   'x'], y=obs_collection.loc[name, 'y'],
                              meta=obs_collection.loc[name, 'obs'].meta)
            mo_list.append(mo)

        mobs_df = pd.DataFrame([mo.to_collection_dict() for mo in mo_list],
                               columns=mo_list[0].to_collection_dict().keys())
        mobs_df.set_index('name', inplace=True)

        return cls(mobs_df)

    @classmethod
    def from_list(cls, obs_list, name='', verbose=False):
        """read observations from a list of obs objects

        Parameters
        ----------
        obs_list : list of observation.Obs
            list of observations
        name : str, optional
            name of the observation collection
        verbose : boolean, optional
            Print additional information to the screen (default is False).


        """

        obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                              columns=obs_list[0].to_collection_dict().keys())
        obs_df.set_index('name', inplace=True)

        return cls(obs_df, name=name)

    @classmethod
    def from_pastas_project(cls, pr, ObsClass=obs.GroundwaterObs,
                            name=None, pr_meta=None):

        if name is None:
            name = pr.name

        if pr_meta is None:
            pr_meta = pr.file_info

        obs_list = []
        for index, row in pr.oseries.iterrows():
            metadata = row.to_dict()
            s = pd.DataFrame(metadata.pop('series').series_original)
            s.rename(columns={index: 'Stand_m_tov_NAP'}, inplace=True)

            o = ObsClass(s, meta=metadata, name=metadata['name'],
                         x=metadata['x'], y=metadata['y'],
                         locatie=metadata['locatie'],
                         filternr=metadata['filternummer'],
                         metadata_available=metadata['metadata_available'],
                         maaiveld=metadata['maaiveld'],
                         meetpunt=metadata['meetpunt'],
                         bovenkant_filter=metadata['bovenkant_filter'],
                         onderkant_filter=metadata['onderkant_filter'])
            obs_list.append(o)

        return cls(cls.from_list(obs_list, name=name), meta=pr_meta)

    @classmethod
    def from_wiski(cls, dirname, ObsClass=obs.GroundwaterObs, suffix='.csv',
                   unpackdir=None, force_unpack=False, preserve_datetime=False,
                   verbose=False, keep_all_obs=True, **kwargs):

        from . import io_wiski

        meta = {'dirname': dirname,
                'type': ObsClass,
                'suffix': suffix,
                'unpackdir': unpackdir,
                'force_unpack': force_unpack,
                'preserver_datetime': preserve_datetime,
                'verbose': verbose,
                'keep_all_obs': keep_all_obs
                }

        name = "wiski_import"

        obs_df = io_wiski.read_wiski_dir(dirname,
                                         ObsClass=ObsClass,
                                         suffix=suffix,
                                         unpackdir=unpackdir,
                                         force_unpack=force_unpack,
                                         preserve_datetime=preserve_datetime,
                                         verbose=verbose,
                                         keep_all_obs=keep_all_obs,
                                         **kwargs)

        return cls(obs_df, name=name, meta=meta)

    def data_frequency_plot(self, column_name='Stand_m_tov_NAP', intervals=None,
                            ignore=['seconde', 'minuut', '14-daags'],
                            normtype='log', cmap='viridis_r', set_yticks=False,
                            figsize=(10, 8), **kwargs):
        """plot data availability and frequency of observation collection

        Parameters
        ----------
        column_name : str, optional
            column name of the timeseries data
        intervals: dict, optional
            A dict with frequencies as keys and number of seconds as values
        ignore: list, optional
            A list with frequencies in intervals to ignore
        normtype: str, optional
            Determines the type of color normalisations, default is 'log'
        cmap: str, optional
            A reference to a matplotlib colormap
        set_yticks: bool, optional
            Setthe names of the series as ytciks
        figsize: tuple, optional
            The size of the new figure in inches (h,v)
        **kwargs: dict, optional
            Extra arguments are passed to matplotlib.pyplot.subplots()



        """

        art = util._import_art_tools()
        series_list = []
        for o in self.obs.values:
            series = o[column_name]
            series.name = o.name
            series_list.append(series)

        ax = art.data_availablity(series_list, intervals=intervals, ignore=ignore,
                                  normtype=normtype, cmap=cmap, set_yticks=set_yticks,
                                  figsize=figsize, **kwargs)

        return ax

    def get_bounding_box(self, xcol='x', ycol='y'):
        """returns the bounding box of all observations

        Parameters
        ----------
        xcol : str, optional
            column name with x values
        ycol : str, optional
            column name with y values

        Returns
        -------
        xmin, ymin, xmax, ymax
            coordinates of bouding box

        """

        xmin = self[xcol].min()
        xmax = self[xcol].max()
        ymin = self[ycol].min()
        ymax = self[ycol].max()

        return (xmin, ymin, xmax, ymax)

    def get_filternr(self, radius=1, xcol='x', ycol='y', overwrite=False):
        # ken filternummers toe aan peilbuizen die dicht bij elkaar staan
        straal_pb = 5

        if 'filternr' in self.columns and overwrite:
            raise RuntimeError(
                "the column 'filternr' already exist, set overwrite=True to replace the current column")

        self['filternr'] = 0
        for name in self.index:
            x = self.loc[name, xcol]
            y = self.loc[name, ycol]
            distance_to_other_filters = np.sqrt(
                (self[xcol]-x)*(self[xcol]-x) + (self[ycol]-y)*(self[ycol]-y))
            dup_x = self.loc[distance_to_other_filters < straal_pb]
            if dup_x.shape[0] == 1:
                self.loc[name, 'filternr'] = 1
            else:
                dup_x2 = dup_x.sort_values('onderkant_filter', ascending=False)
                for i, pb_dub in enumerate(dup_x2.index):
                    self.loc[pb_dub, 'filternummer'] = i+1

    def get_first_last_obs_date(self):
        """adds two columns to the ObsCollection with the date of the first
        and the last measurement


        """

        self['date_first_measurement'] = [o.index.min()
                                          for o in self.obs.values]
        self['date_last_measurement'] = [o.index.max()
                                         for o in self.obs.values]

    def get_min_max(self, obs_column='Stand_m_tov_NAP'):
        """adds two columns to the Obscollection with the minimum and the
        maximum of a column (defined by obs_column)

        returns the absolute minimum and maximum of all observations

        Parameters
        ----------
        obs_column: str, optional
            column name for which min and max is calculated

        Returns
        -------
        min, max: float
            the minimum and maximum of the obs_column in all observations

        """

        self['max'] = [o[obs_column].max() for o in self.obs.values]
        self['min'] = [o[obs_column].min() for o in self.obs.values]

        return(self['min'].min(), self['max'].max())

    def to_fieldlogger(self, fname, additional_collection=None, otype=None,
                       infer_otype=True,
                       use_default_otype_settings=True,
                       heading=None,
                       name=None, subname=None, inputfield=None,
                       properties=None, group=None, group_color='blue',
                       verbose=False):
        """Write a csv file that can be read by fieldlogger

        Notes
        -----
        1. option to change input fields for measurements is not yet added.
        2. option to add min/max values for measurements is not yet added.

        Parameters
        ----------
        fname : str
            name of the fieldlogger csv file
        additional_collections : list, optional
            additional ObsCollections that are added to the fieldlogger csv
        otype : str, optional
            observation type, 'groundwater' and 'surfacewater' are supported
        infer_otype : boolean, optional
            infer mtype from the obs column in self.obs, default is True
        use_default_otype_settings : boolean, optional
            use the default settings for this measurement type
        heading : str, optional
            heading of the csv, if None use default heading
        name : series or str, optional
            1st column in fieldlogger csv, if None use index of ObsCollection
        subname : series or str, optional
            2nd column in fieldlogger csv, if None use name + _sub
        inputfield : series or str, optional
            3th column in fieldlogger csv, if None use 'Stand|Opmerking'
        properties : series or str, optional
            4th column in fieldlogger csv, if None use empty string
        group : series or str, optional
            5th column in fieldlogger csv, if None do not add group to csv
        group_color : series or str, optional
            color of the group
        """
        from . import fieldlogger

        if infer_otype:
            otype = self._infer_otype(verbose)
            if type(otype) != type:  # check of dit ook echt een type is
                raise TypeError(
                    'could not infer observation type, use infer_otype=False')

        f_df = fieldlogger.df_to_fieldlogger_csv(self, fname, otype,
                                                 use_default_otype_settings,
                                                 heading,
                                                 name, subname, inputfield,
                                                 properties, group, group_color,
                                                 verbose)

        return f_df

    def to_interactive_plots(self, savedir,
                             tmin=None, tmax=None,
                             verbose=False, **kwargs):
        """Create interactive plots of the observations using bokeh

        Parameters
        ----------
        savedir : str
            directory used for the folium map and bokeh plots
        tmin : dt.datetime, optional
            start date for timeseries plot
        tmax : dt.datetime, optional
            end date for timeseries plot

        verbose : boolean, optional
            Print additional information to the screen (default is False).
        **kwargs :
            will be passed to the Obs.to_interactive_plot method

            plot_columns : list of str
                name of the column in the obs df that will be plotted with bokeh
            hoover_names : list of str, optional
                names will be displayed together with the plot_column value
                when hoovering over plot
            plot_freq : str, optional
                bokeh plot is resampled with this frequency to reduce the size
            plot_legend_name : str, optional
                legend in bokeh plot
            ylabel : str, optional
                label on the y-axis
            color : str, optional
                color of the lines on the plots
            add_filter_to_legend : boolean, optional
                if True the attributes bovenkant_filter and onderkant_filter
                are added to the graph legend

        Returns
        -------

        """

        _color_cycle = ('blue', 'olive', 'lime', 'red', 'orange', 'yellow', 'purple', 'silver', 'powderblue', 'salmon', 'tan')
        _same_loc_list = []
        for o in self.obs.values:
            # check for multiple observations at the same location (usually multiple filters)
            if (self.x == o.x).sum() == 1:
                # one observation point (filter) at this location
                try:
                    o.iplot_fname = o.to_interactive_plot(savedir=savedir,
                                                          p=None,
                                                          tmin=tmin, tmax=tmax,
                                                          return_filename=True,
                                                          **kwargs)
                except ValueError:
                    if verbose:
                        print('{} has no data between {} and {}'.format(
                            o.name, tmin, tmax))
                    o.iplot_fname = None

            # check if observation was not added to another plot yet
            elif o.name != _same_loc_list:
                # plot multiple filters in same graph
                try:
                    p = o.to_interactive_plot(savedir=None,
                                              p=None,
                                              tmin=tmin, tmax=tmax,
                                              return_filename=False,
                                              **kwargs)
                except ValueError:
                    if verbose:
                        print('{} has no data between {} and {}'.format(
                            o.name, tmin, tmax))
                    o.iplot_fname = None
                    continue

                same_pb = self[(self.x == o.x) & (self.index != o.name)]

                for i, o_pb in enumerate(same_pb.obs.values):
                    if i == 10:
                        raise NotImplementedError(
                            'cannot add more than 4 lines to a single plot')
                    try:
                        o.iplot_fname = o_pb.to_interactive_plot(savedir=savedir,
                                                                 p=p,
                                                                 tmin=tmin, tmax=tmax,
                                                                 colors=[
                                                                     _color_cycle[i+1]],
                                                                 return_filename=True,
                                                                 **kwargs)
                    except ValueError:
                        if verbose:
                            print('{} has no data between {} and {}'.format(
                                o.name, tmin, tmax))
                        o.iplot_fname = None
                    _same_loc_list.append(o_pb.name)

            else:
                o.iplot_fname = None

    def to_interactive_map(self, plot_dir, m=None,
                           tiles='OpenStreetMap', fname=None,
                           color='blue', legend_name=None,
                           add_legend=True,
                           map_label='', map_label_size=20,
                           col_name_lat='lat', col_name_lon='lon',
                           zoom_start=13,
                           create_interactive_plots=True,
                           verbose=False, **kwargs):
        """ create an interactive map with interactive plots using folium
        and bokeh.

        Notes
        -----
        - if you want to have multiple obs collections on one folium map, only
        the last one should have add_legend = True to create a correct legend
        - the color of the observation point on the map is now the same color
        as the line of the observation measurements. Also a built-in color
        cycle is used for different measurements on the same location.


        Parameters
        ----------
        plot_dir : str
            directory used for the folium map and bokeh plots
        m : folium.Map, str, optional
            current map to add observations too, if None a new map is created
        tiles : str, optional
            background tiles, default is openstreetmap
        fname : str, optional
            name of the folium map
        color : str, optional
            color of the observation points on the map
        legend_name : str, optional
            the name of the observation points shown in the map legend
        add_legend : boolean, optional
            add a legend to a plot
        map_label : str, optional
            add a label to the obs locations on the map, this label is
            picked from the meta attribute of the obs points.
        map_label_size : int, optional
            label size of the map_label in pt.
        col_name_lat : str, optional
            name of the column in the obs_collection dic with the lat values
            of the observation points
        col_name_lon : str, optional
            see col_name_lat
        zoom_start : int, optional
            start zoom level of the folium ma
        create_interactive_plots : boolean, optinal
            if True interactive plots will be created, if False the iplot_fname
            attribute of the observations is used.
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        **kwargs :
            will be passed to the to_interactive_plots method options are:

            plot_columns : list str, optional
                name of the column in the obs df that will be plotted with bokeh
            hoover_names : list of str, optional
                names will be displayed together with the plot_column value
                when hoovering over plot
            plot_legend_name : str, optional
                the name of the observation points in the time series plot
            ylabel : str, optional
                label on the y-axis
            add_filter_to_legend : boolean, optional
                if True the attributes bovenkant_filter and onderkant_filter
                are added to the graph title
            plot_freq : str, optional
                bokeh plot is resampled with this frequency to reduce the size
                of the complete folium map
            tmin : dt.datetime, optional
                start date for timeseries plot
            tmax : dt.datetime, optional
                end date for timeseries plot

        Returns
        -------
        m : folium.Map
            the folium map
        """

        import branca
        import folium
        from folium.features import DivIcon

        # create interactive bokeh plots
        if create_interactive_plots:
            self.to_interactive_plots(savedir=plot_dir,
                                      verbose=verbose,
                                      **kwargs)

        # determine start location of map
        northing = np.mean(
            (self[col_name_lat].min(), self[col_name_lat].max()))
        easting = np.mean((self[col_name_lon].min(), self[col_name_lon].max()))

        # create map if no map is given
        if m is None:
            m = folium.Map([northing, easting], zoom_start=zoom_start)

        # add the point observations with plots to the map
        group_name = '<span style=\\"color: {};\\">{}</span>'.format(
            color, legend_name)
        group = folium.FeatureGroup(name=group_name)
        for o in self.obs.values:
            if o.iplot_fname is not None:
                with open(o.iplot_fname, 'r') as f:
                    bokeh_html = f.read()

                iframe = branca.element.IFrame(
                    html=bokeh_html, width=620, height=420)
                popup = folium.Popup(iframe, max_width=620)

                folium.CircleMarker([o.meta[col_name_lat], o.meta[col_name_lon]],
                                    icon=folium.Icon(icon='signal'), fill=True,
                                    color=color,
                                    popup=popup).add_to(group)

                if map_label is not '':
                    folium.map.Marker([o.meta[col_name_lat], o.meta[col_name_lon]],
                                      icon=DivIcon(icon_size=(150, 36),
                                                   icon_anchor=(0, 0),
                                                   html='<div style="font-size: %ipt">%s</div>' % (map_label_size, o.meta[map_label]))).add_to(group)
            elif verbose:
                print('no iplot available for {}'.format(o.name))

        group.add_to(m)

        # add legend
        if add_legend:
            folium.map.LayerControl('topright', collapsed=False).add_to(m)

        # save map
        #filename and path
        if fname is None:
            fname = self.name+'.html'
        else:
            if not fname.endswith('.html'):
                fname = fname + '.html'
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)
            m.save(os.path.join(plot_dir, fname))

        return m

    def to_gdf(self, xcol='x', ycol='y'):
        """convert ObsCollection to GeoDataFrame

        Parameters
        ----------
        xcol : str
            column name with x values
        ycol : str
            column name with y values

        Returns
        -------
        gdf : geopandas.GeoDataFrame

        """
        art = util._import_art_tools()
        return art.util.df2gdf(self, xcol, ycol)

    def to_map(self, ax=None, figsize=(15, 15), label='gws',
               edgecolor='black', facecolor='green',
               marker='o', markersize=100, xlim=None, ylim=None, add_topo=True):
        """plot observation points on a map

        Parameters
        ----------
        ax : matplotlib.axes
            reference to axes object
        figsize : tuple
            figure size
        label : str
            label used in legend
        edgecolor : str

        facecolor : str

        marker : str

        markersize : int

        xlim : tuple

        ylim : tuple

        add_topo : boolean, optional
            topo is added using art tools

        Returns
        -------
        ax : matplotlib.axes
            reference to axes object

        """

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)

        ax.scatter(self.x, self.y, label=label, edgecolor=edgecolor,
                   marker=marker, facecolor=facecolor,
                   s=markersize)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if add_topo:
            import art_tools as art
            art.OpenTopo(ax=ax).plot(output=False, alpha=0.5)

        ax.legend()

        return ax

    def to_mapgraphs(self, graph=None, plots_per_map=10, figsize=(16, 10),
                     extent=None, plot_column='Stand_m_tov_NAP',
                     plot_xlim=None, plot_ylim=(None,None), min_dy=2.0, 
                     vlines=[],
                     savefig=None, map_gdf=None, map_gdf_kwargs={}, 
                     verbose=True):
        """make mapgraph plots of obs collection data

        Parameters
        ----------
        graph : np.ndarray, optional
            passed to the art.MapGraph funcion
        plots_per_map : int, optional
            number of plots per mapgraph, should be consistent with graph argument
        figsize : tuple, optional
            figure size
        extent : list or tuple, optional
            extent of the map [xmin, xmax, ymin, ymax]
        plot_column : str, optional
            column name of data to be plotted in graphs
        plot_xlim : datetime, optional
            xlim of the graphs
        plot_ylim : tuple or str, optional
            the ylim of all the graphs or 
            'max_dy' -> all graphs get the maximum ylim
            'min_dy' -> all graphs get the same dy unless the series has a dy
                        bigger than min_dy
        min_dy : float, optional
            only used if plot_ylim is 'min_dy'
        vlines : list, optional
            list with x values on which a vertical line is plotted
        savefig : str, optional
            path to save the figure. {plot_number}.png will be append to this path
        map_gdf : GeoDataFrame, optional
            if not None the GeoDataFrame is also plotted on the map
        map_gdf_kwargs : dic, optional
            kwargs that will be passed to the map_gdf.plot() method


        Returns
        -------
        mg_list : list of art_tools.plots.MapGraph
            list of mapgraph objects

        """

        import art_tools as art

        if graph is None:
            graph = np.full((4, 4), True)
            graph[1:, 1:-1] = False
            
        if plot_ylim == 'max_dy':
            plot_ylim = self.get_min_max(obs_column=plot_column)
        
        mg_list = []
        for k, xyo in self.groupby(np.arange(self.shape[0])//plots_per_map):
        
            mg = art.MapGraph(xy=xyo[['x', 'y']].values, graph=graph, 
                              figsize=figsize, extent=extent)
        
            plt.sca(mg.mapax)
            plt.yticks(rotation=90, va="center")
            art.OpenTopo(ax=mg.mapax, verbose=verbose).plot(alpha=0.75, 
                        verbose=verbose)
            if map_gdf is not None:
                map_gdf.plot(ax=mg.mapax, **map_gdf_kwargs)
        
            for i, o in enumerate(xyo.obs.values):
                ax = mg.ax[i]
                pb = o.name
                try:
                    o[plot_column][plot_xlim[0]:plot_xlim[1]].plot(ax=ax, lw=.5, 
                     marker='.', markersize=1.,label=pb)
                except TypeError:
                    o[plot_column].plot(ax=ax, lw=.5, marker='.', 
                                         markersize=1.,label=pb)
                    ax.set_xlim(plot_xlim)
                if plot_ylim == 'min_dy':
                    ax.autoscale(axis='y', tight=True)
                    ylim = ax.get_ylim()
                    dy = ylim[1]-ylim[0]
                    if dy<min_dy:
                        ylim=(ylim[0]-(min_dy-dy)/2, ylim[1]+(min_dy-dy)/2)
                    ax.set_ylim(ylim)
                else:
                    ax.set_ylim(plot_ylim)
                    
                for line in vlines:
                    ylim = ax.get_ylim()
                    ax.vlines(line, ylim[0]-100, ylim[1]+100, ls='--', 
                              lw=2.0, color='r')

        
                ax.legend(loc='best')
 
    
            mg.figure.tight_layout()
            if savefig is not None:
                mg.figure.savefig(savefig+"{0}.png".format(k),
                                  dpi=300, bbox_inches="tight")
            mg_list.append(mg)

        return mg_list

    def to_report_table(self, columns=['locatie', 'filternr',
                                       'Van', 'Tot', '# metingen']):

        if 'Van' in columns:
            self['Van'] = self.obs.apply(lambda x: x.index[0])
        if 'Tot' in columns:
            self['Tot'] = self.obs.apply(lambda x: x.index[-1])
        if '# metingen' in columns:
            self['# metingen'] = self.obs.apply(lambda x: x.shape[0])

        return self[columns]

    def to_pastas_project(self, pr=None, project_name='',
                          obs_column='Stand_m_tov_NAP',
                          kind='oseries', verbose=False):
        """add observations to a new or existing pastas project

        Parameters
        ----------
        pr : pastas.project, optional
            Existing pastas project, if None a new project is created
        project_name : str, optional
            Name of the pastas project only used if pr is None
        obs_column : str, optional
            Name of the column in the Obs dataframe to be used
        kind : str, optional
            The kind of series that is added to the pastas project
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        pr : pastas.project
            the pastas project with the series from the ObsCollection
        """

        import pastas as ps

        if pr is None:
            pr = ps.Project(project_name)

        for o in self.obs.values:
            if verbose:
                print('add to pastas project ->{}'.format(o.name))
            series = ps.TimeSeries(o[obs_column], name=o.name, metadata=o.meta)
            pr.add_series(series, kind=kind)

        return pr

    def to_shapefile(self, fname, xcol='x', ycol='y'):
        """save ObsCollection as shapefile

        Parameters
        ----------
        fname : str
            filename of shapefile, ends with .shp
        xcol : str
            column name with x values
        ycol : str
            column name with y values

        """
        art = util._import_art_tools()
        gdf = art.util.df2gdf(self, xcol, ycol)
        gdf.to_file(fname)

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326',
                    add_to_meta=True, add_to_df=True):
        """get lattitude and longitude from x and y attributes

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon
        add_to_meta : boolean, optional
            if True the new coordinates are added to the meta dictionary
        add_to_df : boolean, optional
            if True lon and lat columns are added to the ObsCollection


        """

        for o in self.obs.values:
            o.get_lat_lon(in_epsg, out_epsg, add_to_meta)

        if add_to_df:
            self.add_meta_to_df('lat')
            self.add_meta_to_df('lon')

    def get_nearest_point(self, obs_collection2=None, gdf2=None,
                          xcol_obs1='x', ycol_obs1='y',
                          xcol_obs2='x', ycol_obs2='y', verbose=False):
        """get nearest point of another obs collection for each
        point in the current obs collection.

        Parameters
        ----------
        obs_collection2 : ObsCollection, optional
            collection of observations of which the nearest point is found
        gdf2 : GeoDataFrame, optional
            dataframe to look for nearest observation point
        xcol_obs1 : str, optional
            x column in self used to get geometry
        ycol_obs1 : str, optional
            y column in self used to get geometry
        xcol_obs2 : str, optional
            x column in obs_collection2 used to get geometry
        ycol_obs2 : str, optional
            y column in self used to get geometry
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        pandas.DataFrame
            with columns 'nearest point' and 'distance nearest point'
        """

        from shapely.ops import nearest_points

        gdf1 = self.to_gdf(xcol=xcol_obs1, ycol=ycol_obs1)
        if obs_collection2 is not None:
            gdf2 = obs_collection2.to_gdf(xcol=xcol_obs2, ycol=ycol_obs2)
        elif gdf2 is None:
            raise ValueError('obs_collecction2 or gdf2 should be defined')

        pts_gdf2 = gdf2.geometry.unary_union

        def nearest_point(point_gdf1, pts=pts_gdf2):
            # find the nearest point and return the corresponding Place value
            nearest_point_gdf2 = nearest_points(point_gdf1, pts_gdf2)[1]
            nearest = gdf2[gdf2.geometry == nearest_point_gdf2]
            return nearest.index.values[0]

        def distance_nearest_point(point_gdf1, pts=pts_gdf2):
            # find the nearest point and return the corresponding Place value
            nearest_point_gdf2 = nearest_points(point_gdf1, pts_gdf2)[1]
            distance = point_gdf1.distance(nearest_point_gdf2)
            return distance

        gdf1['nearest point'] = gdf1.apply(
            lambda row: nearest_point(row.geometry), axis=1)
        gdf1['distance nearest point'] = gdf1.apply(
            lambda row: distance_nearest_point(row.geometry), axis=1)

        return gdf1[['nearest point', 'distance nearest point']]

    def get_pb_modellayers(self, dis, zgr=None, verbose=False):
        """Get the modellayers from the dis file

        Parameters
        ----------
        dis : flopy.modflow.ModflowDis
            object containing DIS of model (grid information)
        zgr : np.3darray, optional
            array containing model layer elevation
            information (if None , this information is obtained
            from the dis object)
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """

        for o in self.obs.values:
            o.get_pb_modellayer(dis, zgr, verbose)
            if verbose:
                print(o.name)

        self.add_meta_to_df('modellayer')

    def within_extent(self, extent, inplace=False):
        """Slice dataframe by extent

        Parameters
        ----------
        extent : tuple
            format (xmin, xmax, ymin, ymax), default dis.sr.get_extent() format

        Returns
        -------
        new_oc : obs_collection.ObsCollection
            dataframe with collection of observations within extent
        """

        new_oc = self[(self.x > extent[0]) & (self.x < extent[1])
                      & (self.y > extent[2]) & (self.y < extent[3])]
        if inplace:
            self._update_inplace(new_oc)
        else:
            return new_oc

    def add_meta_to_df(self, name):
        """Add data from the metadata dictionary as a column

        to the ObsCollection

        Parameters
        ----------
        name : str
            variable name in metadata dictionary


        """

        self[name] = [o.meta[name] for o in self.obs.values]

    @property
    def n_observations(self):
        return self.obs.apply(lambda o: o.shape[0])

    @property
    def dates_first_obs(self):
        return self.obs.apply(lambda o: o.index[0])

    @property
    def dates_last_obs(self):
        return self.obs.apply(lambda o: o.index[-1])

    @property
    def obs_periods(self):
        return self.dates_last_obs - self.dates_first_obs

    def obs_per_year(self, col="Stand_m_tov_NAP"):
        pblist = {o.name: o.obs_per_year(col=col) for o in self.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def consecutive_obs_years(self, min_obs=12, col="Stand_m_tov_NAP"):
        pblist = {o.name: o.consecutive_obs_years(
            min_obs=min_obs, col=col) for o in self.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def mean_in_period(self, tmin=None, tmax=None, col="Stand_m_tov_NAP"):
        if tmin is None:
            tmin = self.dates_first_obs.min()
        if tmax is None:
            tmax = self.dates_last_obs.max()

        return self.obs.apply(lambda o: o.loc[tmin:tmax, col].mean())

    def get_series(self, tmin=None, tmax=None, col="Stand_m_tov_NAP"):
        if tmin is None:
            tmin = self.dates_first_obs.min()
        if tmax is None:
            tmax = self.dates_last_obs.max()
        return self.obs.apply(lambda o: o.loc[tmin:tmax, col])
