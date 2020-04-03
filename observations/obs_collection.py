"""
module with ObsCollection class for a collection of observations.

The ObsCollection class is a subclass of a pandas DataFrame with
additional attributes and methods.

More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

"""
import os
import warnings

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
        # self.plots = CollectionPlots(self)

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

    def _set_metadata_value(self, iname, att_name, value, add_to_meta=False,
                            verbose=False):
        """ Set a value on three different levels at once:
            1. the value in an ObsCollection DataFrame
            2. the attribute of the observation
            3. the value in the meta dictionary of an observation

        Parameters
        ----------
        iname : str, int, float, ...
            observation name. Must be same type as self.index.
            e.g. B52D0111_3
        att_name : str, int, float, ...
            name of the column in self.columns and attribute
            of the observation. e.g. 'x'
        value : str, int, float, ...
            value of the the att_name. e.g. 116234
        add_to_meta : bool, optional
            if True the att_name, value pair is added to the meta dictionary
            of an observation. The default is False.
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Raises
        ------
        ValueError
            if the iname is not in self.index the value cannot be set.

        Returns
        -------
        None.

        """
        if iname not in self.index:
            raise ValueError(f"{iname}  not in index")

        self.loc[iname, att_name] = value
        if verbose:
            print(f'set {iname}, {att_name} to {value}')

        o = self.loc[iname, 'obs']
        if att_name in o._metadata:
            setattr(o, att_name, value)
            if verbose:
                print(f'set attribute {att_name} of {iname} to {value}')
        if add_to_meta:
            if verbose:
                print(f'set attribute {att_name} of {iname} to {value}')
            o.meta.update({att_name: value})
            
    @classmethod
    def from_arctic(cls, connstr, libname, ObsClass=obs.GroundwaterObs,
                    verbose=False):
        """Load ObsCollection from MongoDB using arctic

        Parameters
        ----------
        connstr : str
            database connection string
        libname : str
            library name
        ObsClass : class, optional
            observation class to store single timeseries, by
            default obs.GroundwaterObs
        verbose : bool, optional
            show progress bar and database read summary, by default False

        Returns
        -------
        ObsCollection
            ObsCollection DataFrame containing all the obs

        """
        from .io.io_arctic import read_arctic

        meta = {'type': obs.GroundwaterObs}
        obs_list = read_arctic(connstr, libname, ObsClass, verbose=verbose)
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, meta=meta)

    @classmethod
    def from_dino(cls, dirname=None,
                  extent=None, bbox=None,
                  ObsClass=obs.GroundwaterObs,
                  subdir='Grondwaterstanden_Put',
                  suffix='1.csv',
                  unpackdir=None,
                  force_unpack=False,
                  preserve_datetime=False,
                  keep_all_obs=True,
                  name=None,
                  verbose=False,
                  **kwargs):
        """ Read dino data within an extent from the server or from a 
        directory with downloaded files.
        
        Parameters
        ----------
        dirname : str, optional
            directory name, can be a .zip file or the parent directory of subdir
        extent : list, tuple or numpy-array (user must specify extent or bbox)
            get dinodata online within this extent [xmin, xmax, ymin, ymax]
        bbox : list, tuple or numpy-array (user must specify extent or bbox)
            The bounding box, in RD-coordinates, for which you want to retreive locations
            [xmin, ymin, xmax, ymax]
        ObsClass : type
            class of the observations, so far only GroundwaterObs is supported
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
            add all observation points to the collection, even without data or
            metadata
        name : str, optional
            the name of the observation collection
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs:
            kwargs are passed to the io_dino.download_dino_within_extent() or
            the io_dino.read_dino_dir() function

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """
        from .io.io_dino import read_dino_dir, download_dino_within_extent

        if dirname is not None:
            # read dino directory
            if name is None:
                name = subdir

            meta = {'dirname': dirname,
                    'type': ObsClass,
                    'suffix': suffix,
                    'unpackdir': unpackdir,
                    'force_unpack': force_unpack,
                    'preserve_datetime': preserve_datetime,
                    'verbose': verbose,
                    'keep_all_obs': keep_all_obs
                    }

            obs_list = read_dino_dir(dirname,
                                    ObsClass,
                                    subdir,
                                    suffix,
                                    unpackdir,
                                    force_unpack,
                                    preserve_datetime,
                                    verbose,
                                    keep_all_obs,
                                    **kwargs)

        elif extent is not None or bbox is not None:
            # read dino data within extent
            if ObsClass == obs.GroundwaterObs:
                layer = 'grondwatermonitoring'
            else:
                raise NotImplementedError(
                    'cannot download {} from Dino'.format(ObsClass))

            if name is None:
                name = '{} from DINO'.format(layer)

            meta = kwargs.copy()
            meta.update({'verbose': verbose,
                         'extent': extent,
                         'bbox': bbox,
                         'layer': layer,
                         'keep_all_obs': keep_all_obs,
                         'verbose': verbose})

            obs_list = download_dino_within_extent(
                extent=extent, bbox=bbox, ObsClass=ObsClass, layer=layer,
                keep_all_obs=keep_all_obs, verbose=verbose, **kwargs)

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, bbox=bbox, meta=meta)

    @classmethod
    def from_dino_server(cls, extent=None, bbox=None,
                          ObsClass=obs.GroundwaterObs,
                          name=None, keep_all_obs=True,
                          verbose=False, **kwargs
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
        keep_all_obs : boolean, optional
            add all observation points to the collection, even without data or
            metadata
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs:
            kwargs are passed to the io_dino.download_dino_within_extent() function

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations

        """

        warnings.warn("this method will be removed in future versions, use from_dino instead", DeprecationWarning)

        from .io.io_dino import download_dino_within_extent

        if ObsClass == obs.GroundwaterObs:
            layer = 'grondwatermonitoring'
        else:
            raise NotImplementedError(
                'cannot download {} from Dino'.format(ObsClass))

        if name is None:
            name = '{} from DINO'.format(layer)

        meta = kwargs.copy()
        meta.update({'verbose': verbose})

        obs_list = download_dino_within_extent(
            extent=extent, bbox=bbox, ObsClass=ObsClass, layer=layer,
            keep_all_obs=keep_all_obs, verbose=verbose, **kwargs)

        obs_df = util._obslist_to_frame(obs_list)

        if bbox is None:
            bbox = [extent[0], extent[2], extent[1], extent[3]]

        return cls(obs_df, name=name, bbox=bbox, meta=meta)

    @classmethod
    def from_dino_dir(
            cls,
            dirname=None,
            ObsClass=obs.GroundwaterObs,
            subdir='Grondwaterstanden_Put',
            suffix='1.csv',
            unpackdir=None,
            force_unpack=False,
            preserve_datetime=False,
            keep_all_obs=True,
            name=None,
            verbose=False,
            **kwargs):
        """ Read a dino directory

        Parameters
        ----------
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
            add all observation points to the collection, even without data or
            metadata
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

        warnings.warn("this method will be removed in future versions, use from_dino instead", DeprecationWarning)

        from .io.io_dino import read_dino_dir

        if name is None:
            name = subdir

        meta = {'dirname': dirname,
                'type': ObsClass,
                'suffix': suffix,
                'unpackdir': unpackdir,
                'force_unpack': force_unpack,
                'preserve_datetime': preserve_datetime,
                'verbose': verbose,
                'keep_all_obs': keep_all_obs
                }

        obs_list = read_dino_dir(
            dirname,
            ObsClass,
            subdir,
            suffix,
            unpackdir,
            force_unpack,
            preserve_datetime,
            verbose,
            keep_all_obs,
            **kwargs)

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_artdino_dir(
            cls,
            dirname=None,
            ObsClass=obs.GroundwaterObs,
            subdir='csv',
            suffix='.csv',
            unpackdir=None,
            force_unpack=False,
            preserve_datetime=False,
            keep_all_obs=True,
            name=None,
            verbose=False,
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
            add all observation points to the collection, even without data or
            metadata
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

        from .io.io_dino import read_artdino_dir

        if name is None:
            name = subdir

        meta = {'dirname': dirname,
                'type': ObsClass,
                'suffix': suffix,
                'unpackdir': unpackdir,
                'force_unpack': force_unpack,
                'preserve_datetime': preserve_datetime,
                'verbose': verbose,
                'keep_all_obs': keep_all_obs
                }

        obs_list = read_artdino_dir(
            dirname,
            ObsClass,
            subdir,
            suffix,
            unpackdir,
            force_unpack,
            preserve_datetime,
            verbose,
            keep_all_obs,
            **kwargs)

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews(cls, file_or_dir, ObsClass=obs.GroundwaterObs, name='fews',
                  translate_dic={'locationId': 'locatie'}, locations=None,
                  to_mnap=True, remove_nan=True, low_memory=True,
                  unpackdir=None, force_unpack=False,
                  preserve_datetime=False, verbose=False):
        """Read one or several XML-files with measurements from FEWS.

        Parameters
        ----------
        file_or_dir :  str
            zip, xml or directory with zips or xml files to read
        ObsClass : type
            class of the observations, e.g. GroundwaterObs or WaterlvlObs
        name : str, optional
            name of the observation collection, 'fews' by default
        translate_dic : dict
            translate name of attribute by passing key: value pairs in
            dictionary
        locations : list of str, optional
            list of locationId's to read from XML file, others are skipped.
            If None (default) all locations are read. Only supported by
            low_memory=True method!
        to_mnap : boolean, optional
            if True a column with 'stand_m_tov_nap' is added to the dataframe
        remove_nan : boolean, optional
            remove nan values from measurements, flag information about the
            nan values is also lost
        low_memory : bool, optional
            whether to use xml-parsing method with lower memory footprint,
            default is True
        unpackdir : str
            destination directory to unzip file if fname is a .zip
        force_unpack : boolean, optional
            force unpack if dst already exists
        preserve_datetime : boolean, optional
            whether to preserve datetime from zip archive
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations

        """
        from .io.io_xml import parse_xml_filelist

        # get files
        dirname, unzip_fnames = util.get_files(file_or_dir, ext=".xml",
                                               force_unpack=force_unpack,
                                               preserve_datetime=preserve_datetime)
        meta = {'filename': dirname,
                'type': ObsClass,
                'verbose': verbose
                }

        obs_list = parse_xml_filelist(unzip_fnames,
                                      ObsClass,
                                      directory=dirname,
                                      translate_dic=translate_dic,
                                      locations=locations,
                                      to_mnap=to_mnap,
                                      remove_nan=remove_nan,
                                      low_memory=low_memory,
                                      verbose=verbose)
        obs_df = util._obslist_to_frame(obs_list)
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

        from .io.io_fieldlogger import fieldlogger_csv_to_obs_list

        obs_list, fieldlogger_meta = fieldlogger_csv_to_obs_list(
            fname, ObsClass=obs.GroundwaterObs)
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, meta=fieldlogger_meta)
    
    @classmethod
    def from_imod(cls, obs_collection, ml, runfile, mtime, model_ws,
                  modelname='', nlay=None, exclude_layers=0, verbose=False):
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
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """
        from .io.io_modflow import read_imod_results
        mo_list = read_imod_results(obs_collection, ml, runfile,
                                    mtime, model_ws,
                                    modelname=modelname, nlay=nlay,
                                    exclude_layers=exclude_layers,
                                    verbose=verbose)
        obs_df = util._obslist_to_frame(mo_list)
        return cls(obs_df, name=modelname)
    
    @classmethod
    def from_knmi(cls, locations=None, stns=None, 
                  xmid=None, ymid=None,
                  meteo_vars=("RD"), name='', 
                  start=[None, None], end=[None, None],
                  ObsClass=obs.KnmiObs,
                  **kwargs
                  ):
        """ get knmi observations from a list of locations or a list of
        stations
        

        Parameters
        ----------
        locations : pd.DataFrame or None
            dataframe with x and y coordinates. The default is None
        stns : list of str or None
            list of knmi stations. The default is None
        xmid : np.array, optional
            x coördinates of the cell centers of your grid shape(ncol)
        ymid : np.array, optional
            y coördinates of the cell centers of your grid shape(nrow)
        meteo_vars : list or tuple of str
            meteo variables e.g. ["RD", "EV24"]. The default is ("RD")
        name : str, optional
            name of the obscollection. The default is ''.
        start : list of str, datetime or None]
            start date of observations per meteo variable. The default is 
            [None, None]
        end : list of str, datetime or None]
            end date of observations per meteo variable. The default is 
            [None, None]
        ObsClass : type or None
            class of the observations, only KnmiObs is supported for now. The 
            default is None
        **kwargs : 
            kwargs are passed to the io_knmi.get_knmi_obslist function
        """
        
        from .io.io_knmi import get_knmi_obslist
        
        meta = {}
        meta['start'] = start
        meta['end'] = end
        meta['name'] = name
        meta['ObsClass'] = ObsClass
        meta['meteo_vars'] = meteo_vars
        
        obs_list = get_knmi_obslist(locations, stns, 
                                    xmid, ymid,
                                    meteo_vars,
                                    ObsClass=ObsClass,
                                    start=start,
                                    end=end, **kwargs)
        
        obs_df = util._obslist_to_frame(obs_list)
        
        return cls(obs_df, name=name, meta=meta)
    
    
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
        obs_df = util._obslist_to_frame(obs_list)
        return cls(obs_df, name=name)

    @classmethod
    def from_menyanthes(cls, fname, name='', ObsClass=obs.GroundwaterObs,
                        verbose=False):

        from .io.io_menyanthes import read_file

        menyanthes_meta = {'filename': fname,
                           'type': ObsClass,
                           'verbose': verbose}

        obs_list = read_file(fname, ObsClass, verbose)
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, meta=menyanthes_meta)

    

    @classmethod
    def from_modflow(cls, obs_collection, ml, hds_arr, mtime,
                     modelname='', nlay=None, exclude_layers=0, verbose=False):
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
        exclude_layers : int, optional
            exclude the observations up to these modellayers
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """
        from .io.io_modflow import read_modflow_results
        mo_list = read_modflow_results(obs_collection, ml, hds_arr,
                                       mtime, modelname=modelname,
                                       nlay=nlay,
                                       exclude_layers=exclude_layers,
                                       verbose=verbose)
        obs_df = util._obslist_to_frame(mo_list)

        return cls(obs_df)

    @classmethod
    def from_pastas_project(cls, pr, ObsClass=obs.GroundwaterObs,
                            name=None, pr_meta=None, rename_dic={}):

        from .io.io_pastas import read_project

        if pr_meta is None:
            pr_meta = pr.file_info

        if name is None:
            name = pr.name

        obs_list = read_project(pr, obs.GroundwaterObs,
                                rename_dic=rename_dic)
        obs_df = util._obslist_to_frame(obs_list)
        return cls(obs_df, meta=pr_meta, name=name)


    @classmethod
    def from_pystore(cls, storename, pystore_path,
                     ObsClass=obs.GroundwaterObs,
                     extent=None, collection_names=None,
                     item_names=None, nameby="item",
                     read_series=True, verbose=True, progressbar=False):
        """Create ObsCollection from pystore store

        Parameters
        ----------
        storename : str
            name of the store
        pystore_path : str
            path in which stores are stored
        ObsClass : type of Obs, optional
            the Observation type used for reading in the individual series,
            by default obs.GroundwaterObs
        extent : list or tuple, optional
            if not None only Observations within this extent are read
            [xmin, xmax, ymin, ymax]
        collection_names : list of str
            collection names that will be extracted from the store. Default is
            None which reads all collections.
        item_names : list of str
            item (Observation) names that will be extracted from the store.
            The other items (Observations) will be ignored. if None all items
            are read.
        read_series : bool, optional
            if False, read only metadata, default is True which
            loads the full dataset
        nameby : str
            pick whether Obs are named by 'item' or 'collection'
        progressbar : bool, optional
            whether to show progress bar, default is False, for
            best result set verbose to False!

        Returns
        -------
        ObsCollection
            Collection of observations

        """

        from .io.io_pystore import read_pystore

        meta = {'fname': os.path.join(pystore_path, storename),
                'verbose': verbose}

        obs_list = read_pystore(storename, pystore_path,
                                obs.GroundwaterObs,
                                extent=extent,
                                collection_names=collection_names,
                                item_names=item_names, nameby="item",
                                read_series=read_series, verbose=verbose,
                                progressbar=progressbar)
        # if read series is False, returns dataframe with only metadata
        if isinstance(obs_list, list):
            obs_df = util._obslist_to_frame(obs_list)
            meta['type'] = obs.GroundwaterObs
        else:
            # only metadata in df
            obs_df = obs_list

        return cls(obs_df, name=storename, meta=meta)


    @classmethod
    def from_waterinfo(cls, file_or_dir, name="", ObsClass=obs.WaterlvlObs,
                       progressbar=True):
        """Read waterinfo file or directory.

        Parameters
        ----------
        file_or_dir : str
            path to file or directory. Files can be .csv or .zip
        name : str, optional
            name of the collection, by default ""
        ObsClass : Obs, optional
            type of Obs to read data as, by default obs.WaterlvlObs
        progressbar : bool, optional
            show progressbar, by default True

        Returns
        -------
        ObsCollection
            ObsCollection containing data

        """
        from .io import io_waterinfo

        meta = {
            'name': name,
            'type': ObsClass,
            'filename': file_or_dir
        }

        obs_list = io_waterinfo.read_waterinfo_obs(
            file_or_dir, ObsClass, progressbar=progressbar)
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)
    
    
    @classmethod
    def from_wiski(cls, dirname, ObsClass=obs.GroundwaterObs, suffix='.csv',
                   unpackdir=None, force_unpack=False, preserve_datetime=False,
                   verbose=False, keep_all_obs=True, **kwargs):

        from .io.io_wiski import read_wiski_dir

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
        obs_list = read_wiski_dir(dirname,
                                  ObsClass=ObsClass,
                                  suffix=suffix,
                                  unpackdir=unpackdir,
                                  force_unpack=force_unpack,
                                  preserve_datetime=preserve_datetime,
                                  verbose=verbose,
                                  keep_all_obs=keep_all_obs,
                                  **kwargs)
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

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
        from .io import io_fieldlogger

        if infer_otype:
            otype = self._infer_otype(verbose)
            if not isinstance(
                    otype, type):  # check of dit ook echt een type is
                raise TypeError(
                    'could not infer observation type, use infer_otype=False')

        f_df = io_fieldlogger.df_to_fieldlogger_csv(
            self,
            fname,
            otype,
            use_default_otype_settings,
            heading,
            name,
            subname,
            inputfield,
            properties,
            group,
            group_color,
            verbose)

        return f_df

    def to_pi_xml(self, fname, timezone="", version="1.24"):
        from .io import io_xml
        io_xml.write_pi_xml(self, fname, timezone=timezone, version=version)

    def to_pystore(self, store_name, pystore_path, groupby, item_name=None,
                   overwrite=False):
        """Write timeseries and metadata to Pystore format. Series are
        grouped by 'groupby'. Each group is a Collection, each series within
        that group an Item.

        Parameters
        ----------
        store_name : str
            name of the store
        pystore_path : str
            path where store should be saved
        groupby : str
            column name to groupby, (for example, group by location,
            which would create a collection per location, and write
            all timeseries at that location as Items into that Collection).
        item_name : str, optional
            name of column to use as item name, default is None, Item then
            takes name from obs.name
        overwrite : bool, optional
            if True overwrite current data in store, by default False
        """
        import pystore
        pystore.set_path(pystore_path)
        store = pystore.store(store_name)
        for name, group in self.groupby(by=groupby):
            # Access a collection (create it if not exist)
            collection = store.collection(name, overwrite=overwrite)
            for i, o in enumerate(group.obs):
                imeta = o.meta.copy()
                if 'datastore' in imeta.keys():
                    imeta['datastore'] = str(imeta['datastore'])
                # add extra columns to item metadata
                for icol in group.columns:
                    if icol != "obs" and icol != 'meta':
                        # check if type is numpy integer
                        # numpy integers are not json serializable
                        if isinstance(group.iloc[i].loc[icol], np.integer):
                            imeta[icol] = int(group.iloc[i].loc[icol])
                        else:
                            imeta[icol] = group.iloc[i].loc[icol]
                if item_name is None:
                    name = o.name
                else:
                    name = o.meta[item_name]
                collection.write(name, o, metadata=imeta,
                                 overwrite=overwrite)

    def to_arctic(self, connstr, libname, verbose=False):
        """Write ObsCollection to MongoDB using Arctic

        Parameters
        ----------
        connstr : str
            database connection string
        libname : str
            name of the library to store data
        verbose : bool, optional
            show progress bar, by default False

        """
        import arctic
        from tqdm import tqdm
        arc = arctic.Arctic(connstr)

        # get library
        try:
            lib = arc.get_library(libname)
        except arctic.exceptions.LibraryNotFoundException:
            print("Library '{}' not found, initializing library...".format(
                libname))
            arc.initialize_library(libname)
            lib = arc[libname]

        # write data
        for o in (tqdm(self.obs) if verbose else self.obs):
            metadata = o.meta
            lib.write(o.name, o, metadata=metadata)

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
        return util.df2gdf(self, xcol, ycol)

    def to_report_table(self, columns=['locatie', 'filternr',
                                       'Van', 'Tot', '# metingen']):

        if 'Van' in columns:
            self['Van'] = self.obs.apply(lambda x: x.index[0])
        if 'Tot' in columns:
            self['Tot'] = self.obs.apply(lambda x: x.index[-1])
        if '# metingen' in columns:
            self['# metingen'] = self.obs.apply(lambda x: x.shape[0])

        return self[columns]
    
    
    def to_pastastore(self, pstore=None, pstore_name='',
                      obs_column='stand_m_tov_nap',
                      kind='oseries', add_metadata=True,
                      verbose=False):
        """add observations to a new or existing pastastore

        Parameters
        ----------
        oc : observation.ObsCollection
            collection of observations
        pstore : pastastore.PastasProject, optional
            Existing pastastore, if None a new project is created
        pstore_name : str, optional
            Name of the pastastore only used if pstore is None
        obs_column : str, optional
            Name of the column in the Obs dataframe to be used
        kind : str, optional
            The kind of series that is added to the pastas project
        add_metadata : boolean, optional
            If True metadata from the observations added to the project.
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        project : pastastore.PastasProject
            the pastas project with the series from the ObsCollection
        """
        from .io.io_pastas import create_pastastore
        
        pstore = create_pastastore(self, pstore, pstore_name, 
                                   add_metadata=add_metadata,
                                   kind=kind,
                                   obs_column=obs_column)
    
        
        

        return pstore

    def to_pastas_project(self, pr=None, project_name='',
                          obs_column='stand_m_tov_nap',
                          kind='oseries', add_metadata=True,
                          verbose=False, **kwargs):
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
        add_metadata : boolean, optional
            If True metadata from the observations added to the project.
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs
            arguments are passed to the create_pastas_project funtion

        Returns
        -------
        pr : pastas.project
            the pastas project with the series from the ObsCollection
        """
        warnings.warn("this method will be removed in future versions, use to_pastastore instead", DeprecationWarning)

        from .io.io_pastas import create_pastas_project
        
        project = create_pastas_project(self, pr=pr, 
                                        project_name=project_name,
                                        obs_column=obs_column,
                                        kind=kind, add_metadata=add_metadata,
                                        verbose=verbose, **kwargs)
    
        
        

        return project
        

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
        gdf = util.df2gdf(self, xcol, ycol)

        # remove obs column
        if 'obs' in gdf.columns:
            gdf.drop(columns='obs', inplace=True)

        # cast boolean columns to int
        for colname, coltype in gdf.dtypes.items():
            if coltype == bool:
                gdf[colname] = gdf[colname].astype(int)
        gdf.to_file(fname)

    def add_meta_to_df(self, key):
        """Get the values from the meta dictionary of each observation object
        and add these to the ObsCollection as a column.


        to the ObsCollection

        Parameters
        ----------
        key : str
            key in meta dictionary of observation object


        """

        self[key] = [o.meta[key] for o in self.obs.values]

    def get_series(self, tmin=None, tmax=None, col="stand_m_tov_nap"):
        if tmin is None:
            tmin = self.dates_first_obs.min()
        if tmax is None:
            tmax = self.dates_last_obs.max()
        return self.obs.apply(lambda o: o.loc[tmin:tmax, col])
