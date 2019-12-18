'''
module with ObsCollection class for a collection of observations.

The ObsCollection class is a subclass of a pandas DataFrame with
additional attributes and methods.

More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

'''
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
            
    def _set_metadata_value(self, iname, att_name, value, add_to_meta=True,
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
            of an observation. The default is True.
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
    def from_dino_server(cls, extent=None, bbox=None,
                         ObsClass=obs.GroundwaterObs,
                         name=None, get_metadata=True,
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
        get_metadata : boolean, optional

        verbose : boolean, optional
            Print additional information to the screen (default is False).
        kwargs:
            kwargs are passed to the io_dino.download_dino_within_extent() function


        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """

        from .io import io_dino

        if ObsClass == obs.GroundwaterObs:
            layer = 'grondwatermonitoring'
        else:
            raise NotImplementedError('cannot download {} from Dino'.format(ObsClass))

        obs_df = io_dino.download_dino_within_extent(extent, bbox,
                                                     ObsClass, layer=layer,
                                                     get_metadata=get_metadata,
                                                     verbose=verbose, **kwargs)

        if bbox is None:
            bbox = [extent[0], extent[2], extent[1], extent[3]]

        if name is None:
            name = '{} from DINO'.format(layer)

        meta = kwargs
        meta.update({'verbose': verbose})

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

        from .io import io_dino

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

        obs_df = io_dino.read_dino_dir(
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

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews(cls, file_or_dir, ObsClass=obs.GroundwaterObs, name='fews',
                  to_mnap=True, remove_nan=True,
                  unpackdir=None, force_unpack=False,
                  preserve_datetime=False, verbose=False):
        """ read one or several XML-files with measurements from FEWS

        Parameters
        ----------
        file_or_dir :  str
            zip, xml or directory with zips or xml files to read
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
        from .io import io_xml

        # get files
        dirname, unzip_fnames = util.get_files(file_or_dir, ext=".xml",
                                               force_unpack=force_unpack,
                                               preserve_datetime=preserve_datetime)
        meta = {'filename': dirname,
                'type': ObsClass,
                'verbose': verbose
                }

        obs_list = []
        nfiles = len(unzip_fnames)
        for j, ixml in enumerate(unzip_fnames):
            if verbose:
                print("{0}/{1} read {2}".format(j + 1, nfiles, ixml))
            fullpath = os.path.join(dirname, ixml)
            olist = io_xml.read_xml(
                fullpath,
                ObsClass=ObsClass,
                to_mnap=to_mnap,
                remove_nan=remove_nan,
                verbose=False)
            obs_list += olist

        obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                              columns=obs_list[0].to_collection_dict().keys())

        obs_df.set_index('name', inplace=True)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews2(cls, file_or_dir, ObsClass=obs.GroundwaterObs, name='fews',
                   locations=None, to_mnap=True, remove_nan=True,
                   unpackdir=None, force_unpack=False,
                   preserve_datetime=False,
                   verbose=False):
        """ read a XML-file with measurements from FEWS

        Parameters
        ----------
        file_or_dir : str
            zip, xml file or directory to read
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
        from .io import io_xml

        # get files
        dirname, unzip_fnames = util.get_files(file_or_dir, ext=".xml",
                                               unpackdir=unpackdir,
                                               force_unpack=force_unpack,
                                               preserve_datetime=preserve_datetime)
        meta = {'filename': dirname,
                'type': ObsClass,
                'verbose': verbose
                }

        obs_list = []
        nfiles = len(unzip_fnames)
        for j, ixml in enumerate(unzip_fnames):
            if verbose:
                print("{0}/{1} read {2}".format(j + 1, nfiles, ixml))
            fullpath = os.path.join(dirname, ixml)
            _, olist = io_xml.iterparse_pi_xml(fullpath, ObsClass,
                                               locationIds=locations,
                                               verbose=verbose)
            obs_list += olist

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

        from .io import io_fieldlogger

        obs_list, fieldlogger_meta = io_fieldlogger.fieldlogger_csv_to_obs_list(
            fname, ObsClass=obs.GroundwaterObs)

        return cls(cls.from_list(obs_list, name=name), meta=fieldlogger_meta)

    @classmethod
    def from_menyanthes(cls, fname, name='', ObsClass=obs.GroundwaterObs,
                        verbose=False):

        from .io import io_menyanthes

        obs_list = io_menyanthes.read_file(fname, ObsClass, verbose)

        menyanthes_meta = {'filename': fname,
                           'type': ObsClass,
                           'verbose': verbose}

        return cls(cls.from_list(obs_list, name=name), meta=menyanthes_meta)

    @classmethod
    def from_imod(cls, obs_collection, ml, runfile, mtime, model_ws,
                  modelname='', nlay=None, exclude_layers=0, verbose=False):
        """read 'observations' from an imod model

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

        import imod

        if ml.modelgrid.xoffset == 0 or ml.modelgrid.yoffset == 0:
            warnings.warn(
                'you probably want to set the xll and/or yll attributes of ml.modelgrid')

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
                if verbose:
                    print('read {}'.format(fname))
                ihds, _attrs = imod.idf.read(fname)
                hm = util.interpolate(ihds, vtx, wts)
                hm_ts[mask, t] = hm[mask]

        mo_list = []
        for i, name in enumerate(obs_collection.index):
            mo = obs.ModelObs(index=mtime,
                              data=hm_ts[i],
                              name=name,
                              model=modelname,
                              x=obs_collection.loc[name,
                                                   'x'],
                              y=obs_collection.loc[name,
                                                   'y'],
                              meta=obs_collection.loc[name,
                                                      'obs'].meta)
            mo_list.append(mo)

        mobs_df = pd.DataFrame([mo.to_collection_dict() for mo in mo_list],
                               columns=mo_list[0].to_collection_dict().keys())
        mobs_df.set_index('name', inplace=True)
        if verbose:
            print(mobs_df)

        return cls(mobs_df, name=modelname)

    @classmethod
    def from_modflow(cls, obs_collection, ml, hds_arr, mtime,
                     modelname='', nlay=None, exclude_layers=0, verbose=False):
        """ Read moflow groundwater heads at the location of the points in
        obs_collection.

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

        if ml.modelgrid.xoffset == 0 or ml.modelgrid.yoffset == 0:
            warnings.warn(
                'you probably want to set the xll and/or yll attributes of dis.sr')

        if nlay is None:
            nlay = ml.modelgrid.nlay

        if modelname == '':
            modelname = ml.name

        xmid, ymid, _ = ml.modelgrid.xyzcellcenters

        xy = np.array([xmid.ravel(), ymid.ravel()]).T
        uv = obs_collection.loc[:, ("x", "y")].dropna(how="any", axis=0).values
        vtx, wts = util.interp_weights(xy, uv)

        # get interpolated timeseries from hds_arr
        hm_ts = np.zeros((obs_collection.shape[0], hds_arr.shape[0]))

        # loop over layers
        for m in range(nlay):
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
            mo = obs.ModelObs(index=mtime,
                              data=hm_ts[i],
                              name=name,
                              model=modelname,
                              x=obs_collection.loc[name,
                                                   'x'],
                              y=obs_collection.loc[name,
                                                   'y'],
                              meta=obs_collection.loc[name,
                                                      'obs'].meta)
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
                            name=None, pr_meta=None, rename_dic={}):

        if name is None:
            name = pr.name

        if pr_meta is None:
            pr_meta = pr.file_info

        obs_list = []
        for index, row in pr.oseries.iterrows():
            metadata = row.to_dict()
            for key in rename_dic.keys():
                if key in metadata.keys():
                    metadata[rename_dic[key]] = metadata.pop(key)

            s = pd.DataFrame(metadata.pop('series').series_original)
            s.rename(columns={index: 'Stand_m_tov_NAP'}, inplace=True)

            keys_o = ['name', 'x', 'y', 'locatie', 'filternr',
                      'metadata_available', 'maaiveld', 'meetpunt',
                      'bovenkant_filter', 'onderkant_filter']
            meta_o = {k: metadata[k] for k in keys_o if k in metadata}

            o = ObsClass(s, meta=metadata, **meta_o)
            obs_list.append(o)

        return cls(cls.from_list(obs_list, name=name), meta=pr_meta)

    @classmethod
    def from_wiski(cls, dirname, ObsClass=obs.GroundwaterObs, suffix='.csv',
                   unpackdir=None, force_unpack=False, preserve_datetime=False,
                   verbose=False, keep_all_obs=True, **kwargs):

        from .io import io_wiski

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

        from .io import io_pystore
        io_pystore.set_pystore_path(pystore_path)
        if not os.path.isdir(os.path.join(pystore_path, storename)):
            raise FileNotFoundError(
                "pystore -> '{}' does not exist".format(storename))
        # obtain item names within extent
        if extent is not None:
            meta_list = io_pystore.read_store_metadata(
                storename, items="first")
            obs_df = pd.DataFrame(meta_list)
            obs_df.set_index('item_name', inplace=True)
            obs_df['x'] = pd.to_numeric(obs_df.x, errors='coerce')
            obs_df['y'] = pd.to_numeric(obs_df.y, errors='coerce')
            item_names = obs_df[(obs_df.x > extent[0]) & (obs_df.x < extent[1]) & (
                obs_df.y > extent[2]) & (obs_df.y < extent[3])].index

        if read_series:
            obs_list = io_pystore.store_to_obslist(
                storename,
                ObsClass=ObsClass,
                collection_names=collection_names,
                item_names=item_names,
                nameby=nameby,
                verbose=verbose,
                progressbar=progressbar)
            columns = []
            coldict = []
            for o in obs_list:
                d = o.to_collection_dict()
                coldict.append(d)
                columns |= d.keys()
            obs_df = pd.DataFrame(coldict, columns=columns)
            obs_df.set_index('name', inplace=True)
            meta = {'fname': obs_list[0].meta["datastore"],
                    'type': obs.GroundwaterObs,
                    'verbose': True}
        else:
            if item_names is None:
                item_names = "all"
            meta_list = io_pystore.read_store_metadata(storename,
                                                       items=item_names,
                                                       verbose=verbose)
            meta = {}
            obs_df = pd.DataFrame(meta_list)
            if nameby == "collection":
                obs_df.set_index('collection_name', inplace=True)
            elif nameby == "item":
                obs_df.set_index('item_name', inplace=True)
            elif nameby == "both":
                obs_df["name"] = obs_df["collection_name"] + "__" + \
                    obs_df["item_name"]
            else:
                raise ValueError("'{}' is not a valid option "
                                 "for 'nameby'".format(nameby))

        return cls(obs_df, name=storename, meta=meta)

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

        from tqdm import tqdm
        import arctic
        from timeit import default_timer

        arc = arctic.Arctic(connstr)
        lib = arc.get_library(libname)

        start = default_timer()
        obs_list = []
        rows_read = 0
        for sym in (tqdm(lib.list_symbols())
                    if verbose else lib.list_symbols()):
            item = lib.read(sym)
            o = ObsClass(item.data, name=sym, meta=item.metadata)
            obs_list.append(o)
            if verbose:
                rows_read += len(item.data.index)

        end = default_timer()
        if verbose:
            print("Symbols: {0:.0f}  "
                  "Rows: {1:,.0f}  "
                  "Time: {2:.2f}s  "
                  "Rows/s: {3:,.1f}".format(len(lib.list_symbols()),
                                            rows_read,
                                            (end - start),
                                            rows_read / (end - start)))
        # create dataframe
        columns = []
        coldict = []
        for o in obs_list:
            d = o.to_collection_dict()
            coldict.append(d)
            columns |= d.keys()
        obs_df = pd.DataFrame(coldict, columns=columns)
        obs_df.set_index('name', inplace=True)
        meta = {'type': obs.GroundwaterObs}

        return cls(obs_df, meta=meta)

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
            collection = store.collection(name)
            for i, o in enumerate(group.obs):
                imeta = o.meta.copy()
                if 'datastore' in imeta.keys():
                    imeta['datastore'] = str(imeta['datastore'])
                # add extra columns to item metadata
                for icol in group.columns:
                    if icol not in imeta.keys() and icol != "obs":
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
        art = util._import_art_tools()
        return art.util.df2gdf(self, xcol, ycol)

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
                print('add to pastas project -> {}'.format(o.name))
            # clean up metadata so there is no dict in dict
            meta = dict()
            for k, v in o.meta.items():
                if isinstance(v, dict):
                    meta.update(v)
                else:
                    meta[k] = v
            series = ps.TimeSeries(o[obs_column], name=o.name, metadata=meta)
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

        # remove obs column
        if 'obs' in gdf.columns:
            gdf.drop(columns='obs', inplace=True)

        # cast boolean columns to int
        for colname, coltype in gdf.dtypes.items():
            if coltype == bool:
                gdf[colname] = gdf[colname].astype(int)
        gdf.to_file(fname)

    def add_meta_to_df(self, name):
        """Add data from the metadata dictionary as a column

        to the ObsCollection

        Parameters
        ----------
        name : str
            variable name in metadata dictionary


        """

        self[name] = [o.meta[name] for o in self.obs.values]

    def get_series(self, tmin=None, tmax=None, col="Stand_m_tov_NAP"):
        if tmin is None:
            tmin = self.dates_first_obs.min()
        if tmax is None:
            tmax = self.dates_last_obs.max()
        return self.obs.apply(lambda o: o.loc[tmin:tmax, col])
