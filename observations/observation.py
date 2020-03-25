"""
Module with a number of observation classes.

The Obs class is a subclass of a pandas DataFrame with
additional attributes and methods. The specific classes (GroundwaterObs,
WaterlvlObs, ...) are subclasses of the Obs class.

The subclasses of a dataframe can have additional attributes and methods.
Additional attributes have to be defined in the '_metadata' attribute. In order
to keep the subclass methods and attributes when selecting or slicing an object
you need the '_constructor' method.


More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

"""

import warnings
import numpy as np
from pandas import DataFrame

from .io import io_dino


class Obs(DataFrame):
    """class for point observations.

    An Obs object is a subclass of a pandas.DataFrame and allows for additional
    attributes and methods.
    pandas can be found here:
    http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

    Parameters
    ----------
    x : int or float
        x coordinate of observation point
    y : int or float
        y coordinate of observation point
    name : str
        name
    meta : dictionary
        metadata
    filename : str
        filename with data of observation point

    """
    # temporary properties
    _internal_names = DataFrame._internal_names + ['none']
    _internal_names_set = set(_internal_names)

    # normal properties
    _metadata = ['x', 'y', 'name', 'meta', 'filename']

    def __init__(self, *args, **kwargs):
        """ constructor of Obs class

        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.
        """
        self.x = kwargs.pop('x', np.nan)
        self.y = kwargs.pop('y', np.nan)
        self.name = kwargs.pop('name', '')
        self.meta = kwargs.pop('meta', {})
        self.filename = kwargs.pop('filename', '')

        super(Obs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return Obs

    def to_collection_dict(self, include_meta=False):
        """get dictionary with registered attributes and their values
        of an Obs object.

        This method can be used to create a dataframe from a collection
        of Obs objects.
        
        Parameters
        ----------
        include_meta : boolean, optional
            include the meta dictionary in the collection dictionary,
            default is false

        Returns
        -------
        d : dictionary
            dictionary with Obs information
        """
        
        attrs = self._metadata.copy()
        if not include_meta:
            attrs.remove('meta')
        
        d = {}
        for att in attrs:
            d[att] = getattr(self, att)

        d['obs'] = self

        return d


class GroundwaterObs(Obs):
    """class for groundwater quantity point observations

    Subclass of the Obs class. Can have the following attributes:
        - locatie: 2 filters at one piezometer should have the same 'locatie'
        - filternr: 2 filters at one piezometer should have a different 'filternr'.
        a higher filter number is preferably deeper than a lower filter number.
        - bovenkant_filter: top op the filter in m NAP
        - onderkant_filter: bottom of the filter in m NAP
        - maaiveld: surface level in m NAP
        - meetpunt: ? in m NAP
        - metadata_available: boolean indicating if metadata is available for
        the measurement point.

    """

    _metadata = Obs._metadata + \
        ['locatie', 'filternr',
         'bovenkant_filter', 'onderkant_filter',
         'maaiveld', 'meetpunt', 'metadata_available'
         ]

    def __init__(self, *args, **kwargs):
        """
        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.

        if the pandas.DataFrame has a column 'stand_m_tov_nap' a lot of
        plotting and other methods will work automatically without changing
        the default arguments.
        """
        self.locatie = kwargs.pop('locatie', '')
        self.filternr = kwargs.pop('filternr', '')
        self.maaiveld = kwargs.pop('maaiveld', np.nan)
        self.meetpunt = kwargs.pop('meetpunt', np.nan)
        self.bovenkant_filter = kwargs.pop('bovenkant_filter', np.nan)
        self.onderkant_filter = kwargs.pop('onderkant_filter', np.nan)
        self.metadata_available = kwargs.pop('metadata_available', np.nan)

        super(GroundwaterObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return GroundwaterObs

    @classmethod
    def from_dino(cls, fname=None, location=None, filternr=1.,
                  tmin="1900-01-01", tmax="2040-01-01",
                  **kwargs):
        """download dino data from the server.

        Parameters
        ----------
        fname : str, optional
            dino csv filename
        location : str, optional
            location of the peilbuis, i.e. B57F0077
        filternr : float, optional
            filter_nr of the peilbuis, i.e. 1.
        tmin : str
            start date in format YYYY-MM-DD
        tmax : str
            end date in format YYYY-MM-DD
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_csv if
            fname is not None and otherwise to io_dino.findMeetreeks


        """
        if fname is not None:
            # read dino csv file
            measurements, meta = io_dino.read_dino_groundwater_csv(
                fname, **kwargs)

            return cls(measurements, meta=meta, **meta)

        elif location is not None:
            measurements, meta = io_dino.download_dino_groundwater(location,
                                                                   filternr,
                                                                   tmin, tmax,
                                                                   **kwargs)

            if meta['metadata_available']:
                return cls(measurements, meta=meta,
                           x=meta.pop('x'), y=meta.pop('y'),
                           onderkant_filter=meta.pop('onderkant_filter'),
                           bovenkant_filter=meta.pop('bovenkant_filter'),
                           name=meta.pop('name'),
                           metadata_available=meta.pop('metadata_available'),
                           locatie=meta.pop('locatie'),
                           maaiveld=meta.pop('maaiveld'),
                           filternr=meta.pop('filternr'))
            else:
                return cls(measurements, meta=meta)
        else:
            raise ValueError(
                'specify fname or location to obtain groundwater heads')

    @classmethod
    def from_dino_server(cls, location, filternr=1.,
                         tmin="1900-01-01", tmax="2040-01-01",
                         **kwargs):
        """download dino data from the server.

        Parameters
        ----------
        location : str
            location of the peilbuis, i.e. B57F0077
        filternr : float
            filter_nr of the peilbuis, i.e. 1.0
        tmin : str
            start date in format YYYY-MM-DD
        tmax : str
            end date in format YYYY-MM-DD
        kwargs : key-word arguments
            these arguments are passed to dino.findMeetreeks functie
        """

        warnings.warn(
            "this method will be removed in future versions, use from_dino instead", DeprecationWarning)

        measurements, meta = io_dino.download_dino_groundwater(location,
                                                               filternr,
                                                               tmin, tmax,
                                                               **kwargs)

        if meta['metadata_available']:
            return cls(measurements, meta=meta,
                       x=meta.pop('x'), y=meta.pop('y'),
                       onderkant_filter=meta.pop('onderkant_filter'),
                       bovenkant_filter=meta.pop('bovenkant_filter'),
                       name=meta.pop('name'),
                       locatie=meta.pop('locatie'),
                       maaiveld=meta.pop('maaiveld'),
                       filternr=meta.pop('filternr'))
        else:
            return cls(measurements, meta=meta)

    @classmethod
    def from_dino_file(cls, fname=None, **kwargs):
        """read a dino csv file.

        Parameters
        ----------
        name : str, optional
            name of the peilbuis, i.e. B57F0077
        fname : str, optional
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_csv
        """

        warnings.warn(
            "this method will be removed in future versions, use from_dino instead", DeprecationWarning)

        if fname is not None:
            # read dino csv file
            measurements, meta = io_dino.read_dino_groundwater_csv(
                fname, **kwargs)

            return cls(measurements, meta=meta, **meta)
        else:
            raise ValueError(
                'specify either the name or the filename of the measurement point')

    @classmethod
    def from_artdino_file(cls, fname=None, **kwargs):
        """read a dino csv file (artdiver style).

        Parameters
        ----------
        name : str, optional
            name of the peilbuis, i.e. B57F0077
        fname : str, optional
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_csv

        """

        if fname is not None:
            # read dino csv file
            measurements, meta = io_dino.read_artdino_groundwater_csv(
                fname, **kwargs)

            return cls(measurements, meta=meta, **meta)
        else:
            raise ValueError('specify either the name or the filename of the '
                             'measurement point!')

    @classmethod
    def from_wiski(cls, fname, **kwargs):
        """[summary]

        Parameters
        ----------
        fname : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        from .io import io_wiski
        data, metadata = io_wiski.read_wiski_file(fname, **kwargs)

        return cls(data, meta=metadata, **metadata)

    @classmethod
    def from_pystore_item(cls, item):
        """Create GroundwaterObs DataFrame from Pystore item

        Parameters
        ----------
        item : pystore.item.Item
            Pystore item

        Returns
        -------
        GroundwaterObs
            GroundwaterObs DataFrame

        """

        df = item.to_pandas()
        try:
            x = item.metadata["x"]
            y = item.metadata["y"]
        except KeyError:
            x = np.nan
            y = np.nan
        item.metadata["datastore"] = item.datastore
        return cls(df, x=x, y=y, meta=item.metadata)


class GroundwaterQualityObs(Obs):
    """class for groundwater quality (grondwatersamenstelling)
    point observations.

    Subclass of the Obs class

    """

    _metadata = Obs._metadata + \
        ['locatie', 'filternr', 'maaiveld', 'metadata_available']

    def __init__(self, *args, **kwargs):

        self.locatie = kwargs.pop('locatie', '')
        self.filternr = kwargs.pop('filternr', '')
        self.maaiveld = kwargs.pop('maaiveld', np.nan)
        self.metadata_available = kwargs.pop('metadata_available', np.nan)

        super(GroundwaterQualityObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return GroundwaterQualityObs

    @classmethod
    def from_dino(cls, fname, **kwargs):
        """read ad dino file with groundwater quality data

        Parameters
        ----------
        fname : str
            dino txt filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_quality_txt
        """

        measurements, meta = io_dino.read_dino_groundwater_quality_txt(
            fname, **kwargs)

        return cls(measurements, meta=meta, **meta)


class WaterlvlObs(Obs):
    """class for water level point observations.

    Subclass of the Obs class

    """

    _metadata = Obs._metadata + ['locatie', 'metadata_available']

    def __init__(self, *args, **kwargs):

        self.locatie = kwargs.pop('locatie', '')
        self.metadata_available = kwargs.pop('metadata_available', np.nan)

        super(WaterlvlObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return WaterlvlObs

    @classmethod
    def from_dino(cls, fname, **kwargs):
        """read a dino file with waterlvl data

        Parameters
        ----------
        fname : str
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_waterlvl_csv
        """

        measurements, meta = io_dino.read_dino_waterlvl_csv(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_waterinfo(cls, fname, **kwargs):
        """
        Read data from waterinfo csv-file or zip.

        Parameters
        ----------
        fname : str
            path to file (file can zip or csv)

        Returns
        -------
        df : WaterlvlObs
            WaterlvlObs object

        Raises
        ------
        ValueError
            if file contains data for more than one location

        """
        from .io import io_waterinfo
        df, metadata = io_waterinfo.read_waterinfo_file(fname,
                                                        return_metadata=True)
        return cls(df, meta=metadata, **metadata)


class ModelObs(Obs):
    """class for model point results.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + ['model']

    def __init__(self, *args, **kwargs):

        self.model = kwargs.pop('model', '')

        super(ModelObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return ModelObs


class KnmiObs(Obs):
    """class for KNMI timeseries.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + ['station']

    def __init__(self, *args, **kwargs):

        self.station = kwargs.pop('station', np.nan)

        super(KnmiObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return KnmiObs

    @classmethod
    def from_knmi(cls, stn, variable, startdate=None, enddate=None,
                  fill_missing_obs=True, verbose=False):
        from .io import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_stn(stn, variable,
                                                   startdate, enddate,
                                                   fill_missing_obs,
                                                   verbose=verbose)
        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])

    @classmethod
    def from_nearest_xy(cls, x, y, variable, startdate=None, enddate=None,
                        fill_missing_obs=True, verbose=False):
        from .io import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_xy(x, y, variable,
                                                  startdate, enddate,
                                                  fill_missing_obs,
                                                  verbose=verbose)

        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])

    @classmethod
    def from_obs(cls, obs, variable, startdate=None, enddate=None,
                 fill_missing_obs=True, verbose=False):

        from .io import io_knmi

        x = obs.x
        y = obs.y

        if startdate is None:
            startdate = obs.index[0]
        if enddate is None:
            enddate = obs.index[-1]

        ts, meta = io_knmi.get_knmi_timeseries_xy(x, y, variable,
                                                  startdate, enddate,
                                                  fill_missing_obs,
                                                  verbose=verbose)

        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])
