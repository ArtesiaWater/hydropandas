"""Module with a number of observation classes.

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

from _io import StringIO
from pandas._config import get_option
from pandas.io.formats import console, format as fmt
from html import escape


class Obs(DataFrame):
    """class for point observations.

    An Obs object is a subclass of a pandas.DataFrame and allows for additional
    attributes and methods.
    More information about subclassing pandas DataFrames pandas can be
    found here:
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
    _metadata = ['name', 'x', 'y', 'meta', 'filename']

    def __init__(self, *args, **kwargs):
        """constructor of Obs class.

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

    def __repr__(self) -> str:
        """
        Return a string representation for a particular Observation.
        """
        buf = StringIO("")

        # write metadata properties
        buf.write('-----metadata------\n')
        for att in self._metadata:
            buf.write(f'{att} : {getattr(self, att)} \n')
        buf.write('\n')

        if self._info_repr():
            self.info(buf=buf)
            return buf.getvalue()

        max_rows = get_option("display.max_rows")
        min_rows = get_option("display.min_rows")
        max_cols = get_option("display.max_columns")
        max_colwidth = get_option("display.max_colwidth")
        show_dimensions = get_option("display.show_dimensions")
        if get_option("display.expand_frame_repr"):
            width, _ = console.get_console_size()
        else:
            width = None
        self.to_string(
            buf=buf,
            max_rows=max_rows,
            min_rows=min_rows,
            max_cols=max_cols,
            line_width=width,
            max_colwidth=max_colwidth,
            show_dimensions=show_dimensions,
        )

        return buf.getvalue()

    def _repr_html_(self, metadata_out=True):
        """
        Return a html representation for a particular Observation.

        Mainly for IPython notebook.
        """

        if self._info_repr():
            buf = StringIO("")
            self.info(buf=buf)
            # need to escape the <class>, should be the first line.
            val = buf.getvalue().replace("<", r"&lt;", 1)
            val = val.replace(">", r"&gt;", 1)
            return "<pre>" + val + "</pre>"

        if get_option("display.notebook_repr_html"):
            max_rows = get_option("display.max_rows")
            min_rows = get_option("display.min_rows")
            max_cols = get_option("display.max_columns")
            show_dimensions = get_option("display.show_dimensions")

            formatter = fmt.DataFrameFormatter(
                self,
                columns=None,
                col_space=None,
                na_rep="NaN",
                formatters=None,
                float_format=None,
                sparsify=None,
                justify=None,
                index_names=True,
                header=True,
                index=True,
                bold_rows=True,
                escape=True,
                max_rows=max_rows,
                min_rows=min_rows,
                max_cols=max_cols,
                show_dimensions=show_dimensions,
                decimal=".",
            )

            df_out = fmt.DataFrameRenderer(formatter).to_html(notebook=True)

            if metadata_out:
                obj_type = "hydropandas.{}".format(type(self).__name__)
                header = f"<div class='xr-header'><div class='xr-obj-type'>{escape(obj_type)}</div></div><hr>"
                meta_dic = {key: getattr(self, key) for key in self._metadata}
                meta_dic.pop('meta')
                meta_out = DataFrame(meta_dic, index=['metadata']).T
                formatter2 = fmt.DataFrameFormatter(
                    meta_out,
                    columns=None,
                    col_space=None,
                    na_rep="NaN",
                    formatters=None,
                    float_format=None,
                    sparsify=None,
                    justify=None,
                    index_names=True,
                    header=True,
                    index=True,
                    bold_rows=True,
                    escape=True,
                    max_rows=len(self._metadata) + 1,
                    min_rows=0,
                    max_cols=2,
                    show_dimensions=show_dimensions,
                    decimal=".",)
                df_meta_out = fmt.DataFrameRenderer(formatter2).to_html(notebook=True)
                return header + df_meta_out + df_out

            else:
                return df_out
        else:
            return None

    @property
    def _constructor(self):
        return Obs

    def to_collection_dict(self, include_meta=False):
        """get dictionary with registered attributes and their values of an Obs
        object.

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
    """Class for groundwater quantity observations.

    Subclass of the Obs class. Has the following attributes:

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
                  split_cluster=True,
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
        split_cluster : bool
            if False and the piezometer belongs to a cluster, the combined
            time series of the cluster is used. if True the indvidual time
            series of each piezometer is used. Default is True
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_csv if
            fname is not None and otherwise to io_dino.findMeetreeks
        """
        from .io import io_dino
        if fname is not None:
            # read dino csv file
            measurements, meta = io_dino.read_dino_groundwater_csv(
                fname, **kwargs)

            return cls(measurements, meta=meta, **meta)

        elif location is not None:
            measurements, meta = io_dino.download_dino_groundwater(location,
                                                                   filternr,
                                                                   tmin, tmax,
                                                                   split_cluster,
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

        from .io import io_dino
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

        from .io import io_dino
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
        from .io import io_dino
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
        """Create GroundwaterObs DataFrame from Pystore item.

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
    """class for groundwater quality (grondwatersamenstelling) point
    observations.

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
        """read ad dino file with groundwater quality data.

        Parameters
        ----------
        fname : str
            dino txt filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_groundwater_quality_txt
        """
        from .io import io_dino
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
        """read a dino file with waterlvl data.

        Parameters
        ----------
        fname : str
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_waterlvl_csv
        """
        from .io import io_dino
        measurements, meta = io_dino.read_dino_waterlvl_csv(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_waterinfo(cls, fname, **kwargs):
        """Read data from waterinfo csv-file or zip.

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
                                                        return_metadata=True,
                                                        **kwargs)
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
                  fill_missing_obs=True, interval='daily', inseason=False,
                  use_api=True, raise_exceptions=True):
        """Get a KnmiObs object.

        Parameters
        ----------
        stn : int or str
            measurement station e.g. 829.
        variable : str, optional
            observation type e.g. "RD" or "EV24". See list with all variables
            below.
        startdate : str, datetime or None, optional
            start date of observations. The default is None.
        enddate : str, datetime or None, optional
            end date of observations. The default is None.
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is True.
        interval : str, optional
            desired time interval for observations. The default is 'daily'.
        inseason : boolean, optional
            flag to obtain inseason data. The default is False
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary folder and the data is read from there.
            Default is True since the api is back online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        Returns
        -------
        KnmiObs object with time series and attributes

        List of possible variables:
            DDVEC     = Vectorgemiddelde windrichting in graden (360=noord, 90=oost, 180=zuid, 270=west, 0=windstil/variabel). Zie http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken / Vector mean wind direction in degrees (360=north, 90=east, 180=south, 270=west, 0=calm/variable)
            FHVEC     = Vectorgemiddelde windsnelheid (in 0.1 m/s). Zie http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken / Vector mean windspeed (in 0.1 m/s)
            FG        = Etmaalgemiddelde windsnelheid (in 0.1 m/s) / Daily mean windspeed (in 0.1 m/s)
            FHX       = Hoogste uurgemiddelde windsnelheid (in 0.1 m/s) / Maximum hourly mean windspeed (in 0.1 m/s)
            FHXH      = Uurvak waarin FHX is gemeten / Hourly division in which FHX was measured
            FHN       = Laagste uurgemiddelde windsnelheid (in 0.1 m/s) / Minimum hourly mean windspeed (in 0.1 m/s)
            FHNH      = Uurvak waarin FHN is gemeten / Hourly division in which FHN was measured
            FXX       = Hoogste windstoot (in 0.1 m/s) / Maximum wind gust (in 0.1 m/s)
            FXXH      = Uurvak waarin FXX is gemeten / Hourly division in which FXX was measured
            TG        = Etmaalgemiddelde temperatuur (in 0.1 graden Celsius) / Daily mean temperature in (0.1 degrees Celsius)
            TN        = Minimum temperatuur (in 0.1 graden Celsius) / Minimum temperature (in 0.1 degrees Celsius)
            TNH       = Uurvak waarin TN is gemeten / Hourly division in which TN was measured
            TX        = Maximum temperatuur (in 0.1 graden Celsius) / Maximum temperature (in 0.1 degrees Celsius)
            TXH       = Uurvak waarin TX is gemeten / Hourly division in which TX was measured
            T10N      = Minimum temperatuur op 10 cm hoogte (in 0.1 graden Celsius) / Minimum temperature at 10 cm above surface (in 0.1 degrees Celsius)
            T10NH     = 6-uurs tijdvak waarin T10N is gemeten / 6-hourly division in which T10N was measured; 6=0-6 UT, 12=6-12 UT, 18=12-18 UT, 24=18-24 UT
            SQ        = Zonneschijnduur (in 0.1 uur) berekend uit de globale straling (-1 voor <0.05 uur) / Sunshine duration (in 0.1 hour) calculated from global radiation (-1 for <0.05 hour)
            SP        = Percentage van de langst mogelijke zonneschijnduur / Percentage of maximum potential sunshine duration
            Q         = Globale straling (in J/cm2) / Global radiation (in J/cm2)
            DR        = Duur van de neerslag (in 0.1 uur) / Precipitation duration (in 0.1 hour)
            RH        = Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) / Daily precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
            RHX       = Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) / Maximum hourly precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
            RHXH      = Uurvak waarin RHX is gemeten / Hourly division in which RHX was measured
            PG        = Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1 hPa) berekend uit 24 uurwaarden / Daily mean sea level pressure (in 0.1 hPa) calculated from 24 hourly values
            PX        = Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa) / Maximum hourly sea level pressure (in 0.1 hPa)
            PXH       = Uurvak waarin PX is gemeten / Hourly division in which PX was measured
            PN        = Laagste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa) / Minimum hourly sea level pressure (in 0.1 hPa)
            PNH       = Uurvak waarin PN is gemeten / Hourly division in which PN was measured
            VVN       = Minimum opgetreden zicht / Minimum visibility; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)
            VVNH      = Uurvak waarin VVN is gemeten / Hourly division in which VVN was measured
            VVX       = Maximum opgetreden zicht / Maximum visibility; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)
            VVXH      = Uurvak waarin VVX is gemeten / Hourly division in which VVX was measured
            NG        = Etmaalgemiddelde bewolking (bedekkingsgraad van de bovenlucht in achtsten, 9=bovenlucht onzichtbaar) / Mean daily cloud cover (in octants, 9=sky invisible)
            UG        = Etmaalgemiddelde relatieve vochtigheid (in procenten) / Daily mean relative atmospheric humidity (in percents)
            UX        = Maximale relatieve vochtigheid (in procenten) / Maximum relative atmospheric humidity (in percents)
            UXH       = Uurvak waarin UX is gemeten / Hourly division in which UX was measured
            UN        = Minimale relatieve vochtigheid (in procenten) / Minimum relative atmospheric humidity (in percents)
            UNH       = Uurvak waarin UN is gemeten / Hourly division in which UN was measured
            EV24      = Referentiegewasverdamping (Makkink) (in 0.1 mm) / Potential evapotranspiration (Makkink) (in 0.1 mm)
        """
        from .io import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_stn(
            stn, variable, startdate, enddate, fill_missing_obs,
            interval=interval, inseason=inseason,
            use_api=use_api,
            raise_exceptions=raise_exceptions)

        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])

    @classmethod
    def from_nearest_xy(cls, x, y, variable, startdate=None, enddate=None,
                        fill_missing_obs=True, interval='daily',
                        inseason=False, raise_exceptions=False):
        """Get KnmiObs object with measurements from station closest to the
        given (x,y) coördinates.

        Parameters
        ----------
        x : int or float
            x coördinate in m RD.
        y : int or float
            y coördinate in m RD.
        variable : str
            e.g. 'EV24'.
        startdate : str, datetime or None, optional
            start date of observations. The default is None.
        enddate : str, datetime or None, optional
            end date of observations. The default is None.
        fill_missing_obs : bool
            if True missing observations are filled with values of next closest
            KNMI station
        interval : str, optional
            desired time interval for observations. The default is 'daily'.
        inseason : boolean, optional
            flag to obtain inseason data. The default is False
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        Returns
        -------
        KnmiObs object with time series and attributes
        """
        from .io import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_xy(
            x, y, variable, startdate, enddate, fill_missing_obs,
            interval=interval, inseason=inseason,
            raise_exceptions=raise_exceptions)

        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])

    @classmethod
    def from_obs(cls, obs, variable, startdate=None, enddate=None,
                 fill_missing_obs=True, interval='daily', inseason=False,
                 raise_exceptions=False):

        from .io import io_knmi

        x = obs.x
        y = obs.y

        if startdate is None:
            startdate = obs.index[0]
        if enddate is None:
            enddate = obs.index[-1]

        ts, meta = io_knmi.get_knmi_timeseries_xy(
            x, y, variable, startdate, enddate, fill_missing_obs,
            interval=interval, inseason=inseason,
            raise_exceptions=raise_exceptions)

        return cls(ts, meta=meta, station=meta['station'], x=meta['x'],
                   y=meta['y'], name=meta['name'])
