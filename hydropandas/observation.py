"""Module with observation classes.

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

import logging

import numpy as np
import pandas as pd
from _io import StringIO
from pandas._config import get_option
from pandas.api.types import is_numeric_dtype
from pandas.io.formats import console

logger = logging.getLogger(__name__)


class Obs(pd.DataFrame):
    """generic class for a time series with measurements at a certain location.

    Unless specified explicitly the first numeric column in the observation is
    used for analysis and plotting.

    Parameters
    ----------
    name : str
        name
    x : int or float
        x coordinate of observation point
    y : int or float
        y coordinate of observation point
    meta : dictionary
        metadata
    filename : str
        filename with data of observation point
    source : str
        source of the observation e.g. BRO or KNMI
    unit : str
        unit of the first numerical column in the observation
    """

    # temporary properties
    _internal_names = pd.DataFrame._internal_names + ["none"]
    _internal_names_set = set(_internal_names)

    # normal properties
    _metadata = ["name", "x", "y", "meta", "filename", "source", "unit"]

    def __init__(self, *args, **kwargs):
        """constructor of Obs class.

        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.
        """
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in Obs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        self.x = kwargs.pop("x", np.nan)
        self.y = kwargs.pop("y", np.nan)
        self.name = kwargs.pop("name", "")
        self.meta = kwargs.pop("meta", {})
        self.filename = kwargs.pop("filename", "")
        self.source = kwargs.pop("source", "")
        self.unit = kwargs.pop("unit", "")

        super(Obs, self).__init__(*args, **kwargs)

    def __repr__(self) -> str:
        """Return a string representation for a particular Observation."""
        buf = StringIO("")

        buf.write(f"{type(self).__name__} {self.name}\n")

        # write metadata properties
        buf.write("-----metadata------\n")
        for att in self._metadata:
            if not att == "meta":
                buf.write(f"{att} : {getattr(self, att)} \n")
        buf.write("\n")

        if self._info_repr():
            self.info(buf=buf)
            return buf.getvalue()

        buf.write("-----time series------\n")
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

    @property
    def _constructor(self):
        return Obs

    def _get_first_numeric_col_name(self):
        """get the first numeric column name of the observations

        Returns
        -------
        col : str or int
            column name.

        """
        if self.empty:
            return None

        for col in self.columns:
            if is_numeric_dtype(self[col]):
                break

        return col

    def copy(self, deep=True):
        """create a copy of the observation.

        When ``deep=True`` (default), a new object will be created with a
        copy of the calling object's data and indices. Modifications to
        the data or indices of the copy will not be reflected in the
        original object (see notes below).

        When ``deep=False``, a new object will be created without copying
        the calling object's data or index (only references to the data
        and index are copied). Any changes to the data of the original
        will be reflected in the shallow copy (and vice versa).

        Parameters
        ----------
        deep : bool, default True
            Make a deep copy, including a copy of the data and the indices.
            With ``deep=False`` neither the indices nor the data are copied.

        Returns
        -------
        o : TYPE
            DESCRIPTION.
        """

        if not deep:
            return super().copy(deep=deep)

        o = super().copy(deep=deep)

        # creat a copy of all mutable attributes
        for att in self._metadata:
            val = getattr(self, att)
            if isinstance(val, (list, dict)):
                setattr(o, att, val.copy())

        return o

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
            attrs.remove("meta")

        d = {}
        for att in attrs:
            d[att] = getattr(self, att)

        d["obs"] = self

        return d

    def merge_metadata(self, right, overlap="error"):
        """Merge the metadata of an Obs object with new metadata.

        Parameters
        ----------
        right : dict
            dictionary with the metadata of an Obs object.
        overlap : str, optional
            How to deal with overlapping metadata with different values.
            Options are:
                error : Raise a ValueError
                use_left : Use the metadata from self
                use_right : Use the given metadata
            Default is 'error'.

        Raises
        ------
        ValueError
            if the metadata differs and overlap='error'.

        Returns
        -------
        new_metadata : dict
            metadata after merge.
        """
        new_metadata = {}
        same_metadata = True
        for key in self._metadata:
            v1 = getattr(self, key)
            v2 = right[key]
            try:
                if v1 != v2:
                    same_metadata = False
                    if overlap == "error":
                        raise ValueError(
                            f"existing observation {key} differs from new observation"
                        )
                    elif overlap == "use_left":
                        logger.info(
                            f"existing observation {key} differs from new "
                            "observation, use existing"
                        )
                        new_metadata[key] = v1
                    elif overlap == "use_right":
                        logger.info(
                            f"existing observation {key} differs from new "
                            "observation, use new"
                        )
                        new_metadata[key] = v2
                else:
                    new_metadata[key] = v1
            except TypeError:
                same_metadata = False
                if overlap == "error":
                    raise ValueError(
                        f"existing observation {key} differs from new observation"
                    )
                elif overlap == "use_left":
                    logger.info(
                        f"existing observation {key} differs from new "
                        "observation, use existing"
                    )
                    new_metadata[key] = v1
                elif overlap == "use_right":
                    logger.info(
                        f"existing observation {key} differs from new "
                        "observation, use new"
                    )
                    new_metadata[key] = v2
        if same_metadata:
            logger.info("new and existing observation have the same metadata")

        return new_metadata

    def _merge_timeseries(self, right, overlap="error"):
        """merge two timeseries.

        Parameters
        ----------
        right : hpd.observation.Obs
            Observation object.
        overlap : str, optional
            How to deal with overlapping timeseries with different values.
            Options are:
                error : Raise a ValueError
                use_left : use the part of the overlapping timeseries from the
                existing observation
                use_right : use the part of the overlapping timeseries from the
                new observation
            Default is 'error'.

        Raises
        ------
        ValueError
            when the time series have different values on the same timestamp.

        Returns
        -------
        Observation object.
        """
        o = self.iloc[0:0, 0:0]

        # check if time series are the same
        if self.equals(right):
            logger.info("new and existing observation have the same time series")
            return self

        logger.info("new observation has a different time series")
        logger.info("merge time series")

        # merge overlapping columns
        overlap_cols = list(set(self.columns) & set(right.columns))

        # get timeseries for overlapping columns and indices
        dup_ind_o = right.loc[right.index.isin(self.index)]
        if not dup_ind_o.empty:
            # check if overlapping timeseries have different values
            if not self.loc[dup_ind_o.index, overlap_cols].equals(
                dup_ind_o[overlap_cols]
            ):
                logger.warning(
                    f"timeseries of observation {right.name} overlap with "
                    "different values"
                )
                if overlap == "error":
                    raise ValueError(
                        "observations have different values for same time steps"
                    )
                elif overlap == "use_left":
                    dup_o = self.loc[dup_ind_o.index, overlap_cols]
                elif overlap == "use_right":
                    dup_o = dup_ind_o[overlap_cols]
            else:
                dup_o = self.loc[dup_ind_o.index, overlap_cols]
        else:
            dup_o = dup_ind_o[overlap_cols]

        # get unique observations from overlapping columns
        unique_o_left = self.loc[
            ~self.index.isin(right.index)
        ]  # get unique observations left
        unique_o_right = right.loc[
            ~right.index.isin(self.index)
        ]  # get unique observations right

        # merge unique observations from overlapping columns with duplicate observations
        # from overlapping columns
        dfcol = pd.concat(
            [dup_o, unique_o_left[overlap_cols], unique_o_right[overlap_cols]]
        )
        o = pd.concat([o, dfcol], axis=1)

        # merge non overlapping columns
        for col in right.columns:
            if col not in o.columns:
                o = pd.concat([o, right[[col]]], axis=1)
        for col in self.columns:
            if col not in o.columns:
                o = pd.concat([o, self[[col]]], axis=1)

        o.sort_index(inplace=True)

        return o

    def merge_observation(self, right, overlap="error", merge_metadata=True):
        """Merge with another observation of the same type.

        Parameters
        ----------
        right : hpd.observation.Obs
            Observation object.
        overlap : str, optional
            How to deal with overlapping timeseries or metadata with different
            values.
            Options are:
                error : Raise a ValueError
                use_left : use the part of the overlapping timeseries from self
                use_right : use the part of the overlapping timeseries from
                right
            Default is 'error'.
        merge_metadata : bool, optional
            If True the metadata of the two objects are merged. If there are
            any differences the overlap parameter is used to determine which
            metadata is used. If merge_metadata is False, the metadata of self
            is always used for the merged observation. The default is True.

        Raises
        ------
        TypeError
            when the observation types are not the same.
        ValueError
            when the time series have different values on the same date or
            different values for the same metadata.

        Returns
        -------
        Observation object.
        """

        if overlap not in ["error", "use_left", "use_right"]:
            raise ValueError(
                "invalid value for overlap, choose between error, use_left and"
                "use_right"
            )

        # check observation type
        if not isinstance(right, type(self)):
            raise TypeError(
                f"existing observation has a different type {type(self)} than"
                f"new observation {type(right)}"
            )

        # merge timeseries
        o = self._merge_timeseries(right, overlap=overlap)

        # merge metadata
        if merge_metadata:
            metadata = {key: getattr(right, key) for key in right._metadata}
            new_metadata = self.merge_metadata(metadata, overlap=overlap)
        else:
            new_metadata = {key: getattr(self, key) for key in self._metadata}

        for key, item in new_metadata.items():
            setattr(o, key, item)

        return o


class GroundwaterObs(Obs):
    """Class for groundwater quantity observations.

    Subclass of the Obs class. Has the following attributes:

    - monitoring_well: 2 tubes at one piezometer should have the same 'monitoring_well'
    - tube_nr: 2 tubes at one piezometer should have a different 'tube_nr'.
    - screen_top: top op the filter in m above date (NAP)
    - screen_bottom: bottom of the filter in m above date (NAP)
    - ground_level: surface level in m above date (NAP) (maaiveld in Dutch)
    - tube_top: top of the tube in m above date (NAP)
    - metadata_available: boolean indicating if metadata is available for
      the measurement point.
    """

    _metadata = Obs._metadata + [
        "monitoring_well",
        "tube_nr",
        "screen_top",
        "screen_bottom",
        "ground_level",
        "tube_top",
        "metadata_available",
    ]

    def __init__(self, *args, **kwargs):
        """
        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.
        """
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in GroundwaterObs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        self.monitoring_well = kwargs.pop("monitoring_well", "")
        self.tube_nr = kwargs.pop("tube_nr", "")
        self.ground_level = kwargs.pop("ground_level", np.nan)
        self.tube_top = kwargs.pop("tube_top", np.nan)
        self.screen_top = kwargs.pop("screen_top", np.nan)
        self.screen_bottom = kwargs.pop("screen_bottom", np.nan)
        self.metadata_available = kwargs.pop("metadata_available", np.nan)

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return GroundwaterObs

    @classmethod
    def from_bro(
        cls,
        bro_id,
        tube_nr=None,
        tmin="1900-01-01",
        tmax="2040-01-01",
        to_wintertime=True,
        drop_duplicate_times=True,
        only_metadata=False,
    ):
        """download BRO groundwater observations from the server.


        Parameters
        ----------
        bro_id : str
            can be a GLD id or GMW id. If a GMW id is given a tube number is
            required as well. e.g. 'GLD000000012893'.
        tube_nr : str or None, optional
            if the bro_id is a GMW object the tube number should be given.
        tmin : str or None, optional
            start date in format YYYY-MM-DD
        tmax : str or None, optional
            end date in format YYYY-MM-DD
        to_wintertime : bool, optional
            if True the time index is converted to Dutch winter time. The default
            is True.
        drop_duplicate_times : bool, optional
            if True rows with a duplicate time stamp are removed keeping only the
            first row. The default is True.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        from .io import bro

        measurements, meta = bro.get_bro_groundwater(
            bro_id,
            tube_nr,
            tmin=tmin,
            tmax=tmax,
            to_wintertime=to_wintertime,
            drop_duplicate_times=drop_duplicate_times,
            only_metadata=only_metadata,
        )

        return cls(
            measurements,
            meta=meta,
            name=meta.pop("name"),
            x=meta.pop("x"),
            y=meta.pop("y"),
            source=meta.pop("source"),
            unit=meta.pop("unit"),
            screen_bottom=meta.pop("screen_bottom"),
            screen_top=meta.pop("screen_top"),
            ground_level=meta.pop("ground_level"),
            metadata_available=meta.pop("metadata_available"),
            monitoring_well=meta.pop("monitoring_well"),
            tube_nr=meta.pop("tube_nr"),
            tube_top=meta.pop("tube_top"),
        )

    @classmethod
    def from_bronhouderportaal_bro(
        cls,
        fn_xml,
        tube_nr,
        full_meta=False,
    ):
        """load BRO groundwater metadata from XML file. Mind that
        bro_id is applicable, because file is not yet imported in BRO


        Parameters
        ----------
        fn_xml : str
            filename of XML file.
        tube_nr : int
            tube number.
        full_meta : bool
            process not only the standard metadata to ObsCollection.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        from .io import bronhouderportaal_bro

        meta = bronhouderportaal_bro.get_metadata_from_gmw(
            fn_xml, tube_nr, full_meta=full_meta
        )

        empty_df = pd.DataFrame()

        return cls(
            empty_df,
            name=meta.pop("name"),
            x=meta.pop("x"),
            y=meta.pop("y"),
            filename=meta.pop("filename"),
            source=meta.pop("source"),
            unit=meta.pop("unit"),
            screen_bottom=meta.pop("screen_bottom"),
            screen_top=meta.pop("screen_top"),
            ground_level=meta.pop("ground_level"),
            metadata_available=meta.pop("metadata_available"),
            monitoring_well=meta.pop("monitoring_well"),
            tube_nr=meta.pop("tube_nr"),
            tube_top=meta.pop("tube_top"),
            meta=meta,
        )

    @classmethod
    def from_dino(
        cls,
        fname=None,
        **kwargs,
    ):
        """download dino data from the server.

        Parameters
        ----------
        fname : str, optional
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to hydropandas.io.dino.read_dino_groundwater_csv
            if fname is not None and otherwise to hydropandas.io.dino.findMeetreeks
        """
        from .io import dino

        # read dino csv file
        measurements, meta = dino.read_dino_groundwater_csv(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_artdino_file(cls, fname=None, **kwargs):
        """read a dino csv file (artdiver style).

        Parameters
        ----------
        fname : str, optional
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to hydropandas.io._dino.read_dino_groundwater_csv
        """
        from .io import dino

        # read artdino csv file
        measurements, meta = dino.read_artdino_groundwater_csv(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_wiski(cls, fname, **kwargs):
        """
        Read data from a WISKI file.

        Parameters:
        -----------
        fname : str
            The name of the file to be read.
        sep : str, optional (default=";")
            The delimiter used to separate fields in the file.
        header_sep : str, optional (default=None)
            The delimiter used to separate fields in the header. If None, the
            function will try to automatically detect the separator.
        header_identifier : str, optional (default="#")
            The character used to identify header lines.
        read_series : bool, optional (default=True)
            Whether to read the time series data from the file.
        infer_datetime_format : bool, optional (default=True)
            Whether to infer the datetime format of the timestamp column.
        translate_dic : dict, optional (default=None)
            A dictionary mapping header field names to the desired output names.
        tz_localize : bool, optional (default=True)
            Whether to localize the datetime index to the machine's timezone.
        unit : str, optional (default="")
            The unit of measurement of the data.
        **kwargs : keyword arguments
            Additional arguments to pass to the pandas `read_csv` function.
        """
        from .io import wiski

        data, metadata = wiski.read_wiski_file(fname, **kwargs)

        return cls(data, meta=metadata, **metadata)


class WaterQualityObs(Obs):
    """class for water quality ((grond)watersamenstelling) point
    observations.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + [
        "monitoring_well",
        "tube_nr",
        "ground_level",
        "metadata_available",
    ]

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in WaterQualityObs._metadata) and (
                        key not in kwargs.keys()
                    ):
                        kwargs[key] = getattr(args[0], key)

        self.monitoring_well = kwargs.pop("monitoring_well", "")
        self.tube_nr = kwargs.pop("tube_nr", "")
        self.ground_level = kwargs.pop("ground_level", np.nan)
        self.metadata_available = kwargs.pop("metadata_available", np.nan)

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return WaterQualityObs

    @classmethod
    def from_dino(cls, fname, **kwargs):
        """read dino file with groundwater quality data.

        Parameters
        ----------
        fname : str
            dino txt filename
        kwargs : key-word arguments
            these arguments are passed to
            hydropandas.io.dino.read_dino_groundwater_quality_txt
        """
        from .io import dino

        measurements, meta = dino.read_dino_groundwater_quality_txt(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)


class WaterlvlObs(Obs):
    """class for water level point observations.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + ["monitoring_well", "metadata_available"]

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in WaterlvlObs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        self.monitoring_well = kwargs.pop("monitoring_well", "")
        self.metadata_available = kwargs.pop("metadata_available", np.nan)

        super().__init__(*args, **kwargs)

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
            these arguments are passed to hydropandas.io.dino.read_dino_waterlvl_csv
        """
        from .io import dino

        measurements, meta = dino.read_dino_waterlvl_csv(fname, **kwargs)

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
        from .io import waterinfo

        df, metadata = waterinfo.read_waterinfo_file(
            fname, return_metadata=True, **kwargs
        )
        return cls(df, meta=metadata, **metadata)


class ModelObs(Obs):
    """class for model point results.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + ["model"]

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in ModelObs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        self.model = kwargs.pop("model", "")

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return ModelObs


class MeteoObs(Obs):
    """class for meteorological timeseries.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + ["station", "meteo_var"]

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in MeteoObs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        self.station = kwargs.pop("station", np.nan)
        self.meteo_var = kwargs.pop("meteo_var", "")

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return MeteoObs

    @classmethod
    def from_knmi(
        cls,
        stn,
        meteo_var="RH",
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        interval="daily",
        inseason=False,
        use_api=True,
        raise_exceptions=True,
    ):
        """Get a MeteoObs timeseries from the KNMI meteo data.

        Parameters
        ----------
        stn : int or str
            measurement station e.g. 829.
        meteo_var : str,
            meteo variable e.g. "RH" or "EV24". See list with al options below.
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
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        List of possible variables:
            neerslagstations:
            RD    = de 24-uurs neerslagsom, gemeten van 0800 utc op de
            voorafgaande dag tot 0800 utc op de vermelde datum

            meteostations:
            DDVEC = Vectorgemiddelde windrichting in graden (360=noord,
            90=oost, 180=zuid, 270=west, 0=windstil/variabel). Zie
            http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
            / Vector mean wind direction in degrees (360=north, 90=east,
            180=south, 270=west, 0=calm/variable)
            FHVEC = Vectorgemiddelde windsnelheid (in 0.1 m/s). Zie
            http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
            / Vector mean windspeed (in 0.1 m/s)
            FG    = Etmaalgemiddelde windsnelheid (in 0.1 m/s) / Daily mean
            windspeed (in 0.1 m/s)
            FHX   = Hoogste uurgemiddelde windsnelheid (in 0.1 m/s) / Maximum
            hourly mean windspeed (in 0.1 m/s)
            FHXH  = Uurvak waarin FHX is gemeten / Hourly division in which
            FHX was measured
            FHN   = Laagste uurgemiddelde windsnelheid (in 0.1 m/s) / Minimum
            hourly mean windspeed (in 0.1 m/s)
            FHNH  = Uurvak waarin FHN is gemeten / Hourly division in which
            FHN was measured
            FXX   = Hoogste windstoot (in 0.1 m/s) / Maximum wind gust (in
            0.1 m/s)
            FXXH  = Uurvak waarin FXX is gemeten / Hourly division in which
            FXX was measured
            TG    = Etmaalgemiddelde temperatuur (in 0.1 graden Celsius) /
            Daily mean temperature in (0.1 degrees Celsius)
            TN    = Minimum temperatuur (in 0.1 graden Celsius) / Minimum
            temperature (in 0.1 degrees Celsius)
            TNH   = Uurvak waarin TN is gemeten / Hourly division in which TN
            was measured
            TX    = Maximum temperatuur (in 0.1 graden Celsius) / Maximum
            temperature (in 0.1 degrees Celsius)
            TXH   = Uurvak waarin TX is gemeten / Hourly division in which TX
            was measured
            T10N  = Minimum temperatuur op 10 cm hoogte (in 0.1 graden
            Celsius) / Minimum temperature at 10 cm above surface (in 0.1
            degrees Celsius)
            T10NH = 6-uurs tijdvak waarin T10N is gemeten / 6-hourly division
            in which T10N was measured; 6=0-6 UT, 12=6-12 UT, 18=12-18 UT,
            24=18-24 UT
            SQ    = Zonneschijnduur (in 0.1 uur) berekend uit de globale
            straling (-1 voor <0.05 uur) / Sunshine duration (in 0.1 hour)
            calculated from global radiation (-1 for <0.05 hour)
            SP    = Percentage van de langst mogelijke zonneschijnduur /
            Percentage of maximum potential sunshine duration
            Q     = Globale straling (in J/cm2) / Global radiation (in J/cm2)
            DR    = Duur van de neerslag (in 0.1 uur) / Precipitation duration
            (in 0.1 hour)
            RH    = Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) /
            Daily precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
            RHX   = Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05
            mm) / Maximum hourly precipitation amount (in 0.1 mm) (-1 for
            <0.05 mm)
            RHXH  = Uurvak waarin RHX is gemeten / Hourly division in which
            RHX was measured
            PG    = Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1
            hPa) berekend uit 24 uurwaarden / Daily mean sea level pressure
            (in 0.1 hPa) calculated from 24 hourly values
            PX    = Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau
            (in 0.1 hPa) / Maximum hourly sea level pressure (in 0.1 hPa)
            PXH   = Uurvak waarin PX is gemeten / Hourly division in which PX
            was measured
            PN    = Laagste uurwaarde van de luchtdruk herleid tot zeeniveau
            (in 0.1 hPa) / Minimum hourly sea level pressure (in 0.1 hPa)
            PNH   = Uurvak waarin PN is gemeten / Hourly division in which PN
            was measured
            VVN   = Minimum opgetreden zicht / Minimum visibility; 0: <100 m,
            1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
            56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
            ..., 89: >70 km)
            VVNH  = Uurvak waarin VVN is gemeten / Hourly division in which
            VVN was measured
            VVX   = Maximum opgetreden zicht / Maximum visibility; 0: <100 m,
            1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
            56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
            ..., 89: >70 km)
            VVXH  = Uurvak waarin VVX is gemeten / Hourly division in which
            VVX was measured
            NG    = Etmaalgemiddelde bewolking (bedekkingsgraad van de
            bovenlucht in achtsten, 9=bovenlucht onzichtbaar) / Mean daily
            cloud cover (in octants, 9=sky invisible)
            UG    = Etmaalgemiddelde relatieve vochtigheid (in procenten) /
            Daily mean relative atmospheric humidity (in percents)
            UX    = Maximale relatieve vochtigheid (in procenten) / Maximum
            relative atmospheric humidity (in percents)
            UXH   = Uurvak waarin UX is gemeten / Hourly division in which UX
            was measured
            UN    = Minimale relatieve vochtigheid (in procenten) / Minimum
            relative atmospheric humidity (in percents)
            UNH   = Uurvak waarin UN is gemeten / Hourly division in which UN
            was measured
            EV24  = Referentiegewasverdamping (Makkink) (in 0.1 mm) /
            Potential evapotranspiration (Makkink) (in 0.1 mm)

        Returns
        -------
        MeteoObs object with meteorological observations
        """
        from .io import knmi

        if fill_missing_obs and interval == "hourly":
            fill_missing_obs = False

        if interval == "hourly" and meteo_var == "RD":
            raise NotImplementedError(
                "hourly values are not available for precipitation stations"
            )

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
        }

        ts, meta = knmi.get_knmi_timeseries_stn(
            stn, meteo_var, startdate, enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=meteo_var,
        )

    @classmethod
    def from_nearest_xy(
        cls,
        xy,
        meteo_var,
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        interval="daily",
        inseason=False,
        use_precipitation_stn=True,
        use_api=True,
        raise_exceptions=False,
    ):
        """Get a MeteoObs object from the KNMI station closest to the given
        (x,y) coordinates.

        Parameters
        ----------
        xy : list. tuple or numpy array of int or floats
            xy coordinates. e.g. [10.1,25.05]
        meteo_var : str,
            meteo variable e.g. "RH" or "EV24". See list with al options below.
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
        use_precipitation_stn : bool, optional
            if True a combination of neerslagstations and meteo stations are used.
            If False only meteo stations are used to obtain precipitation data.
            Default is True.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        Returns
        -------
        MeteoObs object with meteorological observations
        """
        from .io import knmi

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
            "use_precipitation_stn": use_precipitation_stn,
        }

        ts, meta = knmi.get_knmi_timeseries_xy(
            xy, meteo_var, startdate, enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=meteo_var,
        )

    @classmethod
    def from_obs(
        cls,
        obs,
        meteo_var,
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        interval="daily",
        inseason=False,
        use_precipitation_stn=True,
        use_api=True,
        raise_exceptions=False,
    ):
        """Get a MeteoObs object with measurements from the KNMI station
        closest to the given observation. Uses the x and y coordinates of the
        observation to obtain the nearest KNMI station. Uses the start- and
        enddate of the observation as start- and enddate of the meteo time
        series (unless startdate and enddatet are specified explicitly).

        Parameters
        ----------
        obs : hydropandas.Obs
            Observation object.
        meteo_var : str,
            meteo variable e.g. "RH" or "EV24". See list with al options below.
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
        use_precipitation_stn : bool, optional
            if True a combination of neerslagstations and meteo stations are used.
            If False only meteo stations are used to obtain precipitation data.
            Default is True.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        List of possible variables:
            neerslagstations:
            RD    = de 24-uurs neerslagsom, gemeten van 0800 utc op de
            voorafgaande dag tot 0800 utc op de vermelde datum

            meteostations:
            DDVEC = Vectorgemiddelde windrichting in graden (360=noord,
            90=oost, 180=zuid, 270=west, 0=windstil/variabel). Zie
            http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
            / Vector mean wind direction in degrees (360=north, 90=east,
            180=south, 270=west, 0=calm/variable)
            FHVEC = Vectorgemiddelde windsnelheid (in 0.1 m/s). Zie
            http://www.knmi.nl/kennis-en-datacentrum/achtergrond/klimatologische-brochures-en-boeken
            / Vector mean windspeed (in 0.1 m/s)
            FG    = Etmaalgemiddelde windsnelheid (in 0.1 m/s) / Daily mean
            windspeed (in 0.1 m/s)
            FHX   = Hoogste uurgemiddelde windsnelheid (in 0.1 m/s) / Maximum
            hourly mean windspeed (in 0.1 m/s)
            FHXH  = Uurvak waarin FHX is gemeten / Hourly division in which
            FHX was measured
            FHN   = Laagste uurgemiddelde windsnelheid (in 0.1 m/s) / Minimum
            hourly mean windspeed (in 0.1 m/s)
            FHNH  = Uurvak waarin FHN is gemeten / Hourly division in which
            FHN was measured
            FXX   = Hoogste windstoot (in 0.1 m/s) / Maximum wind gust (in
            0.1 m/s)
            FXXH  = Uurvak waarin FXX is gemeten / Hourly division in which
            FXX was measured
            TG    = Etmaalgemiddelde temperatuur (in 0.1 graden Celsius) /
            Daily mean temperature in (0.1 degrees Celsius)
            TN    = Minimum temperatuur (in 0.1 graden Celsius) / Minimum
            temperature (in 0.1 degrees Celsius)
            TNH   = Uurvak waarin TN is gemeten / Hourly division in which TN
            was measured
            TX    = Maximum temperatuur (in 0.1 graden Celsius) / Maximum
            temperature (in 0.1 degrees Celsius)
            TXH   = Uurvak waarin TX is gemeten / Hourly division in which TX
            was measured
            T10N  = Minimum temperatuur op 10 cm hoogte (in 0.1 graden
            Celsius) / Minimum temperature at 10 cm above surface (in 0.1
            degrees Celsius)
            T10NH = 6-uurs tijdvak waarin T10N is gemeten / 6-hourly division
            in which T10N was measured; 6=0-6 UT, 12=6-12 UT, 18=12-18 UT,
            24=18-24 UT
            SQ    = Zonneschijnduur (in 0.1 uur) berekend uit de globale
            straling (-1 voor <0.05 uur) / Sunshine duration (in 0.1 hour)
            calculated from global radiation (-1 for <0.05 hour)
            SP    = Percentage van de langst mogelijke zonneschijnduur /
            Percentage of maximum potential sunshine duration
            Q     = Globale straling (in J/cm2) / Global radiation (in J/cm2)
            DR    = Duur van de neerslag (in 0.1 uur) / Precipitation duration
            (in 0.1 hour)
            RH    = Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) /
            Daily precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
            RHX   = Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05
            mm) / Maximum hourly precipitation amount (in 0.1 mm) (-1 for
            <0.05 mm)
            RHXH  = Uurvak waarin RHX is gemeten / Hourly division in which
            RHX was measured
            PG    = Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1
            hPa) berekend uit 24 uurwaarden / Daily mean sea level pressure
            (in 0.1 hPa) calculated from 24 hourly values
            PX    = Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau
            (in 0.1 hPa) / Maximum hourly sea level pressure (in 0.1 hPa)
            PXH   = Uurvak waarin PX is gemeten / Hourly division in which PX
            was measured
            PN    = Laagste uurwaarde van de luchtdruk herleid tot zeeniveau
            (in 0.1 hPa) / Minimum hourly sea level pressure (in 0.1 hPa)
            PNH   = Uurvak waarin PN is gemeten / Hourly division in which PN
            was measured
            VVN   = Minimum opgetreden zicht / Minimum visibility; 0: <100 m,
            1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
            56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
            ..., 89: >70 km)
            VVNH  = Uurvak waarin VVN is gemeten / Hourly division in which
            VVN was measured
            VVX   = Maximum opgetreden zicht / Maximum visibility; 0: <100 m,
            1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km,
            56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,
            ..., 89: >70 km)
            VVXH  = Uurvak waarin VVX is gemeten / Hourly division in which
            VVX was measured
            NG    = Etmaalgemiddelde bewolking (bedekkingsgraad van de
            bovenlucht in achtsten, 9=bovenlucht onzichtbaar) / Mean daily
            cloud cover (in octants, 9=sky invisible)
            UG    = Etmaalgemiddelde relatieve vochtigheid (in procenten) /
            Daily mean relative atmospheric humidity (in percents)
            UX    = Maximale relatieve vochtigheid (in procenten) / Maximum
            relative atmospheric humidity (in percents)
            UXH   = Uurvak waarin UX is gemeten / Hourly division in which UX
            was measured
            UN    = Minimale relatieve vochtigheid (in procenten) / Minimum
            relative atmospheric humidity (in percents)
            UNH   = Uurvak waarin UN is gemeten / Hourly division in which UN
            was measured
            EV24  = Referentiegewasverdamping (Makkink) (in 0.1 mm) /
            Potential evapotranspiration (Makkink) (in 0.1 mm)

        Returns
        -------
        MeteoObs object with meteorological observations
        """
        from .io import knmi

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
            "use_precipitation_stn": use_precipitation_stn,
        }

        xy = (obs.x, obs.y)

        if startdate is None:
            startdate = obs.index[0]
        if enddate is None:
            enddate = obs.index[-1]

        ts, meta = knmi.get_knmi_timeseries_xy(
            xy, meteo_var, startdate, enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=meteo_var,
        )

    @classmethod
    def from_knmi_file(cls, fname, meteo_var="RH", startdate=None, enddate=None):
        """Get a MeteoObs timeseries from the KNMI meteo data.

        Parameters
        ----------
        fname : str
            full path of a knmi .txt file
        meteo_var : str,
            meteo variable e.g. "RH" or "EV24".
        startdate : str, datetime or None, optional
            start date of observations. The default is None.
        enddate : str, datetime or None, optional
            end date of observations. The default is None.

        Returns
        -------
        MeteoObs object with meteorological observations
        """
        from .io import knmi

        if not fname.endswith(".txt"):
            fname += ".txt"

        knmi_df, meta = knmi.read_knmi_timeseries_file(
            fname, meteo_var, startdate, enddate
        )

        return cls(
            knmi_df,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=meteo_var,
        )


class EvaporationObs(MeteoObs):
    """class for evaporation timeseries.

    Subclass of the MeteoObs class
    """

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in EvaporationObs._metadata) and (key not in kwargs.keys()):
                        kwargs[key] = getattr(args[0], key)

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return EvaporationObs

    @classmethod
    def from_knmi(
        cls,
        stn,
        et_type="EV24",
        meteo_var=None,
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        interval="daily",
        inseason=False,
        use_api=True,
        raise_exceptions=False,
    ):
        """Get an EvaporationObs timeseries from the KNMI evaporation in m.

        Parameters
        ----------
        stn : int or str
            measurement station e.g. 829.
        et_type: str
            type of evapotranspiration to get from KNMI. Choice between
            'EV24', 'penman', 'makkink' or 'hargraves'. Defaults to
            'EV24' which collects the KNMI Makkink EV24 evaporation.
        meteo_var : None
            not used in this method.
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
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).


        Returns
        -------
        EvaporationObs object with an evaporation time series and attributes
        """
        from .io import knmi

        if interval == "hourly":
            fill_missing_obs = False

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
        }

        ts, meta = knmi.get_evaporation(
            stn, et_type, start=startdate, end=enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=et_type,
        )

    @classmethod
    def from_nearest_xy(
        cls,
        xy,
        et_type="EV24",
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        interval="daily",
        inseason=False,
        use_api=True,
        raise_exceptions=False,
    ):
        """Get an EvaporationObs object with evaporation measurements from the
        KNMI station closest to the given (x,y) coordinates.

        Parameters
        ----------
        xy : list. tuple or numpy array of int or floats
            xy coordinates. e.g. (10.1,25.05)
        et_type: str
            type of evapotranspiration to get from KNMI. Choice between
            'EV24', 'penman', 'makkink' or 'hargraves'. Defaults to
            'EV24' which collects the KNMI Makkink EV24 evaporation.
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
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.


        Returns
        -------
        EvaporationObs object with an evaporation time series and attributes
        """
        from .io import knmi

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
        }

        stn = knmi.get_n_nearest_stations_xy(xy, "EV24")[0]
        ts, meta = knmi.get_evaporation(
            stn, et_type, start=startdate, end=enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var=et_type,
        )

    @classmethod
    def from_obs(
        cls,
        obs,
        et_type="EV24",
        startdate=None,
        enddate=None,
        fill_missing_obs=True,
        use_api=True,
        interval="daily",
        inseason=False,
        raise_exceptions=False,
    ):
        """Get an EvaporationObs object with evaporation measurements from the
        KNMI station closest to the given observation. Uses the x and y
        coordinates of the observation to obtain the nearest KNMI evaporation
        time series. Uses the start- and enddate of the observation as start-
        and enddate of the evaporation time series (unless startdate and
        enddatet are specified explicitly).

        Parameters
        ----------
        obs : hydropandas.Obs
            Observation object.
        et_type: str
            type of evapotranspiration to get from KNMI. Choice between
            'EV24', 'penman', 'makkink' or 'hargraves'. Defaults to
            'EV24' which collects the KNMI Makkink EV24 evaporation.
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
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary directory and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        Returns
        -------
        EvaporationObs object with an evaporation time series and attributes
        """
        from .io import knmi

        settings = {
            "fill_missing_obs": fill_missing_obs,
            "interval": interval,
            "inseason": inseason,
            "use_api": use_api,
            "raise_exceptions": raise_exceptions,
        }

        xy = (obs.x, obs.y)

        stn = knmi.get_n_nearest_stations_xy(xy, "EV24")[0]
        ts, meta = knmi.get_evaporation(
            stn, et_type, start=startdate, end=enddate, settings=settings
        )

        return cls(
            ts,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            meteo_var=et_type,
            source=meta["source"],
            unit=meta["unit"],
        )

    @classmethod
    def from_knmi_file(cls, fname, startdate=None, enddate=None):
        """Get a EvaporationObs timeseries from the KNMI meteo data.

        Parameters
        ----------
        fname : str
            full path of a knmi .txt file
        startdate : str, datetime or None, optional
            start date of observations. The default is None.
        enddate : str, datetime or None, optional
            end date of observations. The default is None.

        Returns
        -------
        MeteoObs object with meteorological observations
        """
        from .io import knmi

        if not fname.endswith(".txt"):
            fname += ".txt"

        knmi_df, meta = knmi.read_knmi_timeseries_file(
            fname, "EV24", startdate, enddate
        )

        return cls(
            knmi_df,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var="EV24",
        )


class PrecipitationObs(MeteoObs):
    """class for precipitation timeseries.

    Subclass of the MeteoObs class
    """

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            if isinstance(args[0], Obs):
                for key in args[0]._metadata:
                    if (key in PrecipitationObs._metadata) and (
                        key not in kwargs.keys()
                    ):
                        kwargs[key] = getattr(args[0], key)

        super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return PrecipitationObs

    @classmethod
    def from_knmi(cls, stn, stn_type="meteo", meteo_var=None, **kwargs):
        """Get a PrecipitationObs timeseries from the KNMI precipitation. The
        precipitation is the Daily precipitation amount (in 0.1 mm) (-1 for.

        <0.05 mm).

        There are 3 different ways to obtain precipitation data from the knmi:
            1. Daily data from precipitation (neerslag) stations
            2. Daily data from meteo stations
            3. Hourly data from meteo stations

        1. If you want to get data from a neerslagstation stn should be the
        station number with '_neerslag_station' at the end e.g. for De Bilt:
            stn = '550_neerslag_station'.
        2. If you want to get daily data from a meteo station stn should be the
        number of the meteo station, can be string or integer
        3. If you want to get hourly data from a meteo station, stn should be
        the number of the meteo station and interval should 'hourly'.

        More information about the differences between neerslag and meteo
        stations can be found in the 02_knmi_observations notebook inside the
        examples directory on github:
            https://github.com/ArtesiaWater/hydropandas

        Parameters
        ----------
        stn : int or str
            measurement station e.g. 829.
        stn_type : str, optional
            type of measurements station. Can be 'meteo' or 'precipitation'.
            Default is 'meteo'.
        meteo_var : None
            not used.
        **kwargs:
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
                if False a text file is downloaded into a temporary directory and
                the data is read from there. Default is True since the api is
                back online (July 2021).
            raise_exceptions : bool, optional
                if True you get errors when no data is returned. The default is False.

        Returns
        -------
        PrecipitationObs object with a precipitation time series and attributes
        """
        if stn_type == "meteo":
            meteo_var = "RH"
            kwargs.pop("meteo_var", None)
        elif stn_type == "precipitation":
            meteo_var = "RD"
            kwargs.pop("meteo_var", None)
        else:
            raise ValueError(f"invalid measurement station type -> {stn_type}")

        return super().from_knmi(stn, meteo_var=meteo_var, **kwargs)

    @classmethod
    def from_nearest_xy(cls, xy, stn_type="meteo", **kwargs):
        """Get a PrecipitationObs object with precipitation measurements from
        the meteo or precipitation station closest to the given (x,y)
        coordinates.

        Parameters
        ----------
        xy : list. tuple or numpy array of ints or floats
            xy coordinates. e.g. [10.1,25.05]
        stn_type : str, optional
            type of measurements station. Can be 'meteo' or 'precipitation'.
            Default is 'meteo'.
        **kwargs:
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
            use_api : bool, optional
                if True the api is used to obtain the data, API documentation is here:
                    https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
                if False a text file is downloaded into a temporary directory and
                the data is read from there. Default is True since the api is
                back online (July 2021).
            raise_exceptions : bool, optional
                if True you get errors when no data is returned. The default is False.

        Returns
        -------
        PrecipitationObs object with a precipitation time series and attributes
        """
        if stn_type == "meteo":
            meteo_var = "RH"
        elif stn_type == "precipitation":
            meteo_var = "RD"
        else:
            raise ValueError(f"invalid measurement station type -> {stn_type}")

        return super().from_nearest_xy(xy, meteo_var=meteo_var, **kwargs)

    @classmethod
    def from_obs(cls, obs, stn_type="meteo", **kwargs):
        """Get a PrecipitationObs object with precipitation measurements from the
        the meteo or precipitation station closest to the given observation.
        Uses the x and y coordinates of the observation to obtain the nearest
        KNMI precipitation time series. Uses the start- and enddate of the
        observation as start- and enddate of the time series (unless startdate
        and enddate are specified explicitly).

        Parameters
        ----------
        obs : hydropandas.Obs
            Observation object.
        stn_type : str, optional
            type of measurements station. Can be 'meteo' or 'precipitation'.
            Default is 'meteo'.
        **kwargs
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
            use_api : bool, optional
                if True the api is used to obtain the data, API documentation is here:
                    https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
                if False a text file is downloaded into a temporary directory and
                the data is read from there. Default is True since the api is
                back online (July 2021).
            raise_exceptions : bool, optional
                if True you get errors when no data is returned. The default is False.

        Returns
        -------
        PrecipitationObs object with a precipitation time series and attributes
        """

        if stn_type == "meteo":
            meteo_var = "RH"
        elif stn_type == "precipitation":
            meteo_var = "RD"
        else:
            raise ValueError(f"invalid measurement station type -> {stn_type}")

        return super().from_obs(obs, meteo_var=meteo_var, **kwargs)

    @classmethod
    def from_knmi_file(cls, fname, startdate=None, enddate=None):
        """Get a PrecipitationObs timeseries from the KNMI meteo data.

        Parameters
        ----------
        fname : str
            full path of a knmi .txt file
        startdate : str, datetime or None, optional
            start date of observations. The default is None.
        enddate : str, datetime or None, optional
            end date of observations. The default is None.

        Returns
        -------
        PrecipitationObs object with precipitation observations
        """
        from .io import knmi

        if not fname.endswith(".txt"):
            fname += ".txt"

        knmi_df, meta = knmi.read_knmi_timeseries_file(fname, "RD", startdate, enddate)

        return cls(
            knmi_df,
            meta=meta,
            station=meta["station"],
            x=meta["x"],
            y=meta["y"],
            name=meta["name"],
            source=meta["source"],
            unit=meta["unit"],
            meteo_var="RD",
        )
