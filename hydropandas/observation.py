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
import numbers
import os
import warnings
from _io import StringIO
from typing import List, Optional

import numpy as np
import pandas as pd
from pandas._config import get_option
from pandas.api.types import is_numeric_dtype
from pandas.io.formats import console

logger = logging.getLogger(__name__)


class Obs(pd.DataFrame):
    """Generic class for a time series with measurements at a certain location.

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
        """Constructor of Obs class.

        *args must be input for the pandas.DataFrame constructor, **kwargs can be one of
        the attributes listed in _metadata or keyword arguments for the constructor of a
        pandas.DataFrame.
        """
        if (len(args) > 0) and isinstance(args[0], Obs):
            for key in args[0]._metadata:
                if (key in Obs._metadata) and (key not in kwargs):
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

    def _repr_html_(self, collapse=False):
        """Uses the pandas DataFrame html representation with the metadata prepended."""
        obs_type = f'<p style="color:#808080";>hydropandas.{type(self).__name__}</p>\n'

        metadata_dic = {key: getattr(self, key) for key in self._metadata}
        metadata_dic.pop("meta")
        metadata_df = pd.DataFrame(
            columns=[metadata_dic.pop("name")],
            index=metadata_dic.keys(),
            data=metadata_dic.values(),
        )
        metadata = metadata_df._repr_html_()

        observations = super()._repr_html_()
        if collapse:
            collapse_button_meta = (
                '<button type="button" class="collapsible '
                'active" id="meta"><i class="arrow right">'
                "</i> Metadata</button>\n"
            )
            collapse_button_obs = (
                '<button type="button" class="collapsible '
                'active" id="obs"><i class="arrow right">'
                "</i> Observations</button>\n"
            )

            with open(
                os.path.join(os.path.dirname(__file__), "static/style.css"), "r"
            ) as fo:
                css_arrow = fo.read()

            metadata = metadata.replace(
                "<div>\n<style scoped>", '<div  style="display: none;">\n<style scoped>'
            )
            observations = observations.replace(
                "<div>\n<style scoped>", '<div  style="display: none;">\n<style scoped>'
            )

            with open(
                os.path.join(os.path.dirname(__file__), "static/js_collapse.html"), "r"
            ) as fo:
                js_collapse_button = fo.read()

            return (
                obs_type
                + css_arrow
                + collapse_button_meta
                + metadata
                + "<br>\n"
                + collapse_button_obs
                + observations
                + js_collapse_button
            )
        return obs_type + metadata + "<br>\n" + observations

    @property
    def _constructor(self):
        return Obs

    def _get_first_numeric_col_name(self):
        """Get the first numeric column name of the observations.

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
        """Create a copy of the observation.

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
        """Get dictionary with registered attributes and their values of an Obs object.

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
        """Merge the metadata of an Obs object with metadata from another Obs object.

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
                    # check if both are nan
                    if isinstance(v1, numbers.Number) and isinstance(
                        v2, numbers.Number
                    ):
                        if np.isnan(v1) and np.isnan(v2):
                            continue
                    same_metadata = False
                    if overlap == "error":
                        raise ValueError(
                            f"left observation {key} differs from right observation"
                        )
                    elif overlap == "use_left":
                        logger.info(
                            f"left observation {key} differs from right "
                            "observation, use left"
                        )
                        new_metadata[key] = v1
                    elif overlap == "use_right":
                        logger.info(
                            f"left observation {key} differs from right "
                            "observation, use right"
                        )
                        new_metadata[key] = v2
                else:
                    new_metadata[key] = v1
            except TypeError:
                same_metadata = False
                if overlap == "error":
                    raise ValueError(
                        f"left observation {key} differs from right observation"
                    )
                elif overlap == "use_left":
                    logger.info(
                        f"left observation {key} differs from right "
                        "observation, use left"
                    )
                    new_metadata[key] = v1
                elif overlap == "use_right":
                    logger.info(
                        f"left observation {key} differs from right "
                        "observation, use right"
                    )
                    new_metadata[key] = v2
        if same_metadata:
            logger.info("left and right observation have the same metadata")

        return new_metadata

    def _merge_timeseries(self, right, overlap="error"):
        """Merge two timeseries.

        Parameters
        ----------
        right : hpd.observation.Obs
            Observation object.
        overlap : str, optional
            How to deal with overlapping timeseries with different values.
            Options are:
                error : Raise a ValueError
                use_left : use the part of the overlapping timeseries from self
                use_right : use the part of the overlapping timeseries from
                right
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
            logger.info("left and right observation have the same time series")
            return self

        logger.info("right observation has a different time series")
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
                "invalid value for overlap, choose between error, use_left anduse_right"
            )

        # check observation type
        if not isinstance(right, type(self)):
            raise TypeError(
                f"observation left has a different type {type(self)} than"
                f"observation right {type(right)}"
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
        """Constructor for ObsCollection.

        Parameters
        ----------
        *args must be input for the pandas.DataFrame constructor
        **kwargs can be one of the attributes listed in _metadata or keyword arguments
        for the constructor of a pandas.DataFrame.
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
        """Download BRO groundwater observations from the server.

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
        only_metadata : bool, optional
            if True only metadata is returned and no time series data. The
            default is False

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
    def from_lizard(
        cls,
        code,
        tube_nr=None,
        tmin=None,
        tmax=None,
        type_timeseries="merge",
        only_metadata=False,
    ):
        """Extracts the metadata and timeseries of a observation well from a LIZARD-API
        based on the code of a monitoring well.

        Parameters
        ----------
        code : str
            code of the measuring well
        tube_nr : int, optional
            select specific tube top
            Default selects tube_nr = 1
        tmin : str YYYY-m-d, optional
            start of the observations, by default the entire serie is returned
        tmax : Ttr YYYY-m-d, optional
            end of the observations, by default the entire serie is returned
        type_timeseries : str, optional
            hand: returns only hand measurements
            diver: returns only diver measurements
            merge: the hand and diver measurements into one time series (default)
            combine: keeps hand and diver measurements separated
        only_metadata : bool, optional
            if True only metadata is returned and no time series data. The
            default is False.

        Returns
        -------
        ObsCollection
            Returns a DataFrame with metadata and timeseries
        """

        from .io import lizard

        measurements, meta = lizard.get_lizard_groundwater(
            code,
            tube_nr,
            tmin,
            tmax,
            type_timeseries,
            only_metadata=only_metadata,
        )
        return cls(
            measurements,
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
            meta=meta,
        )

    @classmethod
    def from_bronhouderportaal_bro(
        cls,
        path,
        tube_nr,
        full_meta=False,
    ):
        """Load BRO groundwater metadata from XML file. Mind that bro_id is applicable,
        because file is not yet imported in BRO.

        Parameters
        ----------
        path : str
            filepath of XML file.
        tube_nr : int
            tube number.
        full_meta : bool
            process not only the standard metadata to ObsCollection.

        Returns
        -------
        ObsCollection
            ObsCollection containing observations from XML file.
        """

        from .io import bronhouderportaal_bro

        meta = bronhouderportaal_bro.get_metadata_from_gmw(
            path, tube_nr, full_meta=full_meta
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
        path=None,
        **kwargs,
    ):
        """Download dino data from the server.

        Parameters
        ----------
        path : str, optional
            path of dino csv file
        kwargs : key-word arguments
            these arguments are passed to hydropandas.io.dino.read_dino_groundwater_csv
            if path is not None and otherwise to hydropandas.io.dino.findMeetreeks
        """
        from .io import dino

        # read dino csv file
        measurements, meta = dino.read_dino_groundwater_csv(path, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_artdino_file(cls, path=None, **kwargs):
        """Read a dino csv file (artdiver style).

        Parameters
        ----------
        path : str, optional
            path of dino csv filename
        kwargs : key-word arguments
            these arguments are passed to hydropandas.io._dino.read_dino_groundwater_csv
        """
        from .io import dino

        # read artdino csv file
        measurements, meta = dino.read_artdino_groundwater_csv(path, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_wiski(cls, path, **kwargs):
        """Read data from a WISKI file.

        Parameters:
        -----------
        path : str
            The path of the file to be read.
        sep : str, optional (default=";")
            The delimiter used to separate fields in the file.
        header_sep : str, optional (default=None)
            The delimiter used to separate fields in the header. If None, the
            function will try to automatically detect the separator.
        header_identifier : str, optional (default="#")
            The character used to identify header lines.
        read_series : bool, optional (default=True)
            Whether to read the time series data from the file.
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

        data, metadata = wiski.read_wiski_file(path, **kwargs)

        return cls(data, meta=metadata, **metadata)

    @classmethod
    def from_pastastore(cls, pstore, libname, name, metadata_mapping=None):
        """Read item from pastastore library.

        Parameters
        ----------
        pstore : pastastore.PastaStore
            pastastore object
        libname : str
            name of library containinig item
        name : str
            name of item
        metadata_mapping : dict, optional
            dictionary containing map between metadata field names in pastastore (keys)
            and metadata field names expected by hydropandas (values), by default None.
        """
        from .io import pastas

        data, metadata = pastas.read_pastastore_item(pstore, libname, name)

        if metadata_mapping is not None:
            for pstore_name, oc_name in metadata_mapping.items():
                metadata[oc_name] = metadata.get(pstore_name, None)

        metadata["source"] = "pastastore"
        kwargs = {}
        for key, value in metadata.items():
            if key in cls._metadata:
                kwargs[key] = value

        return cls(data, meta=metadata, **kwargs)

    @classmethod
    def from_solinst(
        cls,
        path,
        transform_coords=True,
        screen_bottom=None,
        screen_top=None,
        ground_level=None,
        tube_nr=None,
        tube_top=None,
    ):
        """Read data from Solinst xle file.

        Parameters
        ----------
        path : str
            path to file (file can zip or xle)

        """
        from .io import solinst

        df, meta = solinst.read_solinst_file(path, transform_coords=transform_coords)

        return cls(
            df,
            meta=meta,
            name=meta.pop("name"),
            x=meta.pop("x"),
            y=meta.pop("y"),
            filename=meta.pop("filename"),
            source=meta.pop("source"),
            unit=meta.pop("unit"),
            screen_bottom=screen_bottom,
            screen_top=screen_top,
            ground_level=ground_level,
            metadata_available=meta.pop("metadata_available"),
            monitoring_well=meta.pop("monitoring_well"),
            tube_nr=tube_nr,
            tube_top=tube_top,
        )


class WaterQualityObs(Obs):
    """Class for water quality ((grond)watersamenstelling) point observations.

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
    def from_dino(cls, path, **kwargs):
        """Read dino file with groundwater quality data.

        Parameters
        ----------
        path : str
            path of dino txt filename
        kwargs : key-word arguments
            these arguments are passed to
            hydropandas.io.dino.read_dino_groundwater_quality_txt
        """
        from .io import dino

        measurements, meta = dino.read_dino_groundwater_quality_txt(path, **kwargs)

        return cls(measurements, meta=meta, **meta)


class WaterlvlObs(Obs):
    """Class for water level point observations.

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
    def from_dino(cls, path, **kwargs):
        """Read a dino file with waterlvl data.

        Parameters
        ----------
        path : str
            path of dino csv filename
        kwargs : key-word arguments
            these arguments are passed to hydropandas.io.dino.read_dino_waterlvl_csv
        """
        from .io import dino

        measurements, meta = dino.read_dino_waterlvl_csv(path, **kwargs)

        return cls(measurements, meta=meta, **meta)

    @classmethod
    def from_waterinfo(cls, path, **kwargs):
        """Read data from waterinfo csv-file or zip.

        Parameters
        ----------
        path : str
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
            path, return_metadata=True, **kwargs
        )
        return cls(df, meta=metadata, **metadata)


class ModelObs(Obs):
    """Class for model point results.

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
    """Class for meteorological timeseries.

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
        meteo_var,
        stn=None,
        fname=None,
        xy=None,
        start=None,
        end=None,
        fill_missing_obs=False,
        interval="daily",
        use_api=True,
        raise_exceptions=True,
        startdate=None,
        enddate=None,
    ):
        """Get a MeteoObs timeseries from the KNMI meteo data.

        Parameters
        ----------
        meteo_var : str
            meteo variable e.g. "RH" or "EV24". For a list of possible
            variables see the hydropandas.read_knmi function.
        stn : int, str or None, optional
            measurement station e.g. 829. The default is None.
        fname : str, path object, file-like object or None, optional
            filename of a knmi file. The default is None.
        xy : list, tuple or None, optional
            RD coördinates of a location in the Netherlands. The station nearest
            to this location used. The Default is None.
        start : str, datetime or None, optional
            start date of observations. The default is None.
        end : str, datetime or None, optional
            end date of observations. The default is None.
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False. Note: if the given stn has no data between start and
            end the data from nearby stations is used. In this case the metadata of the
            Observation is the metadata from the nearest station that has any
            measurement in the given period.
        interval : str, optional
            desired time interval for observations. Options are 'daily' and
            'hourly'. The default is 'daily'.
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

        if startdate is not None:
            warnings.warn(
                "please use 'start' instead of startdate'", DeprecationWarning
            )
            start = startdate
        if enddate is not None:
            warnings.warn("please use 'end' instead of enddate'", DeprecationWarning)
            end = enddate

        ts, meta = knmi.get_knmi_obs(
            stn=stn,
            fname=fname,
            xy=xy,
            meteo_var=meteo_var,
            start=start,
            end=end,
            fill_missing_obs=fill_missing_obs,
            interval=interval,
            use_api=use_api,
            raise_exceptions=raise_exceptions,
        )

        return cls(
            ts,
            meta=meta,
            station=meta.pop("station"),
            x=meta.pop("x"),
            y=meta.pop("y"),
            name=meta.pop("name"),
            source=meta.pop("source"),
            unit=meta.pop("unit") if "unit" in meta else "",
            meteo_var=meteo_var,
        )

    @classmethod
    def from_wow(
        cls,
        meteo_var: str,
        stn: str = None,
        xy: List[float] = None,
        start: Optional[pd.Timestamp] = None,
        end: Optional[pd.Timestamp] = None,
    ):
        """Get a MeteoObs timeseries from a wow.knmi.nl station.

        Parameters
        ----------
        meteo_var : str
            wow meteo variable
        stn : Optional[int, str], optional
            station name
        xy : Optinal[List[float]], optinial
            longitude latitude of location [lon, lat] eg: [4.85, 51.95]
        start : Optional[pd.Timestamp], optional
            start date of observations, by default None
        end : Optional[pd.Timestamp], optional
            start date of observations, by default None

        Returns
        -------
        MeteoObs
        """
        from .io import wow

        wow_df, meta = wow.get_wow(
            stn=stn, xy=xy, meteo_var=meteo_var, start=start, end=end
        )

        return cls(
            wow_df,
            meta=meta,
            station=meta["site"]["id"],
            x=meta["site"]["geo"]["coordinates"][1],
            y=meta["site"]["geo"]["coordinates"][0],
            name=meta["site"]["name"],
            source="wow.knmi.nl",
            unit=wow_df.columns[0].split(" ")[-1],
            meteo_var=meteo_var,
        )


class EvaporationObs(MeteoObs):
    """Class for evaporation timeseries.

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
        meteo_var="EV24",
        stn=None,
        fname=None,
        xy=None,
        start=None,
        end=None,
        fill_missing_obs=False,
        interval="daily",
        use_api=True,
        raise_exceptions=True,
        startdate=None,
        enddate=None,
    ):
        """Get an EvaporationObs timeseries from the KNMI evaporation in m.

        Parameters
        ----------
        meteo_var : str, optional
            meteo variable should be "EV24".
        stn : int, str or None, optional
            measurement station e.g. 829. The default is None.
        fname : str, path object, file-like object or None, optional
            filename of a knmi file. The default is None.
        xy : list, tuple or None, optional
            RD coördinates of a location in the Netherlands. The station nearest
            to this location used. The Default is None.
        start : str, datetime or None, optional
            start date of observations. The default is None.
        end : str, datetime or None, optional
            end date of observations. The default is None.
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False.
        interval : str, optional
            desired time interval for observations. Options are 'daily' and
            'hourly'. The default is 'daily'.
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

        return super().from_knmi(
            meteo_var,
            stn=stn,
            fname=fname,
            xy=xy,
            start=start,
            end=end,
            fill_missing_obs=fill_missing_obs,
            interval=interval,
            use_api=use_api,
            raise_exceptions=raise_exceptions,
            startdate=startdate,
            enddate=enddate,
        )


class PrecipitationObs(MeteoObs):
    """Class for precipitation timeseries.

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
    def from_knmi(
        cls,
        meteo_var="RH",
        stn=None,
        fname=None,
        xy=None,
        start=None,
        end=None,
        fill_missing_obs=False,
        interval="daily",
        use_api=True,
        raise_exceptions=True,
        startdate=None,
        enddate=None,
    ):
        """Get a PrecipitationObs timeseries from the KNMI precipitation. The
        precipitation is the Daily precipitation amount (in 0.1 mm) (-1 for.

        <0.05 mm).

        There are 3 different ways to obtain precipitation data from the knmi:
            1. Daily data from precipitation (neerslag) stations
            2. Daily data from meteo stations
            3. Hourly data from meteo stations

        1. For daily data from a precipitation station (neerslagstation) meteo_var
        should be 'RD'.
        2. For daily data from a meteo station meteo_var should be 'RH' and
        interval should be 'daily' (default).
        3. For hourly data from a meteo station meteo_var should be 'RH' and
        interval should be 'hourly'.

        More information about the differences between neerslag and meteo
        stations can be found in the hydropandas documentation ->
        02_knmi_observations notebook.

        Parameters
        ----------
        meteo_var : str, optional
            meteo variable can be "RH" or "RD". "RD" if you want data from a
            precipitation station (neerslagstation). "RH" if you want data from
            a meteo station. The default is "RH".
        stn : int, str or None, optional
            measurement station e.g. 829. The default is None.
        fname : str, path object, file-like object or None, optional
            filename of a knmi file. The default is None.
        xy : list, tuple or None, optional
            RD coördinates of a location in the Netherlands. The station nearest
            to this location used. The Default is None.
        start : str, datetime or None, optional
            start date of observations. The default is None.
        end : str, datetime or None, optional
            end date of observations. The default is None.
        fill_missing_obs : bool, optional
            if True nan values in time series are filled with nearby time series.
            The default is False.
        interval : str, optional
            desired time interval for observations. Options are 'daily' and
            'hourly'. The default is 'daily'.
        use_api : bool, optional
            if True the api is used to obtain the data, API documentation is here:
                https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
            if False a text file is downloaded into a temporary folder and the
            data is read from there. Default is True since the api is back
            online (July 2021).
        raise_exceptions : bool, optional
            if True you get errors when no data is returned. The default is False.

        Returns
        -------
        PrecipitationObs object with a precipitation time series and attributes
        """

        assert meteo_var in [
            "RD",
            "RH",
        ], 'meteo_var for precipitation should be "RD" or "RH"'

        return super().from_knmi(
            meteo_var,
            stn=stn,
            fname=fname,
            xy=xy,
            start=start,
            end=end,
            fill_missing_obs=fill_missing_obs,
            interval=interval,
            use_api=use_api,
            raise_exceptions=raise_exceptions,
            startdate=startdate,
            enddate=enddate,
        )

    @classmethod
    def from_wow(
        cls,
        stn: str = None,
        xy: List[float] = None,
        start: Optional[pd.Timestamp] = None,
        end: Optional[pd.Timestamp] = None,
    ):
        """Get a PrecipitationObs timeseries from a wow.knmi.nl station.

        Parameters
        ----------
        stn : Optional[int, str], optional
            station name
        xy : Optinal[List[float]], optinial
            longitude latitude of location [lon, lat] eg: [4.85, 51.95]
        start : Optional[pd.Timestamp], optional
            start date of observations, by default None
        end : Optional[pd.Timestamp], optional
            start date of observations, by default None

        Returns
        -------
        PrecipitationObs
        """
        return super().from_wow(
            meteo_var="rain_rate", stn=stn, xy=xy, start=start, end=end
        )
