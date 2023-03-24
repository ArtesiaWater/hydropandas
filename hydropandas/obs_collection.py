"""module with ObsCollection class for a collection of observations.

The ObsCollection class is a subclass of a pandas DataFrame with
additional attributes and methods.

More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas
"""

import warnings
import numbers
import numpy as np
import pandas as pd

from . import observation as obs
from . import util

import logging

logger = logging.getLogger(__name__)


def read_bro(
    extent=None,
    bro_id=None,
    name="",
    tmin=None,
    tmax=None,
    only_metadata=False,
    keep_all_obs=False,
    epsg=28992,
):
    """ get all the observations within an extent or within a
    groundwatermonitoring net.
    

    Parameters
    ----------
    extent : list, tuple, numpy-array or None, optional
        get groundwater monitoring wells within this extent
        [xmin, xmax, ymin, ymax]
    bro_id : str or None, optional
        starts with 'GMN'.
    name : str, optional
        name of the observation collection
    tmin : str or None, optional
        start time of observations. The default is None.
    tmax : str or None, optional
        end time of observations. The default is None.
    only_metadata : bool, optional
        if True download only metadata, significantly faster. The default
        is False.
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without
        measurements
    epsg : int, optional
        epsg code of the extent. The default is 28992 (RD).

    Returns
    -------
    ObsCollection
        ObsCollection DataFrame with the 'obs' column

    """

    oc = ObsCollection.from_bro(
        extent=extent,
        bro_id=bro_id,
        name=name,
        tmin=tmin,
        tmax=tmax,
        only_metadata=only_metadata,
        keep_all_obs=keep_all_obs,
        epsg=epsg,
    )

    return oc


def read_dino(
    dirname=None,
    extent=None,
    bbox=None,
    locations=None,
    ObsClass=obs.GroundwaterObs,
    subdir="Grondwaterstanden_Put",
    suffix="1.csv",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    name=None,
    **kwargs,
):
    """Read dino observations within an extent from the server or from a
    directory with downloaded files.

    Parameters
    ----------
    dirname : str, optional
        directory name, can be a .zip file or the parent directory
        of subdir
    extent : list, tuple or numpy-array (user must specify extent or bbox)
        get dinodata online within this extent [xmin, xmax, ymin, ymax]
    bbox : list, tuple or numpy-array (user must specify extent or bbox)
        The bounding box, in RD-coordinates, for which you want to
        retrieve locations [xmin, ymin, xmax, ymax]
    locations : list of str, optional
        list of names with location and filter number, separated by
        'filtersep'
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
        add all observation points to the collection, even the points
        without measurements or metadata
    name : str, optional
        the name of the observation collection
    kwargs:
        kwargs are passed to the io_dino.download_dino_within_extent() or
        the io_dino.read_dino_dir() function

    Returns
    -------
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_dino(
        dirname=dirname,
        ObsClass=ObsClass,
        subdir=subdir,
        suffix=suffix,
        unpackdir=unpackdir,
        force_unpack=force_unpack,
        preserve_datetime=preserve_datetime,
        keep_all_obs=keep_all_obs,
        name=name,
        **kwargs,
    )

    return oc


def read_fews(
    file_or_dir=None,
    xmlstring=None,
    ObsClass=obs.GroundwaterObs,
    name="fews",
    translate_dic=None,
    filterdict=None,
    locations=None,
    remove_nan=True,
    low_memory=True,
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    **kwargs,
):
    """Read one or several FEWS PI-XML files.

    Parameters
    ----------
    file_or_dir :  str
        zip, xml or directory with zips or xml files to read
    xmlstring : str or None
        string with xml data, only used if file_or_dir is None. Default is
        None
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    name : str, optional
        name of the observation collection, 'fews' by default
    translate_dic : dic or None, optional
        translate names from fews. If None this default dictionary is used:
        {'locationId': 'locatie'}.
    filterdict : dict, optional
        dictionary with tag name to apply filter to as keys, and list of
        accepted names as dictionary values to keep in final result,
        i.e. {"locationId": ["B001", "B002"]}
    locations : list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read. Only supported by
        low_memory=True method!
    low_memory : bool, optional
        whether to use xml-parsing method with lower memory footprint,
        default is True
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost, only used if low_memory=False
    unpackdir : str
        destination directory to unzip file if fname is a .zip
    force_unpack : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        whether to preserve datetime from zip archive

    Returns
    -------
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_fews_xml(
        file_or_dir=file_or_dir,
        xmlstring=xmlstring,
        ObsClass=ObsClass,
        name=name,
        translate_dic=translate_dic,
        filterdict=filterdict,
        locations=locations,
        remove_nan=remove_nan,
        low_memory=low_memory,
        unpackdir=unpackdir,
        force_unpack=force_unpack,
        preserve_datetime=preserve_datetime,
        **kwargs,
    )

    return oc


def read_imod(
    obs_collection,
    ml,
    runfile,
    mtime,
    model_ws,
    modelname="",
    nlay=None,
    exclude_layers=0,
):
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
        
    Returns
    -------
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_imod(
        obs_collection=obs_collection,
        ml=ml,
        runfile=runfile,
        mtime=mtime,
        model_ws=model_ws,
        modelname=modelname,
        nlay=nlay,
        exclude_layers=exclude_layers,
    )

    return oc


def read_knmi(
    locations=None,
    stns=None,
    xy=None,
    meteo_vars=("RH",),
    name="",
    starts=None,
    ends=None,
    ObsClasses=None,
    method="nearest",
    **kwargs,
):
    """Get knmi observations from a list of locations or a list of
    stations.

    Parameters
    ----------
    locations : pandas DataFrame or None
        dataframe with columns 'x' and 'y' as coordinates. The
        default is None
    stns : list of str or None
        list of knmi stations. The default is None
    xy : list or numpy array, optional
        xy coördinates of the locations. e.g. [[10,25], [5,25]]
    meteo_vars : list or tuple of str
        meteo variables e.g. ["RH", "EV24"]. The default is ("RH").
        See list of all possible variables below
    name : str, optional
        name of the obscollection. The default is ''
    starts : None, str, datetime or list, optional
        start date of observations per meteo variable. The start date is
        included in the time series.
        If start is None the start date will be January 1st of the
        previous year.
        If start is str it will be converted to datetime.
        If start is a list it should be the same length as meteo_vars and
        the start time for each variable. The default is None
    ends : list of str, datetime or None
        end date of observations per meteo variable. The end date is
        included in the time series.
        If end is None the start date will be January 1st of the
        previous year.
        If end is a str it will be converted to datetime.
        If end is a list it should be the same length as meteo_vars and
        the end time for each meteo variable. The default is None
    ObsClasses : list of type or None
        class of the observations, can be PrecipitationObs, EvaporationObs
        or MeteoObs. If None the type of observations is derived from the
        meteo_vars.
    method : str, optional
        specify whether EvaporationObs should be collected from the nearest
        meteo station (fast) or interpolated using thin plate spline (slow).
        Choiche betweeen 'nearest' or 'interpolation'
    **kwargs :
        kwargs are passed to the io_knmi.get_knmi_obslist function

    List of possible variables:
        neerslagstations:
        RD    = de 24-uurs neerslagsom, gemeten van 0800 utc op de
        voorafgaande dag tot 0800 utc op de vermelde datum meteostations:
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
        P     = Luchtdruk (in 0.1 hPa) herleid tot zeeniveau, op het moment
        van meten / Air pressure (in 0.1 hPa) reduced to mean sea level, at
        the time of observation 
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
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_knmi(
        locations=locations,
        stns=stns,
        xy=xy,
        meteo_vars=meteo_vars,
        name=name,
        starts=starts,
        ends=ends,
        ObsClasses=ObsClasses,
        method=method,
        **kwargs,
    )

    return oc


def read_menyanthes(
    fname, name="", ObsClass=obs.Obs, load_oseries=True, load_stresses=True
):
    """ read a Menyanthes file    

    Parameters
    ----------
    fname : str
        full path of the .men file.
    name : str, optional
        name of the observation collection. The default is "".
    ObsClass : type, optional
        class of the observations, e.g. GroundwaterObs. The default is
        obs.Obs.
    load_oseries : bool, optional
        if True the observations are read. The default is True.
    load_stresses : bool, optional
        if True the stresses are read. The default is True.

    Returns
    -------
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_menyanthes(
        fname=fname,
        name=name,
        ObsClass=ObsClass,
        load_oseries=load_oseries,
        load_stresses=load_stresses,
    )

    return oc


def read_modflow(
    obs_collection,
    ml,
    hds_arr,
    mtime,
    modelname="",
    nlay=None,
    exclude_layers=None,
    method="linear",
):
    """Read modflow groundwater heads at locations in obs_collection.

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
    exclude_layers : list of int, optional
        exclude the observations in these model layers
    method : str, optional
        interpolation method, either 'linear' or 'nearest',
        default is linear

    Returns
    -------
    ObsCollection
        collection of multiple point observations
    """

    oc = ObsCollection.from_modflow(
        obs_collection=obs_collection,
        ml=ml,
        hds_arr=hds_arr,
        mtime=mtime,
        modelname=modelname,
        nlay=nlay,
        exclude_layers=exclude_layers,
        method=method,
    )

    return oc


def read_waterinfo(
    file_or_dir, name="", ObsClass=obs.WaterlvlObs, progressbar=True, **kwargs
):
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
    oc = ObsCollection.from_waterinfo(
        file_or_dir=file_or_dir,
        name=name,
        ObsClass=ObsClass,
        progressbar=progressbar,
        **kwargs,
    )

    return oc


def read_wiski(
    dirname,
    ObsClass=obs.GroundwaterObs,
    suffix=".csv",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    **kwargs,
):
    """

    Parameters
    ----------
    dirname : str
        path of the zipfile with wiski data.
    ObsClass : type, optional
        type of Obs. The default is obs.GroundwaterObs.
    suffix : str, optional
        extension of filenames to read. The default is ".csv".
    unpackdir : str or None, optional
        directory to unpack zipped directory. The default is None.
    force_unpack : bool, optional
        force unzip, by default False.
    preserve_datetime : bool, optional
        preserve datetime of unzipped files, by default False
        (useful for checking whether data has changed)
    keep_all_obs : bool, optional
        If True keep all observations even those without metadata. The default
        is True.
    **kwargs

    Returns
    -------
    ObsCollection
        ObsCollection containing observation data
    """

    oc = ObsCollection.from_wiski(
        dirname=dirname,
        ObsClass=ObsClass,
        suffix=suffix,
        unpackdir=unpackdir,
        force_unpack=force_unpack,
        preserve_datetime=preserve_datetime,
        keep_all_obs=keep_all_obs,
        **kwargs,
    )

    return oc


class ObsCollection(pd.DataFrame):
    """class for a collection of point observations.

    An ObsCollection object is a subclass of a pandas.DataFrame and allows for
    additional attributes and methods. Additional attributes are
    defined in the '_metadata' attribute.

    Parameters
    ----------
    name : str
        name of the observation collection
    meta : dic
        metadata of the observation collection
    """

    # temporary properties
    _internal_names = pd.DataFrame._internal_names + ["none"]
    _internal_names_set = set(_internal_names)

    # normal properties
    _metadata = [
        "name",
        "meta",
    ]

    def __init__(self, *args, **kwargs):
        """constructor of the ObsCollection.

        *args must be input for the pandas.DataFrame constructor,
        **kwargs can be one of the attributes listed in _metadata or
        keyword arguments for the constructor of a pandas.DataFrame.
        """
        self.name = kwargs.pop("name", "")
        self.meta = kwargs.pop("meta", {})

        super(ObsCollection, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return ObsCollection

    def _infer_otype(self):
        """Infer observation type from the obs column.

        Parameters
        ----------

        Returns
        -------
        otype, otypes
            type or list of types of the observation objects
        """
        otypes = self.obs.apply(lambda x: type(x)).unique()
        if otypes.shape[0] == 1:
            logger.info("inferred observation type: {}".format(otypes[0]))
            return otypes[0]
        elif otypes.shape[0] > 1:
            logger.info("inferred multiple otypes, types: {}".format(otypes))
            return otypes
        else:
            raise TypeError("could not infer observation type")

    def _set_metadata_value(self, iname, att_name, value, add_to_meta=False):
        """Set a value on three different levels at once:
            1. the value in an ObsCollection DataFrame
            2. the attribute of the observation
            3. the value in the meta dictionary of an observation (optional)

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

        o = self.loc[iname, "obs"]
        if att_name in o._metadata:
            setattr(o, att_name, value)
            logger.debug(f"set attribute {att_name} of {iname} to {value}")

        if att_name == "name":
            # name is the index of the ObsCollection dataframe
            self.rename(index={iname: value}, inplace=True)
        else:
            self.loc[iname, att_name] = value
        logger.debug(f"set {iname}, {att_name} to {value} in obscollection")

        if add_to_meta:
            o.meta.update({att_name: value})
            logger.debug(f"add {att_name} of {iname} with value {value} to meta")

    def _is_consistent(self, check_individual_obs=True):
        """check if an observation collection is consistent. An observation
        collection is consistent if:
            1. all observations have a unique name
            2. there are no nan values in the obs column
            3. (optional) the metadata of each observation has the same type
            and value as the corresponding row in the observation collection
            dataframe. Only checked if check_individual_obs is True.


        Parameters
        ----------
        check_individual_obs : bool, optional
            If True the third condition in the list above is checked. The
            default is True.

        Returns
        -------
        bool
            True -> consistent
            False -> inconsistent.

        """
        # check unique index
        if not self.index.is_unique:
            logger.warning(f"index of observation collection -> {self.name} not unique. non unique indices are:")
            logger.warning(" ".join(self.index[self.index.duplicated()]))
            return False

        # check nan values in observations
        if self.obs.isna().any():
            logger.warning(f"missing observation object in collection -> {self.name} ")
            return False

        # check oc data with individual object attributes
        if check_individual_obs:
            for o in self.obs.values:
                for att in o._metadata:
                    if att not in ["name", "meta"]:
                        v1 = self.loc[o.name, att]
                        v2 = getattr(o, att)
                        # check if values are equal
                        try:
                            if v1 != v2:
                                # check if both are nan
                                if isinstance(v1, numbers.Number) and isinstance(
                                    v2, numbers.Number
                                ):
                                    if np.isnan(v1) and np.isnan(v2):
                                        continue

                                # otherwise return Nan
                                logger.warning(
                                    f"observation collection -> {self.name} not consistent with observation -> {o.name} {att} value"
                                )
                                return False
                        except TypeError:
                            logger.warning(
                                f"observation collection -> {self.name} not consistent with observation -> {o.name} {att} value"
                            )
                            return False
                    elif att == "name":
                        if o.name not in self.index:
                            logger.warning(
                                f"observation collection -> {self.name} not consistent with observation -> {o.name} name"
                            )
                            return False

        return True

    def add_observation(self, o, check_consistency=True, **kwargs):
        """add an observation to an existing observation collection. If the
        observation exists the two observations are merged.

        Parameters
        ----------
        o : hpd.observation.Obs
            Observation object.
        check_consistency : bool, optional
            If True the consistency of the collection is first checked. The
            default is True.
        **kwargs passed to Obs.merge_observation:
            merge_metadata : bool, optional
                If True and observations are merged the metadata of the two
                objects are merged. If there are any differences the overlap
                parameter is used to determine which metadata is used. If
                merge_metadata is False, the metadata of the original
                observation is always used for the merged observation. The
                default is True.
            overlap : str, optional
                How to deal with overlapping timeseries with different values.
                Options are:
                - error : Raise a ValueError
                - use_left : use the overlapping part from the existing
                observations
                - use_right : use the overlapping part from the new observation
                Default is 'error'.

        Raises
        ------
        RuntimeError
            when the observation collection is inconsistent.
        TypeError
            when the observation type is wrong.

        Returns
        -------
        None.

        """
        if check_consistency:
            if not self._is_consistent():
                raise RuntimeError("inconsistent observation collection")

        if not isinstance(o, obs.Obs):
            raise TypeError("Observation should be of type hydropandas.observation.Obs")

        # add new observation to collection
        if o.name not in self.index:
            logger.info(f"adding {o.name} to collection")
            self.loc[o.name] = o.to_collection_dict()
        else:
            logger.info(
                f"observation name {o.name} already in collection, merging observations"
            )

            o1 = self.loc[o.name, "obs"]
            omerged = o1.merge_observation(o, **kwargs)

            # overwrite observation in collection
            self.loc[o.name] = omerged.to_collection_dict()

    def add_obs_collection(
        self, obs_collection, check_consistency=True, inplace=False, **kwargs
    ):
        """add one observation collection to another observation
        collection. See add_observation method for more details

        Parameters
        ----------
        obs_collection : hpd.ObsCollection
            ObsCollection object.
        check_consistency : bool, optional
            If True the consistency of both collections is first checked. The
            default is True.
        inplace : bool, optional
            If True, modifies the ObsCollection in place (do not create a new
            object). The default is False.
        **kwargs passed to Obs.merge_observation:
            merge_metadata : bool, optional
                If True and observations are merged the metadata of the two
                objects are merged. If there are any differences the overlap
                parameter is used to determine which metadata is used. If
                merge_metadata is False, the metadata of the original
                observation is always used for the merged observation. The
                default is True.
            overlap : str, optional
                How to deal with overlapping timeseries with different values.
                Options are:
                - error : Raise a ValueError
                - use_left : use the overlapping part from the existing
                observations
                - use_right : use the overlapping part from the new observation
                Default is 'error'.

        Raises
        ------
        RuntimeError
            when the observation collection is inconsistent.

        Returns
        -------
        ObsCollection or None
            merged ObsCollection if ``inplace=True``.

        """
        if check_consistency:
            if not self._is_consistent():
                raise RuntimeError(
                    f"inconsistent observation collection -> {self.name}"
                )

            if not obs_collection._is_consistent():
                raise RuntimeError(
                    f"inconsistent observation collection -> {obs_collection.name}"
                )

        if inplace:
            for o in obs_collection.obs.values:
                self.add_observation(o, check_consistency=False, **kwargs)

        else:
            oc = self.copy()
            for o in obs_collection.obs.values:
                oc.add_observation(o, check_consistency=False, **kwargs)

            return oc

    @classmethod
    def from_bro(
        cls,
        extent=None,
        bro_id=None,
        name="",
        tmin=None,
        tmax=None,
        only_metadata=False,
        keep_all_obs=True,
        epsg=28992,
        ignore_max_obs=False,
    ):
        """ get all the observations within an extent or within a
        groundwatermonitoring net.
        

        Parameters
        ----------
        extent : list, tuple, numpy-array or None, optional
            get groundwater monitoring wells within this extent
            [xmin, xmax, ymin, ymax]
        bro_id : str or None, optional
            starts with 'GMN'.
        name : str, optional
            name of the observation collection
        tmin : str or None, optional
            start time of observations. The default is None.
        tmax : str or None, optional
            end time of observations. The default is None.
        only_metadata : bool, optional
            if True download only metadata, significantly faster. The default
            is False.
        keep_all_obs : boolean, optional
            add all observation points to the collection, even without
            measurements
        epsg : int, optional
            epsg code of the extent. The default is 28992 (RD).
        ignore_max_obs : bool, optional
            by default you get a prompt if you want to download over a 1000
            observations at once. if ignore_max_obs is True you won't get the
            prompt. The default is False

        Returns
        -------
        ObsCollection
            ObsCollection DataFrame with the 'obs' column

        """

        from .io.io_bro import get_obs_list_from_gmn, get_obs_list_from_extent

        if bro_id is None and (extent is not None):
            obs_list = get_obs_list_from_extent(
                extent,
                obs.GroundwaterObs,
                tmin=tmin,
                tmax=tmax,
                only_metadata=only_metadata,
                keep_all_obs=keep_all_obs,
                epsg=epsg,
                ignore_max_obs=ignore_max_obs,
            )
            meta = {}
        elif bro_id is not None:
            obs_list, meta = get_obs_list_from_gmn(
                bro_id,
                obs.GroundwaterObs,
                only_metadata=only_metadata,
                keep_all_obs=keep_all_obs,
            )
            name = meta.pop("name")
        else:
            raise ValueError("specify bro_id or extent")

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_dataframe(cls, df, obs_list=None, ObsClass=obs.GroundwaterObs):
        """Create an observation collection from a DataFrame by adding a column
        with empty observations.

        Parameters
        ----------
        df : pandas DataFrame
            input dataframe. If this dataframe has a column named 'obs' the
            column is replaced with empty observation objects.
        obs_list : list of observation.Obs, optional
            list of observations. Default is None
        ObsClass : class, optional
            observation class used to create empty obs object, by
            default obs.GroundwaterObs

        Returns
        -------
        ObsCollection
            ObsCollection DataFrame with the 'obs' column
        """
        meta = {"type": obs.GroundwaterObs}
        if isinstance(df, pd.DataFrame):
            if obs_list is None:
                obs_list = [ObsClass() for i in range(len(df))]
            df["obs"] = obs_list
        else:
            raise TypeError(f"df should be type pandas.DataFrame not {type(df)}")

        return cls(df, meta=meta)

    @classmethod
    def from_dino(
        cls,
        dirname=None,
        ObsClass=obs.GroundwaterObs,
        subdir="Grondwaterstanden_Put",
        suffix="1.csv",
        unpackdir=None,
        force_unpack=False,
        preserve_datetime=False,
        keep_all_obs=True,
        name=None,
        **kwargs,
    ):
        """Read dino data within an extent from the server or from a directory
        with downloaded files.

        Parameters
        ----------
        dirname : str, optional
            directory name, can be a .zip file or the parent directory
            of subdir
        extent : list, tuple or numpy-array (user must specify extent or bbox)
            get dinodata online within this extent [xmin, xmax, ymin, ymax]
        bbox : list, tuple or numpy-array (user must specify extent or bbox)
            The bounding box, in RD-coordinates, for which you want to
            retrieve locations [xmin, ymin, xmax, ymax]
        locations : list of str, optional
            list of names with location and filter number, separated by
            'filtersep'
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
            add all observation points to the collection, even the points
            without measurements or metadata
        name : str, optional
            the name of the observation collection
        kwargs:
            kwargs are passed to the io_dino.read_dino_dir() function

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """
        from .io.io_dino import read_dino_dir

        # read dino directory
        if name is None:
            name = subdir

        meta = {
            "dirname": dirname,
            "type": ObsClass,
            "suffix": suffix,
            "unpackdir": unpackdir,
            "force_unpack": force_unpack,
            "preserve_datetime": preserve_datetime,
            "keep_all_obs": keep_all_obs,
        }

        obs_list = read_dino_dir(
            dirname,
            ObsClass,
            subdir,
            suffix,
            unpackdir,
            force_unpack,
            preserve_datetime,
            keep_all_obs,
            **kwargs,
        )

        obs_df = util._obslist_to_frame(obs_list)
        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_artdino_dir(
        cls,
        dirname=None,
        ObsClass=obs.GroundwaterObs,
        subdir="csv",
        suffix=".csv",
        unpackdir=None,
        force_unpack=False,
        preserve_datetime=False,
        keep_all_obs=True,
        name=None,
        **kwargs,
    ):
        """Read a dino directory.

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

        meta = {
            "dirname": dirname,
            "type": ObsClass,
            "suffix": suffix,
            "unpackdir": unpackdir,
            "force_unpack": force_unpack,
            "preserve_datetime": preserve_datetime,
            "keep_all_obs": keep_all_obs,
        }

        obs_list = read_artdino_dir(
            dirname,
            ObsClass,
            subdir,
            suffix,
            unpackdir,
            force_unpack,
            preserve_datetime,
            keep_all_obs,
            **kwargs,
        )

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_fews_xml(
        cls,
        file_or_dir=None,
        xmlstring=None,
        ObsClass=obs.GroundwaterObs,
        name="fews",
        translate_dic=None,
        filterdict=None,
        locations=None,
        remove_nan=True,
        low_memory=True,
        unpackdir=None,
        force_unpack=False,
        preserve_datetime=False,
        **kwargs,
    ):
        """Read one or several FEWS PI-XML files.

        Parameters
        ----------
        file_or_dir :  str
            zip, xml or directory with zips or xml files to read
        xmlstring : str or None
            string with xml data, only used if file_or_dir is None. Default is
            None
        ObsClass : type
            class of the observations, e.g. GroundwaterObs or WaterlvlObs
        name : str, optional
            name of the observation collection, 'fews' by default
        translate_dic : dic or None, optional
            translate names from fews. If None this default dictionary is used:
            {'locationId': 'locatie'}.
        filterdict : dict, optional
            dictionary with tag name to apply filter to as keys, and list of
            accepted names as dictionary values to keep in final result,
            i.e. {"locationId": ["B001", "B002"]}
        locations : list of str, optional
            list of locationId's to read from XML file, others are skipped.
            If None (default) all locations are read. Only supported by
            low_memory=True method!
        low_memory : bool, optional
            whether to use xml-parsing method with lower memory footprint,
            default is True
        remove_nan : boolean, optional
            remove nan values from measurements, flag information about the
            nan values is also lost, only used if low_memory=False
        unpackdir : str
            destination directory to unzip file if fname is a .zip
        force_unpack : boolean, optional
            force unpack if dst already exists
        preserve_datetime : boolean, optional
            whether to preserve datetime from zip archive

        Returns
        -------
        cls(obs_df) : ObsCollection
            collection of multiple point observations
        """
        from .io.io_fews import read_xml_filelist, read_xmlstring

        if translate_dic is None:
            translate_dic = {"locationId": "monitoring_well"}

        meta = {"type": ObsClass}

        if file_or_dir is not None:
            # get files
            dirname, unzip_fnames = util.get_files(
                file_or_dir,
                ext=".xml",
                unpackdir=unpackdir,
                force_unpack=force_unpack,
                preserve_datetime=preserve_datetime,
            )
            meta.update({"filename": dirname})

            obs_list = read_xml_filelist(
                unzip_fnames,
                ObsClass,
                directory=dirname,
                translate_dic=translate_dic,
                filterdict=filterdict,
                locations=locations,
                remove_nan=remove_nan,
                low_memory=low_memory,
                **kwargs,
            )

            obs_df = util._obslist_to_frame(obs_list)
            return cls(obs_df, name=name, meta=meta)

        elif (file_or_dir is None) and (xmlstring is not None):
            obs_list = read_xmlstring(
                xmlstring,
                ObsClass,
                translate_dic=translate_dic,
                filterdict=filterdict,
                locationIds=locations,
                low_memory=low_memory,
                remove_nan=remove_nan,
                **kwargs,
            )
            obs_df = util._obslist_to_frame(obs_list)
            return cls(obs_df, name=name, meta=meta)

        else:
            raise ValueError("either specify variables file_or_dir or xmlstring")

    @classmethod
    def from_imod(
        cls,
        obs_collection,
        ml,
        runfile,
        mtime,
        model_ws,
        modelname="",
        nlay=None,
        exclude_layers=0,
    ):
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
        """
        from .io.io_modflow import read_imod_results

        mo_list = read_imod_results(
            obs_collection,
            ml,
            runfile,
            mtime,
            model_ws,
            modelname=modelname,
            nlay=nlay,
            exclude_layers=exclude_layers,
        )
        obs_df = util._obslist_to_frame(mo_list)
        return cls(obs_df, name=modelname)

    @classmethod
    def from_knmi(
        cls,
        locations=None,
        stns=None,
        xy=None,
        meteo_vars=("RH",),
        name="",
        starts=None,
        ends=None,
        ObsClasses=None,
        method="nearest",
        **kwargs,
    ):
        """Get knmi observations from a list of locations or a list of
        stations.

        Parameters
        ----------
        locations : pandas DataFrame or None
            dataframe with columns 'x' and 'y' as coordinates. The
            default is None
        stns : list of str or None
            list of knmi stations. The default is None
        xy : list or numpy array, optional
            xy coördinates of the locations. e.g. [[10,25], [5,25]]
        meteo_vars : list or tuple of str
            meteo variables e.g. ["RH", "EV24"]. The default is ("RH").
            See list of all possible variables below
        name : str, optional
            name of the obscollection. The default is ''
        starts : None, str, datetime or list, optional
            start date of observations per meteo variable. The start date is
            included in the time series.
            If start is None the start date will be January 1st of the
            previous year.
            If start is str it will be converted to datetime.
            If start is a list it should be the same length as meteo_vars and
            the start time for each variable. The default is None
        ends : list of str, datetime or None
            end date of observations per meteo variable. The end date is
            included in the time series.
            If end is None the start date will be January 1st of the
            previous year.
            If end is a str it will be converted to datetime.
            If end is a list it should be the same length as meteo_vars and
            the end time for each meteo variable. The default is None
        ObsClasses : list of type or None
            class of the observations, can be PrecipitationObs, EvaporationObs
            or MeteoObs. If None the type of observations is derived from the
            meteo_vars.
        method : str, optional
            specify whether EvaporationObs should be collected from the nearest
            meteo station (fast) or interpolated using thin plate spline (slow).
            Choiche betweeen 'nearest' or 'interpolation'
        **kwargs :
            kwargs are passed to the io_knmi.get_knmi_obslist function

        List of possible variables:
            neerslagstations:
            RD    = de 24-uurs neerslagsom, gemeten van 0800 utc op de
            voorafgaande dag tot 0800 utc op de vermelde datum meteostations:
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
        """

        from .io.io_knmi import get_knmi_obslist

        # obtain ObsClass
        if ObsClasses is None:
            ObsClasses = []
            for meteovar in meteo_vars:
                if meteovar in ("RH", "RD"):
                    ObsClasses.append(obs.PrecipitationObs)
                elif meteovar == "EV24":
                    ObsClasses.append(obs.EvaporationObs)
                else:
                    ObsClasses.append(obs.MeteoObs)

        elif isinstance(ObsClasses, type):
            if issubclass(
                ObsClasses, (obs.PrecipitationObs, obs.EvaporationObs, obs.MeteoObs)
            ):
                ObsClasses = [ObsClasses] * len(meteo_vars)
            else:
                TypeError(
                    "must be None, PrecipitationObs, EvaporationObs, MeteoObs, list or tuple"
                )
        elif isinstance(ObsClasses, (list, tuple)):
            pass
        else:
            TypeError(
                "must be None, PrecipitationObs, EvaporationObs, MeteoObs, list or tuple"
            )

        meta = {}
        meta["starts"] = starts
        meta["ends"] = ends
        meta["name"] = name
        meta["ObsClasses"] = ObsClasses
        meta["meteo_vars"] = meteo_vars

        obs_list = get_knmi_obslist(
            locations,
            stns,
            xy,
            meteo_vars,
            ObsClasses=ObsClasses,
            starts=starts,
            ends=ends,
            method=method,
            **kwargs,
        )

        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_list(cls, obs_list, name=""):
        """read observations from a list of obs objects.

        Parameters
        ----------
        obs_list : list of observation.Obs
            list of observations
        name : str, optional
            name of the observation collection
        """
        obs_df = util._obslist_to_frame(obs_list)
        return cls(obs_df, name=name)

    @classmethod
    def from_menyanthes(
        cls, fname, name="", ObsClass=obs.Obs, load_oseries=True, load_stresses=True
    ):

        from .io.io_menyanthes import read_file

        menyanthes_meta = {"filename": fname, "type": ObsClass}

        obs_list = read_file(
            fname, ObsClass, load_oseries=load_oseries, load_stresses=load_stresses
        )
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, meta=menyanthes_meta, name=name)

    @classmethod
    def from_modflow(
        cls,
        obs_collection,
        ml,
        hds_arr,
        mtime,
        modelname="",
        nlay=None,
        exclude_layers=None,
        method="linear",
    ):
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
        exclude_layers : list of int, optional
            exclude the observations in these model layers
        method : str, optional
            interpolation method, either 'linear' or 'nearest',
            default is linear
        """
        from .io.io_modflow import read_modflow_results

        mo_list = read_modflow_results(
            obs_collection,
            ml,
            hds_arr,
            mtime,
            modelname=modelname,
            nlay=nlay,
            method=method,
            exclude_layers=exclude_layers,
        )
        obs_df = util._obslist_to_frame(mo_list)

        return cls(obs_df)

    @classmethod
    def from_waterinfo(
        cls, file_or_dir, name="", ObsClass=obs.WaterlvlObs, progressbar=True, **kwargs
    ):
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

        meta = {"name": name, "type": ObsClass, "filename": file_or_dir}

        obs_list = io_waterinfo.read_waterinfo_obs(
            file_or_dir, ObsClass, progressbar=progressbar, **kwargs
        )
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    @classmethod
    def from_wiski(
        cls,
        dirname,
        ObsClass=obs.GroundwaterObs,
        suffix=".csv",
        unpackdir=None,
        force_unpack=False,
        preserve_datetime=False,
        keep_all_obs=True,
        **kwargs,
    ):

        from .io.io_wiski import read_wiski_dir

        meta = {
            "dirname": dirname,
            "type": ObsClass,
            "suffix": suffix,
            "unpackdir": unpackdir,
            "force_unpack": force_unpack,
            "preserver_datetime": preserve_datetime,
            "keep_all_obs": keep_all_obs,
        }

        name = "wiski_import"
        obs_list = read_wiski_dir(
            dirname,
            ObsClass=ObsClass,
            suffix=suffix,
            unpackdir=unpackdir,
            force_unpack=force_unpack,
            preserve_datetime=preserve_datetime,
            keep_all_obs=keep_all_obs,
            **kwargs,
        )
        obs_df = util._obslist_to_frame(obs_list)

        return cls(obs_df, name=name, meta=meta)

    def to_excel(self, path, main_sheet_name=None):
        """write an ObsCollection to an excel, the first sheet in the
        excel contains the metadata, the other tabs are the timeseries of each
        observation.


        Parameters
        ----------
        path : str
            full path of xlsx file.
        main_sheet_name : str or None, optional
            sheetname with metadata, if None the name of the ObsCollection is
            used. The default is None.

        Raises
        ------
        RuntimeError
            If the ObsCollection is inconsistent.

        Returns
        -------
        None.

        """

        if main_sheet_name is None:
            main_sheet_name = self.name

        if not self._is_consistent():
            raise RuntimeError("inconsistent observation collection")

        with pd.ExcelWriter(path) as writer:
            # write ObsCollection dataframe without observations to first sheet
            super(ObsCollection, self.drop(columns="obs")).to_excel(
                writer, main_sheet_name
            )

            # write each observation time series to next sheets
            for o in self.obs.values:
                sheetname = o.name
                for ch in ["[", "]", ":", "*", "?", "/", "\\"]:
                    sheetname = sheetname.replace(ch, "_")
                o.to_excel(writer, sheetname)

    def to_pi_xml(self, fname, timezone="", version="1.24"):
        from .io import io_fews

        io_fews.write_pi_xml(self, fname, timezone=timezone, version=version)

    def to_gdf(self, xcol="x", ycol="y"):
        """convert ObsCollection to GeoDataFrame.

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

    def to_report_table(
        self,
        columns=(
            "monitoring_well",
            "tube_nr",
            "startdate",
            "enddate",
            "# measurements",
        ),
    ):

        if "startdate" in columns:
            self["startdate"] = self.obs.apply(lambda x: x.index[0])
        if "enddate" in columns:
            self["enddate"] = self.obs.apply(lambda x: x.index[-1])
        if "# measurements" in columns:
            self["# measurements"] = self.obs.apply(lambda x: x.shape[0])

        return self[columns]

    def to_pastastore(
        self,
        pstore=None,
        pstore_name="",
        col=None,
        kind="oseries",
        add_metadata=True,
        conn=None,
        overwrite=False,
    ):
        """add observations to a new or existing pastastore.

        Parameters
        ----------
        pstore : pastastore.PastaStore, optional
            Existing pastastore, if None a new pastastore is created
        pstore_name : str, optional
            Name of the pastastore only used if pstore is None
        col : str, optional
            Name of the column in the Obs dataframe to be used. If None the
            first numeric column in the Obs Dataframe is used.
        kind : str, optional
            The kind of series that is added to the pastastore
        add_metadata : boolean, optional
            If True metadata from the observations added to the pastastore
        conn : pastastore.connectors or None, optional
            type of connector, if None the DictConnector is used. Default is
            None.
        overwrite : boolean, optional
            if True, overwrite existing series in pastastore, default is False

        Returns
        -------
        pstore : pastastore.PastaStore
            the pastastore with the series from the ObsCollection
        """
        from .io.io_pastas import create_pastastore

        pstore = create_pastastore(
            self,
            pstore,
            pstore_name,
            add_metadata=add_metadata,
            kind=kind,
            col=col,
            conn=conn,
            overwrite=overwrite,
        )

        return pstore

    def to_shapefile(self, fname, xcol="x", ycol="y"):
        """save ObsCollection as shapefile.

        Parameters
        ----------
        fname : str
            filename of shapefile, ends with .shp
        xcol : str
            column name with x values
        ycol : str
            column name with y values
        """
        from geopandas.array import GeometryDtype

        gdf = util.df2gdf(self, xcol, ycol)

        # remove obs column
        if "obs" in gdf.columns:
            gdf.drop(columns="obs", inplace=True)

        # change dtypes that are not accepted for shapefiles
        for colname, coltype in gdf.dtypes.items():
            # ommit geometry dtype
            if isinstance(coltype, GeometryDtype):
                pass
            # cast boolean columns to int
            elif coltype == bool:
                gdf[colname] = gdf[colname].astype(int)
            # cast datetime columns to str
            elif np.issubdtype(coltype, np.datetime64):
                gdf[colname] = gdf[colname].astype(str)

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

        self[key] = [
            o.meta[key] if key in o.meta.keys() else None for o in self.obs.values
        ]

    def get_series(self, tmin=None, tmax=None, col=None):
        """        

        Parameters
        ----------
        tmin : datetime, optional
            start time for series. The default is None.
        tmax : datetime, optional
            end time for series. The default is None.
        col : str or None, optional
            the column of the obs dataframe to get measurements from. The
            first numeric column is used if col is None, by default None.

        Returns
        -------
        series of Series
            series of a series of observations within a time frame.

        """

        if tmin is None:
            tmin = self.stats.dates_first_obs.min()
        if tmax is None:
            tmax = self.stats.dates_last_obs.max()

        def get_s(o, tmin=tmin, tmax=tmax, col=col):
            if col is None:
                col = o._get_first_numeric_col_name()
            return o.loc[tmin:tmax, col]

        return self.obs.apply(lambda o: o.loc[tmin:tmax, col])
