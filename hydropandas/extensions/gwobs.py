import logging

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = None

from . import accessor

logger = logging.getLogger(__name__)


def get_model_layer_z(z, zvec, left=-999, right=999):
    """Get index of model layer based on elevation.

    Assumptions:

    - the highest value in zvec is the top of model layer 0.
    - if z is equal to the bottom of a layer, the model layer above that
      layer is assigned.

    Parameters
    ----------
    z : int, float
        elevation.
    zvec : list, np.array
        elevations of model layers. shape is nlay + 1
    left : int, optional
        if z is below the lowest value in zvec, this value is returned.
        The default is -999.
    right : TYPE, optional
        if z is above the highest value in zvec, this value is returned.
        The default is 999.

    Returns
    -------
    int : int
        model layer

    Examples
    --------
    >>> zvec = [0, -10, -20, -30]
    >>> get_model_layer_z(-5, zvec)
    0

    >>> get_model_layer_z(-25, zvec)
    2

    >>> get_model_layer_z(-50, zvec)
    -999

    >>> get_model_layer_z(100, zvec)
    999

    >>> get_model_layer_z(-20, zvec)
    1
    """
    # make sure zvec is descending
    zvec = np.sort(zvec)[::-1]
    if z in zvec:
        z += 1e-10
    lay = int(
        np.interp(z, zvec[::-1], np.arange(len(zvec))[::-1], left=left, right=right)
    )

    return lay


def check_if_var_is_invalid(var):
    if var is None:
        return True
    elif np.isnan(var):
        return True

    return False


def get_modellayer_from_screen_depth(ftop, fbot, zvec, left=-999, right=999):
    """


    Parameters
    ----------
    ftop : int or float
        top of screen.
    fbot : int or float
        bottom of screen, has to be lower than ftop.
    zvec : list or np.array
        elevations of the modellayers at the location of the tube.
    left : int, optional
        value to return if tube screen is below the modellayers.
        The default is -999.
    right : int, optional
        value to return if tube screen is above the modellayers.
        The default is-999.

    Raises
    ------
    ValueError
        raised if something unexpected happens.

    Returns
    -------
    int or np.nan
        modellayer.

    examples
    --------
    >>> zvec = [0, -10, -20, -30, -40]
    >>> get_modellayer_from_screen_depth(-5, -7, zvec)
    0

    >>> get_modellayer_from_screen_depth(-25, -27, zvec)
    2

    >>> get_modellayer_from_screen_depth(-15, -27, zvec)
    2

    >>> get_modellayer_from_screen_depth(-5, -27, zvec)
    1

    >>> get_modellayer_from_screen_depth(-5, -37, zvec)
    1

    >>> get_modellayer_from_screen_depth(15, -5, zvec)
    0

    >>> get_modellayer_from_screen_depth(15, 5, zvec)
    999

    >>> get_modellayer_from_screen_depth(-55, -65, zvec)
    -999

    >>> get_modellayer_from_screen_depth(15, -65, zvec)
    nan

    >>> get_modellayer_from_screen_depth(None, -7, zvec)
    0

    >>> get_modellayer_from_screen_depth(None, None, zvec)
    nan

    """
    if isinstance(zvec, np.ndarray):
        if np.isnan(zvec).all():
            logger.warning("vertical datum invalid returning nan")
            return np.nan

    zvec = np.sort(zvec)[::-1]
    ftop_invalid = check_if_var_is_invalid(ftop)
    fbot_invalid = check_if_var_is_invalid(fbot)
    if ftop_invalid and fbot_invalid:
        logger.error("- values for ftop and fbot are invalid!")
        return np.nan
    elif fbot_invalid:
        lay_ftop = get_model_layer_z(ftop, zvec, left=left, right=right)
        logger.error(f"- fbot invalid. selected based on ftop: {lay_ftop}")
        return lay_ftop
    elif ftop_invalid:
        lay_fbot = get_model_layer_z(fbot, zvec, left=left, right=right)
        logger.error(f"- ftop invalid. selected based on fbot: {lay_fbot}")
        return lay_fbot

    if ftop < fbot:
        logger.warning("- tube screen top below tube screen bot, switching top and bot")
        fbot, ftop = ftop, fbot

    lay_ftop = get_model_layer_z(ftop, zvec, left=left, right=right)
    lay_fbot = get_model_layer_z(fbot, zvec, left=left, right=right)

    # Piezometer in layer in which majority of screen is located
    if lay_fbot == lay_ftop:
        logger.debug(f"- selected layer {lay_fbot}")
        return lay_fbot

    else:
        if lay_fbot == left and lay_ftop == right:
            logger.debug("- tube screen spans all layers. " "return nan")
            return np.nan
        elif lay_ftop == right:
            logger.debug(
                "- tube screen top higher than top layer. " f"selected layer {lay_fbot}"
            )
            return lay_fbot

        elif lay_fbot == left:
            logger.debug(
                "- tube screen bot lower than bottom layer. "
                f"selected layer {lay_ftop}"
            )
            return lay_ftop

        logger.debug(
            "- tube screen crosses layer boundary:\n"
            f"  - layers: {lay_ftop}, {lay_fbot}"
        )

        logger.debug(
            f"- tube screen spans {lay_fbot - lay_ftop +1} layers."
            " checking length per layer\n"
            "  - length per layer:"
        )

        # check which layer has the biggest length of the tube screen
        length_layers = np.zeros(lay_fbot - lay_ftop + 1)
        for i in range(len(length_layers)):
            if i == 0:
                length_layers[i] = ftop - zvec[lay_ftop + 1]
            elif (i + 1) == len(length_layers):
                length_layers[i] = zvec[lay_ftop + i] - fbot
            else:
                length_layers[i] = zvec[lay_ftop + i] - zvec[lay_ftop + 1 + i]

            logger.debug(f"    - lay {lay_ftop+i}: {length_layers[i]:.2f}")

        # choose layer with biggest length
        rel_layer = np.argmax(length_layers)
        lay_out = lay_ftop + rel_layer
        logger.debug(f"  - selected layer: {lay_out}")
        return lay_out


def get_zvec(x, y, gwf=None, ds=None):
    """get a list with the vertical layer boundaries at a point in the model.

    Parameters
    ----------
    x : int or float
        x coordinate.
    y : int or float
        y coordinate.
    gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
        modflow model with top and bottoms
    ds : xarray.Dataset
        xarray Dataset typically created in nlmod. Must have
        variables 'top' and 'botm'.

    Raises
    ------
    NotImplementedError
        not all grid types are supported yet.

    Returns
    -------
    zvec : list
        list of vertical layer boundaries. length is nlay + 1.
    """
    import flopy
    from shapely.geometry import Point

    if gwf and not ds:
        ix = flopy.utils.GridIntersect(gwf.modelgrid)
        if gwf.modelgrid.grid_type == "structured":
            res = ix.intersect(Point(x, y))
            if len(res) > 0:
                r, c = res["cellids"][0]
                zvec = np.array(
                    [gwf.modelgrid.top[r, c]]
                    + [gwf.modelgrid.botm[i, r, c] for i in range(gwf.modelgrid.nlay)]
                )
            else:
                print("Point is not within model extent! Returning NaN.")
                zvec = np.nan
        elif gwf.modelgrid.grid_type == "vertex":
            res = ix.intersect([(x, y)], shapetype="point")
            idx = res["cellids"].astype(int)[0]
            if len(res) > 0:
                zvec = [gwf.modelgrid.top[idx]] + [
                    gwf.modelgrid.botm[i, idx] for i in range(gwf.modelgrid.nlay)
                ]
            else:
                print("Point is not within model extent! Returning NaN.")
                zvec = np.nan
        else:
            raise NotImplementedError(
                f"gridtype {gwf.modelgrid.grid_type} not (yet) implemented"
            )
    elif ds and not gwf:
        # assuming modflow type definition with 1 top and N botms
        # first check if point is within the extent
        if "angrot" in ds.attrs and ds.attrs["angrot"] != 0.0:
            import nlmod

            pol = nlmod.dims.get_extent_polygon(ds)
            if not Point(x, y).within(pol):
                return np.nan
        else:
            xmin, xmax, ymin, ymax = ds.attrs["extent"]
            if (x < xmin) or (x > xmax) or (y < ymin) or (y > ymax):
                return np.nan

        if ds.gridtype == "vertex":
            import nlmod

            cid = nlmod.dims.xy_to_icell2d((x, y), ds)
            sel = ds.sel(icell2d=cid)
            zvec = np.concatenate(([sel["top"].values], sel["botm"].values))
            mask = np.isnan(zvec)
            idx = np.where(~mask, np.arange(mask.size), 0)
            np.maximum.accumulate(idx, axis=0, out=idx)
            zvec[mask] = zvec[idx[mask]]
        else:
            sel = ds.sel(x=x, y=y, method="nearest")
            first_notna = np.nonzero(np.isfinite(np.atleast_1d(sel["top"].values)))[0][
                0
            ]
            if sel["top"].values.shape == tuple():
                top = np.atleast_1d(sel["top"].values)
            else:
                top = np.atleast_1d(sel["top"].values[[first_notna]])
            zvec = np.concatenate([top, sel["botm"].values])
            mask = np.isnan(zvec)
            idx = np.where(~mask, np.arange(mask.size), 0)
            np.maximum.accumulate(idx, axis=0, out=idx)
            zvec[mask] = zvec[idx[mask]]
    else:
        raise ValueError("Pass one of 'gwf' or 'ds'!")

    return zvec


@accessor.register_obscollection_accessor("gwobs")
class GwObsAccessor:
    def __init__(self, oc_obj):
        self._obj = oc_obj

    def set_tube_nr(
        self, radius=1, xcol="x", ycol="y", if_exists="error", add_to_meta=False
    ):
        """This method computes the tube numbers based on the location of the
        observations.

        Then it sets the value of the tube number:

        - in the ObsCollection dataframe
        - as the attribute of an Obs object
        - in the meta dictionary of the Obs object (only if add_to_meta is
          True)


        This method is useful for groundwater observations. If two or more
        observation points are close to each other they will be seen as one
        monitoring_well with multiple tubes. The tube_nr is based on the
        'screen_bottom' attribute of the observations in such a way that
        the deepest tube has the highest tube number.

        Parameters
        ----------
        radius : int, optional
            max distance between two observations to be seen as one location,
            by default 1
        xcol : str, optional
            column name with x coordinates, by default 'x'
        ycol : str, optional
            column name with y coordinates, by default 'y'
        if_exists : str, optional
            what to do if an observation point already has a tube_nr, options:
            'error', 'replace' or 'keep', by default 'error'
        add_to_meta : bool, optional
            if True the tube_nr is added to the meta dictionary
            of an observation. The default is False.

        Raises
        ------
        RuntimeError
            if the column tube_nr exists and if_exists='error' an error is
            raised
        """

        # check if column exists in obscollection
        if "tube_nr" in self._obj.columns:
            # set type to numeric
            if not is_numeric_dtype(self._obj["tube_nr"]):
                self._obj["tube_nr"] = pd.to_numeric(
                    self._obj["tube_nr"], errors="coerce"
                )

            # check if name should be replaced
            if if_exists == "error":
                raise RuntimeError(
                    "the column 'tube_nr' already exist, set"
                    "if_exists='replace' to replace the current values"
                )
            elif if_exists == "replace":
                self._obj["tube_nr"] = np.nan
        else:
            self._obj["tube_nr"] = np.nan

        # apply tube numbers to tubes that are close to eachother
        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, "tube_nr"]):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_tubes = np.sqrt(
                    (self._obj[xcol] - x) ** 2 + (self._obj[ycol] - y) ** 2
                )
                dup_x = self._obj.loc[distance_to_other_tubes < radius]
                if dup_x.shape[0] == 1:
                    self._obj._set_metadata_value(
                        name, "tube_nr", 1, add_to_meta=add_to_meta
                    )
                else:
                    dup_x2 = dup_x.sort_values("screen_bottom", ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        self._obj._set_metadata_value(
                            pb_dub, "tube_nr", i + 1, add_to_meta=add_to_meta
                        )

    def set_tube_nr_monitoring_well(
        self,
        loc_col,
        radius=1,
        xcol="x",
        ycol="y",
        if_exists="error",
        add_to_meta=False,
    ):
        """This method sets the tube_nr and monitoring_well name of an observation
        point based on the location of the observations.

        When two or more tubes are close to another, as defined by radius,
        they are set to the same `monitoring_well` and an increasing `tube_nr` based
        on depth.

        The value of the tube_nr and the monitoring_well are set:

        - in the ObsCollection dataframe
        - as the attribute of an Obs object
        - in the meta dictionary of the Obs object (only if add_to_meta is
          True)

        This method is useful for groundwater observations. If two or more
        observation points are close to each other they will be seen as one
        monitoring_well with multiple tubes. The tube_nr is based on the
        'screen_bottom' attribute of the observations in such a way that
        the deepest tube has the highest tube number. The monitoring_well is
        based on the named of the loc_col of the screen with the lowest
        tube_nr.

        Parameters
        ----------
        loc_col : str
            the column name with the names to use for the monitoring_well
        radius : int, optional
            max distance between two observations to be seen as one location,
            by default 1
        xcol : str, optional
            column name with x coordinates, by default 'x'
        ycol : str, optional
            column name with y coordinates, by default 'y'
        if_exists : str, optional
            what to do if an observation point already has a tube_nr, options:
            'error', 'replace' or 'keep', by default 'error'
        add_to_meta : bool, optional
            if True the tube_nr and location are added to the meta dictionary
            of an observation. The default is False.

        Raises
        ------
        RuntimeError
            if the column tube_nr exists and if_exists='error' an error
            is raised
        """

        # check if columns exists in obscollection
        if "tube_nr" in self._obj.columns or "monitoring_well" in self._obj.columns:
            if if_exists == "error":
                raise RuntimeError(
                    "the column 'tube_nr or monitoring_well' already exist, set"
                    "if_exists='replace' to replace the current values"
                )
            elif if_exists == "replace":
                self._obj["tube_nr"] = np.nan
                self._obj["monitoring_well"] = np.nan
        else:
            self._obj["tube_nr"] = np.nan
            self._obj["monitoring_well"] = np.nan

        # apply tube numbers to tubes that are close to eachother
        for name in self._obj.index:
            if np.isnan(self._obj.loc[name, "tube_nr"]):
                x = self._obj.loc[name, xcol]
                y = self._obj.loc[name, ycol]
                distance_to_other_tubes = np.sqrt(
                    (self._obj[xcol] - x) ** 2 + (self._obj[ycol] - y) ** 2
                )
                dup_x = self._obj.loc[distance_to_other_tubes < radius]
                if dup_x.shape[0] == 1:
                    self._obj._set_metadata_value(
                        name, "tube_nr", 1, add_to_meta=add_to_meta
                    )
                    monitoring_well = self._obj.loc[name, loc_col]
                    self._obj._set_metadata_value(
                        name,
                        "monitoring_well",
                        monitoring_well,
                        add_to_meta=add_to_meta,
                    )
                else:
                    dup_x2 = dup_x.sort_values("screen_bottom", ascending=False)
                    for i, pb_dub in enumerate(dup_x2.index):
                        if i == 0:
                            monitoring_well = self._obj.loc[pb_dub, loc_col]
                        self._obj._set_metadata_value(
                            pb_dub, "tube_nr", i + 1, add_to_meta=add_to_meta
                        )
                        self._obj._set_metadata_value(
                            pb_dub,
                            "monitoring_well",
                            monitoring_well,
                            add_to_meta=add_to_meta,
                        )

    def get_modellayers(self, gwf=None, ds=None):
        """Get the modellayer per observation. The layers can be obtained from
        the modflow model or can be defined in zgr.

        Parameters
        ----------
        gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
            modflow model
        ds : xarray.Dataset
            xarray Dataset with with top and bottoms, must have
            dimensions 'x' and 'y' and variables 'top' and 'bot'

        Returns
        -------
        pd.Series with the modellayers of each observation
        """
        modellayers = []
        for o in self._obj.obs.values:
            logger.debug("-" * 10 + f"\n {o.name}:")

            modellayers.append(o.gwobs.get_modellayer_modflow(gwf=gwf, ds=ds))

        modellayers = pd.Series(
            index=self._obj.index, data=modellayers, name="modellayer"
        )

        return modellayers

    def get_regis_layers(self):
        """Get the regis layer per observation.

        Parameters
        ----------


        Returns
        -------
        pd.Series with the names of the regis layer of each observation
        """
        if tqdm is not None:
            tqdm.pandas()
            return self._obj.obs.progress_apply(lambda o: o.gwobs.get_regis_layer())
        else:
            return self._obj.obs.apply(lambda o: o.gwobs.get_regis_layer())


@accessor.register_obs_accessor("gwobs")
class GeoAccessorObs:
    def __init__(self, obs):
        self._obj = obs

    def get_modellayer_modflow(self, gwf=None, ds=None, left=-999, right=999):
        """Add modellayer to meta dictionary.

        Parameters
        ----------
        gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
            modflow model
        ds : xarray.Dataset
            xarray Dataset with with top and bottoms, must have
            dimensions 'x' and 'y' and variables 'top' and 'bot'

        Returns
        -------
        int
            modellayer
        """

        zvec = get_zvec(self._obj.x, self._obj.y, gwf=gwf, ds=ds)

        if np.all(np.isnan(zvec)):
            return np.nan
        else:
            modellayer = get_modellayer_from_screen_depth(
                self._obj.screen_top,
                self._obj.screen_bottom,
                zvec,
                left=left,
                right=right,
            )
            return modellayer

    def get_regis_layer(self):
        """find the name of the REGIS layer based on the tube screen depth.

        Parameters
        ----------


        Returns
        -------
        str
            name of REGIS layer
        """
        import xarray as xr

        if np.isnan(self._obj.screen_top) or np.isnan(self._obj.screen_bottom):
            return "nan"

        # connect to regis netcdf
        regis_url = r"http://www.dinodata.nl:80/opendap/REGIS/REGIS.nc"
        regis_ds = xr.open_dataset(regis_url, decode_times=False)

        # rename layer in regis netcdf
        regis_ds = regis_ds.rename({"layer": "layer_old"})
        regis_ds.coords["layer"] = regis_ds.layer_old.astype(str)
        regis_ds = regis_ds.swap_dims({"layer_old": "layer"})

        # get zvec regis netcdf
        z = (
            regis_ds["bottom"]
            .sel(x=self._obj.x, y=self._obj.y, method="nearest")
            .dropna(dim="layer")
        )
        ztop = regis_ds["top"].sel(x=self._obj.x, y=self._obj.y, method="nearest").max()
        zvec = np.concatenate([[ztop.values], z.values])

        # get index of regis model layer
        layer_i = get_modellayer_from_screen_depth(
            self._obj.screen_top,
            self._obj.screen_bottom,
            zvec,
            left=-999,
            right=999,
        )

        if layer_i == 999:
            return "above"
        elif layer_i == -999:
            return "below"

        return str(z.layer[layer_i].values)
