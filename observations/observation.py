'''
module with a number of observation classes.

The Obs class is a subclass of a pandas DataFrame with
additional attributes and methods. The specific classes (GroundwaterObs,
WaterlvlObs, ...) are subclasses of the Obs class.

The subclasses of a dataframe can have additional attributes and methods.
Additional attributes have to be defined in the '_metadata' attribute. In order
to keep the subclass methods and attributes when selecting or slicing an object
you need the '_constructor' method.


More information about subclassing pandas DataFrames can be found here:
http://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

'''

import copy
import os

import numpy as np
from pandas import DataFrame, Series, datetime

from . import io_dino


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
    _metadata = ['x', 'y', 'name',
                 'meta',
                 'filename']

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

    def to_collection_dict(self):
        """get dictionary of an Obs object.

        Use this method to create a dataframe from a collection of Obs objects.

        Returns
        -------
        d : dictionary
            dictionary with Obs information
        """
        # meta is mutable so only use a copy
        d = copy.deepcopy(self.meta)

        d.update({'obs': self,
                  'name': self.name,
                  'x': self.x,
                  'y': self.y})

        return d

    def to_interactive_plot(self, savedir=None,
                            plot_columns=['Stand_m_tov_NAP'],
                            markers=['line'],
                            p=None,
                            plot_legend_names=[''],
                            plot_freq=[None], tmin=None, tmax=None,
                            hoover_names=['Peil'],
                            hoover_date_format="%Y-%m-%d",
                            ylabel='m NAP', colors=['blue'],
                            add_filter_to_legend=False,
                            return_filename=False):
        """Create an interactive plot of the observations using bokeh

        To-Do
        -----
        add options for hoovers, markers, linestyle

        Parameters
        ----------
        savedir : str, optional
            directory used for the folium map and bokeh plots
        plot_columns : list of str, optional
            name of the column in the obs df that will be plotted with bokeh
        markers : list of str, optional
            type of markers that can be used for plot, 'line' and 'circle' are
            supported
        p : bokeh.plotting.figure, optional
            reference to existing figure, if p is None a new figure is created
        plot_legend_names : list of str, optional
            legend in bokeh plot
        plot_freq : list of str, optional
            bokeh plot is resampled with this frequency to reduce the size
        tmin : dt.datetime, optional
            start date for timeseries plot
        tmax : dt.datetime, optional
            end date for timeseries plot
        hoover_names : list of str, optional
            names will be displayed together with the plot_column value
            when hoovering over plot
        hoover_date_format : str, optional
            date format to use when hoovering over a plot
        ylabel : str, optional
            label on the y-axis
        colors : list of str, optional
            colors used for the plots
        add_filter_to_legend : boolean, optional
            if True the attributes bovenkant_filter and onderkant_filter
            are added to the legend name
        return_filename : boolean, optional
            if True filename will be returned

        Returns
        -------
        fname_plot : str or bokeh plot
            filename of the bokeh plot or reference to bokeh plot
        """

        from bokeh.plotting import figure
        from bokeh.models import ColumnDataSource, HoverTool
        from bokeh.plotting import save
        from bokeh.resources import CDN

        # create plot dataframe
        plot_df = self[tmin:tmax].copy()
        plot_df['date'] = plot_df.index.strftime(hoover_date_format)
        if plot_df.empty or plot_df[plot_columns].isna().all().all():
            raise ValueError(
                '{} has no data between {} and {}'.format(self.name, tmin, tmax))

        # create plot
        if p is None:
            p = figure(plot_width=600, plot_height=400, x_axis_type='datetime',
                       title='')
            p.yaxis.axis_label = ylabel

        # get x axis
        xcol = self.index.name
        if xcol is None:
            xcol = 'index'

        # plot multiple columns
        for i, column in enumerate(plot_columns):
            # legend name
            if add_filter_to_legend:
                lname = '{} {} (NAP {:.2f} - {:.2f})'.format(plot_legend_names[i], self.name,
                                                             self.onderkant_filter,
                                                             self.bovenkant_filter)
            else:
                lname = '{} {}'.format(plot_legend_names[i], self.name)

            # resample data
            if plot_freq[i] is not None:
                source = ColumnDataSource(
                    plot_df[[column, 'date']].resample(plot_freq[i]).nearest())
            else:
                source = ColumnDataSource(plot_df[[column, 'date']])

            # plot data
            if markers[i] == 'line':
                p.line(xcol, column, source=source, color=colors[i], legend=lname,
                       alpha=0.8, muted_alpha=0.2)
            elif markers[i] == 'circle':
                p.circle(xcol, column, source=source, color=colors[i], legend=lname,
                         alpha=0.8, muted_alpha=0.2)
            else:
                raise NotImplementedError("marker '{}' invalid. Only line and"
                                          "circle are currently available".format(markers[i]))

        # hoover options
        tooltips = []
        for i, column in enumerate(plot_columns):
            tooltips.append((hoover_names[i], "@{}".format(column)))
        tooltips.append(('date', "@date"))
        hover = HoverTool(tooltips=tooltips)

        p.add_tools(hover)

        p.legend.location = "top_left"
        p.legend.click_policy = "mute"

        # save plot
        if savedir is not None:
            self.iplot_fname = os.path.join(savedir, self.name + '.html')
            save(p, self.iplot_fname, resources=CDN, title=self.name)

        if return_filename:
            return self.iplot_fname
        else:
            return p

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326',
                    add_to_meta=True):
        """get lattitude and longitude from x and y attributes

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon
        add_to_meta : boolean, optional
            if True the new coordinates are added to the meta dictionary

        Returns
        -------
        lon, lat : longitude and lattitude of x, y coordinates

        """

        from pyproj import Proj, transform

        inProj = Proj(init=in_epsg)
        outProj = Proj(init=out_epsg)

        if np.isnan(self.x) or np.isnan(self.y):
            lon, lat = np.nan, np.nan
        else:
            lon, lat = transform(inProj, outProj, self.x, self.y)

        if add_to_meta:
            self.meta['lon'] = lon
            self.meta['lat'] = lat

        return lat, lon

    def get_seasonal_stat(self, column_name='Stand_m_tov_NAP', stat='mean',
                          winter_months=[1, 2, 3, 4, 11, 12],
                          summer_months=[5, 6, 7, 8, 9, 10]):
        """get statistics per season

        Parameters
        ----------
        column_name : str, optional
            column name of the  observation data you want stats for
        stat : str, optional
            type of statistics, all statisics from df.describe() are available
        winter_months : list of int, optional
            month number of winter months
        summer_months : list of int, optional
            month number of summer months


        Returns
        -------
        winter_stats, summer_stats
            two lists with the statistics for the summer and the winter.

        """

        if self.empty:
            df = DataFrame(index=[self.name], data={'winter_{}'.format(stat): [np.nan],
                                                    'summer_{}'.format(stat): [np.nan]})
        else:
            winter_stat = self.loc[self.index.month.isin(
                winter_months)].describe().loc[stat, column_name]
            summer_stat = self.loc[self.index.month.isin(
                summer_months)].describe().loc[stat, column_name]
            df = DataFrame(index=[self.name], data={'winter_{}'.format(stat): [winter_stat],
                                                    'summer_{}'.format(stat): [summer_stat]})

        return df

    def obs_per_year(self, col):
        if self.empty:
            return Series()
        else:
            return self.groupby(self.index.year).count()[col]

    def consecutive_obs_years(self, col, min_obs=12):

        obs_per_year = self.obs_per_year(col=col)

        # Add missing years
        if obs_per_year.empty:
            # no obs, set series to current year with 0 obs
            obs_per_year_all = Series(index=[datetime.now().year], data=0)
        else:
            obs_per_year_all = Series(index=range(obs_per_year.index[0],
                                                  obs_per_year.index[-1] + 1))
            obs_per_year_all.loc[obs_per_year.index] = obs_per_year

        mask_obs_per_year = obs_per_year_all >= min_obs
        mask_obs_per_year.loc[obs_per_year_all.isna()] = np.nan
        mask_obs_per_year.loc[mask_obs_per_year == 0] = np.nan
        cumsum = mask_obs_per_year.cumsum().fillna(method="pad")
        reset = -cumsum.loc[mask_obs_per_year.isnull()].diff().fillna(cumsum)
        result = mask_obs_per_year.where(
            mask_obs_per_year.notnull(), reset).cumsum()

        return result


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

        if the pandas.DataFrame has a column 'Stand_m_tov_NAP' a lot of
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
    def from_dino_server(cls, name, filternr=1.,
                         tmin="2000-01-01", tmax="2010-01-01",
                         x=np.nan, y=np.nan,
                         get_metadata=True,
                         **kwargs):
        """download dino data from the server.

        Parameters
        ----------
        name : str, optional
            name of the peilbuis, i.e. B57F0077
        filternr : float, optional
            filter_nr of the peilbuis, i.e. 1.
        tmin : str
            start date in format YYYY-MM-DD
        tmax : str
            end date in format YYYY-MM-DD
        x : int, float, optional
            the x coördinate of the measurement point (not read from server)
        y : int, float, optional
            the y coördinate of the measurement point (not read from server)
        get_metadata : boolean, optional
            download extra metadata from the server (see Notes)
        kwargs : key-word arguments
            these arguments are passed to dino.findMeetreeks functie

        Notes
        -----
        For now only the maaiveld is used from the extra metadata, this method
        should be improved to add more metadata from the server if get_metadata
        is True
        """

        measurements, meta = io_dino.download_dino_groundwater(name,
                                                               filternr,
                                                               tmin, tmax,
                                                               x, y,
                                                               **kwargs)

        if get_metadata:
            import art_tools as art
            raw_meta = art.dino_wfs.get_dino_piezometer_metadata(
                [meta['locatie']])
            if raw_meta != []:
                meta_extra = raw_meta[0]

                for piezometer in meta_extra.pop('levels'):
                    if piezometer['piezometerNr'] == format(filternr, '03'):
                        meta_extra.update(piezometer)
                        break
                maaiveld = meta_extra.pop('surfaceElevation')
                # the meta_extra dictionary has a lot
                # of information, some in nested dictionaries and not always
                # with the same keys. This will give issues when I automatically
                # update the meta dictionary. Therefore this is not implemented.
                # For now only maaiveld is used. If you want to use more info
                # you have to modify the code here.
                # meta.update(meta_extra)
            else:
                maaiveld = np.nan
        else:
            maaiveld = np.nan
        meta['maaiveld'] = maaiveld

        return cls(measurements, meta=meta, x=x, y=y,
                   onderkant_filter=meta['onderkant_filter'],
                   bovenkant_filter=meta['bovenkant_filter'],
                   name=meta['name'], locatie=meta['locatie'],
                   maaiveld=maaiveld,
                   meetpunt=meta['meetpunt'], filternr=meta['filternr'])

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

        if fname is not None:
            # read dino csv file

            measurements, meta = io_dino.read_dino_groundwater_csv(
                fname, **kwargs)

            return cls(measurements, meta=meta, **meta)
        else:
            raise ValueError(
                'specify either the name or the filename of the measurement point')

    @classmethod
    def from_wiski(cls, fname, **kwargs):

        from . import io_wiski

        header, data = io_wiski.read_wiski_file(fname, **kwargs)
        metadata = {}
        if 'Station Site' in header.keys():
            metadata['locatie'] = header['Station Site']
            header['locatie'] = header['Station Site']

        if 'x' in header.keys():
            metadata['x'] = header['x']
        if "y" in header.keys():
            metadata['y'] = header['y']
        if 'name' in header.keys():
            metadata['name'] = header['name']

        return cls(data, meta=header, **metadata)

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

    def get_pb_modellayer(self, ml, zgr=None, verbose=False):
        """Add modellayer to meta dictionary

        Parameters
        ----------
        ml : flopy.modflow.mf.Modflow
            modflow model
        zgr : np.3darray, optional
            array containing model layer elevation
            information (the default is None, which
            gets this information from the dis object)
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        """
        from .modflow import get_pb_modellayer
        self.meta['modellayer'] = \
            get_pb_modellayer(np.array([self.x]) - ml.modelgrid.xoffset,
                              np.array([self.y]) - ml.modelgrid.yoffset,
                              np.array([self.bovenkant_filter]),
                              np.array([self.onderkant_filter]),
                              ml, zgr, verbose)[0]


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
    def from_dino_file(cls, fname, **kwargs):
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

    _metadata = Obs._metadata + \
        ['locatie', 'metadata_available'
         ]

    def __init__(self, *args, **kwargs):

        self.locatie = kwargs.pop('locatie', '')
        self.metadata_available = kwargs.pop('metadata_available', np.nan)

        super(WaterlvlObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return WaterlvlObs

    @classmethod
    def from_dino_file(cls, fname, **kwargs):
        '''read a dino file with waterlvl data

        Parameters
        ----------
        fname : str
            dino csv filename
        kwargs : key-word arguments
            these arguments are passed to io_dino.read_dino_waterlvl_csv
        '''

        measurements, meta = io_dino.read_dino_waterlvl_csv(fname, **kwargs)

        return cls(measurements, meta=meta, **meta)


class ModelObs(Obs):
    """class for model point results.

    Subclass of the Obs class
    """

    _metadata = Obs._metadata + \
        ['model']

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

    _metadata = Obs._metadata + \
        ['station']

    def __init__(self, *args, **kwargs):

        self.station = kwargs.pop('station', np.nan)

        super(KnmiObs, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return KnmiObs

    @classmethod
    def from_knmi(cls, stn, variable, startdate=None, enddate=None,
                  normalize_index=True, verbose=False):
        from . import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_stn(stn, variable,
                                                   startdate, enddate,
                                                   normalize_index=normalize_index,
                                                   verbose=verbose)
        return cls(ts, meta=meta, **meta)

    @classmethod
    def from_nearest_xy(cls, x, y, variable, startdate=None, enddate=None,
                        normalize_index=True, verbose=False):
        from . import io_knmi

        ts, meta = io_knmi.get_knmi_timeseries_xy(x, y, variable,
                                                  startdate, enddate,
                                                  normalize_index=normalize_index,
                                                  verbose=verbose)

        return cls(ts, meta=meta, **meta)

    @classmethod
    def from_obs(cls, obs, variable, startdate=None, enddate=None,
                 normalize_index=True, verbose=False):

        from . import io_knmi

        x = obs.meta["x"]
        y = obs.meta["y"]

        if startdate is None:
            startdate = obs.index[0]
        if enddate is None:
            enddate = obs.index[-1]

        ts, meta = io_knmi.get_knmi_timeseries_xy(x, y, variable,
                                                  startdate, enddate,
                                                  normalize_index=normalize_index,
                                                  verbose=verbose)

        return cls(ts, meta=meta, **meta)
