import os

import numpy as np

from . import accessor

import logging
logger = logging.getLogger(__name__)


@accessor.register_obscollection_accessor("plots")
class CollectionPlots:

    def __init__(self, oc_obj):
        """Object containing plotting methods for ObsCollections.

        Parameters
        ----------
        oc : ObsCollection
            ObsCollection instance
        """
        self._obj = oc_obj

    def interactive_plots(self, savedir,
                          tmin=None, tmax=None,
                          per_location=True,
                          **kwargs):
        """Create interactive plots of the observations using bokeh.

        Parameters
        ----------
        savedir : str
            directory used for the folium map and bokeh plots
        tmin : dt.datetime, optional
            start date for timeseries plot
        tmax : dt.datetime, optional
            end date for timeseries plot
        per_location : bool, optional
            if True plot multiple filters on the same location in one figure
        **kwargs :
            will be passed to the Obs.to_interactive_plot method, options
            include:

            - plot_columns : list of str
            - hoover_names : list of str
            - plot_freq : str
            - plot_legend_name : str
            - ylabel : str
            - color : str
            - add_filter_to_legend : boolean
        """
        _color_cycle = (
            'blue',
            'olive',
            'lime',
            'red',
            'orange',
            'yellow',
            'purple',
            'silver',
            'powderblue',
            'salmon',
            'tan')

        if per_location:
            plot_names = self._obj.groupby('locatie').count().index
        else:
            plot_names = self._obj.index

        for name in plot_names:

            if per_location:
                oc = self._obj.loc[self._obj.locatie
                                   == name, 'obs'].sort_index()
            else:
                oc = self._obj.loc[[name], 'obs']

            p = None
            for i, o in enumerate(oc.values):
                if i == 10:
                    raise NotImplementedError(
                        'cannot add more than 10 lines to a single plot')
                try:
                    p = o.plots.interactive_plot(savedir=savedir,
                                                 p=p,
                                                 tmin=tmin, tmax=tmax,
                                                 colors=[_color_cycle[i + 1]],
                                                 return_filename=False,
                                                 **kwargs)
                    logger.info(f'created iplot -> {o.name}')
                except ValueError:
                    logger.error(f'{o.name} has no data between {tmin} and {tmax}')
                    o.iplot_fname = None

    def interactive_map(self, plot_dir, m=None,
                        tiles='OpenStreetMap',
                        fname=None,
                        per_location=True,
                        color='blue',
                        legend_name=None,
                        add_legend=True,
                        map_label='',
                        map_label_size=20,
                        col_name_lat='lat',
                        col_name_lon='lon',
                        zoom_start=13,
                        create_interactive_plots=True,
                        **kwargs):
        """Create an interactive map with interactive plots using folium and
        bokeh.

        Notes
        -----
        Some notes on this method:

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
        per_location : bool, optional
            if True plot multiple filters on the same location in one figure
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
        create_interactive_plots : boolean, optional
            if True interactive plots will be created, if False the iplot_fname
            attribute of the observations is used.
        **kwargs :
            will be passed to the to_interactive_plots method options are:

            - plot_columns : list of str
            - hoover_names : list of str
            - plot_legend_name : str
            - ylabel : str
            - add_filter_to_legend : boolean
            - plot_freq : str
            - tmin : dt.datetime
            - tmax : dt.datetime

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
            self._obj.plots.interactive_plots(savedir=plot_dir,
                                              per_location=per_location,
                                              **kwargs)

        # check if observation collection has lat and lon values
        if (not col_name_lat in self._obj.columns) and (not col_name_lon in self._obj.columns):
            self._obj.geo.set_lat_lon()

        # determine start location of map
        northing = np.mean(
            (self._obj[col_name_lat].min(), self._obj[col_name_lat].max()))
        easting = np.mean(
            (self._obj[col_name_lon].min(), self._obj[col_name_lon].max()))

        # create map if no map is given
        if m is None:
            m = folium.Map([northing, easting], zoom_start=zoom_start)

        # add the point observations with plots to the map
        group_name = '<span style=\\"color: {};\\">{}</span>'.format(
            color, legend_name)
        group = folium.FeatureGroup(name=group_name)

        if per_location:
            plot_names = self._obj.groupby('locatie').count().index
        else:
            plot_names = self._obj.index

        for name in plot_names:
            if per_location:
                oc = self._obj.loc[self._obj.locatie
                                   == name, 'obs'].sort_index()
                o = oc.iloc[-1]
                name = o.name
            else:
                o = self._obj.loc[name, 'obs']

            if o.iplot_fname is not None:
                with open(o.iplot_fname, 'r') as f:
                    bokeh_html = f.read()

                iframe = branca.element.IFrame(
                    html=bokeh_html, width=620, height=420)
                popup = folium.Popup(iframe, max_width=620)

                folium.CircleMarker([self._obj.loc[o.name, col_name_lat],
                                     self._obj.loc[o.name, col_name_lon]],
                                    icon=folium.Icon(icon='signal'), fill=True,
                                    color=color,
                                    popup=popup).add_to(group)

                if map_label != '':
                    folium.map.Marker(
                        [self._obj.loc[name, col_name_lat], self._obj.loc[name, col_name_lon]], icon=DivIcon(
                            icon_size=(
                                150, 36), icon_anchor=(
                                0, 0), html='<div style="font-size: %ipt">%s</div>' %
                            (map_label_size, o.meta[map_label]))).add_to(group)
            else:
                logger.info(f'no iplot available for {o.name}')

        group.add_to(m)

        # add legend
        if add_legend:
            folium.map.LayerControl('topright', collapsed=False).add_to(m)

        # save map
        #filename and path
        if fname is not None:
            if not fname.endswith('.html'):
                fname = fname + '.html'
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)
            m.save(os.path.join(plot_dir, fname))

        return m


@accessor.register_obs_accessor("plots")
class ObsPlots:
    def __init__(self, obs):
        self._obj = obs

    def interactive_plot(self,
                         savedir=None,
                         plot_columns=['stand_m_tov_nap'],
                         markers=['line'],
                         p=None,
                         plot_legend_names=[''],
                         plot_freq=[None],
                         tmin=None,
                         tmax=None,
                         hoover_names=['Peil'],
                         hoover_date_format="%Y-%m-%d",
                         ylabel='m NAP',
                         colors=['blue'],
                         add_filter_to_legend=False,
                         return_filename=False):
        """Create an interactive plot of the observations using bokeh.

        Todo:

        - add options for hoovers, markers, linestyle

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
        plot_df = self._obj[tmin:tmax].copy()
        plot_df['date'] = plot_df.index.strftime(hoover_date_format)
        if plot_df.empty or plot_df[plot_columns].isna().all().all():
            raise ValueError(
                '{} has no data between {} and {}'.format(self._obj.name, tmin, tmax))

        # create plot
        if p is None:
            p = figure(plot_width=600, plot_height=400, x_axis_type='datetime',
                       title='')
            p.yaxis.axis_label = ylabel

        # get x axis
        xcol = self._obj.index.name
        if xcol is None:
            xcol = 'index'

        # get color
        if len(colors) < len(plot_columns):
            colors = colors * len(plot_columns)

        # plot multiple columns
        for i, column in enumerate(plot_columns):
            # legend name
            if add_filter_to_legend:
                lname = '{} {} (NAP {:.2f} - {:.2f})'.format(plot_legend_names[i], self._obj.name,
                                                             self._obj.onderkant_filter,
                                                             self._obj.bovenkant_filter)
            else:
                lname = '{} {}'.format(plot_legend_names[i], self._obj.name)

            # resample data
            if plot_freq[i] is None:
                source = ColumnDataSource(plot_df[[column, 'date']])
            else:
                source = ColumnDataSource(
                    plot_df[[column, 'date']].resample(plot_freq[i]).first())

            # plot data
            if markers[i] == 'line':
                p.line(xcol, column, source=source, color=colors[i],
                       legend_label=lname,
                       alpha=0.8, muted_alpha=0.2)
            elif markers[i] == 'circle':
                p.circle(xcol, column, source=source, color=colors[i],
                         legend_label=lname,
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
            if not os.path.isdir(savedir):
                os.makedirs(savedir)
            self._obj.iplot_fname = os.path.join(
                savedir, self._obj.name + '.html')
            save(p, self._obj.iplot_fname, resources=CDN, title=self._obj.name)

        if return_filename:
            return self._obj.iplot_fname
        else:
            return p
