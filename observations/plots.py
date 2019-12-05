import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import accessor, util


@accessor.register_obscollection_accessor("plots")
class CollectionPlots:

    def __init__(self, oc_obj):
        """Object containing plotting methods for ObsCollections

        Parameters
        ----------
        oc : ObsCollection
            ObsCollection instance

        """
        self._obj = oc_obj

    def data_frequency(
            self,
            column_name='Stand_m_tov_NAP',
            intervals=None,
            ignore=[
                'seconde',
                'minuut',
                '14-daags'],
            normtype='log',
            cmap='viridis_r',
            set_yticks=False,
            figsize=(
                10,
                8),
            **kwargs):
        """plot data availability and frequency of observation collection

        Parameters
        ----------
        column_name : str, optional
            column name of the timeseries data
        intervals: dict, optional
            A dict with frequencies as keys and number of seconds as values
        ignore: list, optional
            A list with frequencies in intervals to ignore
        normtype: str, optional
            Determines the type of color normalisations, default is 'log'
        cmap: str, optional
            A reference to a matplotlib colormap
        set_yticks: bool, optional
            Setthe names of the series as ytciks
        figsize: tuple, optional
            The size of the new figure in inches (h,v)
        **kwargs: dict, optional
            Extra arguments are passed to matplotlib.pyplot.subplots()



        """

        art = util._import_art_tools()
        series_list = []
        for o in self._obj.obs.values:
            series = o[column_name]
            series.name = o.name
            series_list.append(series)

        ax = art.data_availablity(
            series_list,
            intervals=intervals,
            ignore=ignore,
            normtype=normtype,
            cmap=cmap,
            set_yticks=set_yticks,
            figsize=figsize,
            **kwargs)

        return ax

    def interactive_plots(self, savedir,
                          tmin=None, tmax=None,
                          per_location=True,
                          verbose=True, **kwargs):
        """Create interactive plots of the observations using bokeh

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


        verbose : boolean, optional
            Print additional information to the screen (default is False).
        **kwargs :
            will be passed to the Obs.to_interactive_plot method

            plot_columns : list of str
                name of the column in the obs df that will be plotted with bokeh
            hoover_names : list of str, optional
                names will be displayed together with the plot_column value
                when hoovering over plot
            plot_freq : str, optional
                bokeh plot is resampled with this frequency to reduce the size
            plot_legend_name : str, optional
                legend in bokeh plot
            ylabel : str, optional
                label on the y-axis
            color : str, optional
                color of the lines on the plots
            add_filter_to_legend : boolean, optional
                if True the attributes bovenkant_filter and onderkant_filter
                are added to the graph legend

        Returns
        -------

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
                oc = self._obj.loc[self._obj.locatie ==
                                   name, 'obs'].sort_index()
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
                    if verbose:
                        print('created iplot -> {}'.format(o.name))

                except ValueError:
                    if verbose:
                        print('{} has no data between {} and {}'.format(
                            o.name, tmin, tmax))
                    o.iplot_fname = None

    def interactive_map(self, plot_dir, m=None,
                        tiles='OpenStreetMap', fname=None,
                        per_location=True,
                        color='blue', legend_name=None,
                        add_legend=True,
                        map_label='', map_label_size=20,
                        col_name_lat='lat', col_name_lon='lon',
                        zoom_start=13,
                        create_interactive_plots=True,
                        verbose=False, **kwargs):
        """ create an interactive map with interactive plots using folium
        and bokeh.

        Notes
        -----
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
        create_interactive_plots : boolean, optinal
            if True interactive plots will be created, if False the iplot_fname
            attribute of the observations is used.
        verbose : boolean, optional
            Print additional information to the screen (default is False).
        **kwargs :
            will be passed to the to_interactive_plots method options are:

            plot_columns : list str, optional
                name of the column in the obs df that will be plotted with bokeh
            hoover_names : list of str, optional
                names will be displayed together with the plot_column value
                when hoovering over plot
            plot_legend_name : str, optional
                the name of the observation points in the time series plot
            ylabel : str, optional
                label on the y-axis
            add_filter_to_legend : boolean, optional
                if True the attributes bovenkant_filter and onderkant_filter
                are added to the graph title
            plot_freq : str, optional
                bokeh plot is resampled with this frequency to reduce the size
                of the complete folium map
            tmin : dt.datetime, optional
                start date for timeseries plot
            tmax : dt.datetime, optional
                end date for timeseries plot

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
                                              verbose=verbose,
                                              per_location=per_location,
                                              **kwargs)

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
                oc = self._obj.loc[self._obj.locatie ==
                                   name, 'obs'].sort_index()
                o = oc.iloc[-1]
            else:
                o = self._obj.loc[name, 'obs']

            if o.iplot_fname is not None:
                with open(o.iplot_fname, 'r') as f:
                    bokeh_html = f.read()

                iframe = branca.element.IFrame(
                    html=bokeh_html, width=620, height=420)
                popup = folium.Popup(iframe, max_width=620)

                folium.CircleMarker([o.meta[col_name_lat], o.meta[col_name_lon]],
                                    icon=folium.Icon(icon='signal'), fill=True,
                                    color=color,
                                    popup=popup).add_to(group)

                if map_label != '':
                    folium.map.Marker(
                        [
                            o.meta[col_name_lat], o.meta[col_name_lon]], icon=DivIcon(
                            icon_size=(
                                150, 36), icon_anchor=(
                                0, 0), html='<div style="font-size: %ipt">%s</div>' %
                            (map_label_size, o.meta[map_label]))).add_to(group)
            elif verbose:
                print('no iplot available for {}'.format(o.name))

        group.add_to(m)

        # add legend
        if add_legend:
            folium.map.LayerControl('topright', collapsed=False).add_to(m)

        # save map
        #filename and path
        if fname is None:
            fname = self._obj.name + '.html'
        else:
            if not fname.endswith('.html'):
                fname = fname + '.html'
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)
            m.save(os.path.join(plot_dir, fname))

        return m

    def map(self, ax=None, figsize=(15, 15), label='gws',
            edgecolor='black', facecolor='green',
            marker='o', markersize=100, xlim=None, ylim=None, add_topo=True,
            verbose=False):
        """plot observation points on a map

        Parameters
        ----------
        ax : matplotlib.axes
            reference to axes object
        figsize : tuple
            figure size
        label : str
            label used in legend
        edgecolor : str

        facecolor : str

        marker : str

        markersize : int

        xlim : tuple

        ylim : tuple

        add_topo : boolean, optional
            topo is added using art tools

        Returns
        -------
        ax : matplotlib.axes
            reference to axes object

        """

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)

        ax.scatter(self._obj.x, self._obj.y, label=label, edgecolor=edgecolor,
                   marker=marker, facecolor=facecolor,
                   s=markersize)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if add_topo:
            # attempt art_tools import
            art = util._import_art_tools()
            art.OpenTopo(ax=ax, verbose=verbose).plot(
                verbose=verbose, alpha=0.5)

        ax.legend()

        return ax

    def mapgraphs(self, graph=None, plots_per_map=10, figsize=(16, 10),
                  extent=None, plot_column='Stand_m_tov_NAP',
                  per_location=True,
                  plot_func=None,
                  plot_xlim=(None, None), plot_ylim=(None, None),
                  min_dy=2.0, vlines=[],
                  savefig=None, map_gdf=None, map_gdf_kwargs={},
                  verbose=True):
        """make mapgraph plots of obs collection data

        Parameters
        ----------
        graph : np.ndarray, optional
            passed to the art.MapGraph funcion
        plots_per_map : int, optional
            number of plots per mapgraph, should be consistent with graph argument
        figsize : tuple, optional
            figure size
        extent : list or tuple, optional
            extent of the map [xmin, xmax, ymin, ymax]
        plot_column : str, optional
            column name of data to be plotted in graphs
        per_location : bool, optional
            if True plot multiple filters on the same location in one figure
        plot_func : function, optional,
            if not None this function is used to make the plots
        plot_xlim : list of datetime, optional
            xlim of the graphs
        plot_ylim : tuple or str, optional
            the ylim of all the graphs or
            'max_dy' -> all graphs get the maximum ylim
            'min_dy' -> all graphs get the same dy unless the series has a dy
                        bigger than min_dy
        min_dy : float, optional
            only used if plot_ylim is 'min_dy'
        vlines : list, optional
            list with x values on which a vertical line is plotted
        savefig : str, optional
            path to save the figure. {plot_number}.png will be append to this path
        map_gdf : GeoDataFrame, optional
            if not None the GeoDataFrame is also plotted on the map
        map_gdf_kwargs : dic, optional
            kwargs that will be passed to the map_gdf.plot() method


        Returns
        -------
        mg_list : list of art_tools.plots.MapGraph
            list of mapgraph objects

        """

        # attempt art_tools import
        art = util._import_art_tools()

        if graph is None:
            graph = np.full((4, 4), True)
            graph[1:, 1:-1] = False

        if plot_ylim == 'max_dy':
            plot_ylim = self._obj.get_min_max(obs_column=plot_column)

        mg_list = []
        if per_location:
            plot_names = self._obj.groupby(
                ['locatie', 'x', 'y'], as_index=False).count()
            plot_names.set_index('locatie', inplace=True)
        else:
            plot_names = self._obj.copy()

        for k, xyo in plot_names.groupby(np.arange(plot_names.shape[0]) // plots_per_map):

            mg = art.MapGraph(xy=xyo[['x', 'y']].values, graph=graph,
                              figsize=figsize, extent=extent)

            plt.sca(mg.mapax)
            plt.yticks(rotation=90, va="center")
            art.OpenTopo(ax=mg.mapax, verbose=verbose).plot(alpha=0.75,
                                                            verbose=verbose)
            if map_gdf is not None:
                map_gdf.plot(ax=mg.mapax, **map_gdf_kwargs)

            for i, name in enumerate(xyo.index):
                if per_location:
                    oc = self._obj.loc[self._obj.locatie == name].sort_index()
                else:
                    oc = self._obj.loc[[name]]

                ax = mg.ax[i]
                if plot_func:
                    ax = plot_func(self._obj, oc, ax,
                                   plot_xlim=plot_xlim,
                                   plot_column=plot_column)
                else:
                    for o in oc.obs.values:
                        try:
                            o[plot_column][plot_xlim[0]:plot_xlim[1]].plot(
                                ax=ax, lw=.5, marker='.', markersize=1., label=o.name)
                        except TypeError:
                            o[plot_column].plot(ax=ax, lw=.5, marker='.',
                                                markersize=1., label=o.name)
                            ax.set_xlim(plot_xlim)
                    if plot_ylim == 'min_dy':
                        ax.autoscale(axis='y', tight=True)
                        ylim = ax.get_ylim()
                        dy = ylim[1] - ylim[0]
                        if dy < min_dy:
                            ylim = (ylim[0] - (min_dy - dy) / 2,
                                    ylim[1] + (min_dy - dy) / 2)
                        ax.set_ylim(ylim)
                    else:
                        ax.set_ylim(plot_ylim)

                    for line in vlines:
                        ylim = ax.get_ylim()
                        ax.vlines(line, ylim[0] - 100, ylim[1] + 100, ls='--',
                                  lw=2.0, color='r')
                    ax.legend(loc='best')

            mg.figure.tight_layout()
            if savefig is not None:
                mg.figure.savefig(savefig + "{0}.png".format(k),
                                  dpi=300, bbox_inches="tight")
            mg_list.append(mg)

        # for some reason you have to set the extent at the end again
        if extent is not None:
            for mg_m in mg_list:
                mg_m.mapax.axis(extent)

        return mg_list


@accessor.register_obs_accessor("plots")
class ObsPlots:
    def __init__(self, obs):
        self._obj = obs

    def interactive_plot(self, savedir=None,
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
            self._obj.iplot_fname = os.path.join(
                savedir, self._obj.name + '.html')
            save(p, self._obj.iplot_fname, resources=CDN, title=self._obj.name)

        if return_filename:
            return self._obj.iplot_fname
        else:
            return p
