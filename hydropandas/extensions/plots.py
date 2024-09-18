import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from tqdm.auto import tqdm

from ..observation import GroundwaterObs
from . import accessor

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

    def interactive_plots(
        self,
        savedir="figures",
        tmin=None,
        tmax=None,
        per_monitoring_well=True,
        **kwargs,
    ):
        """Create interactive plots of the observations using bokeh.

        Parameters
        ----------
        savedir : str
            directory used for the folium map and bokeh plots
        tmin : dt.datetime, optional
            start date for timeseries plot
        tmax : dt.datetime, optional
            end date for timeseries plot
        per_monitoring_well : bool, optional
            if True plot multiple tubes at the same monitoring_well in one
            figure
        **kwargs :
            will be passed to the Obs.interactive_plot method, options
            include:

            - cols : list of str or None
            - hoover_names : list of str
            - plot_freq : list of str
            - plot_legend_names : list of str
            - markers : list of str
            - hoover_names : list of str
            - plot_colors : list of str
            - ylabel : str
            - add_screen_to_legend : boolean
        """
        _color_cycle = (
            "blue",
            "olive",
            "lime",
            "red",
            "orange",
            "yellow",
            "purple",
            "silver",
            "powderblue",
            "salmon",
            "tan",
        )

        # check if observations consist of monitoring wells
        if per_monitoring_well:
            otype = self._obj._infer_otype()
            if isinstance(otype, (list, np.ndarray)):
                per_monitoring_well = False
            elif otype.__name__ == "GroundwaterObs":
                pass
            else:
                per_monitoring_well = False

        if per_monitoring_well:
            plot_names = self._obj.groupby("monitoring_well").count().index
        else:
            plot_names = self._obj.index

        for name in plot_names:
            if per_monitoring_well:
                oc = self._obj.loc[
                    self._obj.monitoring_well == name, "obs"
                ].sort_index()
            else:
                oc = self._obj.loc[[name], "obs"]

            p = None
            for i, o in enumerate(oc.values):
                if i == 10:
                    raise NotImplementedError(
                        "cannot add more than 10 lines to a single plot"
                    )
                try:
                    p = o.plots.interactive_plot(
                        savedir=savedir,
                        p=p,
                        tmin=tmin,
                        tmax=tmax,
                        plot_colors=[_color_cycle[i + 1]],
                        return_filename=False,
                        **kwargs,
                    )
                    logger.debug(f"created iplot -> {o.name}")
                except ValueError:
                    logger.error(f"{o.name} has no data between {tmin} and {tmax}")
                    o.meta["iplot_fname"] = None

    def interactive_map(
        self,
        plot_dir="figures",
        m=None,
        tiles="OpenStreetMap",
        fname=None,
        per_monitoring_well=True,
        color="blue",
        legend_name=None,
        add_legend=True,
        map_label="",
        map_label_size=20,
        col_name_lat="lat",
        col_name_lon="lon",
        zoom_start=13,
        create_interactive_plots=True,
        **kwargs,
    ):
        """Create an interactive map with interactive plots using folium and
        bokeh.

        Notes
        -----
        Some notes on this method:

        - if you want to have multiple obs collections on one folium map, only
          the last one should have add_legend = True to create a correct legend
        - the color of the observation point on the map is now the same color
          as the line of the observation measurements. Also a built-in color
          cycle is used for different measurements at the same monitoring_well.

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
        per_monitoring_well : bool, optional
            if True plot multiple tubes at the same monitoring well in one
            figure
        color : str, optional
            color of the observation points on the map
        legend_name : str, optional
            the name of the observation points shown in the map legend
        add_legend : boolean, optional
            add a legend to a plot
        map_label : str, optional
            add a label to the monitoring wells on the map, the label should be
            1. the attribute of an observation
            2. the key in the meta attribute of the observation
            3. a generic label for each observation in this collection.
            A label is only added if map_label is not ''. The default is ''.
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
            in the meta ditctionary of the observations is used.
        **kwargs :
            will be passed to the interactive_plots method options are:

            - cols : list of str or None
            - hoover_names : list of str
            - plot_legend_names : list of str
            - plot_freq : list of str
            - markers : list of str
            - hoover_names : list of str
            - plot_colors : list of str
            - ylabel : str
            - add_screen_to_legend : boolean
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

        # check for empty observations
        if all([o.empty for o in self._obj.obs.values]):
            logger.warning("all observations in the collection are empty")
            for oname in self._obj.index:
                self._obj._set_metadata_value(
                    oname, "iplot_fname", None, add_to_meta=True
                )
            empty_obs = True
        else:
            empty_obs = False

        # create interactive bokeh plots
        if create_interactive_plots and not empty_obs:
            self._obj.plots.interactive_plots(
                savedir=plot_dir, per_monitoring_well=per_monitoring_well, **kwargs
            )

        # check if observation collection has lat and lon values
        if (col_name_lat not in self._obj.columns) and (
            col_name_lon not in self._obj.columns
        ):
            self._obj.geo.set_lat_lon()

        # determine start location of map
        northing = np.mean(
            (self._obj[col_name_lat].min(), self._obj[col_name_lat].max())
        )
        easting = np.mean(
            (self._obj[col_name_lon].min(), self._obj[col_name_lon].max())
        )

        # create map if no map is given
        if m is None:
            m = folium.Map([northing, easting], zoom_start=zoom_start, tiles=tiles)

        # get oc name if no legend name is given
        if legend_name is None:
            legend_name = self._obj.name

        # add the point observations with plots to the map
        group_name = '<span style=\\"color: {};\\">{}</span>'.format(color, legend_name)
        group = folium.FeatureGroup(name=group_name)

        # check if observations consist of monitoring wells
        if per_monitoring_well:
            otype = self._obj._infer_otype()
            if isinstance(otype, (list, np.ndarray)):
                per_monitoring_well = False
            elif otype.__name__ == "GroundwaterObs":
                pass
            else:
                per_monitoring_well = False

        if per_monitoring_well:
            plot_names = self._obj.groupby("monitoring_well").count().index
        else:
            plot_names = self._obj.index

        for name in plot_names:
            if per_monitoring_well:
                oc = self._obj.loc[
                    self._obj.monitoring_well == name, "obs"
                ].sort_index()
                o = oc.iloc[-1]
                name = o.name
            else:
                o = self._obj.loc[name, "obs"]

            if o.meta["iplot_fname"] is not None:
                with open(o.meta["iplot_fname"], "r") as f:
                    bokeh_html = f.read()

                iframe = branca.element.IFrame(html=bokeh_html, width=620, height=420)
                popup = folium.Popup(iframe, max_width=620)

                folium.CircleMarker(
                    [
                        self._obj.loc[o.name, col_name_lat],
                        self._obj.loc[o.name, col_name_lon],
                    ],
                    icon=folium.Icon(icon="signal"),
                    fill=True,
                    color=color,
                    popup=popup,
                ).add_to(group)

                if map_label != "":
                    if map_label in o._metadata:
                        map_label_val = getattr(o, map_label)
                    elif map_label in o.meta.keys():
                        map_label_val = o.meta[map_label]
                    else:
                        map_label_val = map_label
                    folium.map.Marker(
                        [
                            self._obj.loc[name, col_name_lat],
                            self._obj.loc[name, col_name_lon],
                        ],
                        icon=DivIcon(
                            icon_size=(150, 36),
                            icon_anchor=(0, 0),
                            html='<div style="font-size: %ipt">%s</div>'
                            % (map_label_size, map_label_val),
                        ),
                    ).add_to(group)
            else:
                logger.info(f"no iplot available for {o.name}")
                folium.CircleMarker(
                    [
                        self._obj.loc[o.name, col_name_lat],
                        self._obj.loc[o.name, col_name_lon],
                    ],
                    icon=folium.Icon(icon="signal"),
                    fill=True,
                    color=color,
                ).add_to(group)

        group.add_to(m)

        # add legend
        if add_legend:
            folium.map.LayerControl("topright", collapsed=False).add_to(m)

        # save map
        # filename and path
        if fname is not None:
            if not fname.endswith(".html"):
                fname = fname + ".html"
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)
            m.save(os.path.join(plot_dir, fname))

        return m

    def section_plot(
        self,
        tmin=None,
        tmax=None,
        cols=(None,),
        section_colname_x=None,
        section_label_x=None,
        section_well_layout_color="gray",
        section_markersize=100,
        ylabel="auto",
        fn_save=None,
        check_obs_close_to_screen_bottom=True,
        plot_well_layout_markers=True,
        plot_obs=True,
    ):
        """Create plot with well layout (left) en observations (right).

        Parameters
        ----------
        tmin : dt.datetime, optional
            start date for timeseries plot
        tmax : dt.datetime, optional
            end date for timeseries plot
        cols : tuple of str or None, optional
                the columns of the observation to plot. The first numeric column
                is used if cols is None, by default None.
        section_colname_x : str, optional
            column used for x position on section plot, when None order
            collection is used
        section_label_x : str, optional
            label applied to x-axis in section plot
        section_well_layout_color : str, optional
            color of well layout, default is gray
        section_markersize : int, optional
            size of makers in sectionplot
        ylabel: str or list, optional
            when 'auto' column unit in collection is ylabel, otherwise first
            element of list is label of section plot, second element of observation plot
        fn_save : str, optional
            filename to save plot
        check_obs_close_to_screen_bottom : bool, optional
            plots a horizontal line when minimum observation is close to screen_bottom
        plot_well_layout_markers : bool, optional
            plots ground level, top tube, screen levels and sandtrap via
            makers. Default is True
        plot_obs : bool, optional
            Plots observation. Default is True

        TODO:
            - speficy colors via extra column in ObsCollection
            - addtional visual checks:
                maximum observation is close to or above ground level,
                maximum observation is close to or above tube top
                minimum observation is close to or below tube bottom (sand trap)
            - include some interactive Bokeh fancy
            - apply an offset when two wells are at same location
            - limit y-axis of section plot to observations only
            - remove the checking (if obs are near bottom) from this function
            - moving the legend outside the plot
            - set xlim of observation plot more tight when tmin is not specified
        """

        # prepare column for x location in section plot
        if section_colname_x is None:
            # use order in ObsCollection
            plot_x = np.arange(len(self._obj))
        else:
            plot_x = self._obj[section_colname_x].to_numpy()

        # create figure
        fig = plt.figure(figsize=(15, 5))
        if plot_obs:
            # make figure with section plot on left, and righthand plot for obs
            gs = GridSpec(
                1,
                2,
                width_ratios=[1, 3],
            )
            ax_section = fig.add_subplot(gs[0])
            ax_obs = fig.add_subplot(gs[1])
            axes = [ax_section, ax_obs]
        else:
            # make figure with only one plot on left, observations are not plotted
            ax_section = fig.add_subplot(111)
            axes = [ax_section]

        if plot_well_layout_markers:
            # plot well layout via markers
            ax_section.scatter(
                plot_x,
                self._obj.tube_top.values,
                section_markersize,
                label="tube top",
                marker="*",
                facecolors="none",
                color=section_well_layout_color,
            )
            ax_section.scatter(
                plot_x,
                self._obj.ground_level.values,
                section_markersize,
                label="ground level",
                marker="_",
                color=section_well_layout_color,
            )
            ax_section.scatter(
                plot_x,
                self._obj.screen_top.values,
                section_markersize / 2,
                label="screen top",
                marker="x",
                color=section_well_layout_color,
            )
            ax_section.scatter(
                plot_x,
                self._obj.screen_bottom.values,
                section_markersize,
                label="screen bottom",
                marker="+",
                color=section_well_layout_color,
            )

        # loop over all wells, plot observations and details in section plot
        for counter, name in enumerate(self._obj.index):
            if plot_obs:
                # which column to plot?
                cols = list(cols)
                for i, col in enumerate(cols):
                    if col is None:
                        cols[i] = self._obj.loc[
                            name, "obs"
                        ]._get_first_numeric_col_name()

                # create plot dataframe
                plot_df = self._obj.loc[name, "obs"][tmin:tmax][cols].copy()
                if plot_df.empty or plot_df[cols].isna().all().all():
                    logger.warning(f"{name} has no data between {tmin} and {tmax}")
                    continue

                # PART 1: plot timeseries of observations
                # one or multiple columns to plot?
                if len(cols) == 1:
                    p = ax_obs.plot(plot_df, label=name)
                else:
                    for col in cols:
                        p = ax_obs.plot(plot_df[col], label=f"{name}, {col}")
                plot_color = p[0].get_color()

                if check_obs_close_to_screen_bottom:
                    # add horizonal line to plot when minimum observation in first plot
                    # column is close to bottom of screen
                    offset = 0.1
                    if self._obj.loc[name, "screen_bottom"] > (
                        plot_df[cols[0]].dropna().min() - offset
                    ):
                        ax_obs.axhline(
                            y=self._obj.loc[name, "screen_bottom"],
                            ls="--",
                            lw=4,
                            alpha=0.5,
                            label=(
                                f"screen bottom of {name}\n"
                                "is close to minimum observation"
                            ),
                            color=plot_color,
                        )
            else:
                plot_color = "hotpink"

            # PART 2: fancy section plot with lines along tube

            # highlight filter on section plot
            ax_section.plot(
                [plot_x[counter]] * 2,
                [
                    self._obj.loc[name, "screen_top"],
                    self._obj.loc[name, "screen_bottom"],
                ],
                color="k",
                lw=3,
                ls="-",
            )

            # highlight blind tube on section plot
            ax_section.plot(
                [plot_x[counter]] * 2,
                [self._obj.loc[name, "screen_top"], self._obj.loc[name, "tube_top"]],
                color=plot_color,
                lw=3,
            )

            # add sandtrap when present
            if "tube_bottom" in self._obj.columns:
                ax_section.plot(
                    [plot_x[counter]] * 2,
                    [
                        self._obj.loc[name, "screen_bottom"],
                        self._obj.loc[name, "tube_bottom"],
                    ],
                    color=plot_color,
                    lw=3,
                )

            if plot_obs:
                # PART 3: fancy section plot with bandwith of observations
                ax_section.scatter(
                    plot_x[counter],
                    plot_df.quantile(q=0.95),
                    section_markersize,
                    marker=(3, 0, 0),
                    color="blue",
                    alpha=0.5,
                    label="95% observation",
                    zorder=100,
                )
                ax_section.scatter(
                    plot_x[counter],
                    plot_df.median(),
                    section_markersize,
                    marker=(4, 0, 90),
                    color="green",
                    alpha=0.5,
                    label="median",
                    zorder=100,
                )
                ax_section.scatter(
                    plot_x[counter],
                    plot_df.quantile(q=0.05),
                    section_markersize,
                    marker=(3, 0, 180),
                    color="red",
                    alpha=0.5,
                    label="5% observation",
                    zorder=100,
                )

            logger.info(f"created sectionplot -> {name}")

        # layout
        if section_label_x is None:
            # use name as xtick and rotate
            ax_section.set_xticks(
                plot_x,
                self._obj.index,
                rotation="vertical",
                fontsize="small",
            )
        else:
            ax_section.set_xlabel(section_label_x)

        if plot_obs:
            ax_obs.set_xlim(left=tmin, right=tmax)

            # rotate labels on observation axis
            ax_obs.set_xticks(
                ax_obs.get_xticks(),
                ax_obs.get_xticklabels(),
                rotation="vertical",
                fontsize="small",
            )

        if ylabel == "auto":
            # has collection uniform unit?
            if len(self._obj.unit.unique()) == 1:
                ylabel = self._obj.unit.unique()[0]
            else:
                logger.error(
                    "Collection has more than one unit. Plot has both on one y-axis."
                )
                ylabel = ", ".join(map(str, self._obj.unit.unique()))
            for ax in axes:
                ax.set_ylabel(ylabel)
        else:
            try:
                ax_section.set_ylabel(ylabel[0])
                ax_obs.set_ylabel(ylabel[1])
            except ValueError:
                logger.error(f"Invalid value for ylabel {ylabel}. Plot has no ylabels.")

        # add layout to both plots
        for ax in axes:
            ax.grid()

            # simple legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        fig.suptitle(self._obj.name)

        if fn_save is not None:
            try:
                fig.tight_layout()
                fig.savefig(
                    fn_save,
                    dpi=250,
                )
            except ValueError:
                logger.error(f"Save of figure {name} failed ({fn_save}).")

        return fig, axes

    def series_per_group(
        self,
        plot_column,
        by=None,
        savefig=True,
        outputdir=".",
        naming_method=None,
        units_for_well_screen_depth="mNAP",
    ):
        """Plot time series per group.

        The default groupby is based on identical x, y coordinates, so plots unique
        time series per location.

        Parameters
        ----------
        plot_column : str
            name of column containing time series data
        by : (list of) str or (list of) array-like
            groupby parameters, default is None which sets groupby to
            columns ["x", "y"].
        savefig : bool, optional
            save figures, by default True, if False returns axes handles to plots.
        outputdir : str, optional
            path to output directory, by default the current directory (".")
        naming_method : str, optional
            method to determine file names for plots, default is None, which uses the
            following template:
                "series_by_{groupby columns}_group_{groupby values}.png".
            If a string is passed, it is interpreted as a column name and the filename
            contains the unique values from that column per group:
                "series_group_{unique values in column}.png"
            The final option includes "infer_name_monitoring_well", which attempts to
            obtain the name of the monitoring well from the final name in the group,
            by splitting the string on "-" or "_", and using the first part.
            For example, if the name is "PB001-001", the filename will be "PB001.png".
        units_for_well_screen_depth : str, optional
            if the observations are GroundWaterObs, the units of the screen depth are
            added to the legend. This keyword argument defines the units for the screen
            depth. Default is 'mNAP'.

        Returns
        -------
        axes : list
            returns list of axes handles if savefig=False
        """
        if by is None:
            by = ["x", "y"]
        gr = self._obj.groupby(by=by)

        axes = []

        for groupname, group in tqdm(
            gr, desc="Plotting series per group", total=len(gr)
        ):
            f, ax = plt.subplots(1, 1, figsize=(10, 3))
            for name, row in group.iterrows():
                if isinstance(row.obs, GroundwaterObs):
                    if units_for_well_screen_depth == "mNAP":
                        lbl = (
                            f"{name} (NAP{row['screen_top']:+.1f}"
                            f"-{row['screen_bottom']:+.1f} m)"
                        )
                    else:
                        lbl = (
                            f"{name} ({row['screen_top']:.1f}"
                            f"-{row['screen_bottom']:.1f} "
                            f"{units_for_well_screen_depth})"
                        )
                else:
                    lbl = f"{name}"
                ax.plot(
                    row.obs.index,
                    row.obs[plot_column],
                    label=lbl,
                )
            ax.legend(
                loc=(0, 1),
                frameon=False,
                ncol=min(group.index.size, 3),
                fontsize="x-small",
            )
            ax.set_ylabel(row["unit"])
            ax.grid(True)
            f.tight_layout()

            if savefig:
                if isinstance(by, list):
                    by_name = "-".join(by)
                    groupname = "-".join(groupname)
                else:
                    by_name = by
                if naming_method is None:
                    filename = f"series_by_{by_name}_group_{groupname}.png"
                elif naming_method == "infer_name_monitoring_well":
                    if "-" in name:
                        filename = f"{name.split('-')[0]}.png"
                    elif "_" in name:
                        filename = f"{name.split('_')[0]}.png"
                    else:
                        filename = f"{name}.png"
                elif isinstance(naming_method, str):
                    filename = (
                        f"series_by_{by_name}_group_"
                        f"{'-'.join(group[naming_method].unique().tolist())}.png"
                    )

                f.savefig(
                    os.path.join(outputdir, filename), bbox_inches="tight", dpi=150
                )
                plt.close(f)
            else:
                axes.append(ax)

        if not savefig:
            return axes


@accessor.register_obs_accessor("plots")
class ObsPlots:
    def __init__(self, obs):
        self._obj = obs

    def interactive_plot(
        self,
        savedir=None,
        cols=(None,),
        markers=("line",),
        p=None,
        plot_legend_names=("",),
        plot_freq=(None,),
        tmin=None,
        tmax=None,
        hoover_names=("Peil",),
        hoover_date_format="%Y-%m-%d",
        ylabel=None,
        plot_colors=("blue",),
        add_screen_to_legend=False,
        return_filename=False,
    ):
        """Create an interactive plot of the observations using bokeh.

        Todo:

        - add options for hoovers, markers, linestyle

        Parameters
        ----------
        savedir : str, optional
            directory used for the folium map and bokeh plots
        cols : tuple of str or None, optional
            the columns of the observation to plot. The first numeric column
            is used if cols is None, by default None.
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
            names will be displayed together with the cols values
            when hoovering over plot
        hoover_date_format : str, optional
            date format to use when hoovering over a plot
        ylabel : str or None, optional
            label on the y-axis. If None the unit attribute of the observation
            is used.
        plot_colors : list of str, optional
            plot_colors used for the plots
        add_screen_to_legend : boolean, optional
            if True the attributes screen_top and screen_bottom
            are added to the legend name
        return_filename : boolean, optional
            if True filename will be returned

        Returns
        -------
        fname_plot : str or bokeh plot
            filename of the bokeh plot or reference to bokeh plot
        """

        from bokeh.models import ColumnDataSource, HoverTool
        from bokeh.plotting import figure, save
        from bokeh.resources import CDN

        cols = list(cols)
        for i, col in enumerate(cols):
            if col is None:
                cols[i] = self._obj._get_first_numeric_col_name()

        # create plot dataframe
        plot_df = self._obj[tmin:tmax].copy()
        plot_df["date"] = plot_df.index.strftime(hoover_date_format)
        if plot_df.empty or plot_df[cols].isna().all().all():
            raise ValueError(
                "{} has no data between {} and {}".format(self._obj.name, tmin, tmax)
            )

        # create plot
        if p is None:
            p = figure(width=600, height=400, x_axis_type="datetime", title="")
            if ylabel is None:
                ylabel = self._obj.unit
            p.yaxis.axis_label = ylabel

        # get x axis
        xcol = self._obj.index.name
        if xcol is None:
            xcol = "index"

        # get color
        if len(plot_colors) < len(cols):
            plot_colors = list(plot_colors) * len(cols)

        # get base for hoover tooltips
        plots = []
        tooltips = []
        tooltips.append(("date", "@date"))

        # plot multiple columns
        for i, column in enumerate(cols):
            # legend name
            if add_screen_to_legend:
                lname = "{} {} (NAP {:.2f} - {:.2f})".format(
                    plot_legend_names[i],
                    self._obj.name,
                    self._obj.screen_bottom,
                    self._obj.screen_top,
                )
            else:
                lname = "{} {}".format(plot_legend_names[i], self._obj.name)

            # resample data
            if plot_freq[i] is None:
                source = ColumnDataSource(plot_df[[column, "date"]])
            else:
                source = ColumnDataSource(
                    plot_df[[column, "date"]].resample(plot_freq[i]).first()
                )

            # plot data

            if markers[i] in ["line", "l"]:
                plots.append(
                    p.line(
                        xcol,
                        column,
                        source=source,
                        color=plot_colors[i],
                        legend_label=lname,
                        alpha=0.8,
                        muted_alpha=0.2,
                    )
                )
            elif markers[i] in ["circle", "c"]:
                plots.append(
                    p.circle(
                        xcol,
                        column,
                        source=source,
                        color=plot_colors[i],
                        legend_label=lname,
                        alpha=0.8,
                        muted_alpha=0.2,
                    )
                )
            else:
                raise NotImplementedError(
                    "marker '{}' invalid. Only line and"
                    "circle are currently available".format(markers[i])
                )

            # add columns to hoover tooltips
            tooltips_p = tooltips.copy()
            tooltips_p.append((hoover_names[i], "@{}".format(column)))
            hover = HoverTool(renderers=[plots[i]], tooltips=tooltips_p, mode="vline")
            p.add_tools(hover)

        p.legend.location = "top_left"
        p.legend.click_policy = "mute"

        # save plot
        if savedir is not None:
            if not os.path.isdir(savedir):
                os.makedirs(savedir)
            self._obj.meta["iplot_fname"] = os.path.join(
                savedir, self._obj.name + ".html"
            )
            save(p, self._obj.meta["iplot_fname"], resources=CDN, title=self._obj.name)

        if return_filename:
            return self._obj.meta["iplot_fname"]
        else:
            return p
