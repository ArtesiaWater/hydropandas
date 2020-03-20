# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:46:51 2019

@author: oebbe
"""

import test_to_from as ttf

plot_dir = r".\data\2019-Dino-test\plots"


def test_interactive_plot():
    gw = ttf.test_observation_gw()
    gw.plots.interactive_plot(savedir=plot_dir, plot_columns=['stand_m_tov_nap'],
                              hoover_date_format="{%F}",
                              add_filter_to_legend=True)
    return


def test_obscollection_dino_to_imap():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    dino_gw.geo.set_lat_lon(verbose=True)
    dino_gw.plots.interactive_map(plot_dir, plot_columns=['stand_m_tov_nap'],
                                  fname='imap.html',
                                  legend_name='grondwater DINO',
                                  add_filter_to_legend=True, hoover_names=['gws'],
                                  zoom_start=9,
                                  verbose=True)
    return


def test_obscollection_dino_to_mapgraph():
    try:
        from art_tools import obs_extension
        gw = ttf.test_obscollection_dinozip_gw()
        gw.art.plot_mapgraphs(plot_ylim='min_dy')
        return
    except ModuleNotFoundError as e:
        print(e)
        return


def test_obscollection_to_map():
    try:
        from art_tools import obs_extension
        fews_gw_prod = ttf.test_obscollection_fews_lowmemory()
        ax = fews_gw_prod.art.plot_mapfig()
        return ax
    except ModuleNotFoundError as e:
        print(e)
        return


def test_obscollection_to_imap():
    fname = 'texel_fews.html'
    fews_gw_prod = ttf.test_obscollection_fews_lowmemory()
    # add metadata to obscollection DF
    fews_gw_prod.add_meta_to_df("lat")
    fews_gw_prod.add_meta_to_df("lon")

    fews_gw_prod.gwobs.set_filter_num_location('locatie',
                                               if_exists="replace")

    m = fews_gw_prod.plots.interactive_map(plot_dir,
                                           plot_columns=['value'],
                                           fname=fname,
                                           plot_freq='D',
                                           legend_name='opp water FEWS',
                                           map_label='locationId',
                                           map_label_size=10)
    return m
