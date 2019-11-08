# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:46:51 2019

@author: oebbe
"""

import test_to_from as ttf

plot_dir = r".\data\2019-Dino-test\plots"

def test_interactive_plot():
    gw = ttf.test_observation_gw()
    gw.to_interactive_plot(savedir=plot_dir, plot_columns=['Stand_m_tov_NAP'],
                           hoover_date_format="{%F}",
                           add_filter_to_legend=True)
    return


def test_obscollection_dino_to_map():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    dino_gw.get_lat_lon()
    dino_gw.to_interactive_map(plot_dir, plot_columns=['Stand_m_tov_NAP'],
                               fname='imap.html',
                               legend_name='grondwater DINO',
                               add_filter_to_legend=True, hoover_names=['gws'],
                               zoom_start=9,
                               verbose=True)
    return

def test_obscollection_dino_to_mapgraph():
    gw = ttf.test_obscollection_dinozip_gw()
    gw.to_mapgraphs(plot_ylim='min_dy')

    return

def test_obscollection_to_map():
    fews_gw_prod = ttf.test_obscollection_fews()
    ax = fews_gw_prod.to_map()
    
    return ax

def test_obscollection_to_imap():
    fname = 'texel_fews.html'
    fews_gw_prod = ttf.test_obscollection_fews()
    fews_gw_prod.get_filternr_locatie('locationId')
    m = fews_gw_prod.to_interactive_map(plot_dir,
                                        plot_columns=['value'],
                                        fname=fname,
                                        plot_freq='D',
                                        legend_name='opp water FEWS',
                                        map_label='locationId',
                                        map_label_size=10)
    return m