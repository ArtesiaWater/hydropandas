# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:48:29 2019

@author: oebbe
"""

import test_to_from as ttf


def test_within_extent():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    extent = [210350, 213300, 473300, 474000]
    dino_gw.geo.within_extent(extent, inplace=True)
    assert dino_gw.shape[0] == 4
    return dino_gw


def test_obscollection_consecutive_obs_years():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    coy = gw.stats.consecutive_obs_years()
    return coy


def test_obscollection_get_seasonal_stats():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    st = gw.stats.get_seasonal_stat(stat='mean')
    return st


def test_get_filter_num():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    dino_gw.gwobs.get_filter_num(if_exists='replace')
    return dino_gw


def test_get_filter_num_location():
    fews_gw_prod = ttf.test_obscollection_fews()
    fews_gw_prod.add_meta_to_df('locationId')
    fews_gw_prod.gwobs.get_filter_num_location('locationId',
                                               if_exists='replace')
    return fews_gw_prod


def test_get_nearest_point():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    fl = ttf.test_obscollection_fieldlogger()
    dino_gw[['nearest point', 'distance nearest point']
            ] = dino_gw.geo.get_nearest_point(fl)
    return dino_gw


def test_get_surface_level_oc():
    gw = ttf.test_obscollection_fews()
    zp = gw.geo.get_surface_level()
    return zp


def test_get_surface_level_gwobs():
    gw = ttf.test_observation_gw()
    mv = gw.geo.get_surface_level()
    return mv
