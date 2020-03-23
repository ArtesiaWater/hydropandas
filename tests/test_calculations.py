# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:48:29 2019

@author: oebbe
"""

import test_to_from as ttf
import flopy
import numpy as np


def test_within_extent():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    extent = [210350, 213300, 473300, 474000]
    dino_gw.geo.within_extent(extent, inplace=True)
    assert dino_gw.shape[0] == 4
    return dino_gw

#%% stats

def test_obscollection_consecutive_obs_years():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    coy = gw.stats.consecutive_obs_years()
    return coy

def test_obscollection_get_number_of_obs():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    coy = gw.stats.get_no_of_observations()
    return coy

def test_obscollection_get_first_last_obs_date():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    fl_obs_date = gw.stats.get_first_last_obs_date()
    return fl_obs_date

def test_obscollection_get_seasonal_stats():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    st = gw.stats.get_seasonal_stat(stat='mean')
    return st

def test_obscollection_get_min():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    omin = gw.stats.get_min()
    return omin

def test_obscollection_get_max():
    gw = ttf.test_obscollection_dinozip_gw_keep_all_obs()
    omax = gw.stats.get_max()
    return omax

#%%
def test_set_filter_num():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    dino_gw.gwobs.set_filter_num(if_exists='replace')
    return dino_gw

def test_set_filter_num_pystore():
    obsc = ttf.test_obscollection_from_pystore()
    obsc.gwobs.set_filter_num(if_exists='replace')
    return obsc

def test_set_filter_num_location():
    fews_gw_prod = ttf.test_obscollection_fews_lowmemory()
    fews_gw_prod.gwobs.set_filter_num_location('locatie',
                                               if_exists='replace')
    return fews_gw_prod

def test_get_modellayers():
    modelname = 'tutorial1'
    ml = flopy.modflow.Modflow(modelname, exe_name='mf2005')
    # Model domain and grid definition
    Lx = 1000000.
    Ly = 1000000.
    ztop = 50.
    zbot = -150.
    nlay = 4
    nrow = 10
    ncol = 10
    delr = Lx / ncol
    delc = Ly / nrow
    delv = (ztop - zbot) / nlay
    botm = np.linspace(ztop, zbot, nlay + 1)
    # Create the discretization object
    dis = flopy.modflow.ModflowDis(ml, nlay, nrow, ncol, delr=delr, delc=delc,
                                   top=ztop, botm=botm[1:])

    dino_gw = ttf.test_obscollection_dinozip_gw()
    modellayers = dino_gw.gwobs.get_modellayers(ml)

    return modellayers


def test_get_nearest_point():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    fl = ttf.test_obscollection_fieldlogger()
    dino_gw[['nearest point', 'distance nearest point']
            ] = dino_gw.geo.get_nearest_point(fl)
    return dino_gw


def test_get_surface_level_oc():
    try:
        from art_tools import obs_extension
        gw = ttf.test_obscollection_fews_lowmemory()
        zp = gw.art.geo_get_surface_level()
        return zp
    except ModuleNotFoundError as e:
        print(e)
        return



def test_get_surface_level_gwobs():
    try:
        from art_tools import obs_extension
        gw = ttf.test_observation_gw()
        mv = gw.art.geo_get_surface_level()
        return mv
    except ModuleNotFoundError as e:
        print(e)
        return
