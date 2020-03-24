# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:17:55 2020

@author: oebbe
"""

import test_to_from as ttf
import flopy
import numpy as np


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
