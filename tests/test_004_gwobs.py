# -*- coding: utf-8 -*-
"""Created on Tue Mar 24 09:17:55 2020.

@author: oebbe
"""

import numpy as np
import test_001_to_from as ttf


def test_set_tube_nr():
    dino_gw = ttf.test_obscollection_dinozip_gw()
    dino_gw.gwobs.set_tube_nr(if_exists="replace")
    return dino_gw


def test_set_tube_nr_monitoring_well():
    fews_gw_prod = ttf.test_obscollection_fews_lowmemory()
    fews_gw_prod.gwobs.set_tube_nr_monitoring_well(
        "monitoring_well", if_exists="replace"
    )
    return fews_gw_prod


def test_get_modellayers_mf2005():
    import flopy

    modelname = "test_mf2005"
    ml = flopy.modflow.Modflow(modelname, exe_name="mf2005")
    # Model domain and grid definition
    Lx = 300000.0
    Ly = 400000.0
    ztop = 50.0
    zbot = -150.0
    nlay = 4
    nrow = 40
    ncol = 30
    delr = Lx / ncol
    delc = Ly / nrow
    delv = (ztop - zbot) / nlay
    botm = np.linspace(ztop, zbot, nlay + 1)
    # Create the discretization object
    dis = flopy.modflow.ModflowDis(
        ml,
        nlay,
        nrow,
        ncol,
        delr=delr,
        delc=delc,
        top=ztop,
        botm=botm[1:],
        xul=0,
        yul=700000,
    )

    dino_gw = ttf.test_obscollection_dinozip_gw()
    modellayers = dino_gw.gwobs.get_modellayers(ml)

    return modellayers


def test_get_modellayers_mf6_structured():
    import flopy

    # Create the Flopy simulation object
    model_name = "test_mf6_structured"
    sim = flopy.mf6.MFSimulation(sim_name=model_name, exe_name="mf6", version="mf6")

    # Create the Flopy groundwater flow (gwf) model object
    model_nam_file = "{}.nam".format(model_name)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name, model_nam_file=model_nam_file)
    Lx = 300000.0
    Ly = 400000.0
    ztop = 50.0
    zbot = -150.0
    nlay = 4
    nrow = 40
    ncol = 30
    delr = Lx / ncol
    delc = Ly / nrow
    delv = (ztop - zbot) / nlay
    botm = np.linspace(ztop, zbot, nlay + 1)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        xorigin=0,
        yorigin=300000,
        nlay=4,
        nrow=nrow,
        ncol=ncol,
        delr=Lx / ncol,
        delc=Ly / nrow,
        top=ztop,
        botm=botm[1:],
    )
    dino_gw = ttf.test_obscollection_dinozip_gw()
    modellayers = dino_gw.gwobs.get_modellayers(gwf)

    return modellayers


def test_get_regis_layer():
    dino_gw = ttf.test_obscollection_dinozip_gw()

    return dino_gw.gwobs.get_regis_layers()
