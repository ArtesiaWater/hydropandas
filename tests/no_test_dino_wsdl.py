# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 12:07:37 2021

@author: oebbe
"""
# import os
import numpy as np
import pandas as pd
import pytest
from hydropandas import obs_collection as oc
from hydropandas import observation as obs
from hydropandas.io import io_dino


def test_observation_dino_download():
    # download dino
    location = "B57F0077"
    filternr = 4.0
    gw2 = obs.GroundwaterObs.from_dino(
        location=location, filternr=filternr, tmin="2000-01-01", tmax="2010-01-01"
    )
    return gw2


def test_observation_dino_download2():
    # download dino
    gw2 = obs.GroundwaterObs.from_dino(
        location="B57B0069", filternr=1.0, tmin="2000-01-01", tmax="2030-01-01"
    )
    return gw2


def test_observation_dino_download3():
    # download dino data from pb without extra metadata. For this pb
    # io_dino.get_dino_piezometer_metadata() returns an empty list
    location = "B45G1147"
    filternr = 1.0

    gw3 = obs.GroundwaterObs.from_dino(
        location=location, filternr=filternr, tmin="1900-01-01", tmax="1901-01-01"
    )
    return gw3


def test_observation_dino_download_as_cluster():
    # download dino data from pb that belongs to a cluster

    filternr = 1.0

    gw_cl0 = obs.GroundwaterObs.from_dino(
        location="B45E0458", filternr=filternr, split_cluster=False
    )
    gw_cl1 = obs.GroundwaterObs.from_dino(
        location="B45E0459", filternr=filternr, split_cluster=False
    )

    assert gw_cl0.equals(gw_cl1)

    return gw_cl0, gw_cl1


def test_observation_dino_download_split_cluster():
    # download dino data from pb that belongs to a cluster

    filternr = 1.0

    gw_cl0 = obs.GroundwaterObs.from_dino(location="B45E0458", filternr=filternr)
    gw_cl1 = obs.GroundwaterObs.from_dino(location="B45E0459", filternr=filternr)

    assert not gw_cl0.equals(gw_cl1)

    return gw_cl0, gw_cl1


def test_obscollection_dino_download_extent():
    # download DINO from extent
    extent = [117850, 117980, 439550, 439700]  # Schoonhoven zoomed
    dino_gw_extent = oc.ObsCollection.from_dino(
        extent=extent, ObsClass=obs.GroundwaterObs, verbose=True
    )
    return dino_gw_extent


def test_obscollection_dino_download_extent_cache():
    # download DINO from extent
    extent = [117850, 117980, 439550, 439700]  # Schoonhoven zoomed
    dino_gw_extent = oc.ObsCollection.from_dino(
        extent=extent, ObsClass=obs.GroundwaterObs, cache=True, verbose=True
    )
    return dino_gw_extent


def test_obscollection_dino_download_bbox():
    # download DINO from bbox
    bbox = [117850, 439550, 117980, 439700]  # Schoonhoven zoomed
    bbox = np.array([191608.334, 409880.402, 193072.317, 411477.894])
    dino_gw_bbox = oc.ObsCollection.from_dino(
        bbox=bbox, ObsClass=obs.GroundwaterObs, verbose=True
    )
    return dino_gw_bbox


def test_obscollection_dino_download_bbox_only_metadata():
    # check if the keep_all_obs argument works
    bbox = [120110.8948323, 389471.92587313, 121213.23597266, 390551.29918915]
    dino_gw_bbox = oc.ObsCollection.from_dino(bbox=bbox, verbose=True)

    dino_gw_bbox_empty = oc.ObsCollection.from_dino(
        bbox=bbox, keep_all_obs=False, verbose=True
    )
    assert dino_gw_bbox_empty.empty

    return dino_gw_bbox


def test_obscollection_dino_download_bbox_do_not_keep_all_obs():
    bbox = [120110.8948323, 389471.92587313, 121213.23597266, 390551.29918915]
    dino_gw_bbox = oc.ObsCollection.from_dino(bbox=bbox, verbose=True)
    return dino_gw_bbox


def test_dino_download_single():
    # download single file
    measurements, meta = io_dino.download_dino_groundwater(
        location="B57F0077",
        filternr=4.0,
        tmin="2000-01-01",
        tmax="2010-01-01",
        unit="NAP",
    )
    return measurements, meta


def test_dino_meetreeks():

    dino = io_dino.DinoWSDL()

    # measurements
    measurements = dino.findMeetreeks(
        location="B58A0092", filternr="004", tmin="1900-01-01", tmax="2040-01-01"
    )

    return measurements


def test_dino_download_extent():
    # download extent
    extent = [117850, 117980, 439550, 439700]  # Schoonhoven zoomed

    import time

    start_time = time.time()
    gw_col = io_dino.download_dino_within_extent(
        extent, ObsClass=obs.GroundwaterObs, layer="grondwatermonitoring", verbose=True
    )

    print("--- %s seconds ---" % (time.time() - start_time))
    return gw_col


def test_dino_download_single_empty():

    measurements, meta = io_dino.download_dino_groundwater(
        location="B50E0130",
        filternr=1.0,
        tmin="1900-01-01",
        tmax="2040-01-01",
        unit="NAP",
        verbose=True,
    )
    return
