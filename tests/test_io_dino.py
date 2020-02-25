# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import sys
sys.path.insert(1, "..")
from observations.io import io_dino
from observations import observation as obs
import os

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
sys.path.insert(0, PROJECT_DIR)
os.chdir(TEST_DIR)


def test_dino_download_single():
    # download single file
    measurements, meta = io_dino.download_dino_groundwater(location="B57F0077",
                                                           filternr=4.,
                                                           tmin="2000-01-01",
                                                           tmax="2010-01-01",
                                                           unit="NAP")
    return measurements, meta

def test_dino_metadata():
    # download metadata
    meta = io_dino.get_dino_piezometer_metadata(location='B52C0057',
                                                filternr='001')
    return meta


def test_dino_metadata2():
    # download metadata without sample data
    meta = io_dino.get_dino_piezometer_metadata(location='B57B0069',
                                                filternr='001', 
                                                verbose=True)
    return meta

def test_dino_metadata3():
    # download metadata of a well without metadata
    meta = io_dino.get_dino_piezometer_metadata(location='B45G1147',
                                                filternr='001', 
                                                verbose=True)
    return meta


def test_dino_meetreeks():
    
    dino = io_dino.DinoWSDL()

    # measurements
    measurements = dino.findMeetreeks(location='B58A0092', 
                                      filternr='004', 
                                      tmin="1900-01-01", 
                                      tmax="2040-01-01")
    
    return measurements

def test_dino_download_extent():
    # download extent
    extent = [120300, 120500, 439000, 441000]  # Schoonhoven zoomed

    import time
    start_time = time.time()
    gw_col = io_dino.download_dino_within_extent(extent,
                                                 ObsClass=obs.GroundwaterObs,
                                                 layer='grondwatermonitoring',
                                                 verbose=True)

    print("--- %s seconds ---" % (time.time() - start_time))
    return gw_col
