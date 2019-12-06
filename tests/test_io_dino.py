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
    measurements, meta = io_dino.download_dino_groundwater(name="B57F0077",
                                                           filternr=4.,
                                                           tmin="2000-01-01",
                                                           tmax="2010-01-01",
                                                           unit="NAP")
    return measurements, meta


# def test_dino_download_extent():
#     # download extent
#     extent = [120300, 120500, 439000, 441000]  # Schoonhoven zoomed

#     import time
#     start_time = time.time()
#     gw_col = io_dino.download_dino_within_extent(extent, ObsClass=obs.GroundwaterObs,
#                                                  kind='Grondwaterstand',
#                                                  verbose=True)

#     print("--- %s seconds ---" % (time.time() - start_time))
#     return gw_col
