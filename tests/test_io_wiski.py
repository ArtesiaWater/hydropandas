# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import sys
sys.path.insert(1, "..")
import pandas as pd
from observations import observation as obs
from observations.io import io_wiski
import os

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
sys.path.insert(0, PROJECT_DIR)
os.chdir(TEST_DIR)


def test_read_wiski_csv():
    # download single file

    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv",
                                            sep='\s+', header_sep=':',
                                            header_identifier=':', verbose=True,
                                            parse_dates={"datetime": [0, 1]},
                                            infer_datetime_format=True,
                                            index_col=["datetime"],
                                            translate_dic={'name': 'Station Number',
                                                           'x': 'GlobalX',
                                                           'y': 'GlobalY'})

    return header, data


def test_read_wiski_csv2():
    # download single file

    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\8137_PBF.csv",
                                            sep='\s+', header_sep=':',
                                            header_identifier=':', verbose=True,
                                            parse_dates={"datetime": [0, 1]},
                                            infer_datetime_format=False,
                                            dayfirst=True,
                                            index_col=["datetime"],
                                            translate_dic={'name': 'Station Number',
                                                           'x': 'GlobalX',
                                                           'y': 'GlobalY'})

    return header, data


# %%

def test_read_wiski_zip():
    obs_df = io_wiski.read_wiski_dir(r".\data\2019-WISKI-test\1016_PBF.zip",
                                     ObsClass=obs.GroundwaterObs,
                                     sep='\s+', header_sep=':',
                                     header_identifier=':',
                                     parse_dates={"datetime": [0, 1]},
                                     index_col=["datetime"],
                                     translate_dic={'name': 'Station Number',
                                                    'x': 'GlobalX',
                                                    'y': 'GlobalY'},
                                     verbose=True)

    return obs_df


def test_rijnenijssel_wiski_format():
    o = obs.GroundwaterObs.from_wiski(
        r"./data/2019-WISKI-test/Zwiepse Horstweg Barchem_1024_FT1_WNS9040_MomentaanO.csv",
        header_sep=";", end_header_str="#Timestamp", parse_dates=[0], index_col=[0],
        infer_datetime_format=True, tz_localize=False, to_mnap=False)
    return o
