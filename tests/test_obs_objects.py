# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 13:26:04 2020

@author: oebbe
"""

# import os
from observations import observation as obs
from observations import obs_collection as oc
import numpy as np
import pandas as pd
# import sys
# sys.path.insert(1, "..")


# TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
# sys.path.insert(0, PROJECT_DIR)
# os.chdir(TEST_DIR)


def test_groundwater_obs(name='grondwaterobs_001', filternr=2):
    df = pd.DataFrame(index=pd.date_range('2020-1-1', '2020-1-10'),
                      data={'Stand_m_tov_NAP': np.random.rand(10)})
    maaiveld = np.random.random()
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    gwo = obs.GroundwaterObs(df, name=name,
                             locatie=name.split('_')[0],
                             x=x, y=y,
                             maaiveld=maaiveld,
                             meetpunt=maaiveld - 0.2,
                             onderkant_filter=maaiveld - 10.0,
                             bovenkant_filter=maaiveld - 9.0,
                             metadata_available=True,
                             filternr=filternr,
                             filename='',
                             meta={'info': 'in deze dictionary kan je extra informatie kwijt'})
    return gwo


def test_waterlvl_obs():
    df = pd.DataFrame(index=pd.date_range('2020-1-1', '2020-1-10'),
                      data={'Stand_m_tov_NAP': np.random.rand(10)})
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    wlvl = obs.WaterlvlObs(df, name='waterlvl_obs1', locatie='obs1',
                           x=x, y=y, filename='',
                           meta={'info': 'in deze dictionary kan je extra informatie kwijt'})
    return wlvl


def test_groundwater_quality_obs():
    df = pd.DataFrame(index=pd.date_range(
        '2020-1-1', '2020-1-10'), data={'pH': np.random.rand(10)})
    gwq = obs.WaterlvlObs(df, name='waterquality_obs1', locatie='waterquality',
                          x=3, y=4, filename='',
                          meta={'info': 'in deze dictionary kan je extra informatie kwijt'})
    return gwq


def test_obscollection_from_list():
    o_list = []
    for i in range(10):
        o_list.append(test_groundwater_obs(name=f'grondwaterobs_00{i}',
                                           filternr=i))

    obs_col = oc.ObsCollection.from_list(o_list)

    return obs_col
