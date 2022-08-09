# -*- coding: utf-8 -*-
"""Created on Fri Jan 31 13:26:04 2020.

@author: oebbe
"""

import numpy as np
import pandas as pd
from hydropandas import obs_collection as oc
from hydropandas import observation as obs

# import sys
# sys.path.insert(1, "..")


# TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
# sys.path.insert(0, PROJECT_DIR)
# os.chdir(TEST_DIR)


def test_groundwater_obs(name="grondwaterobs_001", filternr=2):
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"),
        data={"Stand_m_tov_NAP": np.random.rand(10)},
    )
    maaiveld = np.random.random()
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    gwo = obs.GroundwaterObs(
        df,
        name=name,
        locatie=name.split("_")[0],
        x=x,
        y=y,
        maaiveld=maaiveld,
        meetpunt=maaiveld - 0.2,
        onderkant_filter=maaiveld - 10.0,
        bovenkant_filter=maaiveld - 9.0,
        metadata_available=True,
        filternr=filternr,
        filename="",
        meta={"info": "in deze dictionary " "kan je extra informatie kwijt"},
    )
    return gwo


def test_waterlvl_obs():
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"),
        data={"Stand_m_tov_NAP": np.random.rand(10)},
    )
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    wlvl = obs.WaterlvlObs(
        df,
        name="waterlvl_obs1",
        locatie="obs1",
        x=x,
        y=y,
        filename="",
        meta={"info": "in deze dictionary kan je extra informatie kwijt"},
    )
    return wlvl


def test_groundwater_quality_obs():
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"), data={"pH": np.random.rand(10)}
    )
    gwq = obs.WaterlvlObs(
        df,
        name="waterquality_obs1",
        locatie="waterquality",
        x=3,
        y=4,
        filename="",
        meta={"info": "in deze dictionary " "kan je extra informatie kwijt"},
    )
    return gwq


def test_obscollection_from_list():
    o_list = []
    for i in range(10):
        o_list.append(test_groundwater_obs(name=f"grondwaterobs_00{i}", filternr=i))

    obs_col = oc.ObsCollection.from_list(o_list)

    return obs_col


def test_copy_obs():

    o = test_groundwater_obs(name="grondwaterobs_001", filternr=2)
    o2 = o.copy()

    o.meta["hello"] = "world"

    # check deep copy attributes
    assert "hello" not in o2.meta.keys(), "copy method failed"

    o3 = o.copy(deep=False)

    # check shallow copy attributes
    o.meta["answer"] = 42
    assert "answer" in o3.meta.keys(), "copy method failed"


def test_merge_observations_same_timeseries():
    # base
    o = test_groundwater_obs(name="grondwaterobs_010", filternr=10)

    # observation with different metadata, same time series
    o2 = test_groundwater_obs(name="grondwaterobs_010", filternr=10)
    o2.iloc[:, 0] = o.iloc[:, 0]

    omerged = o.merge_observation(o2, merge_metadata=False)
    # check if merged object is identical to first observation
    assert omerged.to_collection_dict() == o.to_collection_dict()


def test_merge_observations_different_timeseries():
    # base
    o = test_groundwater_obs(name="grondwaterobs_010", filternr=10)

    # observation with different time series
    o2 = o.copy()
    o2.index = pd.date_range("2020-1-11", "2020-1-20")
    o2["Stand_m_tov_NAP"] = np.random.rand(10)

    omerged = o.merge_observation(o2)

    assert omerged.shape == (
        20,
        1,
    ), "merged observation should have one column with 20 values"

    return omerged


def test_merge_overlapping():
    # base
    o = test_groundwater_obs(name="grondwaterobs_010", filternr=10)

    # observation with partially overlapping time series and extra columns
    o2 = o.copy()
    o2.index = pd.date_range("2020-1-6", "2020-1-15")
    o2["Stand_m_tov_NAP"].iloc[:5] = o["Stand_m_tov_NAP"].iloc[-5:]
    o2["valid"] = np.random.randint(0, 2, 10)
    o2["Opmerking"] = "test"

    omerged = o.merge_observation(o2)

    assert omerged.shape == (
        15,
        3,
    ), "merged observation should have one column with 20 values"

    return omerged


def test_merge_errors():
    # base
    o = test_groundwater_obs(name="grondwaterobs_010", filternr=10)

    # observation with partially overlapping time series and extra columns
    o2 = test_waterlvl_obs()

    try:
        o.merge_observation(o2)
    except TypeError:
        return 0

    raise RuntimeError("function should raise an error")


def test_add_observation_to_oc():
    oc = test_obscollection_from_list()

    o = test_groundwater_obs(name="grondwaterobs_010", filternr=10)

    oc.add_observation(o)
