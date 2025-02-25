import numpy as np
import pandas as pd
import pytest

import hydropandas as hpd


def _get_groundwater_obs(name="groundwaterobs_001", tube_nr=2):
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"),
        data={"values": np.random.rand(10)},
    )
    ground_level = np.random.random()
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    o = hpd.GroundwaterObs(
        df,
        name=name,
        location=name.split("_")[0],
        x=x,
        y=y,
        source="generated",
        unit="m NAP",
        ground_level=ground_level,
        tube_top=ground_level - 0.2,
        screen_bottom=ground_level - 10.0,
        screen_top=ground_level - 9.0,
        metadata_available=True,
        tube_nr=tube_nr,
        filename="",
        meta={"info": "you can store additional information in this dictionary"},
    )
    return o


def _get_waterlvl_obs():
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"),
        data={"values": np.random.rand(10)},
    )
    x = np.random.randint(0, 10000)
    y = np.random.randint(10000, 20000)
    o = hpd.WaterlvlObs(
        df,
        name="waterlvl_obs1",
        location="obs1",
        x=x,
        y=y,
        filename="",
        meta={"info": "you can store additional information in this dictionary"},
    )
    return o


def _obscollection_from_list():
    o_list = []
    for i in range(10):
        o_list.append(_get_groundwater_obs(name=f"groundwaterobs_00{i}", tube_nr=i))

    oc = hpd.ObsCollection.from_list(o_list)

    return oc


def test_groundwater_quality_obs():
    df = pd.DataFrame(
        index=pd.date_range("2020-1-1", "2020-1-10"), data={"pH": np.random.rand(10)}
    )
    hpd.WaterlvlObs(
        df,
        name="waterquality_obs1",
        location="waterquality",
        x=3,
        y=4,
        filename="",
        meta={"info": "you can store additional information in this dictionary"},
    )


def test_add_meta_to_df():
    oc = _obscollection_from_list()
    oc = oc.add_meta_to_df(key="all")

    assert "info" in oc.columns, "unexpected result for add_meta_to_df"


def test_copy_obs():
    o = _get_groundwater_obs(name="groundwaterobs_001", tube_nr=2)
    o2 = o.copy()

    o.meta["hello"] = "world"

    # check deep copy attributes
    assert "hello" not in o2.meta.keys(), "copy method failed"

    o3 = o.copy(deep=False)

    # check shallow copy attributes
    o.meta["answer"] = 42
    assert "answer" in o3.meta.keys(), "copy method failed"


def test_returns():
    # check if a DataFrame is returned when an ObsCollection is sliced without the
    # 'obs' column
    oc = _obscollection_from_list()

    assert isinstance(oc.loc[:, ["x", "y"]], pd.DataFrame)
    assert not isinstance(oc.loc[:, ["x", "y"]], hpd.ObsCollection)

    assert isinstance(oc.loc[:, ["x", "y", "obs"]], hpd.ObsCollection)


def test_convert_waterlvl_groundwater_obs():
    # create WaterlvlObs
    df = pd.DataFrame(
        {"measurements": np.random.randint(0, 10, 5)},
        index=pd.date_range("2020-1-1", "2020-1-5"),
    )
    o_wl = hpd.WaterlvlObs(
        df,
        name="obs",
        x=54.37326,
        y=-5.57900,
        source="my fantasy",
        location="Weirwood tree",
        meta={"place": "Winterfell"},
    )

    # This is what I want to do, but now I will lose all metadata
    o_gw = hpd.GroundwaterObs(o_wl, ground_level=200)

    assert o_wl.location == o_gw.location, "conversion failed"
    assert o_gw.ground_level == 200, "conversion failed"


def test_merge_observations_same_timeseries():
    # base
    o = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)

    # observation with different metadata, same time series
    o2 = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)
    o2.iloc[:, 0] = o.iloc[:, 0]

    omerged = o.merge_observation(o2, merge_metadata=False)
    # check if merged object is identical to first observation
    assert omerged.to_collection_dict() == o.to_collection_dict()


def test_merge_observations_different_timeseries():
    # base
    o = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)

    # observation with different time series
    o2 = o.copy()
    o2.index = pd.date_range("2020-1-11", "2020-1-20")
    o2["values"] = np.random.rand(10)

    omerged = o.merge_observation(o2)

    assert omerged.shape == (
        20,
        1,
    ), "merged observation should have one column with 20 values"


def test_merge_overlapping():
    # base
    o = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)

    # observation with partially overlapping time series and extra columns
    o2 = o.copy()
    o2.index = pd.date_range("2020-1-6", "2020-1-15")
    o2["values"].iloc[:5] = o["values"].iloc[-5:]
    o2["valid"] = np.random.randint(0, 2, 10)
    o2["qualifier"] = "test"

    omerged = o.merge_observation(o2)

    assert omerged.shape == (
        15,
        3,
    ), "merged observation should have one column with 20 values"


def test_merge_errors():
    # base
    o = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)

    # observation with partially overlapping time series and extra columns
    o2 = _get_waterlvl_obs()

    with pytest.raises(TypeError):
        o.merge_observation(o2)


def test_add_observation_to_oc():
    oc = _obscollection_from_list()

    o = _get_groundwater_obs(name="groundwaterobs_010", tube_nr=10)

    oc.add_observation(o)


def test_interpolate_obscollection():
    oc = _obscollection_from_list()

    xy = [[500, 11000], [9000, 18000]]
    oc.interpolate(xy)


def test_get_obs():
    oc = _obscollection_from_list()

    # by name
    o = oc.get_obs(name="groundwaterobs_001")
    assert isinstance(o, hpd.GroundwaterObs)
    assert o.name == "groundwaterobs_001"

    # by attributes
    o = oc.get_obs(location="groundwaterobs", tube_nr=2)
    assert isinstance(o, hpd.GroundwaterObs)
    assert o.tube_nr == 2

    # multiple observations
    with pytest.raises(ValueError):
        oc.get_obs(location="groundwaterobs")

    # no observations
    with pytest.raises(ValueError):
        oc.get_obs(location="I do not exist")
