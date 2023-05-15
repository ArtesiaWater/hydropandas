# import os
import numpy as np
import pandas as pd
import pytest
from hydropandas import obs_collection as oc
from hydropandas import observation as obs


def test_obscollection_fews_lowmemory():
    fews_gw_prod = oc.ObsCollection.from_fews_xml(
        "./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip",
        locations=None,
        low_memory=True,
    )
    return fews_gw_prod


# %% PYSTORE


def test_obscollection_to_pystore():
    obsc = test_obscollection_fews_lowmemory()
    obsc.to_pystore(
        "test_pystore",
        "./tests/data/2019-Pystore-test",
        groupby="locatie",
        overwrite=True,
    )


def test_obscollection_from_pystore():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore", "./tests/data/2019-Pystore-test"
    )
    return obsc


def test_obscollection_pystore_only_metadata():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore", "./tests/data/2019-Pystore-test", read_series=False
    )
    return obsc


def test_obscollection_pystore_extent():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore",
        "./tests/data/2019-Pystore-test",
        extent=[115534, 115539, 0, 10000000],
    )
    return obsc


def test_obscollection_pystore_item_names():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore", "./tests/data/2019-Pystore-test", item_names=["MPN-N-2"]
    )
    return obsc


def test_obs_from_pystore_item():
    import pystore

    pystore.set_path("./tests/data/2019-Pystore-test")
    store = pystore.store("test_pystore")
    coll = store.collection(store.collections[0])
    item = coll.item(list(coll.list_items())[0])
    o = obs.GroundwaterObs.from_pystore_item(item)
    return o


if __name__ == "__main__":
    obsc = test_obscollection_fews_lowmemory()
    obsc.to_pystore(
        "test_pystore",
        "./tests/data/2019-Pystore-test",
        groupby="locatie",
        overwrite=True,
    )

    import pystore

    pystore.set_path("./tests/data/2019-Pystore-test")
    store = pystore.store("test_pystore")
    coll = store.collection(store.collections[0])
    item = coll.item(list(coll.list_items())[0])
    o = obs.GroundwaterObs.from_pystore_item(item)
