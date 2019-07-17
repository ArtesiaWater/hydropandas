import numpy as np
import pandas as pd
import pystore
from .observation import GroundwaterObs
from .obs_collection import ObsCollection


def set_pystore_path(pystore_path):
    pystore.set_path(pystore_path)


def item_to_obs(item, ObsClass):
    df = item.to_pandas()
    try:
        x = item.metadata["x"]
        y = item.metadata["y"]
    except KeyError:
        x = np.nan
        y = np.nan
    item.metadata["datastore"] = item.datastore
    o = ObsClass(df, x=x, y=y, meta=item.metadata)
    return o


def collection_to_obslist(store, collection, ObsClass=GroundwaterObs):
    collection = store.collection(collection)
    obs_list = []
    for i in collection.items:
        item = collection.item(i)
        o = item_to_obs(item, ObsClass)
        obs_list.append(o)
    return obs_list


def store_to_obslist(store, ObsClass=GroundwaterObs):
    store = pystore.store(store)
    obs_list = []
    for coll in store.collections:
        obs_list += collection_to_obslist(store, coll,
                                          ObsClass=ObsClass)
    return obs_list


def read_store_metadata(store):
    store = pystore.store(store)
    meta_list = []
    for coll in store.collections:
        c = store.collection(coll)
        for i in c.list_items():
            metadata = pystore.utils.read_metadata(c._item_path(i))
            meta_list.append(metadata)
    return meta_list


def pystore_obslist_to_obscollection(obs_list, name="obs_coll"):
    coldict = [o.to_collection_dict() for o in obs_list]
    obs_df = pd.DataFrame(coldict, columns=coldict[0].keys())
    obs_df.set_index('name', inplace=True)
    meta = {'fname': obs_list[0].meta["datastore"],
            'type': GroundwaterObs,
            'verbose': True}
    oc = ObsCollection(obs_df, name=name, meta=meta)
    return oc
