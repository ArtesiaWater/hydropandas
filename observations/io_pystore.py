import numpy as np
import pandas as pd
import pystore
from .observation import GroundwaterObs
from .obs_collection import ObsCollection


def set_pystore_path(pystore_path):
    """Set pystore path

    Parameters
    ----------
    pystore_path : str
        path to location with stores
    """
    pystore.set_path(pystore_path)


def item_to_obs(item, ObsClass):
    """convert pystore Item to ObsClass

    Parameters
    ----------
    item : pystore.item.Item
        pystore Item
    ObsClass : type of Obs
        type of observation DataFrame, e.g. GroundwaterObs

    Returns
    -------
    ObsClass
        DataFrame containing observations
    """
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
    """pystore collection to list of observations

    Parameters
    ----------
    store : pystore.store
        pystore store
    collection : pystore.collection
        pystore collection
    ObsClass : type of Obs
        type of observation, by default GroundwaterObs

    Returns
    -------
    list : list of ObsClass
        list of ObsClass DataFrames
    """
    collection = store.collection(collection)
    obs_list = []
    for i in collection.items:
        item = collection.item(i)
        o = item_to_obs(item, ObsClass)
        obs_list.append(o)
    return obs_list


def store_to_obslist(store, ObsClass=GroundwaterObs):
    """convert pystore to list of ObsClass

    Parameters
    ----------
    store : pystore.store
        pystore store containing data
    ObsClass : type of Obs
        type of observation DataFrames, by default GroundwaterObs

    Returns
    -------
    list : list of obs
        list of ObsClass DataFrames

    """
    store = pystore.store(store)
    obs_list = []
    for coll in store.collections:
        obs_list += collection_to_obslist(store, coll,
                                          ObsClass=ObsClass)
    return obs_list


def read_store_metadata(store):
    """read only metadata from pystore

    Parameters
    ----------
    store : pystore.store
        store containing data

    Returns
    -------
    list : list of dictionaries
        list of dictionaries containing metadata

    """
    store = pystore.store(store)
    meta_list = []
    for coll in store.collections:
        c = store.collection(coll)
        for i in c.list_items():
            metadata = pystore.utils.read_metadata(c._item_path(i))
            meta_list.append(metadata)
    return meta_list


def pystore_obslist_to_obscollection(obs_list, name="obs_coll"):
    """convert list of Obs to ObsCollection

    Parameters
    ----------
    obs_list : list of Obs
        list of Obs DataFrames
    name : str, optional
        name of the collection, by default "obs_coll"

    Returns
    -------
    ObsCollection :
        DataFrame containing all Obs
    """
    coldict = [o.to_collection_dict() for o in obs_list]
    obs_df = pd.DataFrame(coldict, columns=coldict[0].keys())
    obs_df.set_index('name', inplace=True)
    meta = {'fname': obs_list[0].meta["datastore"],
            'type': GroundwaterObs,
            'verbose': True}
    oc = ObsCollection(obs_df, name=name, meta=meta)
    return oc
