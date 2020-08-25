"""A pystore is a datastore for Pandas Dataframes designed to store timeseries.

The functions in this module aim to save an obs_collection to a pystore. The
main advantages of a pystore are:

- smaller file size compared to .csv files
- exchangable format, the pystore format is independent of the pc
  (unlike pickle)

A pystore with an ObsCollection has 3 layers:

1. directory with name of the pystore with all the information of an
   ObsCollection.
2. Inside the pystore directory are directories with the collections of the
   pystore. One pystore collection corresponds to a subcollection of an
   ObsCollection (which is an ObsCollection on its own).

   You can use pystore collections to group observations from an ObsCollection
   by a certain feature, for example the location of the observation. One
   location can have multiple observations.

3. Inside the pystore collection are directories with the items of the
   pystore. An item of a pystore contains all the data from an Observation
   object.
"""

import os

import pandas as pd
from tqdm import tqdm

try:
    import pystore
except ModuleNotFoundError:
    pass

from ..obs_collection import ObsCollection
from ..observation import GroundwaterObs


def set_pystore_path(pystore_path):
    """Set pystore path.

    Parameters
    ----------
    pystore_path : str
        path to location with stores
    """
    pystore.set_path(pystore_path)


def item_to_obs(item, ObsClass, nameby="item"):
    """convert pystore Item to ObsClass.

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
    if len(item.data.index) == 0:
        df = pd.DataFrame(columns=item.data.columns)
    else:
        df = item.to_pandas()
    item.metadata["datastore"] = str(item.datastore)
    if nameby == "item":
        name = item.item
    elif nameby == "collection":
        name = item.collection
    elif nameby == "both":
        name = item.collection + "__" + item.item
    else:
        raise ValueError(
            "'{}' is not a valid option for 'nameby'".format(nameby))

    metadata = item.metadata
    if not "name" in metadata.keys():
        metadata["name"] = name
    else:
        metadata.update({"name": name})

    obs_attr_dic = {}
    for attr in ObsClass._metadata:
        if attr in metadata.keys():
            obs_attr_dic[attr] = metadata[attr]

    o = ObsClass(df, **obs_attr_dic, meta=metadata)
    return o


def collection_to_obslist(store, collection, ObsClass=GroundwaterObs,
                          item_names=None, nameby="item", verbose=True):
    """pystore collection to list of observations.

    Parameters
    ----------
    store : pystore.store
        pystore store
    collection : pystore.collection
        pystore collection
    ObsClass : type of Obs
        type of observation, by default GroundwaterObs
    item_names : list of str
        item (Observation) names that will be extracted from the store
        the other items (Observations) will be ignored.

    Returns
    -------
    list : list of ObsClass
        list of ObsClass DataFrames
    """
    obs_list = []
    if collection in store.collections:
        collection = store.collection(collection)
    else:
        if verbose:
            print("Not in store -> {0}".format(collection))
        return obs_list

    if item_names is None:
        items = collection.items
    else:
        items = set(item_names) & set(collection.items)

    for i in items:
        try:
            item = collection.item(i)
            o = item_to_obs(item, ObsClass, nameby=nameby)
            obs_list.append(o)
        except Exception:
            if verbose:
                print("Skipped -> {0}: {1}".format(collection.collection, i))
    return obs_list


def store_to_obslist(store, ObsClass=GroundwaterObs, collection_names=None,
                     item_names=None, nameby="item", verbose=True,
                     progressbar=False):
    """convert pystore to list of ObsClass.

    Parameters
    ----------
    store : pystore.store
        pystore store containing data
    ObsClass : type of Obs
        type of observation DataFrames, by default GroundwaterObs
    item_names : list of str
        item (Observation) names that will be extracted from the store
        the other items (Observations) will be ignored. if None all items
        are read.
    nameby : str
        pick whether obs are named by collection or item name

    Returns
    -------
    list : list of obs
        list of ObsClass DataFrames
    """
    store = pystore.store(store)
    obs_list = []
    if collection_names is None:
        collections = store.collections
    else:
        collections = collection_names
    for coll in (tqdm(collections) if progressbar else collections):
        obs_list += collection_to_obslist(store, coll,
                                          ObsClass=ObsClass,
                                          item_names=item_names,
                                          nameby=nameby, verbose=verbose)
    return obs_list


def read_store_metadata(store, items='all', verbose=False):
    """read only metadata from pystore.

    Parameters
    ----------
    store : pystore.store
        store containing data
    items : str, list of str, optional
        if 'all' read all items
        if 'first' read first item only
        if list of str, read those items
    verbose : bool, optional
        if True, print progress info

    Returns
    -------
    list : list of dictionaries
        list of dictionaries containing metadata
    """
    store = pystore.store(store)
    meta_list = []
    for coll in store.collections:
        c = store.collection(coll)
        if items == "all":
            item_list = c.list_items()
        elif items == 'first':
            item_list = list(c.list_items())[0:1]
        else:
            item_list = items
        for i in item_list:
            metadata = pystore.utils.read_metadata(c._item_path(i))
            if metadata is None:
                if verbose:
                    print("Cannot read metadata for {0}/{1}".format(coll, i))
                metadata = dict()
                metadata['item_name'] = ""
            else:
                metadata['item_name'] = i
            metadata['collection_name'] = coll
            meta_list.append(metadata)
    return meta_list


def pystore_obslist_to_obscollection(obs_list, name="obs_coll"):
    """convert list of Obs to ObsCollection.

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


def read_pystore(storename, pystore_path,
                 ObsClass, extent=None, collection_names=None,
                 item_names=None, nameby="item",
                 read_series=True, verbose=True, progressbar=False):
    # set path
    set_pystore_path(pystore_path)

    if not os.path.isdir(os.path.join(pystore_path, storename)):
        raise FileNotFoundError("pystore -> '{}' "
                                "does not exist".format(storename))

    # obtain item names within extent
    if extent is not None:
        meta_list = read_store_metadata(storename, items="first")
        obs_df = pd.DataFrame(meta_list)
        obs_df.set_index('item_name', inplace=True)
        obs_df['x'] = pd.to_numeric(obs_df.x, errors='coerce')
        obs_df['y'] = pd.to_numeric(obs_df.y, errors='coerce')
        item_names = obs_df[(obs_df.x > extent[0]) & (obs_df.x < extent[1]) & (
            obs_df.y > extent[2]) & (obs_df.y < extent[3])].index

    if read_series:
        obs_list = store_to_obslist(
            storename,
            ObsClass=ObsClass,
            collection_names=collection_names,
            item_names=item_names,
            nameby=nameby,
            verbose=verbose,
            progressbar=progressbar)

        return obs_list
    else:
        if item_names is None:
            item_names = "all"
        meta_list = read_store_metadata(storename,
                                        items=item_names,
                                        verbose=verbose)
        meta = {}
        obs_df = pd.DataFrame(meta_list)
        if nameby == "collection":
            obs_df.set_index('collection_name', inplace=True)
        elif nameby == "item":
            obs_df.set_index('item_name', inplace=True)
        elif nameby == "both":
            obs_df["name"] = obs_df["collection_name"] + "__" + \
                obs_df["item_name"]
            obs_df.set_index('name', inplace=True)
        else:
            raise ValueError("'{}' is not a valid option "
                             "for 'nameby'".format(nameby))
        return obs_df
