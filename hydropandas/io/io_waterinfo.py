
import os
import zipfile

import numpy as np
import pandas as pd
from pyproj import Proj, transform
from tqdm import tqdm


def read_waterinfo_file(path_to_file, return_metadata=False):
    """Read waterinfo file (CSV or zip)

    Parameters
    ----------
    path_to_file : str
        path to file

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing file content
    metadata : dict, optional
        dict containing metadata, returned if return_metadata is
        True, default is False
    """

    name = os.path.splitext(os.path.basename(path_to_file))[0]

    if path_to_file.endswith(".csv"):
        f = path_to_file
    elif path_to_file.endswith(".zip"):
        zf = zipfile.ZipFile(path_to_file)
        f = zf.open('{}.csv'.format(name))
    else:
        raise NotImplementedError("File type '{}' not supported!".format(
            os.path.splitext(path_to_file)[-1]))

    # read data
    df = pd.read_csv(f,
                     sep=';',
                     decimal=',',
                     encoding="ISO-8859-1",
                     parse_dates=[['WAARNEMINGDATUM', 'WAARNEMINGTIJD']],
                     dayfirst=True,
                     infer_datetime_format=True,
                     index_col='WAARNEMINGDATUM_WAARNEMINGTIJD')

    # do some conversions
    df.loc[df['NUMERIEKEWAARDE'] ==
           999999999, 'NUMERIEKEWAARDE'] = np.NaN
    df['NUMERIEKEWAARDE'] = df['NUMERIEKEWAARDE'] / 100.

    # parse metadata into dict
    if return_metadata:
        if len(df["MEETPUNT_IDENTIFICATIE"].unique()) > 1:
            raise ValueError("File contains data for more than one location!"
                             " Use ObsCollection.from_waterinfo()!")

        metadata = {}
        x, y = transform(Proj('epsg:25831'),
                         Proj('epsg:28992'),
                         df['X'].iloc[-1],
                         df['Y'].iloc[-1])
        metadata["name"] = df["MEETPUNT_IDENTIFICATIE"].iloc[-1]
        metadata["x"] = x
        metadata["y"] = y

        return df, metadata
    else:
        return df


def read_waterinfo_obs(file_or_dir, ObsClass, progressbar=False):
    """Read waterinfo file or directory and extract locations and observations.

    Parameters
    ----------
    file_or_dir : str
        path to file or directory
    ObsClass: Obs type
        type of Obs to store data in
    progressbar : bool, optional
        show progressbar if True, default is False

    Returns
    -------
    obs_collection : list
        list of Obs objects
    """

    # Waterinfo file
    if os.path.isfile(file_or_dir):
        files = [file_or_dir]
    # directory with waterinfo files (zips or csvs)
    elif os.path.isdir(file_or_dir):
        files = [os.path.join(file_or_dir, f)
                 for f in sorted(os.listdir(file_or_dir))]
    else:
        raise NotImplementedError("Provide path to file or directory!")

    # loop over files
    metadata = {}
    obs_collection = []
    for filenm in tqdm(files):
        # read file or zip
        df = read_waterinfo_file(filenm)

        # get location and convert to m RD
        for stn in df["MEETPUNT_IDENTIFICATIE"].unique():
            mask = df['MEETPUNT_IDENTIFICATIE'] == stn
            x, y = transform(Proj(init='epsg:25831'),
                             Proj(init='epsg:28992'),
                             df.loc[mask, 'X'][-1],
                             df.loc[mask, 'Y'][-1])
            metadata = {"name": stn, "x": x, "y": y}

            # add to list
            o = ObsClass(df.loc[mask, :], meta=metadata, **metadata)
            obs_collection.append(o)

    return obs_collection
