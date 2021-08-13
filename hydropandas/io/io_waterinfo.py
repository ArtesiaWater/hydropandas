
import os
import zipfile

import numpy as np
import pandas as pd
from pyproj import Proj, transform
from tqdm import tqdm


def read_waterinfo_file(path_to_file, index_cols=None, return_metadata=False,
                        value_col=None, location_col=None,
                        xcol=None, ycol=None, transform_coords=True):
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

    if index_cols is None:
        index_cols = ['WAARNEMINGDATUM', 'WAARNEMINGTIJD']

    if value_col is None:
        value_col = 'NUMERIEKEWAARDE'

    if location_col is None:
        location_col = "MEETPUNT_IDENTIFICATIE"

    if xcol is None:
        xcol = "X"

    if ycol is None:
        ycol = "Y"

    # read data
    df = pd.read_csv(f,
                     sep=';',
                     decimal=',',
                     encoding="ISO-8859-1",
                     parse_dates=[index_cols],
                     dayfirst=True,
                     infer_datetime_format=True,
                     index_col="_".join(index_cols))

    # do some conversions
    df.loc[df[value_col]
           == 999999999, value_col] = np.NaN
    df[value_col] = df[value_col] / 100.

    # parse metadata into dict
    if return_metadata:
        if len(df[location_col].unique()) > 1:
            raise ValueError("File contains data for more than one location!"
                             " Use ObsCollection.from_waterinfo()!")

        metadata = {}
        if transform_coords:
            x, y = transform(Proj('epsg:25831'),
                             Proj('epsg:28992'),
                             df[xcol].iloc[-1],
                             df[ycol].iloc[-1])
        else:
            x = df[xcol].iloc[-1] / 100.
            y = df[ycol].iloc[-1] / 100.
        metadata["name"] = df[location_col].iloc[-1]
        metadata["x"] = x
        metadata["y"] = y

        return df, metadata
    else:
        return df


def read_waterinfo_obs(file_or_dir, ObsClass, progressbar=False, **kwargs):
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

    location_col = kwargs.pop("location_col", "MEETPUNT_IDENTIFICATIE")

    # loop over files
    metadata = {}
    obs_collection = []
    for filenm in (tqdm(files) if progressbar else files):
        # read file or zip
        df = read_waterinfo_file(filenm, location_col=location_col, **kwargs)

        # get location and convert to m RD
        for stn in df[location_col].unique():
            mask = df[location_col] == stn
            x, y = transform(Proj('epsg:25831'),
                             Proj('epsg:28992'),
                             df.loc[mask, 'X'][-1],
                             df.loc[mask, 'Y'][-1])
            metadata = {"name": stn, "x": x, "y": y}

            # add to list
            o = ObsClass(df.loc[mask, :], meta=metadata, **metadata)
            obs_collection.append(o)

    return obs_collection
