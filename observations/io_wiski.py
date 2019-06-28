import os
import tempfile

import numpy as np
import pandas as pd

from .util import unzip_file


def _read_wiski_header(f):
    line = f.readline()
    header = dict()
    while ":" in line:
        prop, val = line.split(":")
        prop = prop.strip()
        val = val.strip()
        try:
            val = float(val)
        except ValueError:
            pass
        header[prop] = val
        line = f.readline()

    return line, header


def read_wiski_file(fname, read_series=True, verbose=False):
    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))

    # read header
    with open(fname, "r") as f:
        line, header = _read_wiski_header(f)

        # get column names of data
        # columns = split(r'\s{2,}', line)
        columns = line.split("\t")
        columns = [icol.strip("\n") for icol in columns]
        columns = [icol.replace("[", "_[") for icol in columns]

        validator = np.lib._iotools.NameValidator()
        columns = validator(columns)

        if read_series:
            # read data
            data = pd.read_csv(f, sep=r"\s+", parse_dates={"datetime": [0, 1]},
                               header=None, names=columns, index_col=["datetime"])

            # convert Value to float
            data["Value_mNAP"] = pd.to_numeric(
                data["Value_mNAP"], errors="coerce")
        else:
            data = None

    return header, data


def read_wiski_dir(dirname, ObsClass=None, suffix=".csv", verbose=True,
                   unpackdir=None, force_unpack=False, preserve_datetime=False,
                   keep_all_obs=True, **kwargs):

    # unzip dir
    if dirname.endswith('.zip'):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(zipf, dirname, force=force_unpack,
                   preserve_datetime=preserve_datetime)

    # read filenames
    files = os.listdir(os.path.join(dirname))
    if suffix:
        files = [f for f in files if f.endswith(suffix)]

    if not files:
        raise FileNotFoundError('no files were found in {} that end with {}'.format(
            os.path.join(dirname), suffix))

    # gather all obs in list
    obs_list = []
    for csv in files:
        obs = ObsClass.from_wiski(os.path.join(
            dirname, csv), verbose=verbose, **kwargs)
        if obs.metadata_available:
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            if verbose:
                print('not added to collection -> {}'.format(fname))

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df


if __name__ == "__main__":
    fname = r"g:\My Drive\r\01projekt\18038040_WS_BRABANTSE_DELTA_TRA_NNP\03data\BD\1016_PBF.csv"
    header, data = read_wiski_file(fname)

    zipname = r"g:\My Drive\r\01projekt\18038040_WS_BRABANTSE_DELTA_TRA_NNP\03data\BD\1016_PBF.zip"
    import art_tools.observations.obs_collection as obs
    o = obs.ObsCollection.from_wiski(zipname, verbose=True)
