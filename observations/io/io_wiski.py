import os
import tempfile

import numpy as np
import pandas as pd

from ..util import unzip_file, get_files


def _read_wiski_header(f, header_sep=":", header_identifier='#',
                       end_header_str=None):
    line = f.readline()
    header = dict()
    while header_identifier in line:
        prop, val = line.split(header_sep)
        prop = prop.strip().replace(header_identifier, "")
        val = val.strip()
        try:
            val = float(val)
        except ValueError:
            pass
        header[prop] = val
        line = f.readline()
        if end_header_str is not None:
            if end_header_str in line:
                break

    return line, header


def read_wiski_file(fname, sep=";", header_sep=None, header_identifier='#',
                    read_series=True, infer_datetime_format=True,
                    translate_dic={}, verbose=False,
                    tz_localize=True, to_mnap=True, **kwargs):
    if verbose:
        print('reading -> {}'.format(os.path.split(fname)[-1]))

    # manually break header parse at certain point
    if "end_header_str" in kwargs.keys():
        end_header_str = kwargs.pop("end_header_str")
    else:
        end_header_str = None

    # read header
    with open(fname, "r") as f:
        if header_sep is None:
            line, header = _read_wiski_header(f, end_header_str=end_header_str,
                                              header_identifier=header_identifier)
        else:
            line, header = _read_wiski_header(f, header_sep=header_sep,
                                              header_identifier=header_identifier,
                                              end_header_str=end_header_str)

        # specify names
        for key, item in translate_dic.items():
            header[key] = header[item]

        # get column names of data
        # columns = split(r'\s{2,}', line)
        if sep == r"\s+":
            columns = line.split("\t")
        else:
            columns = line.split(sep)
        columns = [icol.strip("\n") for icol in columns]
        columns = [icol.replace("[", "_[") for icol in columns]

        validator = np.lib._iotools.NameValidator()
        columns = validator(columns)

        if read_series:
            # read data
            data = pd.read_csv(f, sep=sep, header=None, names=columns,
                               infer_datetime_format=infer_datetime_format,
                               **kwargs)

            if tz_localize:
                data.index = data.index.tz_localize(None)

            # convert Value to float
            col = [icol for icol in data.columns if
                   icol.lower().startswith("value")][0]

            data[col] = pd.to_numeric(
                data[col], errors="coerce")

            if to_mnap:
                data.rename(columns={col: 'Stand_m_tov_NAP'}, inplace=True)
        else:
            data = None

    return header, data


def read_wiski_dir(dirname, ObsClass=None, suffix=".csv",
                   unpackdir=None, force_unpack=False, preserve_datetime=False,
                   keep_all_obs=True, verbose=True, **kwargs):

    # get files
    dirname, unzip_fnames = get_files(dirname, ext=suffix,
                                      force_unpack=force_unpack,
                                      preserve_datetime=preserve_datetime)

    if not unzip_fnames:
        raise FileNotFoundError("no files were found in '{}' that end with '{}'".format(
            os.path.join(dirname), suffix))

    # gather all obs in list
    obs_list = []
    for i, csv in enumerate(unzip_fnames):
        if verbose:
            print("reading {0}/{1} -> {2}".format(i+1, len(unzip_fnames), csv))
        obs = ObsClass.from_wiski(os.path.join(
            dirname, csv), verbose=False, **kwargs)

        if obs.metadata_available:
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            if verbose:
                print('not added to collection -> {}'.format(csv))

    # create dataframe
    obs_df = pd.DataFrame([o.to_collection_dict() for o in obs_list],
                          columns=obs_list[0].to_collection_dict().keys())
    obs_df.set_index('name', inplace=True)

    return obs_df