import os
import zipfile

import numpy as np
import pandas as pd


def read_soilinst_file(
    path,
    transform_coords=True,
):
    """Read Soilinst file (XLE)

    Parameters
    ----------
    path : str
        path to Soilint file (.xle)
    transform_coords : boolean
        convert coordinates from WGS84 to RD

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing file content
    metadata : dict, optional
        dict containing metadata
    """
    from pyproj import Transformer

    # open file
    name = os.path.splitext(os.path.basename(path))[0]
    if path.endswith(".xle"):
        f = path
    elif path.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        f = zf.open("{}.xle".format(name))
    else:
        raise NotImplementedError(
            "File type '{}' not supported!".format(os.path.splitext(path)[-1])
        )

    # read channel 1 data header
    df_ch1_data_header = pd.read_xml(
        path,
        xpath="/Body_xle/Ch1_data_header"
        )
    series_ch1_data_header = df_ch1_data_header.T.iloc[:, 0]
    colname_ch1 = series_ch1_data_header.Identification.lower() + \
        '_' + series_ch1_data_header.Unit.lower()

    # read channel 2 data header
    df_ch2_data_header = pd.read_xml(
        path,
        xpath="/Body_xle/Ch2_data_header"
        )
    series_ch2_data_header = df_ch2_data_header.T.iloc[:, 0]
    colname_ch2 = series_ch2_data_header.Identification.lower() + \
        '_' + series_ch2_data_header.Unit.lower()

    # read observations
    df = pd.read_xml(
        path,
        xpath="/Body_xle/Data/Log",
        )
    df.rename(
        columns={'ch1': colname_ch1,
                 'ch2': colname_ch2},
        inplace=True
        )
    if 'ms' in df.columns:
        df['date_time'] = pd.to_datetime(
            df['Date'] + ' ' + df['Time']) + \
            pd.to_timedelta(df['ms'], unit='ms')
        drop_cols = ['id', 'Date', 'Time', 'ms']
    else:
        df['date_time'] = pd.to_datetime(
            df['Date'] + ' ' + df['Time'])
        drop_cols = ['id', 'Date', 'Time']
    df.set_index('date_time', inplace=True)

    df.drop(columns=drop_cols, inplace=True)

    # parse metadata into dict, per group in XLE file
    metadata = {}
    # read file info
    df_file_info = pd.read_xml(
        path,
        xpath="/Body_xle/File_info"
        )
    dict_file_info = df_file_info.T.iloc[:, 0].to_dict()

    # read instrument info
    df_instrument_info = pd.read_xml(
        path,
        xpath="/Body_xle/Instrument_info"
        )
    dict_instrument_info = df_instrument_info.T.iloc[:, 0].to_dict()

    # read instrument info
    df_instrument_info_data_header = pd.read_xml(
        path,
        xpath="/Body_xle/Instrument_info_data_header"
        )
    dict_instrument_info_data_header = df_instrument_info_data_header.T.iloc[
        :, 0].to_dict()

    metadata = {**dict_file_info,
                **dict_instrument_info,
                **dict_instrument_info_data_header}

    if transform_coords:
        transformer = Transformer.from_crs("epsg:4326", "epsg:28992")
        x, y = transformer.transform(metadata["Latitude"], metadata["Longtitude"])
        x = np.round(x, 2)
        y = np.round(y, 2)
    else:
        x = metadata["Latitude"]
        y = metadata["Longtitude"]
    metadata["x"] = x
    metadata["y"] = y
    metadata["filename"] = f
    metadata["source"] = metadata["Created_by"]
    metadata["name"] = name
    metadata["monitoring_well"] = name
    metadata["unit"] = series_ch1_data_header.Unit.lower()
    metadata["metadata_available"] = True

    return df, metadata

