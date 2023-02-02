# -*- coding: utf-8 -*-
"""Created on Mon Jun 24 11:43:27 2019.

@author: oebbe
"""

import hydropandas as hpd
from hydropandas.io import io_wiski


def test_read_wiski_csv():
    # download single file

    header, data = io_wiski.read_wiski_file(
        "./tests/data/2019-WISKI-test/1016_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        verbose=True,
        parse_dates={"datetime": [0, 1]},
        infer_datetime_format=True,
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )

    return header, data


def test_read_wiski_csv2():
    # download single file

    header, data = io_wiski.read_wiski_file(
        "./tests/data/2019-WISKI-test/8137_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        verbose=True,
        parse_dates={"datetime": [0, 1]},
        infer_datetime_format=False,
        dayfirst=True,
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )

    return header, data


# %%


def test_read_wiski_zip():
    obs_df = io_wiski.read_wiski_dir(
        "./tests/data/2019-WISKI-test/1016_PBF.zip",
        ObsClass=hpd.GroundwaterObs,
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
        verbose=True,
    )

    return obs_df


def test_rijnenijssel_wiski_format():
    o = hpd.GroundwaterObs.from_wiski(
        (
            "./tests/data/2019-WISKI-test/"
            "Zwiepse Horstweg Barchem_1024_FT1_WNS9040_MomentaanO.csv"
        ),
        header_sep=";",
        end_header_str="#Timestamp",
        parse_dates=[0],
        index_col=[0],
        infer_datetime_format=True,
        tz_localize=False,
    )
    return o
