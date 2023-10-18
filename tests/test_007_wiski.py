from pathlib import Path

import hydropandas as hpd
from hydropandas.io import wiski

wiski_data_path = Path(__file__).parent / "data/2019-WISKI-test"


def test_read_wiski_csv() -> None:
    # download single file

    wiski.read_wiski_file(
        wiski_data_path / "1016_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        verbose=True,
        parse_dates={"datetime": [0, 1]},
        dayfirst=True,
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )


def test_read_wiski_csv2() -> None:
    # download single file

    wiski.read_wiski_file(
        wiski_data_path / "8137_PBF.csv",
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        verbose=True,
        parse_dates={"datetime": [0, 1]},
        dayfirst=True,
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
    )


def test_read_wiski_zip() -> None:
    wiski.read_wiski_dir(
        str(wiski_data_path / "1016_PBF.zip"),
        ObsClass=hpd.GroundwaterObs,
        sep=r"\s+",
        header_sep=":",
        header_identifier=":",
        parse_dates={"datetime": [0, 1]},
        index_col=["datetime"],
        translate_dic={"name": "Station Number", "x": "GlobalX", "y": "GlobalY"},
        verbose=True,
        dayfirst=True,
    )


def test_rijnenijssel_wiski_format() -> None:
    hpd.GroundwaterObs.from_wiski(
        wiski_data_path / "Zwiepse Horstweg Barchem_1024_FT1_WNS9040_MomentaanO.csv",
        header_sep=";",
        end_header_str="#Timestamp",
        parse_dates=[0],
        index_col=[0],
        tz_localize=False,
    )
