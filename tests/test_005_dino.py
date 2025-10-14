import logging
from pathlib import Path

from hydropandas.io import dino

logging.basicConfig(level=logging.DEBUG)

datadir = Path(__file__).parent / "data"
dino2019 = datadir / "2019-Dino-test"
dino2024 = datadir / "2024-Dino-test"


def test_dino_csv_old_style():
    fname = dino2019 / "Grondwaterstanden_Put/B37A0112001_1.csv"
    dino.read_dino_groundwater_csv(fname)


def test_dino_csv_new_style():
    fname = dino2024 / "DINO_Grondwaterstanden" / "B02H0089001.csv"
    dino.read_dino_groundwater_csv(fname)


def test_dino_csv_duplicate_index():
    # contains 1 duplicate index 2019-11-19
    fname = dino2019 / "Grondwaterstanden_Put" / "B22D0155001_1.csv"
    measurements, _ = dino.read_dino_groundwater_csv(fname)

    # check if measurements contains duplicate indices
    assert measurements.index.duplicated().any()

    measurements, _ = dino.read_dino_groundwater_csv(
        fname, remove_duplicates=True, keep_dup="last"
    )

    # check if measurements contains no duplicate indices
    assert ~measurements.index.duplicated().any()
