import logging

from hydropandas.io import dino

logging.basicConfig(level=logging.DEBUG)


def test_dino_csv_old_style():
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/B37A0112001_1.csv"
    dino.read_dino_groundwater_csv(fname)


def test_dino_csv_new_style():
    fname = "./tests/data/2024-Dino-test/DINO_Grondwaterstanden/B02H0089001.csv"
    dino.read_dino_groundwater_csv(fname)


def test_dino_csv_duplicate_index():
    # contains 1 duplicate index 2019-11-19
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/B22D0155001_1.csv"
    measurements, _ = dino.read_dino_groundwater_csv(fname)

    # check if measurements contains duplicate indices
    assert measurements.index.duplicated().any()

    measurements, _ = dino.read_dino_groundwater_csv(
        fname, remove_duplicates=True, keep_dup="last"
    )

    # check if measurements contains no duplicate indices
    assert ~measurements.index.duplicated().any()
