# -*- coding: utf-8 -*-
"""Created on Mon Jun 24 11:43:27 2019.

@author: oebbe
"""
import logging

from hydropandas.io import dino

logging.basicConfig(level=logging.DEBUG)


def test_dino_csv():
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/" "B33F0080001_1.csv"
    dino.read_dino_groundwater_csv(fname)

    return


def test_dino_csv_duplicate_index():
    # contains 1 duplicate index 2019-11-19
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/" "B22D0155001_1.csv"
    measurements, meta = dino.read_dino_groundwater_csv(fname)

    # check if measurements contains duplicate indices
    assert measurements.index.duplicated().any()

    measurements, meta = dino.read_dino_groundwater_csv(
        fname, remove_duplicates=True, keep_dup="last"
    )

    # check if measurements contains no duplicate indices
    assert ~measurements.index.duplicated().any()

    return
