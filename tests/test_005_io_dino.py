# -*- coding: utf-8 -*-
"""Created on Mon Jun 24 11:43:27 2019.

@author: oebbe
"""
from hydropandas.io import io_dino

import logging

logging.basicConfig(level=logging.DEBUG)


def test_dino_csv():
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/" "B33F0080001_1.csv"
    measurements, meta = io_dino.read_dino_groundwater_csv(fname)

    return measurements, meta


def test_dino_csv_duplicate_index():
    # contains 1 duplicate index 2019-11-19
    fname = "./tests/data/2019-Dino-test/Grondwaterstanden_Put/" "B22D0155001_1.csv"
    measurements, meta = io_dino.read_dino_groundwater_csv(fname)

    # check if measurements contains duplicate indices
    assert measurements.index.duplicated().any()

    measurements, meta = io_dino.read_dino_groundwater_csv(
        fname, remove_duplicates=True, keep_dup="last"
    )

    # check if measurements contains no duplicate indices
    assert measurements.index.duplicated().any() == False

    return measurements, meta


def test_dino_metadata():
    # download metadata
    dinorest = io_dino.DinoREST()
    meta = dinorest.get_gwo_metadata(location="B52C0057", tube_nr="001")

    return meta


def test_dino_metadata2():
    # download metadata without sample metadata in json
    dinorest = io_dino.DinoREST()
    meta = dinorest.get_gwo_metadata(location="B57B0069", tube_nr="002")
    assert meta["metadata_available"]
    return meta


def test_dino_metadata3():
    # try to download metadata of a well that does not have metadata
    dinorest = io_dino.DinoREST()
    meta = dinorest.get_gwo_metadata(location="B45G1147", tube_nr="001")

    assert not meta["metadata_available"]

    return meta


def test_get_dino_locations():
    bbox = [130988.58668351, 386465.80657847, 135922.01739107, 391922.4284911]

    gdf = io_dino.get_dino_locations(bbox=bbox, layer="grondwatermonitoring")

    return gdf
