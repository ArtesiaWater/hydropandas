import logging

import hydropandas as hpd
from hydropandas.io import bro

logging.basicConfig(level=logging.DEBUG)


def test_metadata():
    # single observation
    bro_id = "GMW000000036287"
    bro.get_metadata_from_gmw(bro_id, 1)


def test_metadata_full():
    # single observation
    bro_id = "GMW000000036287"
    bro.get_full_metadata_from_gmw(bro_id, 1)


def test_groundwater_observations():
    bro_id = "GLD000000012893"
    bro.get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False)


def test_gld_no_monitoringnet():
    bro_id = "GLD000000013128"
    bro.get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False)


def test_groundwater_observations2():
    bro_id = "GLD000000008061"
    measurements, _ = bro.get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False)

    measurements["values"].plot(color="blue", marker=".")
    measurements.loc[measurements["qualifier"] == "goedgekeurd", "values"].plot(
        color="green", marker="."
    )


def test_get_gld_ids_from_gmw():
    bro_id = "GMW000000036287"
    bro_id = "GMW000000055372"
    bro_id = "GMW000000059186"
    bro_id = "GMW000000063853"  # two gld ids
    bro_id = "GMW000000030953"  # three gld ids
    bro.get_gld_ids_from_gmw(bro_id, tube_nr=1)


def test_obs_list_from_extent():
    # extent = (102395, 103121, 434331, 434750)
    extent = [117800, 118300, 439700, 439800]  # 4 measurements in extent 2025-4-7
    bro.get_obs_list_from_extent(
        extent, hpd.GroundwaterObs, tmin=None, tmax=None, epsg=28992, only_metadata=True
    )
