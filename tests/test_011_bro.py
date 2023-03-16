import logging

import hydropandas as hpd
from hydropandas.io import bro

logging.basicConfig(level=logging.DEBUG)


def test_metadata():
    # single observation
    bro_id = "GMW000000036287"
    bro.get_metadata_from_gmw(bro_id, 1)

    return


def test_metadata_full():
    # single observation
    bro_id = "GMW000000036287"
    bro.get_full_metadata_from_gmw(bro_id, 1)

    return


def test_groundwater_monitoring_net_metadata():
    bro_id = "GMN000000000163"
    bro.get_obs_list_from_gmn(bro_id, hpd.GroundwaterObs, only_metadata=True)

    return


def test_groundwater_observations():
    bro_id = "GLD000000012893"
    bro.get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False)
    return


def test_gld_no_monitoringnet():
    bro_id = "GLD000000013128"
    bro.get_bro_groundwater(bro_id, tube_nr=None, only_metadata=False)
    return


def test_groundwater_observations2():
    bro_id = "GLD000000008061"
    measurements, meta = bro.get_bro_groundwater(
        bro_id, tube_nr=None, only_metadata=False
    )

    measurements["values"].plot(color="blue", marker=".")
    measurements.loc[measurements["qualifier"] == "goedgekeurd", "values"].plot(
        color="green", marker="."
    )

    return


def test_get_gld_id_from_gmw():
    bro_id = "GMW000000036287"
    bro_id = "GMW000000055372"
    bro_id = "GMW000000059186"
    bro.get_gld_id_from_gmw(bro_id, tube_nr=1)

    return


def test_get_gld_id_from_gmw_quality_regime():
    # bro id with two gld's with a different quality regime
    # both gld's have no measurements 13-3-2023
    bro_id = "GMW000000063853"
    gld1 = bro.get_gld_id_from_gmw(bro_id, tube_nr=1, quality_regime="IMBRO/A")

    gld2 = bro.get_gld_id_from_gmw(bro_id, tube_nr=1, quality_regime="IMBRO")
    assert (
        gld1 != gld2
    ), "different quality regimes should return different gld id's for this gmw id"

    return


def test_obs_list_from_extent():
    extent = (102395, 103121, 434331, 434750)
    extent = [116500, 120000, 439000, 442000]
    bro.get_obs_list_from_extent(
        extent, hpd.GroundwaterObs, tmin=None, tmax=None, epsg=28992, only_metadata=True
    )

    return
