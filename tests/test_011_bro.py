import logging

from hydropandas.io import bro

logging.basicConfig(level=logging.DEBUG)


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