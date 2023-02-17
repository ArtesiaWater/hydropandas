import hydropandas as hpd
from hydropandas.io import io_bro
import logging

logging.basicConfig(level=logging.DEBUG)


def test_metadata():
    # single observation
    bro_id = "GMW000000036287"
    meta = io_bro.get_metadata_from_gmw(bro_id, 1)

    return meta


def test_metadata_full():
    # single observation
    bro_id = "GMW000000036287"
    meta = io_bro.get_full_metadata_from_gmw(bro_id, 1)

    return meta


def test_groundwater_monitoring_net_metadata():
    bro_id = "GMN000000000163"
    obs_list, meta = io_bro.get_obs_list_from_gmn(
        bro_id, hpd.GroundwaterObs, only_metadata=True
    )

    return obs_list


def test_groundwater_observations():
    bro_id = "GLD000000012893"
    measurements, meta = io_bro.get_bro_groundwater(
        bro_id, tube_nr=None, only_metadata=False
    )
    return measurements


def test_gld_no_monitoringnet():
    bro_id = "GLD000000013128"
    measurements, meta = io_bro.get_bro_groundwater(
        bro_id, tube_nr=None, only_metadata=False
    )
    return measurements


def test_groundwater_observations2():
    bro_id = "GLD000000008061"
    measurements, meta = io_bro.get_bro_groundwater(
        bro_id, tube_nr=None, only_metadata=False
    )

    ax = measurements["values"].plot(color="blue", marker=".")
    measurements.loc[measurements["qualifier"] == "goedgekeurd", "values"].plot(
        color="green", marker="."
    )

    return measurements


def test_get_gld_id_from_gmw():

    bro_id = "GMW000000036287"
    bro_id = "GMW000000055372"
    bro_id = "GMW000000059186"
    gld = io_bro.get_gld_id_from_gmw(bro_id, tube_nr=1)

    return gld


def test_obs_list_from_extent():

    extent = (102395, 103121, 434331, 434750)
    extent = [116500, 120000, 439000, 442000]
    obs_list = io_bro.get_obs_list_from_extent(
        extent, hpd.GroundwaterObs, tmin=None, tmax=None, epsg=28992
    )

    return obs_list
