import os
import pytest

import hydropandas as hpd


# Original tests when only Vitens was supported
def test_single_observation():
    code = "27BP0003"
    o = hpd.GroundwaterObs.from_lizard(code)
    assert o.tube_nr == 1


def test_extent():
    extent = [201500, 202000, 502000, 502200]
    oc = hpd.read_lizard(extent)
    assert not oc.empty


@pytest.mark.slow
def test_codes():
    oc = hpd.read_lizard(
        codes=["39F-0735", "39F-0736", "39F-0737"], combine_method="merge"
    )
    assert not oc.empty


@pytest.mark.slow
def test_many_tubed_well():
    oc = hpd.read_lizard(codes="EEWP004", tube_nr="all")
    assert not oc.empty


@pytest.mark.slow
def test_complex_well():
    oc = hpd.read_lizard(codes="BUWP014", tube_nr="all")
    assert not oc.empty


def test_combine():
    o = hpd.GroundwaterObs.from_lizard("39F-0736", tube_nr=1, combine_method="combine")
    assert o.tube_nr == 1


# Additional tests for use with the 'Rotterdam' data
api_key_rotterdam = os.environ["LIZARD_ROTTERDAM_API_KEY"]
assert len(api_key_rotterdam) == 41
if api_key_rotterdam is not None:
    auth = ("__key__", api_key_rotterdam)
else:
    auth = None


def test_single_observation_rotterdam():
    code = "GMW000000036819"
    o = hpd.GroundwaterObs.from_lizard(code, organisation="rotterdam", auth=auth)
    assert o.tube_nr == 1


def test_extent_rotterdam():
    extent = [68_500, 69_500, 443_500, 444_500]
    oc = hpd.read_lizard(
        extent=extent,
        which_timeseries=["hand", "diver", "diver_validated"],
        datafilters=None,
        combine_method="merge",
        organisation="rotterdam",
        auth=auth,
    )
    assert not oc.empty


@pytest.mark.slow
def test_codes_rotterdam():
    oc = hpd.read_lizard(
        codes=["GMW000000036819", "GMW000000037933"],
        which_timeseries=["hand", "diver", "diver_validated"],
        datafilters=None,
        combine_method="merge",
        organisation="rotterdam",
        auth=auth,
    )
    assert not oc.empty
