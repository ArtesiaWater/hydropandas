import pytest

import hydropandas as hpd


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
        codes=["39F-0735", "39F-0736", "39F-0737"], type_timeseries="merge"
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
    hpd.GroundwaterObs.from_lizard("39F-0736", tube_nr=1, type_timeseries="combine")
