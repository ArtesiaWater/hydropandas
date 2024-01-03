import hydropandas as hpd


def test_single_observation():
    code = "27BP0003"
    o = hpd.GroundwaterObs.from_lizard(code)
    assert o.tube_nr == 1


def test_extent():
    extent = [201500, 202000, 502000, 502200]
    oc = hpd.read_lizard(extent)
    assert not oc.empty


def test_codes():
    oc = hpd.read_lizard(codes="27BP0003")
    assert not oc.empty
