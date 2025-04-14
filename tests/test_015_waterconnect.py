import hydropandas as hpd


def test_single_well():
    o = hpd.GroundwaterObs.from_waterconnect(119988)

    assert not o.empty
