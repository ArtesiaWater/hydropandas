import pandas as pd

import hydropandas as hpd


def test_single_waterlevel_observation():
    tmin = pd.Timestamp("2020-1-1")
    tmax = pd.Timestamp("2020-1-3")
    o = hpd.WaterlvlObs.from_matroos(
        location="krimpen a/d lek",
        source="observed",
        unit="waterlevel",
        tmin=tmin,
        tmax=tmax,
    )

    assert not o.empty


def test_waterlevel_observation_within_extent():
    extent = [116_500, 118_000, 418_000, 422_000]
    oc = hpd.read_matroos(
        extent=extent, units="waterlevel", sources="observed", keep_all_obs=False
    )
    oc

    assert not oc.empty
    assert all([not o.empty for o in oc.obs.values])
