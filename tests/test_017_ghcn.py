import numpy as np
import pandas as pd
import pytest

import hydropandas as hpd


@pytest.mark.slow
def test_read_ghcn_extent_prcp_rd():
    oc = hpd.read_ghcn(
        extent=[100000, 120000, 450000, 470000],
        crs=28992,
        elements="PRCP",
        tmin="2020",
        tmax="2021",
    )
    assert isinstance(oc, hpd.ObsCollection)

    # test if data from the knmi resturns the same stations and data
    oc_knmi = hpd.read_knmi(
        xy=list(zip(oc.x, oc.y)), meteo_vars=("RD",), starts="2020", ends="2021"
    )
    assert len(oc) == len(oc_knmi)
    for name, name_knmi in zip(oc.index, oc_knmi.index):
        o = oc.loc[name]
        o_knmi = oc_knmi.loc[name_knmi]
        # calculate the distance between the two stations, they should be close to each other (within 2 km)
        d = np.sqrt((o.x - o_knmi.x) ** 2 + (o.y - o_knmi.y) ** 2)
        assert d < 2000, f"distance between GHCN and KNMI station is too large: {d}"

        s_knmi = o_knmi.obs["RD"]
        s_knmi.index = s_knmi.index.normalize()
        # check if data is the same (allowing for small differences)
        df = pd.DataFrame({"ghcn": o.obs["PRCP"], "knmi": s_knmi})
        df = df.dropna()
        if not df.empty:
            corr = df.corr().iloc[0, 1]
            assert corr > 0.99, (
                f"correlation between GHCN and KNMI data is too low: {corr}"
            )
