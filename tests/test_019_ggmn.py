import pytest

import hydropandas as hpd
from hydropandas.io import ggmn


@pytest.mark.slow
def test_read_ggmn_realworld_extent():
    # Small extent around central Netherlands with known GGMN locations.
    oc = hpd.read_ggmn(
        extent=[5.03, 5.08, 52.13, 52.18],
        epsg=4326,
        tmin="2024-01-01",
        tmax="2026-12-31",
        max_locations=1,
        max_pages=5,
        keep_all_obs=True,
        timeout=120,
    )

    assert isinstance(oc, hpd.ObsCollection)
    assert len(oc) >= 1

    o = oc.iloc[0].obs
    assert o.x >= 5.03 and o.x <= 5.08
    assert o.y >= 52.13 and o.y <= 52.18


@pytest.mark.slow
def test_get_level_measurements_known_record_realworld():
    df, unit = ggmn.get_level_measurements(
        695507,
        parameter="Water level elevation a.m.s.l.",
        max_pages=5,
        timeout=120,
    )

    assert not df.empty
    assert df.index.is_monotonic_increasing
    assert unit == "m"
