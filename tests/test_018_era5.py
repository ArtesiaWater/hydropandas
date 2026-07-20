import pandas as pd
import pytest

import hydropandas as hpd


def test_read_era5_daily_precipitation_realworld():
    tmin = "2020-01-01"
    tmax = "2020-01-02"

    oc = hpd.read_era5(
        extent=[5.0, 5.0, 52.0, 52.0],
        epsg=4326,
        variables=("precipitation_sum",),
        tmin=tmin,
        tmax=tmax,
        grid_size=0.25,
        timeout=120,
    )

    assert isinstance(oc, hpd.ObsCollection)
    assert len(oc) == 1

    o = oc.iloc[0].obs
    assert not o.empty
    assert o.unit == "m"
    assert o.index.min() == pd.Timestamp(tmin) + pd.Timedelta(days=1)
    assert o.index.max() == pd.Timestamp(tmax) + pd.Timedelta(days=1)


test_read_era5_daily_precipitation_realworld()


def test_read_era5_daily_precipitation_realworld_xy():
    tmin = "2020-01-01"
    tmax = "2020-01-02"

    oc = hpd.read_era5(
        xy=(5.0, 52.0),
        epsg=4326,
        variables=("precipitation_sum",),
        tmin=tmin,
        tmax=tmax,
        timeout=120,
    )

    assert isinstance(oc, hpd.ObsCollection)
    assert len(oc) == 1

    o = oc.iloc[0].obs
    assert not o.empty
    assert o.unit == "m"
    assert o.x == 5.0
    assert o.y == 52.0


def test_read_era5_land_daily_realworld_xy():
    oc = hpd.read_era5(
        xy=(5.0, 52.0),
        epsg=4326,
        variables=("precipitation_sum",),
        source="era5_land",
        tmin="2020-01-01",
        tmax="2020-01-02",
        keep_all_obs=True,
        timeout=120,
    )

    assert isinstance(oc, hpd.ObsCollection)
    assert len(oc) == 1

    o = oc.iloc[0].obs
    assert o.meta.get("era5_source") == "era5_land"
    assert o.meta.get("interval") == "daily"


def test_read_era5_hourly_source_forces_hourly():
    tmin = "2020-01-01"
    tmax = "2020-01-01"

    oc = hpd.read_era5(
        xy=(5.0, 52.0),
        epsg=4326,
        variables=("precipitation",),
        source="era5_hourly",
        interval="daily",
        tmin=tmin,
        tmax=tmax,
        timeout=120,
    )

    assert isinstance(oc, hpd.ObsCollection)
    assert len(oc) == 1

    o = oc.iloc[0].obs
    assert not o.empty
    assert o.meta.get("era5_source") == "era5_hourly"
    assert o.meta.get("interval") == "hourly"
    assert o.index.min() == pd.Timestamp(tmin) + pd.Timedelta(hours=1)
