import pandas as pd
import pytest

import hydropandas as hpd
from hydropandas.io import wow


def test_wow_strftime() -> None:
    timestamp = pd.Timestamp(year=2010, month=10, day=2, hour=10, minute=34, second=0)
    timestampstr = wow._wow_strftime(timestamp)
    assert timestampstr == "2010-10-02T10%3A30%3A00"


@pytest.mark.skip
def test_get_wow_stn() -> None:
    stn = "423216079"  # Macquarie Island
    start = pd.Timestamp(year=2023, month=6, day=4, hour=1, minute=34, second=0)
    end = pd.Timestamp(year=2023, month=6, day=5, hour=23, minute=54, second=0)
    obs = hpd.MeteoObs.from_wow(stn=stn, meteo_var="temperature", start=start, end=end)
    assert obs.name == "Macquarie Island"


@pytest.mark.skip
def test_get_wow_nearest() -> None:
    lat = -35.4184
    lon = 149.0937
    xy = [lon, lat]
    start = pd.Timestamp(year=2023, month=6, day=4, hour=1, minute=34, second=0)
    end = pd.Timestamp(year=2023, month=6, day=5, hour=23, minute=54, second=0)
    obs = hpd.PrecipitationObs.from_wow(xy=xy, start=start, end=end)
    assert obs.name == "Hillanhome"
