import pandas as pd

import hydropandas as hpd
from hydropandas.io import wow


def test_wow_strftime() -> None:
    timestamp = pd.Timestamp(year=2010, month=10, day=2, hour=10, minute=34, second=0)
    timestampstr = wow._wow_strftime(timestamp)
    assert timestampstr == "2010-10-02T10%3A30%3A00"


def test_get_wow_stn() -> None:
    stn = "423216079"  # Macquarie Island
    startdate = pd.Timestamp(year=2023, month=6, day=4, hour=1, minute=34, second=0)
    enddate = pd.Timestamp(year=2023, month=6, day=5, hour=23, minute=54, second=0)
    obs = hpd.MeteoObs.from_wow(
        stn=stn, meteo_var="temperature", startdate=startdate, enddate=enddate
    )
    assert obs.name == "Macquarie Island"


def test_get_wow_nearest() -> None:
    lat = -35.4184
    lon = 149.0937
    xy = [lon, lat]
    startdate = pd.Timestamp(year=2023, month=6, day=4, hour=1, minute=34, second=0)
    enddate = pd.Timestamp(year=2023, month=6, day=5, hour=23, minute=54, second=0)
    obs = hpd.PrecipitationObs.from_wow_nearest_xy(
        xy=xy, startdate=startdate, enddate=enddate
    )
    assert obs.name == "Hillanhome"
