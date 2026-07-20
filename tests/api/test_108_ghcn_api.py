import pandas as pd


def test_ghcn():
    station_id = "SF003715790" #Utrecht
    url = f"https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/{station_id}.dly"

    colspecs = [
        (0, 11),  # ID
        (11, 15),  # YEAR
        (15, 17),  # MONTH
        (17, 21),  # ELEMENT
    ] + [(21 + i * 8, 26 + i * 8) for i in range(31)]  # VALUE1-31

    colnames = ["id", "year", "month", "element"] + [f"value{i}" for i in range(1, 32)]

    data = pd.read_fwf(url, colspecs=colspecs, names=colnames)
    assert not data.empty
