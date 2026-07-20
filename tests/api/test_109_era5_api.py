import requests


def test_era5():
    ERA5_ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"

    params = {
        "latitude": 5.0,
        "longitude": 52.0,
        "start_date": "2020-01-01",
        "end_date": "2021-01-02",
        "timezone": "UTC",
        "models": "era5_seamless",
        "daily": "precipitation_sum",
    }

    r = requests.get(ERA5_ARCHIVE_URL, params=params, timeout=600)

    r.raise_for_status()
    assert "daily" in r.json()
