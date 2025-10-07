# %%

import requests


def test_waarnemingen():
    URL = "https://noos.matroos.rws.nl/direct/get_series.php?"

    params = {
        "loc": "schoonhoven",
        "source": "observed",
        "unit": "waterlevel",
        "tstart": "202401010000",
        "tstop": "202401020000",
    }

    r = requests.get(URL, params=params, timeout=600)

    r.raise_for_status()
