# %%

import requests


def test_catalog():
    json = {
        "CatalogusFilter": {
            "Eenheden": True,
            "Grootheden": True,
            "Hoedanigheden": True,
            "Groeperingen": True,
            "Parameters": True,
            "Compartimenten": True,
        }
    }
    url = "https://waterwebservices.rijkswaterstaat.nl/METADATASERVICES_DBO/OphalenCatalogus"
    r = requests.post(url, json=json, timeout=300)

    r.raise_for_status()


def test_waarnemingen():
    json = {
        "AquoPlusWaarnemingMetadata": {
            "AquoMetadata": {
                "Compartiment": {"Code": "OW"},
                "Eenheid": {"Code": "cm"},
                "Grootheid": {"Code": "WATHTE"},
                "Hoedanigheid": {"Code": "NAP"},
                "Parameter": {"Code": "NVT"},
                "Groepering": {"Code": "NVT"},
            }
        },
        "Locatie": {"X": 627246.546177034, "Y": 5756327.60323024, "Code": "SCHOONHVN"},
        "Periode": {
            "Begindatumtijd": "2024-01-01T00:00:00.000+00:00",
            "Einddatumtijd": "2024-01-03T00:00:00.000+00:00",
        },
    }
    url = "https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenWaarnemingen"
    r = requests.post(url, json=json, timeout=300)

    r.raise_for_status()
