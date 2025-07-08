import pandas as pd

import hydropandas as hpd


def test_single_well():
    d = {
        "UNIT_NO": 702402090,
        "LAT": -36.7203406,
        "LON": 140.5163628,
        "REF_ELEV": 37.185,
        "GRND_ELEV": 37.13,
        "MAPNUM": "7024-2090",
        "STATUS": "EQP",
        "PH": -9999.0,
        "TDS": 1867,
        "TDSDATE": "2009/12/10",
        "YIELD": 5.05,
        "YIELD_DATE": "1960/04/29",
        "WATER_CUT": -9999.0,
        "DTW": 3.87,
        "SWL": 3.82,
        "SWLDATE": "2025/03/24",
        "GRND_ELEV_": 37.60035706,
        "HGUID": 0,
    }
    meta_series = pd.Series(d)
    o = hpd.GroundwaterObs.from_waterconnect(119988, meta_series=meta_series)

    assert not o.empty
