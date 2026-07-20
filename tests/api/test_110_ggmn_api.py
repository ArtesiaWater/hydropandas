
import requests

def test_wfs():
    GGMN_WFS_URL = "https://ggis.un-igrac.org/geoserver/wfs"
    GGMN_LEVEL_LAYER = "groundwater:GGMN_Levels_Data"

    params = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeNames": GGMN_LEVEL_LAYER,
        "srsName": "EPSG:4326",
        "outputFormat": "application/json",
        "count": 200,
        "bbox": f"{5.03},{52.13},{5.08},{52.18},EPSG:4326",
    }

    r = requests.get(GGMN_WFS_URL, params=params, timeout=120)
    r.raise_for_status()
    result = r.json()

    assert len(result) > 1
    assert 'properties' in result['features'][0]

def test_measurements():

    record_id = '695507'
    url = f"https://ggis.un-igrac.org/groundwater/record/{record_id}/WellLevelMeasurement/list?set=1"
    

    r = requests.get(
            url,
            headers={"X-Requested-With": "XMLHttpRequest"},
            timeout=120,
        )
    r.raise_for_status()
    result = r.json()
    assert 'data' in result