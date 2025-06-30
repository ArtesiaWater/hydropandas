# %%

import json

import requests


def test_get_measurements():
    url = "https://www.waterconnect.sa.gov.au/_layouts/15/dfw.sharepoint.wdd/WDDDMS.ashx/GetWaterLevelDownload?bulkOutput=CSV"
    json_data = {"DHNOs": [119988], "Pumping": True, "Anomalous": True}
    r = requests.post(
        url, verify=True, data={"exportdata": json.dumps(json_data)}, timeout=60
    )
    r.raise_for_status()
