import json
import pathlib
import pyproj
from datetime import date, datetime

import numpy as np
from pandas import Timestamp


class HydropandasEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, (pathlib.Path, pathlib.PurePath)):
            return str(o)
        elif isinstance(o, (datetime, date, Timestamp)):
            return o.isoformat()
        elif isinstance(o, pyproj.CRS):
            return o.to_string()
        elif isinstance(o, type):
            return f"class : {o.__name__}"

        # Add other conversions here
        return super().default(o)
