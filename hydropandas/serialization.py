import json
from datetime import datetime
import pathlib
import numpy as np


class HydropandasEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, pathlib.PurePath):
            return str(obj)
        elif isinstance(obj, datetime):
            return obj.isoformat()

        # Add other conversions here
        return super().default(obj)
