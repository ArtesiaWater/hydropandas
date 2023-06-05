from functools import lru_cache
from typing import List, Optional

import numpy as np
import pandas as pd
import requests
from pyproj import Transformer

from .. import util

URL_WOW_KNMI = "https://wow.knmi.nl/sites"


meteo_vars_wow = (
    "rain_rate",
    "temperature",
    # "weather_code",
    "wind",
    "humidity",
    "pressure_msl",
)

meteo_vars_wow_translate = {
    "neerslagintensiteit (mm/uur)": "precipitation rate [mm/h]",
    "temperatuur (C)": "temperature [C]",
    "windsnelheid (m/s)": "windspeed [m/s]",
    "relatieve vochtigheid  (%)": "relative humidity [%]",
    "luchtdruk (hPa)": "air pressure [hPa]",
}


def get_stations_amateur(
    meteo_var: str = "rain_rate",
    date: Optional[pd.Timestamp] = None,
    bbox: Optional[List[float]] = None,
) -> pd.DataFrame:
    if meteo_var not in meteo_vars_wow:
        raise ValueError(f"{meteo_var} must be one of {meteo_vars_wow}")

    if date is None:
        date = pd.Timestamp.today() - pd.Timedelta(days=0, hours=0, minutes=10)
    datestr = _wow_knmi_strftime(date)

    if bbox is None:
        bbox = [2.5, 50.0, 8.0, 54.0]  # Netherlands extent
    bboxstr = (
        str(bbox[0]) + "%20" + str(bbox[1]) + "," + str(bbox[2]) + "%20" + str(bbox[3])
    )  # lat lon,lat lon

    # filter=wow_observations
    # filter=official_observations
    try:
        r = requests.get(
            f"{URL_WOW_KNMI}?bbox={bboxstr}&layer={meteo_var}&filter=&date={datestr}",
            timeout=120,
        )
        r.raise_for_status()
        sites = r.json()["sites"]
        transformer = Transformer.from_crs(4326, 28992)
        lat_lon = np.array([x["geo"]["coordinates"] for x in sites])
        xy = pd.DataFrame(
            np.column_stack(transformer.transform(lat_lon[:, 1], lat_lon[:, 0])),
            columns=["x", "y"],
        )
        stations = pd.concat([pd.DataFrame(sites), xy], axis=1).set_index(
            "id", drop=False
        )
    except requests.HTTPError as ex:
        raise ex

    return stations


def _wow_knmi_strftime(timestamp: pd.Timestamp) -> str:
    return (
        timestamp.strftime("%Y-%m-%dT%H")
        + "%3A"
        + timestamp.strftime("%M")[0]
        + "0%3A00"
    )


def get_knmi_wow_metadata(stn_id: str, report_id: str) -> dict:
    try:
        r = requests.get(f"{URL_WOW_KNMI}/{stn_id}?report_id={report_id}")
        r.raise_for_status()
        metadata = r.json()
    except requests.HTTPError as ex:
        raise ex
    return metadata


@lru_cache()
def get_knmi_wow_measurements(
    stn_id: str,
    meteo_var: str = "rain_rate",
    start: pd.Timestamp = None,
    end: pd.Timestamp = None,
):
    start, end = util._start_end_to_datetime(start, end)
    if (end - start) > pd.to_timedelta(365, unit="D"):
        t = list(pd.date_range(start, end, freq="365D"))
        if t[-1] != end:
            t.append(end)
        dfs = []
        for i in range(len(t) - 1):
            dfs.append(
                get_knmi_wow_measurements(stn_id, meteo_var, start=t[i], end=t[i + 1])
            )
        df = pd.concat(dfs)
        df = df[~df.index.duplicated()]
        return df

    startstr = _wow_knmi_strftime(start)
    endstr = _wow_knmi_strftime(end)
    # get station measurements
    try:
        meas = pd.read_csv(
            f"{URL_WOW_KNMI}/{stn_id}/export?start={startstr}&end={endstr}&layer={meteo_var}",
            delimiter=";",
            index_col=["datum"],
        )
        meas.index = pd.DatetimeIndex(pd.to_datetime(meas.index.values), name="date")
        measurements = meas.rename(
            columns={
                x: meteo_vars_wow_translate[x.replace(f" [{stn_id}]", "")]
                for x in meas.columns
            }
        )
    except requests.HTTPError as ex:
        raise ex

    if meteo_var == "rain_rate":
        weights = (
            (measurements.index[1:] - measurements.index[:-1])
            / pd.Timedelta(1, unit="H")
        ).values
        measurements = (
            measurements.iloc[0:-1]
            .multiply(weights, axis=0)
            .rename(columns={"precipitation rate [mm/h]": "precipitation [mm]"})
        )

    return measurements
