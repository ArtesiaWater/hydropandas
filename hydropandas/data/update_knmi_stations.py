# %%
# packages
import numpy as np
import pandas as pd
import requests

from hydropandas.io import knmi


class EmptyDataFrameError(Exception):
    pass


# %%
# meteo stations
meteo_df = pd.read_json("knmi_meteostation.json")
meteo_dft = meteo_df.copy()

didx = [pd.Timestamp.today().normalize() - pd.Timedelta(days=1)] + list(
    reversed(
        pd.date_range(pd.Timestamp("1900-01-01"), pd.Timestamp.today(), freq="10YS")
    )
)


for stn in meteo_dft.index:
    df_daily, _ = knmi.get_knmi_daily_meteo_url(stn)
    vars_available_daily = ~df_daily.isna().all().drop(index=["STN"])
    meteo_dft.loc[stn, vars_available_daily.index] = vars_available_daily.values

    for end, start in zip(didx[:-1], didx[1:]):
        try:
            df_hourly, _ = knmi.get_knmi_hourly_meteo_api(stn, start=start, end=end)
            if df_hourly.empty:
                raise EmptyDataFrameError
            vars_available_hourly = ~df_hourly.isna().all().drop(index=["STN"])
            vah = vars_available_hourly[vars_available_hourly]
            meteo_dft.loc[stn, vah.index] = vah.values
        except (
            pd.errors.EmptyDataError,
            requests.ConnectionError,
            EmptyDataFrameError,
        ):
            print(f"Geen data {stn}, {start}, {end}")

meteo_dft = meteo_dft.fillna(False).drop(columns=[""])

# add tmin, tmax
meteo_tminmax = (
    pd.read_html("https://www.knmi.nl/nederland-nu/klimatologie/daggegevens")[0]
    .dropna(how="any")
    .convert_dtypes()
    .set_index("Nummer")
)
meteo_tminmax["tmin"] = [
    pd.to_datetime(x, format="%Y%m%d").strftime("%Y-%m-%d")
    for x in meteo_tminmax["Vanaf"]
]
meteo_tminmax["tmax"] = [
    pd.to_datetime(x, format="%Y%m%d") for x in meteo_tminmax["Tot en met"]
]
meteo_tminmax["tmax"] = [
    None
    if (pd.Timestamp.today() - x) < pd.Timedelta(days=365)
    else x.strftime("%Y-%m-%d")
    for x in meteo_tminmax["tmax"]
]
meteo_dft.loc[meteo_tminmax.index, "tmin"] = meteo_tminmax["tmin"].values
meteo_dft.loc[meteo_tminmax.index, "tmax"] = meteo_tminmax["tmax"].values

meteo_dft.to_json("knmi_meteostation.json")
# %%
# neerslag stations
prec_df = pd.read_json("knmi_neerslagstation.json")
prec_dft = prec_df.copy()
prec_tminmax = (
    pd.concat(
        pd.read_html(
            "https://www.knmi.nl/nederland-nu/klimatologie/monv/reeksen", index_col=0
        )
    )
    .reset_index()
    .set_index("Nr")
)
unknown_location = prec_tminmax.index[~np.isin(prec_tminmax.index, prec_dft.index)]
prec_tminmax = prec_tminmax.drop(unknown_location)
tminmax = [[x[0], x[-1]] for x in prec_tminmax["Periode"].str.split(" ")]
tmin = [pd.to_datetime(x[0], format="%Y%m%d").strftime("%Y-%m-%d") for x in tminmax]
tmax = [pd.to_datetime(x[1], format="%Y%m%d") for x in tminmax]
tmax = [
    None
    if (pd.Timestamp.today() - x) < pd.Timedelta(days=365)
    else x.strftime("%Y-%m-%d")
    for x in tmax
]
prec_dft.loc[prec_tminmax.index, "tmin"] = tmin
prec_dft.loc[prec_tminmax.index, "tmax"] = tmax

# Sometimes there are no measurements for a period after tmin, this is not a problem
# for all station but De Bilt. This is why we manually correct De Bilt.
prec_dft.loc[550, "tmin"] = "1898-01-01"

prec_dft.sort_index().loc[
    :, ["lon", "lat", "name", "x", "y", "altitude", "tmin", "tmax", "RD"]
].to_json("knmi_neerslagstation.json")
