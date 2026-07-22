# %%
# packages
import numpy as np
import pandas as pd
import requests

from hydropandas.io import knmi

from logging import getLogger

logger = getLogger(__name__)

class EmptyDataFrameError(Exception):
    pass

#%%
# read tmin tmax from knmi website
meteo_tminmax_nederland = (
    pd.read_html("https://www.knmi.nl/nederland-nu/klimatologie/daggegevens")[0]
    .dropna(how="any")
    .convert_dtypes()
    .set_index("Nummer")
)
meteo_tminmax_noordzee = pd.read_html(
    "https://www.knmi.nl/nederland-nu/klimatologie/daggegevens_Noordzee"
)[0].dropna(how="any").convert_dtypes().set_index("Nummer")
meteo_tminmax_knmi = pd.concat([meteo_tminmax_nederland, meteo_tminmax_noordzee]).sort_index()
meteo_tminmax_knmi["tmin"] = [
    pd.to_datetime(x, format="%Y%m%d").strftime("%Y-%m-%d")
    for x in meteo_tminmax_knmi["Vanaf"]
]
meteo_tminmax_knmi["tmax"] = [
    pd.to_datetime(x, format="%Y%m%d") if x != "gisteren" else pd.Timestamp.today() for x in meteo_tminmax_knmi["Tot en met"]
]
meteo_tminmax_knmi["tmax"] = [
    "9999-12-31"
    if (pd.Timestamp.today() - x) < pd.Timedelta(days=365)
    else x.strftime("%Y-%m-%d")
    for x in meteo_tminmax_knmi["tmax"]
]

# %%
# read hydropandas knmi meteo stations bookkeeping
meteo_df = pd.read_json("knmi_meteostation.json")
assert meteo_tminmax_knmi.index.difference(meteo_df.index).empty, "KNMI meteo station index does not match KNMI tmin/tmax index"
#%%
# fill meteo_df with tmin, tmax of meteo_tminmax_knmi
meteo_df.loc[meteo_tminmax_knmi.index, "tmin"] = meteo_tminmax_knmi["tmin"]
meteo_df.loc[meteo_tminmax_knmi.index, "tmax"] = meteo_tminmax_knmi["tmax"]

#%%
# some extra metadata file from knmi data platform, but not used in hydropandas
# because our locations are more accurate
# https://dataplatform.knmi.nl/dataset/waarneemstations-csv-1-0
# meteo_meta = pd.read_csv("csv_knmi_nl_20260615.csv", index_col=0)

#%%
# loop through all KNMI meteo stations and check if data is available via url and api
meteo_tminmax_url = pd.DataFrame(index=meteo_tminmax_knmi.index, columns=["tmin", "tmax"])
meteo_tminmax_api = pd.DataFrame(index=meteo_tminmax_knmi.index, columns=["tmin", "tmax"])

didx = [pd.Timestamp.today().normalize() - pd.Timedelta(days=1)] + list(
    reversed(
        pd.date_range(pd.Timestamp("1900-01-01"), pd.Timestamp.today(), freq="10YS")
    )
)
for stn in meteo_df.index:
    logger.info(f"Checking KNMI meteo station {stn}")
    df_daily, _ = knmi.get_knmi_daily_meteo_url(stn)
    vars_available_daily = ~df_daily.isna().all().drop(index=["STN"])
    meteo_df.loc[stn, vars_available_daily.index] = vars_available_daily.values
    meteo_tminmax_url.loc[stn, "tmin"] = df_daily.index[0]
    meteo_tminmax_url.loc[stn, "tmax"] = df_daily.index[-1]

    for end, start in zip(didx[:-1], didx[1:]):
        try:
            df_hourly, _ = knmi.get_hourly_meteo_api(stn, start=start, end=end)
            if df_hourly.empty:
                raise EmptyDataFrameError
            vars_available_hourly = ~df_hourly.isna().all().drop(index=["STN"])
            vah = vars_available_hourly[vars_available_hourly]
            meteo_df.loc[stn, vah.index] = vah.values
        except (
            pd.errors.EmptyDataError,
            requests.ConnectionError,
            requests.HTTPError,
            EmptyDataFrameError,
        ) as e:
            logger.error(f"Geen data {stn}, {start}, {end}, {e}")

    try:
        df_daily, _ = knmi.get_daily_meteo_api(stn)
        meteo_df.loc[stn, "api_available"] = True
        meteo_tminmax_api.loc[stn, "tmin"] = df_daily.index[0]
        meteo_tminmax_api.loc[stn, "tmax"] = df_daily.index[-1]

    except requests.HTTPError as e:
        meteo_df.loc[stn, "api_available"] = False
        logger.error(f"Geen data {stn}, {e}. Setting api_available to False")

#%%
# save meteo data variables to json
meteo_dft = meteo_df.fillna(False).drop(columns=[""])
meteo_dft.to_json("knmi_meteostation.json")

#%%
# compare tminmax api, url and knmi website
meteo_tminmax_compare = pd.concat(
    [
        meteo_tminmax_knmi[["tmin", "tmax"]].rename(columns={"tmin": "tmin_knmi", "tmax": "tmax_knmi"}),
        meteo_tminmax_url.rename(columns={"tmin": "tmin_url", "tmax": "tmax_url"}),
        meteo_tminmax_api.rename(columns={"tmin": "tmin_api", "tmax": "tmax_api"}),
    ],
    axis=1,
)

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
    "9999-12-31"
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
