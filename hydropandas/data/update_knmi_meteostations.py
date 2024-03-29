# %%
import pandas as pd
from requests import ConnectionError

from hydropandas.io import knmi


meteo_df = pd.read_json("knmi_meteostation.json")
meteo_dft = meteo_df.copy()

didx = [pd.Timestamp.today().normalize() - pd.Timedelta(days=1)] + list(
    reversed(
        pd.date_range(
            pd.Timestamp("1900-01-01"), pd.Timestamp.today(), freq="10YS"
        ).to_list()
    )
)

for stn in meteo_dft.index:
    df_daily, _ = knmi.get_knmi_daily_meteo_url(stn)
    vars_available_daily = ~df_daily.isna().all().drop(index=["STN"])
    meteo_dft.loc[stn, vars_available_daily.index] = vars_available_daily.values

    for end, start in zip(didx[:-1], didx[1:]):
        try:
            df_hourly, _ = knmi.get_knmi_hourly_meteo_api(stn, start=start, end=end)
            vars_available_hourly = ~df_hourly.isna().all().drop(index=["STN"])
            vah = vars_available_hourly[vars_available_hourly]
            meteo_dft.loc[stn, vah.index] = vah.values
        except ConnectionError:
            print(f"Geen data {stn}, {start}, {end}")
        except pd.errors.EmptyDataError:
            print(f"Geen data {stn}, {start}, {end}")

meteo_dft = meteo_dft.fillna(False).drop(columns=[""])
meteo_dft.to_json("knmi_meteostation.json")
