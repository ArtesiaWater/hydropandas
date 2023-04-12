from typing import Dict

from hydropandas.observation import (
    EvaporationObs,
    GroundwaterObs,
    WaterQualityObs,
    MeteoObs,
    Obs,
    PrecipitationObs,
    WaterlvlObs,
)


def get_fews_pid(name: str) -> Dict[str, Obs]:
    """Get matching ParameterId's and HydroPandas Observation Classes

    Parameters
    ----------
    name : str
        Waterboard name

    Returns
    -------
    Dict[str, Obs]
        Dictonary with ParameterId and the resulting Observation Class
    """
    return pid[name.lower()]


wsvv_pid = {
    "P.G.d9": PrecipitationObs,  # Gemeten neerslag [mm] - dag (9-9 uur)
    "P.G.h": PrecipitationObs,  # Gemeten neerslag [mm] - uur
    "P.TH.B": PrecipitationObs,  # Doorval berekend [mm]
    "P.f.control": PrecipitationObs,  # Neerslag verwachting [mm] (controle)
    "P.f.determ": PrecipitationObs,  # Neerslag verwachting [mm] (deterministisch)
    "P.f.ensem": PrecipitationObs,  # Neerslag verwachting [mm] (ensemble)
    "P.f.harmonie": PrecipitationObs,  # Neerslag verwachting [mm] (harmonie)
    "P.f.hirlam": PrecipitationObs,  # Neerslag verwachting [mm] (hirlam)
    "P.overschot": PrecipitationObs,  # Neerslagoverschot [mm]
    "P.radar.5": PrecipitationObs,  # Radarneerslag [mm] - 5min
    "P.radar.5.C": PrecipitationObs,  # Radarneerslag [mm] - 5min gecorrigeerd
    "P.radar.5.O": PrecipitationObs,  # Radarneerslag [mm] - 5min ongecorrigeerd
    "P.radar.d": PrecipitationObs,  # Radarneerslag [mm] - dag
    "P.radar.d9": PrecipitationObs,  # Radarneerslag [mm] - dag (9-9 uur)
    "P.radar.h": PrecipitationObs,  # Radarneerslag [mm] - uur
    "Dsnow.f.harmonie": PrecipitationObs,  # Sneeuwhoogte verwachting [m] (harmonie)
    "Pgrau.f.harmonie": PrecipitationObs,  # Neerslag (hagel) verwachting [mm] harmonie
    "Psnow.f.harmonie": PrecipitationObs,  # Neerslag (sneeuw) verwachting [mm] harmonie
    "Eact.d": EvaporationObs,  # Actuele Verdamping [mm] - dag
    "Edef.d": EvaporationObs,  # Verdampingstekort [mm] - dag
    "Eref.d": EvaporationObs,  # Referentiegewasverdamping [mm] - dag
    "EF.GR.f.harmonie": MeteoObs,  # Straling globaal verwacht [J/m2] (harmonie)
    "EF.LH.f.harmonie": MeteoObs,  # Straling latente warmte verwacht [J/m2] (harmonie)
    "EF.LWR.f.harmonie": MeteoObs,  # Straling lange golf verwacht [J/m2] (harmonie)
    "EF.SWR.f.harmonie": MeteoObs,  # Straling korte golf verwacht [J/m2] (harmonie)
    "WGU.f.harmonie": MeteoObs,  # Windstoot verwachting u-component [m/s] (harmonie)
    "WGV.f.harmonie": MeteoObs,  # Windstoot verwachting v-component [m/s] (harmonie)
    "WR.G": MeteoObs,  # Windrichting gemeten [graden]
    "WR.f.control": MeteoObs,  # Windrichting verwacht [graden] (controle)
    "WR.f.determ": MeteoObs,  # Windrichting verwacht [graden] (deterministisch)
    "WR.f.ensem": MeteoObs,  # Windrichting verwacht [graden] (ensemble)
    "WS.G": MeteoObs,  # Windsnelheid gemeten [m/s]
    "WS.f.control": MeteoObs,  # Windsnelheid verwachting [m/s] (controle)
    "WS.f.determ": MeteoObs,  # Windsnelheid verwachting [m/s] (deterministisch)
    "WS.f.ensem": MeteoObs,  # Windsnelheid verwachting [m/s] (ensemble)
    "WSU10.f.harmonie": MeteoObs,  # Windsnelheid verwachting u-component [m/s] harmonie
    "WSU10.f.hirlam": MeteoObs,  # Windsnelheid verwachting u-component [m/s] (hirlam)
    "WSV10.f.harmonie": MeteoObs,  # Windsnelheid verwachting v-component [m/s] harmonie
    "WSV10.f.hirlam": MeteoObs,  # Windsnelheid verwachting v-component [m/s] (hirlam)
    "H.B": WaterlvlObs,  # Waterhoogte berekend [mNAP]
    "H.B.d": WaterlvlObs,  # Waterhoogte berekend [mNAP] - dag
    "H.B.h": WaterlvlObs,  # Waterhoogte berekend [mNAP] - uur
    "H.G": WaterlvlObs,  # Waterhoogte [mNAP]
    "H.G.AA.d": WaterlvlObs,  # Waterhoogte alle [-] - dag
    "H.G.d": WaterlvlObs,  # Waterhoogte [mNAP] - dag
    "H.G.h": WaterlvlObs,  # Waterhoogte [mNAP] - uur
    "H.Max": WaterlvlObs,  # Maximum waterpeil [mNAP]
    "H.S": WaterlvlObs,  # Streefpeil [mNAP]
    "H.S.h": WaterlvlObs,  # Streefpeil [mNAP] - uur
    "H.f": WaterlvlObs,  # Waterhoogte verwachting [mNAP]
    "H.f.determ": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (ecmwf-deterministisch)
    "H.f.ecmwf": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (ecmwf-ensemble)
    "H.f.fr": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (fews rivieren)
    "H.f.matroos": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (matroos)
    "H.f.rws": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (rws)
    "H.f.wdij": WaterlvlObs,  # Waterhoogte verwachting [mNAP] (wdij)
    "H.schuif.G": WaterlvlObs,  # Schuifstand [mNAP]
    "H.wave": WaterlvlObs,  # Golfhoogte [m]
    # "Hk": "Kruinhoogte [mNAP]",
    # "Hk.AA.d": "Kruinhoogte alle [-] - dag",
    # "Hk.d": "Kruinhoogte [mNAP] - dag",
    # "Hk.h": "Kruinhoogte [mNAP] - uur",
    # "Hopv.B": "Opvoerhoogte berekend [mNAP]",
    # "Hpbe.1b": "Peilevaluatie peil eerste bovengrens [mNAP]",
    # "Hpbe.1o": "Peilevaluatie peil eerste ondergrens [mNAP]",
    # "Hpbe.2b": "Peilevaluatie peil tweede bovengrens [mNAP]",
    # "Hpbe.2o": "Peilevaluatie peil tweede ondergrens [mNAP]",
    # "Hu.G": "Uitwateringspeil [mNAP]",
    "Hwk": WaterlvlObs,  # Hoogte waterkolom [m]
    "Hwk.H1": WaterlvlObs,  # Overstortende straal H1 [m]
    "Hwk.H2": WaterlvlObs,  # Hoogteverschil na stuw H2 [m]
    "Hwk.PA.G.Barodiver": WaterlvlObs,  # Hoogteverschil na stuw H2 [m]
    "Hwk.PW.G": WaterlvlObs,  # Gemeten waterdruk [m]
    "Hwk.V.PB": WaterlvlObs,  # Waterstandsafwijking PB [m]
    "Hwk.V.SP": WaterlvlObs,  # Waterstandsafwijking SP [m]
    # "Q.B": WaterlvlObs,  # Debiet FEWS [m3/s]
    # "Q.B.d": WaterlvlObs,  # Debiet FEWS [m3/s] - dag
    # "Q.B.h": WaterlvlObs,  # Debiet FEWS [m3/s] - uur
    # "Q.G": WaterlvlObs,  # Debiet [m3/s]
    # "Q.G.d": WaterlvlObs, # Debiet [m3/s] - dag
    # "Q.G.h": WaterlvlObs,  # Debiet [m3/s] - uur
    # "Q.f": WaterlvlObs,  # Debiet verwachting [m3/s]
    # "Q.f.determ": WaterlvlObs,  # Debiet verwachting [m3/s] (ecmwf-deterministisch)
    # "Q.f.ecmwf": WaterlvlObs,  # Debiet verwachting [m3/s] (ecmwf-ensemble)
    # "Q.f.fr": WaterlvlObs,  # Debiet verwachting [m3/s] (fews rivieren)
    # "Q.f.rws": WaterlvlObs,  # Debiet verwachting [m3/s] (rws)
    # "Q.f.totaal": WaterlvlObs,  # Debiet verwachting [m3/s] (totaal)
    # "Qafval.G": WaterlvlObs,  # Debiet afvalwater [m3/h]
    # "Qafval.G.d": WaterlvlObs,  # Debiet afvalwater [m3/h] - dag
    # "Qafval.G.h": WaterlvlObs,  # Debiet afvalwater [m3/h] - uur
    "T.ow.G": WaterQualityObs,  # Temperatuur oppervlaktewater [C]
    "O2.gehalte": WaterQualityObs,  # Zuurstof gehalte [mg/l]
    "O2.saturatie": WaterQualityObs,  # Zuurstof saturatie [%]
    "PA.G": MeteoObs,  # Luchtdruk gemeten [hPa]
    "PA.G.bar": MeteoObs,  # Luchtdruk gemeten [bar]
    "PA.f.harmonie": MeteoObs,  # Luchtdruk verwachting [hPa] (harmonie)
    "PKEL.V": MeteoObs,  # Drukverschil Keller (P1-P2)
    "RH.G": MeteoObs,  # Luchtvochtigheid [%]
    "RH.f.harmonie": MeteoObs,  # Relatieve luchtvochtigheid [%] (harmonie)
    "T.air.G": MeteoObs,  # Temperatuur lucht [C]
    "T.air.f.control": MeteoObs,  # Temperatuur verwacht (controle)
    "T.air.f.determ": MeteoObs,  # Temperatuur verwacht (deterministisch)
    "T.air.f.ensem": MeteoObs,  # Temperatuur verwacht (ensemble)
    "T.air.f.harmonie": MeteoObs,  # Temperatuur verwacht (harmonie)
    # "PW.G": "Waterspanning [hPa]",
    # "PZ.G": "Zuigspanning [kPa]",
    # "SM.G": "Bodemvochtigheid [m3/m3]",
    # "T.soil.G": "Temperatuur bodem [C]",
    "SS.G": WaterlvlObs,  # Stroomsnelheid gemeten [m/s
    "SS.G.d": WaterlvlObs,  # Stroomsnelheid gemeten [m/s] - dag
    "SS.G.h": WaterlvlObs,  # Stroomsnelheid gemeten [m/s] - uur
    "GW.G": GroundwaterObs,  # Grondwaterstand [mNAP]
    "GW.G.hm": WaterlvlObs,  # Grondwaterstand - handmeting [mNAP]
    "SH": GroundwaterObs,  # Stijghoogte [mNAP]
    "T.gw.G": WaterQualityObs,  # Temperatuur grondwater [C]
}


pid = {"wsvv": wsvv_pid}
