# flake8: noqa
wsvv_pid = {
    "P.G.d9": "PrecipitationObs",  # Gemeten neerslag [mm] - dag (9-9 uur)
    "P.G.h": "PrecipitationObs",  # Gemeten neerslag [mm] - uur
    "P.TH.B": "PrecipitationObs",  # Doorval berekend [mm]
    "P.f.control": "PrecipitationObs",  # Neerslag verwachting [mm] (controle)
    "P.f.determ": "PrecipitationObs",  # Neerslag verwachting [mm] (deterministisch)
    "P.f.ensem": "PrecipitationObs",  # Neerslag verwachting [mm] (ensemble)
    "P.f.harmonie": "PrecipitationObs",  # Neerslag verwachting [mm] (harmonie)
    "P.f.hirlam": "PrecipitationObs",  # Neerslag verwachting [mm] (hirlam)
    "P.overschot": "PrecipitationObs",  # Neerslagoverschot [mm]
    "P.radar.5": "PrecipitationObs",  # Radarneerslag [mm] - 5min
    "P.radar.5.C": "PrecipitationObs",  # Radarneerslag [mm] - 5min gecorrigeerd
    "P.radar.5.O": "PrecipitationObs",  # Radarneerslag [mm] - 5min ongecorrigeerd
    "P.radar.d": "PrecipitationObs",  # Radarneerslag [mm] - dag
    "P.radar.d9": "PrecipitationObs",  # Radarneerslag [mm] - dag (9-9 uur)
    "P.radar.h": "PrecipitationObs",  # Radarneerslag [mm] - uur
    "Dsnow.f.harmonie": "PrecipitationObs",  # Sneeuwhoogte verwachting [m] (harmonie)
    "Pgrau.f.harmonie": "PrecipitationObs",  # Neerslag (hagel) verwachting [mm] harmonie
    "Psnow.f.harmonie": "PrecipitationObs",  # Neerslag (sneeuw) verwachting [mm] harmonie
    "Eact.d": "EvaporationObs",  # Actuele Verdamping [mm] - dag
    "Edef.d": "EvaporationObs",  # Verdampingstekort [mm] - dag
    "Eref.d": "EvaporationObs",  # Referentiegewasverdamping [mm] - dag
    "EF.GR.f.harmonie": "MeteoObs",  # Straling globaal verwacht [J/m2] (harmonie)
    "EF.LH.f.harmonie": "MeteoObs",  # Straling latente warmte verwacht [J/m2] (harmonie)
    "EF.LWR.f.harmonie": "MeteoObs",  # Straling lange golf verwacht [J/m2] (harmonie)
    "EF.SWR.f.harmonie": "MeteoObs",  # Straling korte golf verwacht [J/m2] (harmonie)
    "WGU.f.harmonie": "MeteoObs",  # Windstoot verwachting u-component [m/s] (harmonie)
    "WGV.f.harmonie": "MeteoObs",  # Windstoot verwachting v-component [m/s] (harmonie)
    "WR.G": "MeteoObs",  # Windrichting gemeten [graden]
    "WR.f.control": "MeteoObs",  # Windrichting verwacht [graden] (controle)
    "WR.f.determ": "MeteoObs",  # Windrichting verwacht [graden] (deterministisch)
    "WR.f.ensem": "MeteoObs",  # Windrichting verwacht [graden] (ensemble)
    "WS.G": "MeteoObs",  # Windsnelheid gemeten [m/s]
    "WS.f.control": "MeteoObs",  # Windsnelheid verwachting [m/s] (controle)
    "WS.f.determ": "MeteoObs",  # Windsnelheid verwachting [m/s] (deterministisch)
    "WS.f.ensem": "MeteoObs",  # Windsnelheid verwachting [m/s] (ensemble)
    "WSU10.f.harmonie": "MeteoObs",  # Windsnelheid verwachting u-component [m/s] harmonie
    "WSU10.f.hirlam": "MeteoObs",  # Windsnelheid verwachting u-component [m/s] (hirlam)
    "WSV10.f.harmonie": "MeteoObs",  # Windsnelheid verwachting v-component [m/s] harmonie
    "WSV10.f.hirlam": "MeteoObs",  # Windsnelheid verwachting v-component [m/s] (hirlam)
    "H.B": "WaterlvlObs",  # Waterhoogte berekend [mNAP]
    "H.B.d": "WaterlvlObs",  # Waterhoogte berekend [mNAP] - dag
    "H.B.h": "WaterlvlObs",  # Waterhoogte berekend [mNAP] - uur
    "H.G": "WaterlvlObs",  # Waterhoogte [mNAP]
    "H.G.AA.d": "WaterlvlObs",  # Waterhoogte alle [-] - dag
    "H.G.d": "WaterlvlObs",  # Waterhoogte [mNAP] - dag
    "H.G.h": "WaterlvlObs",  # Waterhoogte [mNAP] - uur
    "H.Max": "WaterlvlObs",  # Maximum waterpeil [mNAP]
    "H.S": "WaterlvlObs",  # Streefpeil [mNAP]
    "H.S.h": "WaterlvlObs",  # Streefpeil [mNAP] - uur
    "H.f": "WaterlvlObs",  # Waterhoogte verwachting [mNAP]
    "H.f.determ": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (ecmwf-deterministisch)
    "H.f.ecmwf": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (ecmwf-ensemble)
    "H.f.fr": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (fews rivieren)
    "H.f.matroos": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (matroos)
    "H.f.rws": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (rws)
    "H.f.wdij": "WaterlvlObs",  # Waterhoogte verwachting [mNAP] (wdij)
    "H.schuif.G": "WaterlvlObs",  # Schuifstand [mNAP]
    "H.wave": "WaterlvlObs",  # Golfhoogte [m]
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
    "Hwk": "WaterlvlObs",  # Hoogte waterkolom [m]
    "Hwk.H1": "WaterlvlObs",  # Overstortende straal H1 [m]
    "Hwk.H2": "WaterlvlObs",  # Hoogteverschil na stuw H2 [m]
    "Hwk.PA.G.Barodiver": "WaterlvlObs",  # Hoogteverschil na stuw H2 [m]
    "Hwk.PW.G": "WaterlvlObs",  # Gemeten waterdruk [m]
    "Hwk.V.PB": "WaterlvlObs",  # Waterstandsafwijking PB [m]
    "Hwk.V.SP": "WaterlvlObs",  # Waterstandsafwijking SP [m]
    # "Q.B": "WaterlvlObs",  # Debiet FEWS [m3/s]
    # "Q.B.d": "WaterlvlObs",  # Debiet FEWS [m3/s] - dag
    # "Q.B.h": "WaterlvlObs",  # Debiet FEWS [m3/s] - uur
    # "Q.G": "WaterlvlObs",  # Debiet [m3/s]
    # "Q.G.d": "WaterlvlObs", # Debiet [m3/s] - dag
    # "Q.G.h": "WaterlvlObs",  # Debiet [m3/s] - uur
    # "Q.f": "WaterlvlObs",  # Debiet verwachting [m3/s]
    # "Q.f.determ": "WaterlvlObs",  # Debiet verwachting [m3/s] (ecmwf-deterministisch)
    # "Q.f.ecmwf": "WaterlvlObs",  # Debiet verwachting [m3/s] (ecmwf-ensemble)
    # "Q.f.fr": "WaterlvlObs",  # Debiet verwachting [m3/s] (fews rivieren)
    # "Q.f.rws": "WaterlvlObs",  # Debiet verwachting [m3/s] (rws)
    # "Q.f.totaal": "WaterlvlObs",  # Debiet verwachting [m3/s] (totaal)
    # "Qafval.G": "WaterlvlObs",  # Debiet afvalwater [m3/h]
    # "Qafval.G.d": "WaterlvlObs",  # Debiet afvalwater [m3/h] - dag
    # "Qafval.G.h": "WaterlvlObs",  # Debiet afvalwater [m3/h] - uur
    "T.ow.G": "WaterQualityObs",  # Temperatuur oppervlaktewater [C]
    "O2.gehalte": "WaterQualityObs",  # Zuurstof gehalte [mg/l]
    "O2.saturatie": "WaterQualityObs",  # Zuurstof saturatie [%]
    "PA.G": "MeteoObs",  # Luchtdruk gemeten [hPa]
    "PA.G.bar": "MeteoObs",  # Luchtdruk gemeten [bar]
    "PA.f.harmonie": "MeteoObs",  # Luchtdruk verwachting [hPa] (harmonie)
    "PKEL.V": "MeteoObs",  # Drukverschil Keller (P1-P2)
    "RH.G": "MeteoObs",  # Luchtvochtigheid [%]
    "RH.f.harmonie": "MeteoObs",  # Relatieve luchtvochtigheid [%] (harmonie)
    "T.air.G": "MeteoObs",  # Temperatuur lucht [C]
    "T.air.f.control": "MeteoObs",  # Temperatuur verwacht (controle)
    "T.air.f.determ": "MeteoObs",  # Temperatuur verwacht (deterministisch)
    "T.air.f.ensem": "MeteoObs",  # Temperatuur verwacht (ensemble)
    "T.air.f.harmonie": "MeteoObs",  # Temperatuur verwacht (harmonie)
    # "PW.G": "Waterspanning [hPa]",
    # "PZ.G": "Zuigspanning [kPa]",
    # "SM.G": "Bodemvochtigheid [m3/m3]",
    # "T.soil.G": "Temperatuur bodem [C]",
    "SS.G": "WaterlvlObs",  # Stroomsnelheid gemeten [m/s]
    "SS.G.d": "WaterlvlObs",  # Stroomsnelheid gemeten [m/s] - dag
    "SS.G.h": "WaterlvlObs",  # Stroomsnelheid gemeten [m/s] - uur
    "GW.G": "GroundwaterObs",  # Grondwaterstand [mNAP]
    "GW.G.hm": "WaterlvlObs",  # Grondwaterstand - handmeting [mNAP]
    "SH": "GroundwaterObs",  # Stijghoogte [mNAP]
    "T.gw.G": "WaterQualityObs",  # Temperatuur grondwater [C]
}

hhnk_pid = {
    "windrichting.meting": "MeteoObs",  # windrichting gemeten [graden]
    "windrichting.rws": "MeteoObs",  # windrichting voorspeld rws [graden]
    "windrichting.voorsp.det": "MeteoObs",  # windrichting voorspeld (det) [graden]
    "windrichting.voorsp.ens": "MeteoObs",  # windrichting voorspeld (ens) [graden]
    "Rns.voorsp": "MeteoObs",  # netto kortgolvige straling (Joules/m2) [J/m2]
    "LH.voorsp": "MeteoObs",  # Latente warmte (Joules/m2) [J/m2]
    "windsnelheid.meting": "MeteoObs",  # windsnelheid gemeten [m/s]
    "Wind.voorsp.u": "MeteoObs",  # voorspelde windsnelheid (u-horiz.) [m/s]
    "Wind.voorsp.v": "MeteoObs",  # voorspelde windsnelheid (v-vert.) [m/s]
    "wind.meting.max": "MeteoObs",  # windsnelheid gemeten max [m/s]
    "windsnelheid.rws": "MeteoObs",  # windsnelheid voorspeld rws [m/s]
    "windsnelheid.voorsp.det": "MeteoObs",  # windsnelheid voorspeld (det) [m/s]
    "windsnelheid.voorsp.ens": "MeteoObs",  # windsnelheid voorspeld (ens) [m/s]
    "WNS9026": "MeteoObs",  # luchtdruk [hPa]",
    "GELDHD.meting": "WaterQualityObs",  # geleidbaarheid (EGV) [mS/m]
    "GELDHD.meting.hand": "WaterQualityObs",  # geleidbaarheid (EGV) handmeting [mS/m]
    "EGV.meting": "WaterQualityObs",  # geleidbaarheid (EGV) [mS/cm]
    "EGVms_cm.meting": "WaterQualityObs",  # geleidbaarheid (EGV mS/cm) [mS/cm]
    "EGVms.meting": "WaterQualityObs",  # geleidbaarheid (EGV mS) [mS]
    "EGVms_m.meting": "WaterQualityObs",  # geleidbaarheid (EGV mS/m) [mS/m]
    "EGVms_m.meting.hand": "WaterQualityObs",  # geleidbaarheid handmeting (EGV mS/m) [mS/m]
    "O2.geh.meting": "WaterQualityObs",  # zuurstofgehalte [mg/l]
    "O2.sat.meting": "WaterQualityObs",  # zuurstofsaturatie [%]
    "CL.berekend": "WaterQualityObs",  # chloridegehalte berekend [mg/l]
    "T.water.meting": "WaterQualityObs",  # watertemperatuur [oC]
    "T.water.meting.5min": "WaterQualityObs",  # watertemperatuur [oC]
    "T.water.meting.hand": "WaterQualityObs",  # watertemperatuur handmeting [oC]
    "GW.meting": "GroundWaterObs",  # grondwaterhoogte gemeten [m]
    "H.meting.ref": "GroundwaterObs",  # cm tov referentiehoogte [cm]
    "H.meting.5min": "GroundWaterObs",  # grondwaterhoogte gemeten [m]
    "GW.meting.1%": "GroundWaterObs",  # 1% percentiel grondwater [m]
    "GW.meting.10%": "GroundWaterObs",  # 10% percentiel grondwater [m]
    "GW.meting.90%": "GroundWaterObs",  # 90% percentiel grondwater [m]
    "GW.meting.25%": "GroundWaterObs",  # 25% percentiel grondwater [m]
    "GW.meting.75%": "GroundWaterObs",  # 75% percentiel grondwater [m]
    "GW.meting.99%": "GroundWaterObs",  # 99% percentiel grondwater [m]
    "H.berekend": "WaterlvlObs",  # waterhoogte berekend [m]
    "H.streef": "WaterlvlObs",  # streefpeil [m]
    "H.streef.boven": "WaterlvlObs",  # bovengrens streefpeil [m]
    "H.streef.onder": "WaterlvlObs",  # ondergrens streefpeil [m]
    "H.meting": "WaterlvlObs",  # waterhoogte [m]
    "H.meting.hand": "WaterlvlObs",  # waterhoogte handmeting [m]
    "H.meting.org": "WaterlvlObs",  # waterhoogte origineel [m]
    "H.meting.int": "WaterlvlObs",  # waterhoogte geinterpoleerd [m]
    "H.meting.gebied": "WaterlvlObs",  # gebiedsgemiddelde waterhoogte [m]
    "H.meting.1%": "WaterlvlObs",  # 1% percentiel waterhoogte [m]
    "H.meting.10%": "WaterlvlObs",  # 10% percentiel waterhoogte [m]
    "H.meting.90%": "WaterlvlObs",  # 90% percentiel waterhoogte [m]
    "H.meting.25%": "WaterlvlObs",  # 25% percentiel waterhoogte [m]
    "H.meting.75%": "WaterlvlObs",  # 75% percentiel waterhoogte [m]
    "H.meting.99%": "WaterlvlObs",  # 99% percentiel waterhoogte [m]
    "H.meting.dag": "WaterlvlObs",  # waterhoogte dagwaarden [m]
    "H.astronomisch": "WaterlvlObs",  # getij astronomisch [m]
    "H.getij.meting": "WaterlvlObs",  # getij gemeten [m]
    "H.getij.meting.opzet": "WaterlvlObs",  # getij opzet gemeten [m]
    "H.getij.hmcn": "WaterlvlObs",  # getij voorspeld hmcn [m]
    "H.getij.rws": "WaterlvlObs",  # getij voorspeld rws [m]
    "H.getij.wdij": "WaterlvlObs",  # getij voorspeld wdij [m]
    "H.getij.opzet.rws": "WaterlvlObs",  # getij opzet voorspeld rws [m]
    "WNS9688": "WaterlvlObs",  # waterhoogte [m]
    # "inslagpeil": "2e inslagpeil [m]
    # "Bov.buis": "bovenkant buis [m]",
    # "Maaiveld": "maaiveld [m]",
    "P.meting.1m": "PrecipitationObs",  # neerslag gemeten (1 min) [mm]
    "P.meting.5m": "PrecipitationObs",  # neerslag gemeten (5 min) [mm]
    "P.meting.10m": "PrecipitationObs",  # neerslag gemeten (10 min) [mm]
    "P.meting.1h": "PrecipitationObs",  # neerslag gemeten (1 uur) [mm]
    "P.meting.24h": "PrecipitationObs",  # neerslag gemeten (24 uur) [mm]
    "P.radar.2m": "PrecipitationObs",  # neerslagradar (2 min) [mm/5min]
    "P.radar.5m": "PrecipitationObs",  # neerslagradar (5 min) [mm/5min]
    "P.radar.hist": "PrecipitationObs",  # neerslag (ongecalibreerde radar) [mm/hr]
    "P.radar.cal": "PrecipitationObs",  # neerslag (gecalibreerde radar) [mm/hr]
    "P.voorsp": "PrecipitationObs",  # neerslag voorspeld [mm/hr]
    "P.voorsp.ens": "PrecipitationObs",  # neerslag voorspeld (ens) [mm/hr]
    "P.voorsp.ens.min": "PrecipitationObs",  # neerslag voorspeld (ens) min [mm/hr]
    "P.voorsp.ens.max": "PrecipitationObs",  # neerslag voorspeld (ens) max [mm/hr]
    "P.voorsp.ens.median": "PrecipitationObs",  # neerslag voorspeld (ens) mediaan [mm/hr]
    "P.voorsp.ens.10": "PrecipitationObs",  # neerslag voorspeld (ens) 10 [mm/hr]
    "P.voorsp.ens.90": "PrecipitationObs",  # neerslag voorspeld (ens) 90 [mm/hr]
    "P.voorsp.ctr": "PrecipitationObs",  # neerslag voorspeld (ctr) [mm/hr]
    "P.voorsp.det": "PrecipitationObs",  # neerslag voorspeld (det) [mm/hr]
    "P.radar.1h": "PrecipitationObs",  # neerslagradar (1 uur) [mm/hr]
    "P.radar.24h": "PrecipitationObs",  # neerslagradar (24 uur) [mm/24hr]
    "E.meting": "EvaporationObs",  # verdamping berekend KNMI (Makkink) [mm]
    "E.ref.Makkink": "EvaporationObs",  # verdamping voorspeld (Makkink gewas ref) [mm]
    "Neerslagoverschot": "MeteoObs",  # Berekende neerslagoverschot (mm/dag) [mm]
    "Neerslagoverschot_cumulatief": "MeteoObs",  # Berekende cumulatieve neerslagoverschot (mm) [mm]
    "T.lucht.meting": "MeteoObs",  # luchttemperatuur [oC]
    "T.meting": "MeteoObs",  # luchttemperatuur [oC]
    "T.voorsp": "MeteoObs",  # temperatuur verwacht [oC]
    "T.voorsp.ens": "MeteoObs",  # temperatuur voorspeld (ens) [oC]
    "T.voorsp.det": "MeteoObs",  # temperatuur voorspeld (det) [oC]
    # "T.golf.meting": "golfperiode [-]",
    # "Hz.meting": "frequentie [Hz]",
    # "H.golf.meting": "golfhoogte [m]",
    # "Stuw.klep.meting": "klephoogte [m]",
    # "H.meting.proc": "vullingsgraad [%]",
    # "H.meting.percentage": "niveau percentage [%]",
    # "Druk.meting": "druk gemeten [kPa]",
    # "P.druk.meting": "druk [mBar]",
    # "WNS2368": "debiet [m3/h]",
    # "WNS2369.sim": "debiet voorspeld [m3/h]",
    # "WNS2369.h": "verpompt volume uur [m3/h]",
    # "WNS2369.h.pred": "verpompt volume uur voorspeld [m3/h]",
    # "WNS2369.h.upp": "verpompt volume uur voorspeld bovengrens [m3/h]",
    # "WNS2369.h.low": "verpompt volume uur voorspeld ondergrens [m3/h]",
    # "Q.meting": "debiet gemeten [m3/s]",
    # "Q.berekend": "debiet berekend [m3/s]",
    # "Q.berekend.totNZK": "debiet totaal NZK [m3/s]",
    # "Q.berekend.totNZK.1h": "debiet totaal NZK (1 uur) [m3/s]",
    # "Q.advies": "debiet advies [m3/s]",
    # "WNS2369.2min": "verpompt volume 2min [m3/2min]",
    # "WNS2369.1min": "verpompt volume min [m3/min]",
    # "Q.berekend.dag": "debiet berekend dag [m3/dag]",
    # "WNS2369.d": "verpompt volume dag [m3/dag]",
    # "AfwijkingInslagpeil": "afwijking peil [min]",
    # "AfwijkingNorm_kwartaal": "tijd niet voldaan aan afnameverplichting [min]",
    # "AfwijkingNorm": "voldoet aan afnameverplichting [min]",
    # "AfwijkingInslagpeil_kwartaal": "tijd boven 2e inslagpeil [min]",
    # "AfwijkingRVW_dag": "Percentage RVW >50% over 24 uur [%]",
    # "AfwijkingRVW_dag_nat": "Percentage RVW >50% over 24 uur nat [%]",
    # "AfwijkingRVW_dag_droog": "Percentage RVW >50% over 24 uur droog [%]",
}

wli_pid = {
    "WATDTE": "GroundwaterObs",  # Waterdiepte [m]
    "DRUK": "GroundwaterObs",  # Drukhoogte [m]
    "STIJGHTE": "GroundwaterObs",  # Stijghoogte [m]
    "STIJGHTE_sim": "GroundwaterObs",  # Stijghoogte gesimuleerd [m]
    "STIJGHTE_res": "GroundwaterObs",  # Resultaat toetsing PASTAS
    "WATHTE": "WaterlvlObs",  # Waterhoogte [m]
    # "HEFHTE": "Hefhoogte [m]",
    # "KRUINHTE": "Kruinhoogte [m]",
    # "STREEFHTE": "Streefhoogte [m]",
    # "MAAIVELDHTE": "Maaiveldhoogte [m]",
    # "KRUINREL": "Kruinoogte [%]"
    # "STROOMSHD": "Stroomsnelheid [m/s]",
    # "Q": "Afvoer [m3/s]",
    # "Q_res" : "Resultaat toetsing PASTAS Afvoer res"
    # "STREEFQ": "Afvoer streefwaarde [m3/s]",
    # "Q_sim": "Afvoer gesimuleerd [m3/s]",
    # "Q_sim.fc": "Afvoer voorspeld is gemiddelde van voorspellingen (sim.fc) [m3/s]",
    # "Q_sim_q05.fc": "Afvoer gesimuleerd 5%-waarde of ondergrens 90% betrouwbaarheidsinterval [m3/s]",
    # "Q_sim_q25.fc": "Afvoer gesimuleerd 25%-waarde of ondergrens 50% betrouwbaarheidsinterval [m3/s]",
    # "Q_sim_q75.fc": "Afvoer gesimuleerd 75%-waarde of bovengrens 50% betrouwbaarheidsinterval [m3/s]",
    # "Q_sim_q95.fc": "Afvoer gesimuleerd 95%-waarde of bovengrens 90% betrouwbaarheidsinterval [m3/s]",
    "GELDHD": "WaterQualityObs",  # Geleidbaarheid (MON) [uS/cm]
    "O2SAT": "WaterQualityObs",  #  Zuurstof saturatie [%]
    "O2": "WaterQualityObs",  # Zuurstofgehalte [mg/l]
    "NEERSG": "PrecipitationObs",  # Neerslag [mm]
    "NEERSG.tekort": "MeteoObs",  # Neerslagtekort [mm]
    "NEERSG.voorsp.det": "PrecipitationObs",  # Neerslag voorspeld deterministisch [mm]
    "NEERSG.voorsp.ens": "PrecipitationObs",  # Neerslag voorspeld ensemble [mm]
    "VERDPG": "EvaporationObs",  # Verdamping [mm]
    "LUCHTDK": "MeteoObs",  # Luchtdruk [m]
    "WINDRTG": "MeteoObs",  # Wind richting [degrees]
    "WINDSHD": "MeteoObs",  # Windhsnelheid [m/s]
    "GLOBSLG": "MeteoObs",  # Globale straling [Joules/cm2]
    "T": "MeteoObs",  # Temperatuur [C]
    "T.voorsp.det": "MeteoObs",  # Temperatuur voorspeld deterministisch [C]
    "T.voorsp.ens": "MeteoObs",  # Temperatuur voorspeld ensemble [C]
    "RELTVLVTHD": "MeteoObs",  # Relatieve luchtvochtigheid [%]
}

pid = {"wsvv": wsvv_pid, "hhnk": hhnk_pid, "wli": wli_pid}
